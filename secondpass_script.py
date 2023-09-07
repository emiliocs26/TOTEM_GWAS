#!/bin/python
import numpy as np
import pandas as pd
import h5py
import itertools
from auxillary import *
import matplotlib.pyplot as plt
from tqdm import tqdm
from itertools import combinations as com
from pathos.multiprocessing import ProcessingPool as Pool
import random
import sys

def Shannon(p) -> float:
    p = p[p>0]                         # microstates with zero probability do not contribute to entropy
    return np.sum(-p*np.log(p))

def constraint_vector(n: int, option: str) -> pd.Series:
    """
    builds the row vector for a given constraint
    
    __option__
    define if a constraint is global or SNP-specific
    "global"   = constraint vector for a global constraint
    "SNP-x"    = constraint, specific for first SNP
    "SNP-y"    = constraint, specific for first SNP
    "both"     = constraint for both SNPs
    
    __n__
    define which moment is being coded 
    0 = normalization/ snp-specific prevalence
    1 = global mean / snp-specific mean
    2 = global variance / snp-specifc variance
    3 = global skewness / snp-sepcific skewness
    4 = global kurtosis / snp-specific kurtosis
    etc.
    
    """
    states = microstates.copy()
                                                     
    if option == 'global': 
        return np.power(((states[phenotype] - global_mean) / global_std), n)  
    
    else:
        constraint = np.zeros(dimS)
      
        if 'SNP' in option: #if SNP-x or SNP-y
            assert option in states.columns, "option must either be \"SNP-x\", \"SNP-y\" or \"both\""
            mask = states[option]
            
        elif option == 'both':
            SNP_names = [col for col in states.columns if 'SNP' in col]
            assert len(SNP_names) == 2, "The microstatespace DataFrame must only contain two SNP columns"
            mask = np.logical_and(states[SNP_names[0]], states[SNP_names[1]])
        else:
            raise RuntimeError("Unknown type of constraint requested")
        
        constraint[mask] = np.power((states.loc[mask, phenotype] - global_mean) / global_std, n)  
        
        return constraint



def combine_SNPs(SNP1:int, SNP2:int) -> pd.DataFrame:
    """
    combines two SNPs from the h5py file along the columns of a dataframe indexed by the accessions

    """
    
    two_SNP_frame = SNP_data[[SNP1, SNP2]]
    two_SNP_frame.columns = "SNP-x", "SNP-y"
    
    return two_SNP_frame

def combine_pheno_genome(two_SNP_frame: pd.DataFrame) -> pd.DataFrame:
    """
    returns one dataframe with counts and relative frequencies of all possible microstates indexed by the phenotype and the two SNPs
    """
    two_SNP_frame = two_SNP_frame.loc[selected_accessions]
    pheno_SNPs = two_SNP_frame.merge(FT_data, left_index=True, right_index=True, how='outer')   
    assert not pheno_SNPs.isnull().values.any(), "pheno-SNP dataframe has NaN (possibly combining pheno with SNP data failed)"

    # merge with reference microstates to get the full dataframe over all possible (SNP-x, SNP-y, phenotype)
    pheno_SNPs = pd.concat([microstates, pheno_SNPs])                                                
    pheno_SNPs = pheno_SNPs.groupby(['SNP-x', 'SNP-y', phenotype]).sum('count')
    assert len(pheno_SNPs) == len(microstates), "merging the pheno-genome dataframe with the microstates failed"

    pheno_SNPs['f'] = pheno_SNPs['count'] / N
    pheno_SNPs.N = N
    
    return pheno_SNPs


def entropy_difference(SNP_pair: tuple) -> dict:
    """
    Calculate the entropy difference between p0 and p1

    """
    
    two_SNP_frame = combine_SNPs(*SNP_pair)

    #joint SNP effects can only be inferred, if all 4 genotype combinations are observed, get rid of the rest
    
    if len(two_SNP_frame.drop_duplicates()) < 4: 
        return 0
    
    pheno_genome = combine_pheno_genome(two_SNP_frame)               # get the microstate counts

    f_vector = pheno_genome['f']                                     # get f_vector 

    # calculate moments

    if "global_moments" not in globals():

        snp_moments = C_SNP@f_vector

        global global_moments                          

        global_moments = C_global@f_vector

        # we also create a start point for the Newton Sovler, increase stability and speed of Newton's method

        global p_start

        p_start =  Newton_solver(C_global, global_moments, conv_tolerance = conv_tolerance)[0]
        
        moments = np.concatenate((global_moments, snp_moments))
        
    else:
     
        snp_moments = C_SNP@f_vector

        moments = np.concatenate((global_moments, snp_moments))
        
        
        
    C1       = C
    moments1 = moments

    C0       = C[:-n_p1_specific_constraints]
    moments0 = moments1[:-n_p1_specific_constraints]

    try:
        p0 = Newton_solver(C0, moments0, p0 = p_start, conv_tolerance = conv_tolerance)[0]       
        p1 = Newton_solver(C1, moments1, p0 = p_start, conv_tolerance = conv_tolerance)[0]       
        
        H0 = Shannon(p0)
        H1 = Shannon(p1)

        return H0 - H1
    
    except ConvergenceError:
        return np.nan



##### data import + formatting #####

annotated_snps = pd.read_csv("./results/annotated_snps.csv", sep = ";")

genotype_data = h5py.File('./data/2029_snps.hdf5', 'r') #genotype_data
accessions = np.array([int(ID.decode("utf-8")) for ID in genotype_data['accessions']])
    
phenotype = "FT"
    
for pheno_data in ["ft10"]:
    
    if pheno_data == "ft10":
        FT_data = pd.read_csv('./data/ft10_final.csv', index_col = "id")
        FT_data.drop("tg_id", axis = 1, inplace = True)
        FT_data = pd.concat([FT_data[col].rename(phenotype) for col in FT_data.columns]).to_frame()      
        FT_data = FT_data.dropna(axis=0) 
    
    if pheno_data == "ft16":
        FT_data = pd.read_csv("./data/ft16_final.csv", index_col = "ID")
        FT_data.columns = [phenotype]
    
    
    selected_accessions = np.array(list(set(accessions).intersection(set(FT_data.index))), dtype=int)

    SNP_data = pd.DataFrame(np.array(genotype_data["snps"][:]), columns = accessions)
    SNP_data = SNP_data[SNP_data.columns[np.isin(accessions, FT_data.index.values)]]

    assert np.isin(selected_accessions, list(SNP_data.columns)).all()
    assert len(SNP_data.columns) == len(selected_accessions)

    SNP_data = SNP_data.T
    print("finished loading genomic and phenotype info")

    FT_data = FT_data.loc[selected_accessions]  # only keep phenotypes from accessions for which genomic info is available            
    FT_data['count'] = 1

    N = FT_data["count"].sum()

    # get global statistics
    global_mean = FT_data[phenotype].mean()
    global_std = FT_data[phenotype].std()
    

    ##### settings #####    


    phenotype_min, phenotype_max = 30, 365                                             
    pheno_states = np.arange(phenotype_min,phenotype_max+1)

    # create microstates for the attributes: phenotype, SNP-x and SNP-y       
    SNP_states = [False, True]                                                
    microstates = pd.DataFrame(itertools.product(SNP_states, SNP_states, pheno_states), columns = ['SNP-x', 'SNP-y', phenotype])   
    microstates['count'] = 0
    dimS = len(microstates)
    
    conv_tolerance = 10**-12


    # build constraint matrix

    C_global = np.array([constraint_vector(n = 0, option = "global"),
                         constraint_vector(n = 1, option = "global")
                         ])

    C_SNP   = np.array([constraint_vector(0, option='SNP-x'),
                        constraint_vector(0, option='SNP-y'),
                        constraint_vector(1, option='SNP-x'),
                        constraint_vector(1, option='SNP-y'),
                        constraint_vector(0, option='both'),
                        constraint_vector(1, option='both')
                        ])

    C = np.vstack((C_global, C_SNP))

    assert C.shape == ((len(C_global) + len(C_SNP)), dimS), "building constraint matrix failed"

    n_p1_specific_constraints = 1

    
    ############# run the code on all selected SNP-pairs############# 

    entropy_difference((218822, 3782584)) #run the code for any random SNP once, so that global_moments is defined in all pathos processes
    
    
    firstpass_results = pd.read_csv(f"./results/{pheno_data}/firstpass/results_{pheno_data}_globalmean.csv")
    firstpass_results.sort_values("entropy_difference", inplace = True, ascending = False)
    
    ranked_genes = annotated_snps.loc[firstpass_results["hdf5_index"]].drop_duplicates("gene_code", keep = "first")

    print("finished loading and formatting first pass results")
    
    snp_selection = ranked_genes.index.values
    
    print("calculating pairs...")
    
    snp_pairs = [x for x in tqdm(com(snp_selection, 2))]
    
    print("shuffeling pairs...")
    
    random.shuffle(snp_pairs)
    
    snp_ones  = [x[0] for x in snp_pairs]
    snp_twos  = [x[1] for x in snp_pairs]
    
    chunksize = 5000000
    
    counter = 1
    
    for i in tqdm(range(0, len(snp_pairs), chunksize)):
    
        with Pool() as pool:
                entropy_differences = pool.map(entropy_difference, snp_pairs[i:i+chunksize])
                pool.clear()
                
    
        results = pd.DataFrame({"SNP1":snp_ones[i:i+chunksize],
                               "SNP2":snp_twos[i:i+chunksize],
                               "byproxy1":ranked_genes.loc[snp_ones[i:i+chunksize], "byproxy"].reset_index(drop = True),
                               "byproxy2":ranked_genes.loc[snp_twos[i:i+chunksize], "byproxy"].reset_index(drop = True),
                               "gene_code1":ranked_genes.loc[snp_ones[i:i+chunksize], "gene_code"].reset_index(drop = True),
                               "gene_code2": ranked_genes.loc[snp_twos[i:i+chunksize], "gene_code"].reset_index(drop = True),
                               "cd1":ranked_genes.loc[snp_ones[i:i+chunksize], "computational_description"].reset_index(drop = True),
                               "cd2": ranked_genes.loc[snp_twos[i:i+chunksize], "computational_description"].reset_index(drop = True),
                               "entropy_difference":entropy_differences
                              })

        results.to_csv(f"./results/{pheno_data}/secondpass/results_{pheno_data}_secondpass_part{counter}.csv", index = False)
    
        counter += 1


