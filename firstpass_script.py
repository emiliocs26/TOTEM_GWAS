import time

start = time.perf_counter()

import numpy as np
import pandas as pd
import h5py
from auxillary import *
import sys
import time
import matplotlib.pyplot as plt
from pathos.multiprocessing import ProcessingPool as Pool

def Shannon(p) -> float:
    p = p[p>0]                         # microstates with zero probability do not contribute to entropy
    return np.sum(-p*np.log(p))

def constraint_vector(option, n):
    """
    builds the row vector for a given constraint

    __option__
    define if it is a SNP-specific constraint or a global constraint
    "global" = constraint vector for a global constraint
    "SNP"    = constraint vector for a SNP specific constraint

    __n__
    define which moment is being coded 
    0 = normalization/snp-specific prevalence
    1 = global mean /snp-specific mean
    2 = global variance / snp-specific variance
    3 = glboal skewness / snp-specific skewness
    4 = global kurtosis / snp-specific kurtosis
    etc.
    """

    if option == "global":
        constraint = np.power((statespace["FT"] - global_mean) / np.sqrt(global_var), n)

    if option == "SNP":
        constraint = statespace["FT"].copy()
        constraint[SNP_off_mask] = 0
        constraint.loc[SNP_on_mask] = np.power((constraint[SNP_on_mask] - global_mean) / np.sqrt(global_var), n)

    return constraint

def snp_moments(SNP):
    """
    builds the f-vector and calculates the snp-specific moments/prevalence

    Since we only need the moments for our further calculation, this function returns only the moments by default
    """

    # we define "off" as 0 and "on" as 1

    # devide the FT_data into 2 dataframes

    FT_data_SNP_off = FT_data.iloc[np.where(SNP_data.loc[SNP] == 0)]
    FT_data_SNP_on = FT_data.iloc[np.where(SNP_data.loc[SNP] == 1)]

    # phenotype value counts

    values_off, counts_off = np.unique(FT_data_SNP_off, return_counts = True) 
    values_on, counts_on  = np.unique(FT_data_SNP_on, return_counts = True)

    # get rid of nan-values

    if len(values_off) != 0: # preventing KeyErrors

        if np.isnan(values_off[-1]):
            values_off, counts_off = values_off[:-1], counts_off[:-1]

    if len(values_on) != 0:
        if np.isnan(values_on[-1]):
            values_on, counts_on   = values_on[:-1], counts_on[:-1]


    # build the f-vector

    f_vector = np.zeros((phenotype_max+1 - phenotype_min)*2)

    f_vector[SNP_off_mask & (np.isin(statespace["FT"], values_off))] = counts_off
    f_vector[SNP_on_mask & (np.isin(statespace["FT"], values_on))]  = counts_on

    f_vector = f_vector/N

    # since the global momemts never change, we can calculate them only once to lower the computational burden

    if "global_moments" not in globals():

        moments = C@f_vector

        global global_moments                          

        global_moments = moments[:n_global_constraints]

        # we also create a start point for the Newton Sovler, increase stability and speed of Newton's method

        global p_start

        p_start =  Newton_solver(C[:n_global_constraints], moments[:n_global_constraints], conv_tolerance = conv_tolerance)[0]


    else:
        snp_moments = C[n_global_constraints:]@f_vector

        moments = np.concatenate((global_moments, snp_moments))


    return moments

def entropy_difference(SNP):
    """
    estimate p0 and p1 and calculate the entropy difference
    """
    
    if macs[SNP] == 0:            #SNPs that don't appear get a defined entropy difference of 0
        return 0

    C1       = C
    moments1 = snp_moments(SNP)

    C0       = C[:-n_p1_specific_constraints]
    moments0 = moments1[:-n_p1_specific_constraints]

    try:
        p0 = Newton_solver(C0, moments0, p0 = p_start, conv_tolerance = conv_tolerance)[0]        
        p1 = Newton_solver(C1, moments1, p0 = p_start, conv_tolerance = conv_tolerance)[0]                          

        H0 = Shannon(p0)
        H1 = Shannon(p1)

        return H0 - H1
    
    except ConvergenceError: # this happens only in rare cases, for SNPs that only appear in one or two accessions
        return np.nan
    

def get_Minor_Allele_Count(SNP):

    row_sum = SNP_data.loc[SNP].sum()
    value = 1 if (row_sum) < (n_accessions/2) else 0

    if value == 1:
        return row_sum
    else:
        return n_accessions-row_sum

##### data import  + settings #####

genotype_data = h5py.File('./data/2029_snps.hdf5', 'r') 
accessions = np.array([int(ID.decode("utf-8")) for ID in genotype_data['accessions']])


for pheno_data in ["ft16","ft10"]:     #run the script for FT16 and FT10

    if pheno_data == "ft10":

        FT_data = pd.read_csv('./data/ft10_final.csv', index_col = "id")
        FT_data.drop("tg_id", axis = 1, inplace = True)
        FT_data = FT_data.loc[accessions[np.isin(accessions, FT_data.index.values)]]

    if pheno_data == "ft16":

        FT_data = pd.read_csv('./data/ft16_final.csv')
        FT_data = FT_data[np.isin(FT_data["ID"] , accessions)]
        FT_data['column_number'] = FT_data.groupby("ID").cumcount()
        FT_data = FT_data.pivot(index='ID', columns='column_number', values='raw-FT16')

        acc_indices = []

        for accession in accessions:
            acc_index = np.where(FT_data.index == accession)[0]

            if len(acc_index) == 0:
                pass
            else:
                acc_indices.append(acc_index[0])

        FT_data = FT_data.iloc[acc_indices]

    SNP_data = pd.DataFrame(np.array(genotype_data["snps"]), columns = accessions)
    SNP_data = SNP_data[SNP_data.columns[np.isin(accessions, FT_data.index.values)]]

    n_snps       = len(SNP_data)
    n_accessions = len(SNP_data.columns)

    assert all(SNP_data.columns == FT_data.index), "accessions are not listed in the same order in SNP_data and FT_data"
    
    ##### get sample size, global moments and minor allele counts #####
    
    N = (~np.isnan(FT_data)).sum().sum()
    global_mean = np.nanmean(FT_data)
    global_var = np.nanvar(FT_data)
    
    with Pool() as pool:
        macs = pool.map(get_Minor_Allele_Count, range(n_snps))
        pool.clear()

    ##### defining a microstatespace #####

    phenotype_min = 30                                                                    # <- option
    phenotype_max = 365                                                                   # <- option

    statespace = pd.DataFrame(itertools.product(np.arange(phenotype_min,phenotype_max+1), ["on", "off"]), columns = ["FT", "SNP"])

    SNP_off_mask = statespace["SNP"] == "off"
    SNP_on_mask  = statespace["SNP"] == "on"


    # building the constraint matrix
    # p1 specific constraints must be on the bottom of the C Matrix

    C_global = np.array([constraint_vector(option = "global", n = 0),
                         constraint_vector(option = "global", n = 1),                   
                         constraint_vector(option = "global", n = 2)                
                        ])
                                                                                       # <- option
    C_SNP    = np.array([constraint_vector(option = "SNP", n = 0),
                         constraint_vector(option = "SNP", n = 1),
                         constraint_vector(option = "SNP", n = 2)
                        ])


    n_p1_specific_constraints = 2                                                       # <- option

    C = np.vstack((C_global, C_SNP))


    n_global_constraints = len(C_global)


    last_global_constraint = ["globalmean", "globalvariance", "globalskewness", "globalkurtosis"]

    last_global_constraint = last_global_constraint[n_global_constraints-2]
    
    
    conv_tolerance = 10**-12 #for Newton's method                                       #<- option
        
    
    ###### running th code on all SNPs ######
    
    entropy_difference(665) # run the code for any random SNP once, so that global_moments is defined in all pathos processes
    
    
    with Pool() as pool:
            entropy_differences = pool.map(entropy_difference, range(n_snps))
            pool.clear()
    
    ###### generate a DataFrame with more information on each SNP ######

    positions_array = np.zeros(len(genotype_data["snps"])).astype(int) #bp position
    chrs_array      = np.zeros(len(genotype_data["snps"])).astype(int) #chromosomal position

    for i in range(5):

        indices_for_chr   = genotype_data["positions"].attrs["chr_regions"][i]

        positions_array[indices_for_chr[0]:indices_for_chr[1]] = genotype_data["positions"][indices_for_chr[0]:indices_for_chr[1]]

        chrs_array[indices_for_chr[0]:indices_for_chr[1]] = i+1

    
    results = pd.DataFrame({"hdf5_index":range(n_snps),
                            "entropy_difference":entropy_differences,
                            "chromosome": chrs_array,
                            "position": positions_array,
                           f"MAC_in_{n_accessions}_accs":macs
              })

    score_indices = pd.read_csv(f"results/{pheno_data}_pvalue_indices.csv")["0"]

    score_data = h5py.File(f"./data/{pheno_data}.hdf5")
    scores      = np.array([]).astype(int) #scores of classical GWAS

    for i in range(1,6):
        scores      = np.concatenate((scores, score_data["pvalues"][f"chr{i}"]["scores"][:]))
    
 
    results.loc[score_indices, "scores"] = scores
    
    #saving dataframe
    
    #results.to_csv(f"./results/{pheno_data}/firstpass/results_{pheno_data}_{last_global_constraint}.csv", index = False)
    results.to_csv(f"./results/{pheno_data}/firstpass/results_{pheno_data}_snp_variance+mean.csv", index = False)
    
    ###### plotting the results ######

    x = range(n_snps)
    y = entropy_differences

    positions = [0 ,2_500_000 , 5_000_000, 7_500_000, 10_000_000]
    labels = ["0", "2.500.000", "5.000.000", "7.500.000", "10.000.000"]
    plt.xticks(positions, labels)


    plt.scatter(x[:2597735], y[:2597735], marker='.', s=0.3, color="#e74c3c", label = "Chromosom 1")
    plt.scatter(x[2597735:4466530], y[2597735:4466530], marker='.', s=0.3, color="#5dade2",label = "Chromosom 2")
    plt.scatter(x[4466530:6660782], y[4466530:6660782], marker='.', s=0.3, color="#52be80",label = "Chromosom 3")
    plt.scatter(x[6660782:8427786], y[6660782:8427786], marker='.', s=0.3, color="#a569bd",label = "Chromosom 4")
    plt.scatter(x[8427786:], y[8427786:], marker='.', s=0.3, color="#f4d03f",label = "Chromosom 5")

    plt.xlabel("SNPs", fontweight="bold")
    plt.ylabel("entropy difference", fontweight='bold')
    
    #saving plot
    
    #plt.savefig(f"./results/{pheno_data}/firstpass/plots/{pheno_data}_{last_global_constraint}.png", dpi = 800)
    plt.savefig(f"./results/{pheno_data}/firstpass/plots/{pheno_data}snp_variance+mean.png", dpi = 800)
    plt.close()
    
    print(sum(np.isnan(results["entropy_difference"])))
    
finish = time.perf_counter()

print(finish-start)