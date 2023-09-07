import pandas as pd
import numpy as np
from tqdm import tqdm
import sys
import h5py
import gffutils as gff
import os
import concurrent.futures

#### importing genotype data ####

snp_data = h5py.File("./data/2029_snps.hdf5")

positions   = np.array(snp_data["positions"])
chr_numbers = np.array([]).astype(int)
indices     = np.arange(len(snp_data["positions"]))

j = 1

for chr_region in snp_data["positions"].attrs["chr_regions"]:
    chr_numbers = np.concatenate((chr_numbers, np.full(((chr_region[1]-chr_region[0]),),j)))
    j+=1
    
#### importing Araport11 data ####

if not "araport11_correct_encoding.gff" in os.listdir("./data"): #correct encoding errors
    
    with open("./data/araport11.gff", 'r', encoding='utf-8', errors='replace') as f:
        data = f.read()

    with open("./data/araport11_correct_encoding.gff", 'w', encoding='utf-8') as f:
        f.write(data)
        
if not "araport11.db" in os.listdir("./data"):
    
    db = gff.create_db("./data/araport11_correct_encoding.gff",\
                       dbfn='./data/araport11.db',\
                       force = True,\
                       merge_strategy='merge',\
                       keep_order = True,\
                       sort_attribute_values=True
                      )
    
else:
    db = gff.FeatureDB("./data/araport11.db")
        
keys = ["gene_code","chromosome","symbol","bp_start","bp_end","gene_length","computational_description","curator_summary"]

annotation_data = []

for gene in db.features_of_type("gene"):
    
    condition1 = (np.isin(gene.id[2], ["1","2","3","4","5"])) #filter out all non-nuclear genes
    try:
        condition2 = (gene.attributes["locus_type"][0] == "protein_coding") #protein_coding genes only
        
    except KeyError:
        if gene.id != "AT1G69572":
            condition2 = True
        else:                         # 30 genes miss the locus_type object. All of them are protein_coding besides AT1G69572
            condition2 = False
        pass
    
    if condition1 & condition2:
        j += 1

        mRNAs = [mRNA for mRNA in db.children(gene, featuretype='mRNA')]

        gene_dict = {key: "-" for key in keys}

        gene_dict["gene_code"]   = gene.id
        gene_dict["chromosome"]  = int(gene.seqid[3])
        gene_dict["bp_start"]    = gene.start
        gene_dict["bp_end"]      = gene.stop
        gene_dict["gene_length"] = gene.stop-gene.start
        
        
        attributes = np.array(["symbol","computational_description", "curator_summary"])
        
        mask1 = np.isin(attributes, list(gene.attributes.keys()))
        
        for attribute in attributes[mask1]:
            gene_dict[attribute] = gene.attributes[attribute]
        
        if len(mRNAs) > 0:
            mask2 = np.isin(attributes, list(mRNAs[0].attributes.keys()))
            
            for attribute in attributes[mask2]:
                try:
                    gene_dict[attribute] = mRNAs[0][attribute][0]
                except IndexError:
                    pass

        annotation_data.append(gene_dict)
        
complete_gene_data = pd.DataFrame.from_records(annotation_data)

complete_gene_data.sort_values(["chromosome", "bp_start"], ascending = [True, True], inplace = True)

#### functions #### 

def divide_and_conquer(df, value):
    """
    df    = Welches Chromosom wir betrachten
    value = Position unseres snps 
    """
    first_index = 0
    last_index  = len(df)-1
    
    if value < df.loc[first_index, "bp_start"]:
        return first_index
    
    if value > df.loc[last_index, "bp_end"]:
        return last_index

    
    while first_index <= last_index:
        
        mid = (last_index + first_index) // 2

        mid_value = df.loc[mid, ["bp_start", "bp_end"]]
        
       # if len(df.loc[first_index:last_index]) <= 50:
       #     display(df.loc[first_index:last_index])

        if value < mid_value[0]:
            if value > df.loc[(mid-1), "bp_end"]:
                return (find_closest_gene(df, value, mid, (mid-1)),"byproxy")
            last_index = mid - 1
            
        elif value > mid_value[1]:
            if value < df.loc[(mid+1), "bp_start"]:
                return (find_closest_gene(df, value, (mid+1), mid), "byproxy")
            first_index = mid + 1
            
        elif (value >= mid_value[0]) & (value <= mid_value[1]):
            return mid 
        
    return "intergenic variant"

def get_annotation_data_for_snp(index):
    
    chromsome_number = chr_numbers[index]
    
    current_data = chromosome_list[(chromsome_number-1)]
    
    row_of_interest = divide_and_conquer(current_data, positions[index])
    
    if isinstance(row_of_interest, tuple):
        row_dict = current_data.loc[row_of_interest[0]].to_dict()
        
        row_dict["byproxy"] = True
        
        return row_dict
     
    row_dict = current_data.loc[row_of_interest].to_dict()
        
    row_dict["byproxy"] = False
        
    return row_dict

def find_closest_gene(chromosome_df ,position, value1, value2):
    
    #Abstand zum Gen, dass danach liegt
    abstand1 = chromosome_df.loc[value1 ,"bp_start"] - position
    
    #Abstand zum Gen, dass davor liegt
    abstand2 = position - chromosome_df.loc[value2 ,"bp_end"]
    
    return value1 if abstand1 < abstand2 else value2 

# Wenn beide Gene gleich weit weg sind, wird by default das davorliegende genommen

#### end of pre-processing and data import ####

chromosome_list = []

for i in range(1,6):
    chromosome_gene_data = complete_gene_data[complete_gene_data["chromosome"] == i].copy()
    chromosome_gene_data.sort_values("bp_start",inplace = True)
    chromosome_gene_data.reset_index(inplace = True, drop = True)
    
    chromosome_list.append(chromosome_gene_data)

with concurrent.futures.ProcessPoolExecutor() as executor:
    results_rows = list(executor.map(get_annotation_data_for_snp, indices))
    
result_df = pd.DataFrame.from_records(results_rows)

# Marker position einfügen
result_df["SNP position"] = positions

# paar columns droppen
result_df.drop(["gene_length", "bp_start", "bp_end"], axis =1,  inplace = True)

# rearrangement of columns

#result_df.reset_index(inplace = True)
#result_df.columns = ['hdf5_index', 'gene_code', 'chromosome', 'symbol','computational_description', 'curator_summary', #'byproxy','SNP position']
#result_df = result_df.loc[:,['gene_code', 'symbol', 'computational_description','curator_summary', 'byproxy', 'chromosome','SNP position', 'hdf5_index']]

result_df.to_csv("./results/annotated_snps.csv", sep = ";", index = False)
