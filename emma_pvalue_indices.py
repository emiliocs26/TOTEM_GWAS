import pandas as pd
import numpy as np
from tqdm import tqdm
import sys
import time
import h5py
from pathos.multiprocessing import ProcessingPool as Pool

snp_data = h5py.File("./data/2029_snps.hdf5")

snp_positions   = np.array(snp_data["positions"])
snp_chr_numbers = np.array([]).astype(int)


j = 1

for chr_region in snp_data["positions"].attrs["chr_regions"]:
    snp_chr_numbers = np.concatenate((snp_chr_numbers, np.full(((chr_region[1]-chr_region[0]),),j)))
    j+=1
    
all_snps = pd.DataFrame({"position":snp_positions,
                         "chromosome":snp_chr_numbers})



for pheno_data in ["ft16","ft10"]:

    data = h5py.File(f"./data/{pheno_data}.hdf5")

    scores      = np.array([]).astype(int)
    positions   = np.array([]).astype(int)
    chromosomes = np.array([]).astype(int)

    for i in range(1,6):
        positions   = np.concatenate((positions, data["pvalues"][f"chr{i}"]["positions"][:]))

        n_positions = len(data["pvalues"][f"chr{i}"]["positions"][:])
        scores      = np.concatenate((scores, data["pvalues"][f"chr{i}"]["scores"][:]))

        chromosomes = np.concatenate((chromosomes, np.full((n_positions,), i)))

    # dividing the result dataframe into 5 dataframes for each chromosome, to increase search speed

    dfs  = []

    for i in range(1,6):

        dfs.append(all_snps[all_snps["chromosome"] == i])

    def get_index(i):

        current_df = dfs[chromosomes[i]-1]

        condition = current_df["position"] == positions[i]

        return current_df.loc[condition].index.item()


    with Pool() as pool:
        indices = pd.Series(pool.map(get_index, tqdm(range(len(positions)))))
        pool.clear()

    indices.to_csv(f"results/{pheno_data}_pvalue_indices.csv", index = False)


