import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from tqdm import tqdm
import sys
from pathos.multiprocessing import ProcessingPool as Pool

def get_position(index_interactors):
    
    condition1 = (np.isin(secondpass_result["gene_code1"], interactors.loc[index_interactors].values))
    condition2 = (np.isin(secondpass_result["gene_code2"], interactors.loc[index_interactors].values))
    
    position = np.where((condition1 & condition2))[0]
    
    if len(position) == 1:
        return position.item()
    
    else:
        return -1
    
### data import and formatting

interactors = pd.read_csv("results/ft_genes_interactions.csv", index_col = 0)
n_results = sum(pd.Series(os.listdir("results/ft10/secondpass/")).str.startswith("results"))

secondpass_result = pd.read_csv(f"results/ft10/secondpass/results_ft10_secondpass_part1.csv")

for i in tqdm(range(2,n_results+1)):
    secondpass_result = pd.concat([secondpass_result, pd.read_csv(f"results/ft10/secondpass/results_ft10_secondpass_part{i}.csv")])
    
secondpass_result.sort_values("entropy_difference", inplace = True, ascending = False)
secondpass_result.reset_index(drop = True, inplace = True)

print(f"evaluating for {len(secondpass_result)} pairs ")

# evaluate overrepresentation and plot it

with Pool() as pool:
    positions_of_interactors = pd.Series(pool.map(get_position, range(len(interactors))))
    pool.clear()
    
n_missing_interactors = sum(positions_of_interactors == -1)

print(f"{n_missing_interactors} are missing")
    
positions_of_interactors = positions_of_interactors[positions_of_interactors != -1].values

print("start building result array")

result_array = np.zeros((len(secondpass_result),))

print("get rid of all missing interactors")

result_array[positions_of_interactors] = 1/ (len(interactors) - n_missing_interactors)

condensed_result_array = []

chunksize = 10000

from tqdm import tqdm

for i in tqdm(range(0,len(result_array),chunksize)):
    
    condensed_result_array.append(np.sum(result_array["0"][i:i+chunksize]))
    

y = [np.sum(condensed_result_array[0:i]) for i in tqdm(range(len(condensed_result_array)))]

x = np.arange(len(condensed_result_array))*chunksize

plt.plot(x,y)
plt.plot([0,314995000], [0,1])
positions = [0 ,60000000 , 180000000, 300000000]
labels = ["0", "60.000.000", "180.000.000", "300.000.000"]
plt.xticks(positions, labels)


plt.savefig("results/ft10/secondpass/secondpass_overrepresentation.png", dpi = 800)

