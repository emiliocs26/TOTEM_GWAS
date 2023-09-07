import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

pheno_datasets = ["ft16", "ft10"]

for pheno_data in pheno_datasets:
    results = pd.read_csv(f"results/{pheno_data}/results_{pheno_data}_globalmean_+_scores.csv")
    results = results[~np.isnan(results["score"])]

    x = results["hdf5_index"]
    y = results["score"]

    positions = [0 ,2_500_000 , 5_000_000, 7_500_000, 10_000_000]
    labels = ["0", "2.500.000", "5.000.000", "7.500.000", "10.000.000"]
    plt.xticks(positions, labels)

    masks  = []
    colors = ["#e74c3c", "#5dade2", "#52be80", "#a569bd", "#f4d03f"]
    lables = ["chromosome 1", "chromosome 2", "chromosome 3", "chromosome 4", "chromosome 5"]

    for i in range(0,5):

        mask = results["chromosome"] == (i+1)

        plt.scatter(x[mask], y[mask],marker='.', s=0.3, color=colors[i], label = lables[i])

    plt.xlabel("SNPs", fontweight="bold")
    plt.ylabel("score", fontweight='bold')
    plt.savefig(f"./results/emma_manhatten_{pheno_data}.png", dpi = 800)
    plt.close()