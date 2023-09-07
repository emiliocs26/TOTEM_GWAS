# make dirs necessary for the scripts to work
import os

os.mkdir("results")

for pheno_data in ["ft10","ft16"]:
    os.mkdir(f"results/{pheno_data}")
    os.mkdir(f"results/{pheno_data}/firstpass")
    os.mkdir(f"results/{pheno_data}/firstpass/plots")
    
os.mkdir(f"results/ft10/secondpass")