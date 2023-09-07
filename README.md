# TOTEM_GWAS

__firstpass_script.py__ -> script to perform the first second pass with TOTEM.

__secondpass_script.py__ -> script to perform the second pass with TOTEM.

__evaluation_firstpass.py__ -> script to plot the overrepresentation of FT genes in the TOTEM GWAS firstpass results.

__auxillary.py__ -> includes the Newton-Raphson-method to solve the non-linear equation system

__evaluation_secondpass.py__ -> script to plot the overrepresentation of FT genes and their interactors in the TOTEM GWAS secondpass results.

__get_interactors.py__ -> extract information about the interactors of FT genes published in FLOR-ID

__gene_annotation_script.py__ -> script to map each SNP to a gene using Araport11 annotation data. The result is saved as a csv file.

__plotting_emma_results.py__ -> script to make Manhatten plots out of the results from classical GWAS.

__emma_pvalue_indices.py__ -> creates a csv file that helps to map the scores from classical GWAS to the corresponding SNPs

------------------------------------------------

The data directory includes phenotype data that is not publicly available (FT for each replicate).

To run firstpass_script.py & evaluation_firstpass.py , the following script must be run first: 
- gene_annotation_script.py 
- emma_pvalue_indices.py 

To run secondpass_script.py & evaluation_secondpass.py , the following script must be run first: \\
- firstpass_script.py 
- evaluation_firstpass.py 
- get_interactors.py
