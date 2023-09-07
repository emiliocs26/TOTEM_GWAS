import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def get_annotated_scoreboard(df: pd.DataFrame, option: str):
    
    """
    Takes the primary list of SNPs with corresponding entropy differences and makes an annotated scoreboard out of it 
    
    options:
    
    "scores"              -> returns list ranked by classical GWAS score
    "entropy_difference" -> returns list ranked by entropy difference
    
    """
    
    df.sort_values(option, ascending = False, inplace = True)
    df = annotated_snps.loc[df["hdf5_index"]]
    df.drop_duplicates("gene_code", keep = "first", inplace =True)
    df.reset_index(drop = True, inplace = True)
    
    df.sorting = option
    
    return df

def get_plot_arrays(ranked_list):
    
    """
    Takes ranked list and returns arrays for the corresponding overrepresentation plot
    """
    
    where_are_they = []

    for gene in FLOR_ID_genes:

        gene_on_scoreboard = ranked_list[ranked_list["gene_code"]== gene]

        if len(gene_on_scoreboard) == 1:

            index = gene_on_scoreboard.index.item()
            where_are_they.append(index)

        else:
            pass


    result_array = np.zeros(len(ranked_list))

    result_array[where_are_they] = (1/len(where_are_they))

    x = range(len(ranked_list))
    y = [np.sum(result_array[0:i]) for i in range(len(result_array))]
    
    return x,y

annotated_snps = pd.read_csv("./results/annotated_snps.csv", sep = ";")

FLOR_ID_genes = ['AT2G13540', 'AT5G48300', 'AT1G75060', 'AT1G19330', 'AT4G24540',
       'AT5G13790', 'AT3G57390', 'AT4G22950', 'AT3G57230', 'AT2G22630',
       'AT2G45650', 'AT1G69120', 'AT2G45430', 'AT3G01090', 'AT4G36920',
       'AT1G79430', 'AT1G18450', 'AT3G33520', 'AT2G37630', 'AT1G51450',
       'AT1G66650', 'AT3G23060', 'AT2G27550', 'AT4G32980', 'AT3G44110',
       'AT4G03090', 'AT1G04870', 'AT2G31650', 'AT1G05830', 'AT5G42400',
       'AT5G14170', 'AT5G44200', 'AT4G38960', 'AT5G62040', 'AT3G18550',
       'AT2G46020', 'AT4G23100', 'AT2G46830', 'AT5G62430', 'AT5G39660',
       'AT3G47500', 'AT2G34140', 'AT1G69570', 'AT2G23380', 'AT5G64960',
       'AT4G34530', 'AT1G26260', 'AT5G48560', 'AT1G10120', 'AT3G60250',
       'AT2G44680', 'AT5G15840', 'AT5G57660', 'AT3G07650', 'AT2G32950',
       'AT1G04400', 'AT4G08920', 'AT1G50700', 'AT2G17290', 'AT4G38680',
       'AT1G71800', 'AT1G17760', 'AT1G26830', 'AT5G46210', 'AT1G01040',
       'AT3G43920', 'AT5G20320', 'AT4G10180', 'AT3G19140', 'AT4G22140',
       'AT2G40080', 'AT2G25930', 'AT2G03500', 'AT1G77300', 'AT5G04240',
       'AT1G79730', 'AT2G06210', 'AT5G11530', 'AT5G51230', 'AT4G15880',
       'AT2G39810', 'AT1G08260', 'AT1G35460', 'AT4G09180', 'AT1G51140',
       'AT2G42280', 'AT4G16280', 'AT5G10140', 'AT1G68050', 'AT3G04610',
       'AT3G10390', 'AT4G35900', 'AT2G17770', 'AT2G33835', 'AT3G20740',
       'AT2G21070', 'AT1G77080', 'AT2G30120', 'AT1G65480', 'AT5G60910',
       'AT2G43410', 'AT5G16320', 'AT4G00650', 'AT5G61920', 'AT1G31814',
       'AT5G24860', 'AT5G06850', 'AT2G19520', 'AT4G25530', 'AT5G13480',
       'AT4G02780', 'AT1G79460', 'AT4G25420', 'AT5G51810', 'AT5G07200',
       'AT1G78440', 'AT1G30040', 'AT2G34555', 'AT1G47990', 'AT1G02400',
       'AT1G50960', 'AT4G21200', 'AT1G22770', 'AT1G15550', 'AT1G80340',
       'AT1G14920', 'AT3G02885', 'AT3G05120', 'AT3G63010', 'AT5G27320',
       'AT5G63960', 'AT2G20570', 'AT5G44190', 'AT3G44680', 'AT5G64610',
       'AT5G09740', 'AT2G21660', 'AT5G63110', 'AT4G40060', 'AT5G61060',
       'AT5G40490', 'AT3G54560', 'AT1G52740', 'AT2G38810', 'AT3G63070',
       'AT3G26744', 'AT5G23150', 'AT2G44950', 'AT1G55250', 'AT2G48160',
       'AT4G29130', 'AT5G44160', 'AT2G18915', 'AT1G62830', 'AT3G13682',
       'AT1G01060', 'AT2G46260', 'AT3G61600', 'AT3G20810', 'AT4G20400',
       'AT4G02560', 'AT3G57300', 'AT2G34880', 'AT3G45880', 'AT5G48890',
       'AT5G61850', 'AT4G00830', 'AT5G64170', 'AT3G54500', 'AT5G65050',
       'AT5G65060', 'AT3G46640', 'AT5G65070', 'AT5G65080', 'AT1G12910',
       'AT3G26640', 'AT3G01460', 'AT4G00450', 'AT1G55325', 'AT4G04920',
       'AT2G22370', 'AT2G25095', 'AT4G30972', 'AT4G31877', 'AT5G10945',
       'AT5G11977', 'AT5G26147', 'AT2G19425', 'AT5G55835', 'AT2G28056',
       'AT5G04275', 'AT3G11435', 'AT4G37280', 'AT1G02740', 'AT4G24680',
       'AT5G58230', 'AT3G28910', 'AT5G12840', 'AT2G34720', 'AT2G38880',
       'AT1G54830', 'AT1G08970', 'AT5G18240', 'AT3G04030', 'AT5G47640',
       'AT4G14540', 'AT5G63470', 'AT1G60220', 'AT1G10570', 'AT2G18790',
       'AT3G22590', 'AT1G25540', 'AT4G24620', 'AT5G51820', 'AT1G72390',
       'AT1G09570', 'AT5G35840', 'AT4G16250', 'AT4G18130', 'AT3G12810',
       'AT5G49020', 'AT2G43010', 'AT3G59060', 'AT1G80070', 'AT5G60100',
       'AT5G24470', 'AT5G02810', 'AT2G46790', 'AT3G48430', 'AT2G01570',
       'AT2G28830', 'AT3G46510', 'AT1G32230', 'AT2G47700', 'AT1G66350',
       'AT3G03450', 'AT1G76710', 'AT5G37055', 'AT5G17490', 'AT1G54440',
       'AT5G35910', 'AT5G23730', 'AT5G37260', 'AT2G44150', 'AT3G54990',
       'AT2G39250', 'AT2G45660', 'AT4G39100', 'AT5G60410', 'AT4G31120',
       'AT2G33810', 'AT1G53160', 'AT3G15270', 'AT3G57920', 'AT2G46340',
       'AT3G15354', 'AT1G53090', 'AT1G30970', 'AT2G42200', 'AT2G35510',
       'AT4G10710', 'AT3G28730', 'AT5G59560', 'AT1G06040', 'AT5G06170',
       'AT3G43190', 'AT2G22540', 'AT5G61380', 'AT5G03840', 'AT1G25560',
       'AT1G68840', 'AT5G17690', 'AT2G23740', 'AT3G22380', 'AT2G28550',
       'AT5G60120', 'AT5G67180', 'AT4G20370', 'AT1G14400', 'AT2G02760',
       'AT1G78580', 'AT1G15750', 'AT5G06600', 'AT3G11910', 'AT3G49600',
       'AT2G30140', 'AT3G49660', 'AT5G61150', 'AT1G61040', 'AT4G30200',
       'AT4G29830', 'AT4G28190', 'AT5G57380', 'AT3G24440', 'AT2G18880',
       'AT4G16845', 'AT1G57820', 'AT1G66050', 'AT5G39550', 'AT1G28520',
       'AT2G42400', 'AT3G18990', 'AT4G11880', 'AT5G57360', 'AT4G26440',
       'AT5G45600']


# plot entropy scoring vs p-value scoring

for pheno_data in ["ft10","ft16"]:
    
    ranked_lists = []
    
    primary_result = pd.read_csv(f"./results/{pheno_data}/firstpass/results_{pheno_data}_globalmean.csv")
    
    ranked_lists.append(get_annotated_scoreboard(primary_result, option = "entropy_difference"))
    ranked_lists.append(get_annotated_scoreboard(primary_result, option = "scores"))
                                   
    for ranked_list in ranked_lists:

        if ranked_list.sorting == "entropy_difference":
            string = "entropy_ranked"
            
        if ranked_list.sorting == "scores":
            string = "score_ranked"

        x,y = get_plot_arrays(ranked_list)

        plt.plot(x,y, label = f"{string}")
        plt.xlabel("top scorers")
        plt.ylabel("relative proportion of FT genes found")
        plt.plot([0, len(ranked_list)], [0,1])
        plt.legend()
        
    
    plt.savefig(f"./results/{pheno_data}/firstpass/plots/overrepresentation_{pheno_data}_entropy_vs_score.png", dpi = 800)
    plt.close()
    
# plot results including differing moments


moments = ["globalmean", "globalvariance", "globalskewness"]

for pheno_data in ["ft10", "ft16"]:
   
    ranked_lists = []
    
    for moment in moments:
        
        primary_result =  pd.read_csv(f"./results/{pheno_data}/firstpass/results_{pheno_data}_{moment}.csv")
        
        ranked_lists.append(get_annotated_scoreboard(primary_result, option = "entropy_difference"))
        
    j = 0 
    
    for ranked_list in ranked_lists:

        string = moments[j]

        x,y = get_plot_arrays(ranked_list)


        plt.plot(x,y, label = f"{string}")
        plt.xlabel("top scorers")
        plt.ylabel("relative proportion of FT genes found")
        plt.plot([0, len(ranked_list)], [0,1])
        plt.legend()
        j += 1

    plt.savefig(f"./results/{pheno_data}/firstpass/plots/overrepresentation_{pheno_data}_moments.png", dpi = 800)
    plt.close()
    
    
    
# snp_specific variance 


for pheno_data in ["ft10", "ft16"]:
   
    ranked_lists = []
        
    primary_result =  pd.read_csv(f"./results/{pheno_data}/firstpass/results_{pheno_data}_snp_variance.csv")
        
    ranked_lists.append(get_annotated_scoreboard(primary_result, option = "entropy_difference"))
        
    j = 0 
    
    for ranked_list in ranked_lists:
        
        x,y = get_plot_arrays(ranked_list)


        plt.plot(x,y, label = f"snp_variance")
        plt.xlabel("top scorers")
        plt.ylabel("relative proportion of FT genes found")
        plt.plot([0, len(ranked_list)], [0,1])
        plt.legend()
        j += 1

    plt.savefig(f"./results/{pheno_data}/firstpass/plots/overrepresentation_{pheno_data}_snp_variance.png", dpi = 800)
    plt.close()
