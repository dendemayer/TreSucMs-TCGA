import os
import pandas as pd
import sys
if "snakemake" in dir():
    sys.stderr = sys.stdout = open(snakemake.log[0], 'w')

    print('# snakemake inputs:')
    [ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

    print('# snakemake output:')
    [ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

    print('# snakemake wildcards:')
    [ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

    print('# snakemake params:')
    [ print(f'{i[0]} = "{i[1]}"') for i in snakemake.params.items()]

    qval_input = snakemake.input.qval_input
    qval_lim = snakemake.output.qval_lim
################################################################################
#                                   testset                                    #
################################################################################
else:
    # snakemake inputs:
    qval_input = "/scr/palinca/gabor/TCGA-pipeline_9_pval_prod/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin_cisplatin,gemcitabine/female/cutoff_8/metilene_qval.0.05.out"
    # snakemake output:
    qval_lim = "/scr/palinca/gabor/TCGA-pipeline_9_pval_prod/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin_cisplatin,gemcitabine/female/cutoff_8/metilene_qval.0.05_limited.out"
    # snakemake wildcards:
    output_path = "/scr/palinca/gabor/TCGA-pipeline_9_pval_prod"
    project = "TCGA-LUSC"
    drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin_cisplatin,gemcitabine"
    gender = "female"
    cutoff = "cutoff_8"
################################################################################
#                                   testset                                    #
################################################################################

DF = pd.read_table(qval_input, header=None).dropna()
if DF.empty:
    DF.to_csv(qval_lim)
    os._exit(0)

DF_mmd_pos = DF[DF[4] > 0]
DF_mmd_pos = DF_mmd_pos.sort_values(3).iloc[0:60, :]
DF_mmd_neg = DF[DF[4] < 0]
DF_mmd_neg = DF_mmd_neg.sort_values(3).iloc[0:60, :]

DF_final = pd.concat([DF_mmd_pos, DF_mmd_neg])
print(f'saving {qval_lim}')
DF_final.to_csv(qval_lim, sep='\t', index=None, header=None)
