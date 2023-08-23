import os
import pandas as pd
from PyPDF2 import PdfMerger
import re
import sys

###############################################################################
#                              snakemake inputs                               #
###############################################################################

sys.stderr = sys.stdout = open(snakemake.log[0], 'w')

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

metilene_lifeline_aggregated = snakemake.input[0]
metilene_lifeline_eval = snakemake.output[0] # metilene_lifeline_eval
metilene_lifeline_eval_pdfs = snakemake.output[1] # metilene_lifeline_eval_pdfs
###############################################################################

# # snakemake inputs:
# metilene_lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_aggregated.tsv.gz"
# # snakemake output:
# metilene_lifeline_eval = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated.tsv.gz"
# metilene_lifeline_eval_pdfs = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated.pdf"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
# project = "TCGA-CESC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_5"
###############################################################################
#                                  test set                                   #
###############################################################################
# # snakemake inputs:
# metilene_lifeline_aggregated = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_4/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene2/metilene2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5/metilene2_lifelines_aggregated.tsv.gz"
# # snakemake output:
# metilene_lifeline_eval = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_4/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene2/metilene2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5/metilene2_lifelines_evaluated.tsv.gz"
# # snakemake wildcards:
# output_path = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_4"
# project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_0"
# end = "tsv.gz"
###############################################################################
#                                  test set                                   #
###############################################################################

DF_eval = pd.read_table(metilene_lifeline_aggregated)
if DF_eval.empty:
    DF_eval.to_csv(metilene_lifeline_eval)
    open(metilene_lifeline_eval_pdfs, 'a').close()
    os._exit(0)

# starting point are the scored base_plots, sort them p_sum wise [and take 20
# best of them]
DF_base_sort_temp = DF_eval.set_index(['plot_type', 'scored']).sort_index().loc[('base_plot', True), :].sort_values('p_sum').reset_index()
DF_base_sort_temp.set_index(['ENSG', 'threshold'], inplace=True)
DF_eval = DF_eval.set_index(DF_base_sort_temp.index.names).loc[DF_base_sort_temp.index,:]
# for each ENSG take the best result:
ENSGs = DF_eval.reset_index()['ENSG'].value_counts().index.to_list()
DF_eval_final = pd.concat([pd.concat([DF_eval.loc[ensg,:].reset_index().iloc[:3,:], pd.DataFrame({'ENSG':[ensg]})], axis=1).fillna(ensg) for ensg in ENSGs]).set_index('ENSG')

DF_eval_final.to_csv(metilene_lifeline_eval, sep='\t')

pdfs_to_merge = [re.sub('tsv.*', 'pdf', i) for i in  DF_eval_final['file_path'].values.tolist()]

merger = PdfMerger()

for pdf in pdfs_to_merge:
    merger.append(pdf)

merger.write(metilene_lifeline_eval_pdfs)
merger.close()

