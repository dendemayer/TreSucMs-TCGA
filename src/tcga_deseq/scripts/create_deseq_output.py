import pandas as pd
import os
import sys
import subprocess

sys.stderr = sys.stdout = open(snakemake.log[0], 'w')

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

R_path = snakemake.input.R_path
summary_table = snakemake.input.summary_table
summary_table_info = snakemake.input.summary_table_info
project = snakemake.wildcards.project
deseq_counts = snakemake.output.deseq_counts[0]
deseq_output = os.path.split(deseq_counts)[0]
gender = snakemake.wildcards.gender

###############################################################################
#                               test input set                                #
###############################################################################
# # snakemake inputs:
# summary_table = "/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC/DESeq2/DESeq2_input_table/carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/summary_for_DESeq2.tsv"
# summary_table_info = "/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC/DESeq2/DESeq2_input_table/carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/summary_for_DESeq2_INFO.tsv"
# R_path = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/tcga_deseq/resources/DESeq2_diffexp_multifactor.R"
# script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/tcga_deseq/scripts/create_deseq_output.py"
# # snakemake output:
# deseq_counts = ["/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC/DESeq2/DESeq2_output/carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/DESeq2_heatmap_log2fINCREASE_nt_counts.tsv" ,"/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC/DESeq2/DESeq2_output/carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/DESeq2_heatmap_log2fINCREASE_raw_counts.tsv" ,"/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC/DESeq2/DESeq2_output/carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/DESeq2_heatmap_log2fINCREASE_norm_counts.tsv" ,"/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC/DESeq2/DESeq2_output/carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/DESeq2_heatmap_log2fINCREASE_vsd_counts.tsv" ,"/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC/DESeq2/DESeq2_output/carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/DESeq2_heatmap_log2fDECREASE_nt_counts.tsv" ,"/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC/DESeq2/DESeq2_output/carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/DESeq2_heatmap_log2fDECREASE_raw_counts.tsv" ,"/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC/DESeq2/DESeq2_output/carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/DESeq2_heatmap_log2fDECREASE_norm_counts.tsv" ,"/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC/DESeq2/DESeq2_output/carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/DESeq2_heatmap_log2fDECREASE_vsd_counts.tsv"][0]
# results = ["/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC/DESeq2/DESeq2_output/carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/DESeq2_log2fINCREASE_result.tsv" ,"/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC/DESeq2/DESeq2_output/carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/DESeq2_log2fDECREASE_result.tsv"]
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_5"
# project = "TCGA-CESC"
# drug_combi = "carboplatin,paclitaxel_cisplatin"
# gender = "female_male"
# cutoff = "cutoff_0"
# deseq_output = os.path.split(deseq_counts)[0]
###############################################################################
#                               test input set                                #
###############################################################################

# if the input table is empty, write empty tables for ever output table which
# shall be created

if pd.read_table(summary_table, sep='\t').empty:
    for out in snakemake.output:
        open(out, 'a').close()
    os._exit(0)


# also if just alive or just dead cases are included, no diff exp analysis is
# possible:
DF_temp = pd.read_table(summary_table_info, sep='\t').T
if DF_temp.iloc[1:,0].nunique() == 1:
    for out in snakemake.output:
        open(out, 'a').close()
    os._exit(0)

# also if female_male is gender, but the summary_table holds just one gender
if len(gender.split('_')) == 2 :
    if DF_temp.iloc[1:,1].nunique() == 1:
        for out in snakemake.output:
            open(out, 'a').close()
        os._exit(0)

sequence = ['Rscript', R_path, summary_table, summary_table_info, deseq_output,
            project]
subprocess.check_call(sequence, stdout=sys.stdout, stderr=sys.stderr)
# subprocess.Popen(sequence, stdout=sys.stdout, stderr=sys.stderr)
