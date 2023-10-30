import os
# from matplotlib.patches import bbox_artist
import pandas as pd
from matplotlib import pyplot as plt
# from pandas.core import apply
# import seaborn as sns
import sys
from pybedtools import BedTool
from natsort import natsort_keygen
import re
from venn import venn
from PyPDF2 import PdfMerger
import subprocess

sys.stdout = sys.stderr = open(snakemake.log[0], 'w')

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

evaluated_tables =  snakemake.input.evaluated_tables
annot = snakemake.input.annot[0]  # single value in list -> due to expand
HM450_annot = snakemake.input.HM450_annot
major_table = snakemake.output.major_table
venn_project =snakemake.output.venn_project
venn_pipeline =snakemake.output.venn_pipeline
md_aggr =snakemake.output.md_aggr
pdf_aggr =snakemake.output.pdf_aggr
pdf_final = snakemake.output.pdf_final
pipeline_str = snakemake.wildcards.pipeline_str
project = snakemake.wildcards.project
drug_str = snakemake.wildcards.drug_str


###############################################################################
#                               test input set                                #
###############################################################################

# # snakemake inputs:
# evaluated_tables = ["/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz"]
# annot = "/scr/palinca/gabor/TCGA-pipeline_7/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# HM450_annot = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../shared/../tcga_metilene/resources/HM450.hg38.manifest.gencode.v36.tsv.gz"
# script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../shared/scripts/majority_vote_all.py"
# # snakemake output:
# major_table = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote.tsv.gz"
# venn_project = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote_venn_projects.pdf"
# venn_pipeline = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote_venn_pipeline.pdf"
# md_aggr = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote_pipeline_project.md"
# pdf_aggr = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote_pipeline_project.pdf"
# pdf_final = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote_pipeline_project_final.pdf"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_7"
# project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
# pipeline_str = "DESeq2_metilene"
# drug_str = "carboplatin_carboplatin,paclitaxel_cisplatin"


# # snakemake inputs:
# evaluated_tables = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz"
# annot = "/scr/palinca/gabor/TCGA-pipeline_7/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# HM450_annot = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../shared/../tcga_metilene/resources/HM450.hg38.manifest.gencode.v36.tsv.gz"
# script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../shared/scripts/majority_vote_all.py"
# # snakemake output:
# major_table = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote.tsv.gz"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_7"
# project = "TCGA-CESC_TCGA-HNSC"
# pipeline_str = "DESeq2_metilene"
# drug_str = "carboplatin_carboplatin,paclitaxel_cisplatin"

# ### test set for TCGA-CESC_TCGA-HNSC_TCGA-LUSC
# # snakemake inputs:
# evaluated_tables = ["/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz"]
# annot = "/scr/palinca/gabor/TCGA-pipeline_7/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# HM450_annot = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../shared/../tcga_metilene/resources/HM450.hg38.manifest.gencode.v36.tsv.gz"
# script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../shared/scripts/majority_vote_all.py"
# # snakemake output:
# major_table = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote.tsv.gz"
# venn_project = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote_venn_projects.pdf"
# venn_pipeline = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote_venn_pipeline.pdf"
# md_aggr = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote_pipeline_project.md"
# pdf_aggr = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote_pipeline_project.pdf"
# pdf_final = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote_pipeline_project_final.pdf"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_7"
# project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
# pipeline_str = "DESeq2_metilene"
# drug_str = "carboplatin_carboplatin,paclitaxel_cisplatin"

# ### test set for single project TCGA-CESC
# snakemake inputs:
# evaluated_tables = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz"
# annot = "/scr/palinca/gabor/TCGA-pipeline_7/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# HM450_annot = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../shared/../tcga_metilene/resources/HM450.hg38.manifest.gencode.v36.tsv.gz"
# script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../shared/scripts/majority_vote_all.py"
# # snakemake output:
# major_table = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote.tsv.gz"
# venn_project = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote_venn_projects.pdf"
# venn_pipeline = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote_venn_pipeline.pdf"
# md_aggr = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote_pipeline_project.md"
# pdf_aggr = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote_pipeline_project.pdf"
# pdf_final = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/DESeq2_metilene/carboplatin_carboplatin,paclitaxel_cisplatin/final_majority_vote_pipeline_project_final.pdf"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_7"
# project = "TCGA-CESC"
# pipeline_str = "DESeq2_metilene"
# drug_str = "carboplatin_carboplatin,paclitaxel_cisplatin"

###############################################################################
#                               test input set                                #
###############################################################################

pipelines = pipeline_str.split('_')
projects = project.split('_')
if len(projects) != 1:
    projects.append(project)


annot_DF = pd.read_table(annot).sort_values(['chr', 'start'], key=natsort_keygen())
final_DF = []
if 'metilene' in pipelines:
    # HM450_annot_DF = pd.read_table(HM450_annot)
    # HM450_annot_DF = HM450_annot_DF.rename({'CpG_chrm': 'chr', 'CpG_beg': 'start', 'CpG_end': 'end'}, axis=1)
    # dropping NA in coordinates, casting start end to int,  sorting coordinate wise,
    # HM450_annot_DF = HM450_annot_DF.dropna(subset=['chr', 'start', 'end']).astype({'start': 'int32', 'end': 'int32'}).sort_values(['chr', 'start'], key=natsort_keygen())
    # met_annot_bed = BedTool.from_dataframe(HM450_annot_DF)
    met_tables = [i for i in evaluated_tables if re.search('metilene', i)]
    met_DF = pd.concat([pd.read_table(i) for i in met_tables])
    # met_DF = met_DF.drop(',Unnamed: 0', axis=1)
    met_DF['Chr'] = met_DF['DMR'].apply(lambda x: x.split('_')[0])
    met_DF['Start'] = met_DF['DMR'].apply(lambda x: x.split('_')[1])
    met_DF['End'] = met_DF['DMR'].apply(lambda x: x.split('_')[2])
    met_DF = pd.concat([met_DF.iloc[:, -3:], met_DF.iloc[:, 0:-3]], axis=1)
    met_DF.sort_values(['Chr', 'Start', 'End'], key=natsort_keygen(), inplace=True)
    met_DF = met_DF.rename({'ENSG': 'ENSG|ENST'}, axis=1)
    met_bed = BedTool.from_dataframe(met_DF)
    annot_bed = BedTool.from_dataframe(annot_DF)
    met_int = met_bed.intersect(annot_bed, loj=True)
    names = list(met_DF.columns) + list(annot_DF.columns)
    met_DF = met_int.to_dataframe(names=names)
    # was mit annot nicht gefunden werden kann wird mit HM450 versucht:
    if not met_DF[met_DF['start'] == -1].empty:
        met_temp = met_DF[met_DF['start'] == -1].iloc[:, :-7]
        met_bed = BedTool.from_dataframe(met_temp)
        HM450 = pd.read_table(HM450_annot).dropna(subset=['CpG_chrm', 'CpG_beg', 'CpG_end']).astype({'CpG_beg': 'int32', 'CpG_end': 'int32'}).loc[:, ['CpG_chrm', 'CpG_beg', 'CpG_end', 'transcriptTypes',  'transcriptIDs', 'geneNames']]
        HM_bed = BedTool.from_dataframe(HM450)
        met_int = met_bed.intersect(HM_bed, loj=True)
        names = list(met_temp.columns) + list(HM450.columns)
        met_temp = met_int.to_dataframe(names=names)
        met_temp = met_temp.drop_duplicates(subset='file_path')
        met_temp.drop(['CpG_chrm', 'CpG_beg', 'CpG_end'], axis=1, inplace=True)
        met_temp['geneNames'] = met_temp['geneNames'].apply(lambda x: x.split(';')[-1])
        met_temp['transcriptTypes'] = met_temp['transcriptTypes'].apply(lambda x: x.split(';')[-1])
        met_temp.rename({'geneNames': 'gene_name'}, axis=1, inplace=True)
        met_temp.rename({'transcriptTypes': 'gene_type'}, axis=1, inplace=True)
    met_DF = met_DF[met_DF['start'] != -1]
    met_DF.drop(['Chr', 'Start', 'End'], axis=1, inplace=True)
    met_DF['#CpGs'] = met_DF['#CpGs'].astype('int32')
    met_DF['pipeline'] = 'metilene'
    final_DF.append(met_DF)
if 'DESeq2' in pipelines:
    deseq_tables = [i for i in evaluated_tables if re.search('DESeq2', i)]
    deseq_DF = pd.concat([pd.read_table(i) for i in deseq_tables])
    # deseq_DF = deseq_DF.drop(',Unnamed: 0', axis=1)
    deseq_DF.dropna(how='all', axis=1, inplace=True)
    deseq_DF = deseq_DF.merge(annot_DF, how='left')
    deseq_DF['pipeline'] = 'DESeq2'
    final_DF.append(deseq_DF)

final_DF = pd.concat(final_DF)
# better namings of ambigous vars:
final_DF = final_DF.rename({'p_value': 'life_p_value', 'pvalue': 'DESeq2_p_value'}, axis=1)


# we want to know which genes can be found over cutoffs, over projects and over
# pipelines, those cols must be created
def apply_project(row):
    for project in projects:
        if re.search(r'\/' + f'{project}' + r'\/', row):
            return project


final_DF['project'] = final_DF['file_path'].apply(apply_project)
final_DF['cutoff'] = final_DF['file_path'].apply(lambda x: re.search(r'cutoff_\d{1,2}', x).group(0))
final_DF['sex'] = final_DF['file_path'].apply(lambda x: re.search(r'\/female\/|\/male\/|\/female_male\/', x).group(0).strip('/'))

print(f'saving {major_table}')
# final_DF.loc[:, ['gene_name', 'project', 'cutoff', 'pipeline']].value_counts().reset_index().value_counts('gene_name')
# gene die in mind 3 versch projecten gefunden wurden:
#       gene_name
# 0
# 3          HOXA9
# 3    B3GALT5-AS1
# 3           GNAS
# 3       RORB-AS1
# 3         GSTA9P
# 3           PAK7
# 3   HOXA10-HOXA9
# 3          HOXA3
# 3          UPK1B
# 3         GABRG3
# 3        MED15P9
# 3  RP11-750B16.1
# 3  RP11-497E19.1
# 3  RP11-497E19.2
# 3  RP11-564D11.3
# 3      LINC01475
# 3       SLC6A10P
# 3    RP11-54A9.1
# 3          FOLR3
# 3          FLRT2
# 3         FAM24B
# final_DF.loc[:, ['gene_name', 'project']].value_counts().reset_index().value_counts('gene_name').reset_index().set_index(0).loc[3]
# gene die in beiden pipelines gefunden wurden:
# (Pdb) final_DF.loc[:, ['gene_name', 'pipeline']].value_counts().reset_index().drop(0, axis=1).value_counts('gene_name').reset_index().set_index(0).loc[2]
#   gene_name
# 0
# 2   SHCBP1L
# 2      TBX5
# 2      HMX2
# 2    PRSS21
# 2      PLD5
# 2     TCP11
# 2       CA4
# 2     FGF19
# 2    SHANK2
# temp = ['SHCBP1L', 'TBX5', 'HMX2', 'PRSS21', 'PLD5', 'TCP11', 'CA4', 'FGF19', 'SHANK2']
# still empty gene_name in metilene, fill them up with the available ENST
# final_DF['gene_name'] = final_DF['gene_name'].replace('.', pd.NA)
# final_DF['gene_name'] = final_DF['gene_name'].fillna(final_DF['ENSG|ENST'])
final_DF.to_csv(major_table, index=None, sep='\t')


set_dict = {}
for project in projects:
    set_dict.update({project: set(final_DF.loc[:, ['gene_name', 'project']].set_index('project').loc[project, :].value_counts().reset_index()['gene_name'].tolist())})
    # set_dict= {project: set(gene_names)}

if len(set_dict) > 1:
    venn(set_dict)
    plt.title(f"Venn Diagram for:\n{', '.join(projects)}\nout of piplines:\n{', '.join(pipelines)}")
    plt.savefig(venn_project)
else:
    open(venn_project, 'w').close

set_dict = {}
for pipeline in pipelines:
    set_dict.update({pipeline: set(final_DF.loc[:, ['gene_name', 'pipeline']].set_index('pipeline').loc[pipeline, :].value_counts().reset_index()['gene_name'].tolist())})
if len(set_dict) > 1:
    venn(set_dict)
    plt.title(f"Venn Diagram for:\n{', '.join(pipelines)}\nout of projects:\n{', '.join(projects)}")
    plt.savefig(venn_pipeline)
else:
    open(venn_pipeline, 'w').close


# final_DF.loc[:, ['gene_name', 'pipeline', 'project']].value_counts().reset_index()

# genes found in more than one project:
# this just works if more than one project and pipeline is applied
if len(projects) > 1 and len(pipelines) > 1:
    gene_project_count_gt_one = final_DF.value_counts(subset=['gene_name', 'project']).reset_index().value_counts('gene_name').to_frame().reset_index().rename({0: 'count'}, axis=1).query('count > 1')
    gene_project_count_gt_one.rename({'count': 'project-count'}, axis=1, inplace=True)
    gene_project_count_list = gene_project_count_gt_one.gene_name.tolist()

    # which of those genes are also found in both pipelines?
    genes_in_both_piplines = final_DF.loc[:, ['gene_name', 'pipeline']].value_counts().reset_index().value_counts(subset=['gene_name']).reset_index().rename({0: 'count'}, axis=1).query("count > 1")
    genes_in_both_piplines.rename({'count': 'pipeline-count'}, axis=1, inplace=True)

    genes_in_proj_pipeline = gene_project_count_gt_one.merge(genes_in_both_piplines)
    gene_list = genes_in_proj_pipeline.gene_name.to_list()
    temp_DF = final_DF.set_index('gene_name').loc[gene_list, :].loc[:, 'project'].reset_index().drop_duplicates()

    for project in projects:
        temp_DF[project] = (temp_DF['project'] == project)
    # temp_DF:
    #     gene_name                        project  TCGA-CESC  TCGA-HNSC  TCGA-LUSC  TCGA-CESC_TCGA-HNSC_TCGA-LUSC
    # 0    TMEM255A                      TCGA-HNSC      False       True      False                          False
    # 6    TMEM255A  TCGA-CESC_TCGA-HNSC_TCGA-LUSC      False      False      False                           True
    # 33   TMEM255A                      TCGA-LUSC      False      False       True                          False

    final_dict = {}
    for gene in gene_list:
        final_dict.update({gene : ' + '.join(temp_DF.set_index('gene_name').loc[gene, 'project'].values.tolist())})

    temp_DF_2 = pd.DataFrame(final_dict, index=range(len(final_dict))).T
    temp_DF_2.index.name = 'gene_name'
    temp_DF_2 = temp_DF_2.reset_index()
    temp_DF_2 = genes_in_proj_pipeline.merge(temp_DF_2)
    temp_DF_2 = temp_DF_2.iloc[:, :3]

    temp_DF_3 = temp_DF.loc[:, ['gene_name'] + projects]
    temp_DF_3 = pd.concat([temp_DF_3.replace(False, pd.NA).set_index('gene_name').T.loc[:, gene].fillna(method='ffill', axis=1).fillna(method='bfill', axis=1) for gene in gene_list], axis=1).T.reset_index().drop_duplicates()
    temp_DF_3 = temp_DF_3.merge(temp_DF_2)

    temp_DF_3 = temp_DF_3.replace(pd.NA, 0).replace(True, 1)
    # temp_DF_3:
    #    gene_name  TCGA-CESC  TCGA-HNSC  TCGA-LUSC  TCGA-CESC_TCGA-HNSC_TCGA-LUSC  project_count  pipeline_count
    # 0   TMEM255A          0          1          1                              1              3               2
    # 1       EBF3          0          1          1                              1              3               2

    with open(md_aggr, 'w') as f:
        f.write(f'# Overview for {" + ".join(projects)}, with therapy {" + ".join(drug_str.split("_"))}\n\n## genes_found in both pipelines and at least 2 projects:\n\n')

    temp_DF_3.rename({'gene_name': 'gene-name'}, axis=1).to_markdown(md_aggr, index=False, mode='a')
    # temp_DF_3.rename({'gene_name':'gene-name'}, axis=1).set_index('gene-name').style.to_latex(md_aggr.replace('.md', '_pipeline.tex'))

    temp_DF_3_annot = temp_DF_3.merge(annot_DF, how='left').loc[:, ['gene_name', 'chr', 'start',  'ENSG', 'gene_type', 'gene_status']]
    temp_DF_3_annot['ENSG'] = temp_DF_3_annot['ENSG'].str.replace('ENSG', 'ENSG-')
    temp_DF_3_annot['gene_type'] = temp_DF_3_annot['gene_type'].str.replace('_', '-')

    with open(md_aggr, 'a') as f:
        f.write('\n\n')
        f.write(r'\newpage')
        f.write('\n\n## Belonging annotations:\n\n')

    temp_DF_3_annot.rename({'gene_name': 'gene-name', 'gene_type': 'gene-type', 'gene_name': 'gene-name', 'gene_status': 'gene-status'}, axis=1).to_markdown(md_aggr, index=False, mode='a')
    # temp_DF_3_annot.rename({'gene_name':'gene-name', 'gene_type':'gene-type', 'gene_name': 'gene-name', 'gene_status': 'gene-status'}, axis=1).set_index('gene-name').style.to_latex(md_aggr.replace('.md', '_pipeline_annot.tex'))

    with open(md_aggr, 'a') as f:
        f.write('\n\n')
        f.write(r'\newpage')
        f.write('\n\n## Venn Diagrams for commonalities between pipelines and between projects :\n\n')

    sequence = ['pandoc', md_aggr, '-t', 'latex', '-o', pdf_aggr]
    # such that pandoc finds the iftex.sty, chdir to this sourcefile, than go one
    # up and into resources/iftex
    print(f"cd into dir: {os.path.join(os.path.dirname(__file__), os.pardir, 'resources', 'iftex')}")
    os.chdir(os.path.join(os.path.dirname(__file__), os.pardir, 'resources', 'iftex'))
    call = subprocess.check_call(sequence, stdout=sys.stdout, stderr=sys.stderr)

    merger = PdfMerger()

    pdfs_to_merge = [pdf_aggr, venn_pipeline, venn_project]

    for pdf in pdfs_to_merge:
        merger.append(pdf)

    merger.write(pdf_final)
    merger.close()

else:
    with open(md_aggr, 'w') as f:
        f.write(f'# Overview for {" + ".join(projects)}, with therapy {" + ".join(drug_str.split("_"))} and the pipelines {" + ".join(pipelines)}\n\n')

    with open(md_aggr, 'a') as f:
        f.write('\n\n')
        f.write(r'\newpage')
        f.write('\n\n## Venn Diagrams for commonalities between pipelines and between projects :\n\n')
    # which of the venn could be created actually?

    def return_venn(venn_list):
        final_venn = []
        for venn_ in venn_list:
            if os.path.getsize(venn_):
                final_venn.append(venn_)
        return final_venn
    final_venn = return_venn([venn_project, venn_pipeline])
    if len(final_venn) == 0:
        open(pdf_final).close
        os._exit(0)
    else:
        sequence = ['pandoc', md_aggr, '-t', 'latex', '-o', pdf_aggr]
        # such that pandoc finds the iftex.sty, chdir to this sourcefile, than go one
        # up and into resources/iftex
        os.chdir(os.path.join(os.path.dirname(__file__), os.pardir, 'resources', 'iftex'))
        call = subprocess.check_call(sequence, stdout=sys.stdout, stderr=sys.stderr)

        merger = PdfMerger()
        pdfs_to_merge = [pdf_aggr] + final_venn

        for pdf in pdfs_to_merge:
            merger.append(pdf)

        merger.write(pdf_final)
        merger.close()
