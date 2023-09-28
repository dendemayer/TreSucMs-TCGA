import sys
import os
import pandas as pd
# from methyl import set_logger
# from lifelines.datasets import load_waltons
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from decimal import Decimal
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter
from lifelines.plotting import add_at_risk_counts
import statistics as st

"""
validate found DMRs
-> need the start position on whichs betavalues the lifeline plot is based on,
    to be found in start_tsv
-> all betavalues, (pre filtered gender and cutoff )
    to be found in summary_for_metilene, and the complement
thats it, consider also here the threshold when dividing in UP and DOWN in the
merged summary for metilene
for validation filtering add cols:
'threshold', 'plot_type', 'fst_life_mean', 'scnd_life_mean'
"""
# #####################
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print('# snakemake inputs:')
[print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

meta_table = snakemake.input[0]
start_tsv = snakemake.input[1]
summary = snakemake.input[2]
summary_complement = snakemake.input[3]
annot_file = snakemake.input[4]
annot_file_2 = snakemake.input[5]

# # start_tsv is the table which is based on the lifeline_DMR_plot, and holds the
# # extract start value
DMR = snakemake.wildcards.DMR
UP_val_plot = snakemake.output.UP_val_plot
UP_val_tsv =snakemake.output.UP_val_tsv
DOWN_val_plot =snakemake.output.DOWN_val_plot
DOWN_val_tsv =snakemake.output.DOWN_val_tsv
threshold = snakemake.wildcards.threshold
drug_combi = snakemake.wildcards.drug_combi
DRUGS = drug_combi.split('_')
drug_combi = snakemake.wildcards.drug_combi.replace('_', ';')
cutoff = snakemake.wildcards.cutoff.split('_')[1]
project = ', '.join(snakemake.wildcards.project.split('_'))
gender = ', '.join(snakemake.wildcards.gender.split('_'))

###############################################################################
#                                 test input                                  #
# ###############################################################################

# meta_table = "/scr/palinca/gabor/TCGA-pipeline_4/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/merged_meta_files/cutoff_5/meta_info_druglist_merged_drugs_combined.tsv"
# start_tsv = "/scr/palinca/gabor/TCGA-pipeline_4/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_5/metilene_intersect_lifeline_plot_chr17_44656158_44657291.tsv"
# summary = "/scr/palinca/gabor/TCGA-pipeline_4/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/summary_for_metilene.tsv"
# summary_complement = "/scr/palinca/gabor/TCGA-pipeline_4/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/summary_for_metilene_complement.tsv"
# annot_file = "/scr/palinca/gabor/TCGA-pipeline_4/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# annot_file_2 = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/tcga_metilene/../shared/resources/annot_from_betafile.tsv.gz"
# script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/tcga_metilene/scripts/create_lifeline_plots_validation.py"
# # snakemake output:
# UP_val_plot = "/scr/palinca/gabor/TCGA-pipeline_4/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_5/metilene_intersect_lifeline_plot_chr17_44656158_44657291_UP_val.pdf"
# UP_val_tsv = "/scr/palinca/gabor/TCGA-pipeline_4/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_5/metilene_intersect_lifeline_plot_chr17_44656158_44657291_UP_val.tsv"
# DOWN_val_plot = "/scr/palinca/gabor/TCGA-pipeline_4/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_5/metilene_intersect_lifeline_plot_chr17_44656158_44657291_DOWN_val.pdf"
# DOWN_val_tsv = "/scr/palinca/gabor/TCGA-pipeline_4/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_5/metilene_intersect_lifeline_plot_chr17_44656158_44657291_DOWN_val.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_4"
# project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "male"
# cutoff = "cutoff_5"
# threshold = "threshold_5"
# DMR = "chr17_44656158_44657291"
# DRUGS = drug_combi.split('_')

# snakemake inputs:
# meta_table = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-LUSC/metilene/merged_meta_files/cutoff_8/meta_info_druglist_merged_drugs_combined.tsv"
# start_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0/metilene_intersect_lifeline_plot_chr7_97012143_97022910.tsv"
# summary = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-LUSC/metilene/metilene_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/summary_for_metilene.tsv"
# summary_complement = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-LUSC/metilene/metilene_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/summary_for_metilene_complement.tsv"
# annot_file = "/scr/palinca/gabor/TCGA-pipeline_2/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# annot_file_2 = "../shared/resources/annot_from_betafile.tsv.gz"
# script_file = "../tcga_metilene/scripts/create_lifeline_plots_validation.py"
# # snakemake output:
# UP_val_plot = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0/metilene_intersect_lifeline_plot_chr7_97012143_97022910_UP_val.pdf"
# UP_val_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0/metilene_intersect_lifeline_plot_chr7_97012143_97022910_UP_val.tsv"
# DOWN_val_plot = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0/metilene_intersect_lifeline_plot_chr7_97012143_97022910_DOWN_val.pdf"
# DOWN_val_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0/metilene_intersect_lifeline_plot_chr7_97012143_97022910_DOWN_val.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
# project = "TCGA-LUSC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female_male"
# cutoff = "cutoff_8"
# threshold = "threshold_0"
# DMR = "chr7_97012143_97022910"

# # snakemake inputs:
# meta_table = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv"
# start_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr6_166530702_166530887.tsv"
# summary = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/summary_for_metilene.tsv"
# summary_complement = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/summary_for_metilene_complement.tsv"
# annot_file = "/scr/palinca/gabor/TCGA-pipeline_2/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# annot_file_2 = "../shared/resources/annot_from_betafile.tsv.gz"
# script_file = "../tcga_metilene/scripts/create_lifeline_plots_validation.py"
# # snakemake output:
# UP_val_plot = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr6_166530702_166530887_UP_val.pdf"
# UP_val_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr6_166530702_166530887_UP_val.tsv"
# DOWN_val_plot = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr6_166530702_166530887_DOWN_val.pdf"
# DOWN_val_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr6_166530702_166530887_DOWN_val.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
# project = "TCGA-CESC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_0"
# threshold = "threshold_0"
# DMR = "chr6_166530702_166530887"

# # snakemake inputs:
# meta_table = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv"
# start_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr14_85529381_85531409.tsv"
# summary = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/summary_for_metilene.tsv"
# summary_complement = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/summary_for_metilene_complement.tsv"
# annot_file = "/scr/palinca/gabor/TCGA-pipeline_2/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# annot_file_2 = "../shared/resources/annot_from_betafile.tsv.gz"
# script_file = "../tcga_metilene/scripts/create_lifeline_plots_validation.py"
# # snakemake output:
# UP_val_plot = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr14_85529381_85531409_UP_val.pdf"
# UP_val_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr14_85529381_85531409_UP_val.tsv"
# DOWN_val_plot = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr14_85529381_85531409_DOWN_val.pdf"
# DOWN_val_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr14_85529381_85531409_DOWN_val.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
# project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_0"
# threshold = "threshold_5"
# DMR = "chr14_85529381_85531409"

###############################################################################
#                                 test input                                  #
###############################################################################

# # snakemake inputs:
# meta_table = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/merged_meta_files/cutoff_5/meta_info_druglist_merged_drugs_combined.tsv"
# start_tsv = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_5/metilene_intersect_lifeline_plot_chr2_38665703_38666377.tsv"
# summary = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/summary_for_metilene.tsv"
# summary_complement = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/summary_for_metilene_complement.tsv"
# annot_file = "/scr/palinca/gabor/TCGA-pipeline/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# script_file = "../tcga_metilene/scripts/create_lifeline_plots_validation.py"
# # snakemake output:
# UP_val_plot = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_5/metilene_intersect_lifeline_plot_chr2_38665703_38666377_UP_val.pdf"
# UP_val_tsv = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_5/metilene_intersect_lifeline_plot_chr2_38665703_38666377_UP_val.tsv"
# DOWN_val_plot = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_5/metilene_intersect_lifeline_plot_chr2_38665703_38666377_DOWN_val.pdf"
# DOWN_val_tsv = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_5/metilene_intersect_lifeline_plot_chr2_38665703_38666377_DOWN_val.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline"
# project = "TCGA-CESC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_5"
# threshold = "threshold_5"
# DMR = "chr2_38665703_38666377"

# # # #####################
# meta_table = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/merged_meta_files/cutoff_5/meta_info_druglist_merged_drugs_combined.tsv"
# start_tsv = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0/metilene_intersect_lifeline_plot_chr7_94654961_94657930.tsv"
# summary = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/summary_for_metilene.tsv"
# summary_complement = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/summary_for_metilene_complement.tsv"
# annot_file = "/scr/palinca/gabor/TCGA-pipeline/metadata_processed/gencode.v36.annotation.gtf_genes.gz"
# script_file = "../tcga_metilene/scripts/create_lifeline_plots_validation.py"
# # snakemake output:
# UP_val_plot = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0/metilene_intersect_lifeline_plot_chr7_94654961_94657930_UP_val.pdf"
# UP_val_tsv = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0/metilene_intersect_lifeline_plot_chr7_94654961_94657930_UP_val.tsv"
# DOWN_val_plot = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0/metilene_intersect_lifeline_plot_chr7_94654961_94657930_DOWN_val.pdf"
# DOWN_val_tsv = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0/metilene_intersect_lifeline_plot_chr7_94654961_94657930_DOWN_val.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline"
# project = "TCGA-CESC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female_male"
# cutoff = "cutoff_5"
# threshold = "threshold_0"
# DMR = "chr7_94654961_94657930"

def write_empty_files():
    for path in [UP_val_tsv, DOWN_val_tsv]:
        pd.DataFrame(columns=['group', 'vital_status', 'drugs', 'gender',
                              'projects', 'beta_values', 'T', 'E',
                              'in_therapy', 'p_value', 'start', 'chr',
                              'threshold', 'fst_life_mean', 'scnd_life_mean',
                              'plot_type', 'mean_median', 'ENSG', 'gene_type',
                              'gene_status',
                              'gene_name']).to_csv(path, sep='\t', index=False)
        print(f'writing empty file {path}')
    for path in [UP_val_plot, DOWN_val_plot]:
        open(path, 'a').close()
        print(f'writing empty file {path}')
    os._exit(0)
# in case the start plot is empty, write empty files:
if pd.read_table(start_tsv).empty:
    write_empty_files()


# input files dependent variables:
start = pd.read_table(start_tsv)['start'].value_counts().index[0]
chr_ = DMR.split('_')[0]
# read in all beta values at the desired start position:
summary_DF = pd.read_table(summary, na_values='-').set_index(['Chromosome', 'Start']).dropna(how='all').sort_index().loc[(chr_,start),:]
summary_complement_DF = pd.read_table(summary_complement, na_values='-').set_index(['Chromosome', 'Start']).dropna(how='all').sort_index().loc[(chr_,start),:]
DF_summary = pd.concat([summary_DF, summary_complement_DF]).to_frame().T
# # make a multiindex of the vital_status;case_id;PROJECT;DRUGS header, s.t. in
col_t = [tuple(x) for x in [i.split(';') for i in DF_summary.columns]]
MI = pd.MultiIndex.from_tuples(col_t, names=('vital_status', 'case_id', 'drugs', 'gender', 'projects'))
# # can be read on that basis by the following pandas methods
DF_summary.columns=MI

# DF_meta path: '/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC/metilene/merged_meta_files/cutoff_2/meta_info_druglist_merged_drugs_combined.tsv'
DF_meta = pd.read_table(meta_table, usecols=['bcr_patient_uuid', 'survivaltime', 'years_to_last_follow_up', 'vital_status', 'pharmaceutical_therapy_drug_name'])
# vital_status                                                                         dead
# case_id                                              b726f18b-b996-4978-9f27-b3a540e06270
# drugs                                                                         carboplatin
# gender                                                                             female
# projects                                                                        TCGA-CESC
# Chromosome Start    End      region
# chr16      51150651 51150652 chr16_51150650_51153014                             0.593852


# mark beta values which shall be excluded:
thresh = float(threshold.split('_')[1])
# set every beta value which is to close to the median to pd.NA:

def apply_thresh(row):
    alive_median = row.loc[('alive', slice(None), slice(None), slice(None), slice(None))].median()
    dead_median = row.loc[('dead', slice(None), slice(None), slice(None), slice(None))].median()
    mean_median = st.mean([alive_median, dead_median])

    limit_val = (thresh/100) * mean_median
    upper_limit = mean_median + limit_val
    lower_limit = mean_median - limit_val
    row_thresh = row.apply(lambda x: pd.NA if float(x) < upper_limit and float(x) > lower_limit else x)
    return row_thresh

DF_summary = DF_summary.apply(apply_thresh, axis=1)
# TODO check whether the dropna is needed here
# DF_summary.dropna(axis=1, inplace=True)
# no, not on the whole DF, when it comes to single plots for every position
# within the range, check on how lifeline deals with na values

# prepare the DF_meta table, do not limit on case_id but make sure we just
# include cases in the DF_meta and DF_metilene for which we have the T value
# available:
DF_meta['T'] = DF_meta['survivaltime'].combine_first(DF_meta['years_to_last_follow_up'])
T = DF_meta['T']
index_to_delete = DF_meta[T.isna()].index
# beta value cases with no T content are dropped:
DF_summary.drop(labels=index_to_delete, axis=1, level=1, inplace=True)
# also the meta is shortened
DF_meta = DF_meta[T.notna()]

DF_temp = DF_summary.T.dropna()
# median = DF_summary.median(axis=1).values[0]
alive_median = DF_summary.loc[:, ('alive', slice(None), slice(None), slice(None), slice(None))].median(axis=1).values[0]
dead_median = DF_summary.loc[:, ('dead', slice(None), slice(None), slice(None), slice(None))].median(axis=1).values[0]
mean_median = st.mean([alive_median, dead_median])
# TypeError: float() argument must be a string or a real number, not 'NAType'
# 1ba4e70c-692c-487a-abc5-36e5f0a03286 carboplatin,docetaxel male TCGA-LUSC <NA>
# '/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/summary_for_metilene.tsv'
DF_summary = DF_summary.dropna(axis=1)
temp_DF = DF_summary.apply(lambda x: 'UP' if float(x) > mean_median else 'DOWN').to_frame()
temp_DF.columns = DF_summary.index
final_DF = pd.concat([temp_DF, DF_temp], axis=1)
final_DF.columns = (final_DF.columns[0], 'beta_values')

final_DF = final_DF.reset_index().set_index('case_id')
DF_meta = DF_meta.set_index('bcr_patient_uuid')
final_DF = pd.concat([final_DF, DF_meta.iloc[:, -1]], axis=1, join='inner')
final_DF['E'] = final_DF['vital_status'].apply(lambda x: True if x=='dead' else False)
final_DF = final_DF.reset_index().rename({'index': 'case_id'}, axis=1)
# we need DOWN and UP values in the col (chr_, start) to be able to plot those
# groups, if the value_counts is not 2, write empty files
if len(final_DF[(chr_, start)].value_counts().index) != 2:
    write_empty_files()
    os._exit(0)

# catch the exception that we have just one representative of a group, that
# would lead to a series, handle that:
try:
    DOWN_DF = final_DF.set_index((chr_, start)).loc['DOWN', :].reset_index().set_index('case_id')
    # for later summaries, the col names must be identical, add the chr_start as
    # extra col:
    DOWN_DF.rename({(chr_, start): 'group'}, axis=1, inplace=True)
except Exception as e:
    print(f'found {e}, handle Series instead of DataFrame')
    DOWN_DF = final_DF.set_index((chr_, start)).loc['DOWN', :].to_frame().T.reset_index().set_index('case_id')
    DOWN_DF.rename({'index': 'group'}, axis=1, inplace=True)
try:
    UP_DF = final_DF.set_index((chr_, start)).loc['UP', :].reset_index().set_index('case_id')
    UP_DF.rename({(chr_, start): 'group'}, axis=1, inplace=True)
except Exception as e:
    print(f'found {e}, handle Series instead of DataFrame')
    UP_DF = final_DF.set_index((chr_, start)).loc['UP', :].to_frame().T.reset_index().set_index('case_id')
    UP_DF.rename({'index': 'group'}, axis=1, inplace=True)


def set_if_therapy(drug_invoked):
    temp_val = False
    for drug in DRUGS:
        if drug == drug_invoked:
            temp_val = True
    return temp_val

DOWN_DF['in_therapy'] = DOWN_DF['drugs'].apply(set_if_therapy)
UP_DF['in_therapy'] = UP_DF['drugs'].apply(set_if_therapy)
# both in_therapy cols must contain True and False values, otherwise no
# comparison in possible:
if len(UP_DF['in_therapy'].value_counts().index) != 2:
    write_empty_files()

if len(DOWN_DF['in_therapy'].value_counts().index) != 2:
    write_empty_files()

# also add the p_value to the tables and the plots:
################################## cox fitter:
# dfA = pd.DataFrame({'E': event_observed_A, 'T': durations_A, 'groupA': 1})
# dfB = pd.DataFrame({'E': event_observed_B, 'T': durations_B, 'groupA': 0})
# df = pd.concat([dfA, dfB])

# cph = CoxPHFitter().fit(df, 'T', 'E')
# cph.print_summary()
#############################################


# *** lifelines.exceptions.ConvergenceError: Convergence halted due to matrix inversion problems. Suspicion is high collinearity. Please see the following tips in the lifelines documentation: https://lifelines
# .readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-modelMatrix is singular.
try:
    cph_UP = CoxPHFitter().fit(UP_DF.loc[:, ['in_therapy', 'T', 'E']], 'T', 'E')
    p_value_UP = cph_UP.summary['p'].values[0]
except Exception as e:
    print(e)
    print('using logrank_test instead of CoxPHFitter:')
    results = logrank_test(UP_DF['T'][UP_DF['in_therapy']], UP_DF['T'][~UP_DF['in_therapy']],UP_DF['E'][UP_DF['in_therapy']], UP_DF['E'][~UP_DF['in_therapy']])
    p_value_UP = results.p_value
    # 0.028546733835376127

try:
    cph_DOWN = CoxPHFitter().fit(DOWN_DF.loc[:, ['in_therapy', 'T', 'E']], 'T', 'E')
    p_value_DOWN = cph_DOWN.summary['p'].values[0]
except Exception as e:
    print(e)
    print('using logrank_test instead of CoxPHFitter:')
    results = logrank_test(DOWN_DF['T'][DOWN_DF['in_therapy']], DOWN_DF['T'][~DOWN_DF['in_therapy']],DOWN_DF['E'][DOWN_DF['in_therapy']], DOWN_DF['E'][~DOWN_DF['in_therapy']])
    p_value_DOWN = results.p_value
# cph_DOWN = CoxPHFitter().fit(DOWN_DF.loc[:, ['in_therapy', 'T', 'E']], 'T', 'E')
# p_value_DOWN = cph_DOWN.summary['p'].values[0]

### plot up and down

T = DOWN_DF['T']
E = DOWN_DF['E']
therapy_bool = DOWN_DF['in_therapy']
fig, ax = plt.subplots(figsize=(8, 6))
kmf_THERAPY = KaplanMeierFitter()
kmf_THERAPY.fit(T[therapy_bool], E[therapy_bool], label='in therapy')
DOWN_in_therapy_life_mean = vars(kmf_THERAPY)['survival_function_'].iloc[:, 0].mean()
ax = kmf_THERAPY.plot_survival_function(ax=ax)

kmf_NO_THERAPY = KaplanMeierFitter()
kmf_NO_THERAPY.fit(T[~therapy_bool], E[~therapy_bool], label='not in therapy')
DOWN_not_in_therapy_life_mean = vars(kmf_NO_THERAPY)['survival_function_'].iloc[:, 0].mean()
kmf_NO_THERAPY.plot_survival_function(ax=ax)

add_at_risk_counts(kmf_THERAPY, kmf_NO_THERAPY, ax=ax)
p_value_DOWN_str = f'p_value_DOWN = {Decimal(str(p_value_DOWN)):.2e}'

ax.set_title(f'{p_value_DOWN_str}, threshold = {round(thresh)}\nDMR: {": ".join(DMR.split("_")[:2]) + "-" + DMR.split("_")[2]}\nStart: {start}\n{project}, {drug_combi}, {gender}, cutoff={cutoff}')

plt.tight_layout()
print(f'saving: {DOWN_val_plot}')
plt.savefig(DOWN_val_plot)
plt.close()


T = UP_DF['T']
E = UP_DF['E']
therapy_bool = UP_DF['in_therapy']
fig, ax = plt.subplots(figsize=(8, 6))
kmf_THERAPY = KaplanMeierFitter()
kmf_THERAPY.fit(T[therapy_bool], E[therapy_bool], label='in therapy')
UP_in_therapy_life_mean = vars(kmf_THERAPY)['survival_function_'].iloc[:, 0].mean()
ax = kmf_THERAPY.plot_survival_function(ax=ax)

kmf_NO_THERAPY = KaplanMeierFitter()
kmf_NO_THERAPY.fit(T[~therapy_bool], E[~therapy_bool], label='not in therapy')
UP_not_in_therapy_life_mean = vars(kmf_NO_THERAPY)['survival_function_'].iloc[:, 0].mean()
kmf_NO_THERAPY.plot_survival_function(ax=ax)

add_at_risk_counts(kmf_THERAPY, kmf_NO_THERAPY, ax=ax)
p_value_UP_str = f'p_value_UP = {Decimal(str(p_value_UP)):.2e}'
ax.set_title(f'{p_value_UP_str}, threshold = {round(thresh)}\nDMR: {": ".join(DMR.split("_")[:2]) + "-" + DMR.split("_")[2]}\nStart: {start}\n{project}, {drug_combi}, {gender}, cutoff={cutoff}')

plt.tight_layout()
print(f'saving: {UP_val_plot}')
plt.savefig(UP_val_plot)
plt.close()

# write the tables:
UP_DF['p_value'] = p_value_UP
UP_DF['start'] = str(start)
UP_DF['chr'] = chr_
UP_DF['threshold'] = thresh
UP_DF['fst_life_mean'] = UP_in_therapy_life_mean
UP_DF['scnd_life_mean'] = UP_not_in_therapy_life_mean
UP_DF['plot_type'] = 'UP_validation'
UP_DF['mean_median'] = mean_median
# # TODO include the ENSGs for the positions found
DF_annot = pd.read_table(annot_file)
# # limit the annotatin DF to the right chromosome:
DF_annot = DF_annot.set_index('chr').loc[chr_,:].reset_index()
# acces the ENSG:
DF_annot = DF_annot[(DF_annot['start'] <= start) & (DF_annot['stop'] > start)].loc[:, ["ENSG", "gene_type", "gene_status", "gene_name"]]
if DF_annot.empty:
    DF_annot_2 = pd.read_table(annot_file_2).set_index('chr').loc[chr_,:].reset_index()
    DF_annot_2 = DF_annot_2[(DF_annot_2['start'] <= start) & (DF_annot_2['stop'] > start)].loc[:, ["ENST", "gene_type", "gene_status", "gene_name"]]
    UP_DF['ENSG'] = DF_annot_2['ENST'].values[0]
    UP_DF['gene_type'] = DF_annot_2['gene_type'].values[0]
    UP_DF['gene_status'] = DF_annot_2['gene_status'].values[0]
    UP_DF['gene_name'] = DF_annot_2['gene_name'].values[0]
else:
    UP_DF['ENSG'] = DF_annot['ENSG'].values[0]
    UP_DF['gene_type'] = DF_annot['gene_type'].values[0]
    UP_DF['gene_status'] = DF_annot['gene_status'].values[0]
    UP_DF['gene_name'] = DF_annot['gene_name'].values[0]
UP_DF.to_csv(UP_val_tsv, sep='\t')
print(f'writing {UP_val_tsv}')

DOWN_DF['p_value'] = p_value_DOWN
DOWN_DF['start'] = str(start)
DOWN_DF['chr'] = chr_
DOWN_DF['threshold'] = thresh
DOWN_DF['fst_life_mean'] = DOWN_in_therapy_life_mean
DOWN_DF['scnd_life_mean'] = DOWN_not_in_therapy_life_mean
DOWN_DF['plot_type'] = 'DOWN_validation'
DOWN_DF['mean_median'] = mean_median

if DF_annot.empty:
    DF_annot_2 = pd.read_table(annot_file_2).set_index('chr').loc[chr_,:].reset_index()
    DF_annot_2 = DF_annot_2[(DF_annot_2['start'] <= start) & (DF_annot_2['stop'] > start)].loc[:, ["ENST", "gene_type", "gene_status", "gene_name"]]
    DOWN_DF['ENSG'] = DF_annot_2['ENST'].values[0]
    DOWN_DF['gene_type'] = DF_annot_2['gene_type'].values[0]
    DOWN_DF['gene_status'] = DF_annot_2['gene_status'].values[0]
    DOWN_DF['gene_name'] = DF_annot_2['gene_name'].values[0]
else:
    DOWN_DF['ENSG'] = DF_annot['ENSG'].values[0]
    DOWN_DF['gene_type'] = DF_annot['gene_type'].values[0]
    DOWN_DF['gene_status'] = DF_annot['gene_status'].values[0]
    DOWN_DF['gene_name'] = DF_annot['gene_name'].values[0]
DOWN_DF.to_csv(DOWN_val_tsv, sep='\t')
print(f'writing {DOWN_val_tsv}')
