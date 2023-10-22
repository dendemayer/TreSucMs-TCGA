import pandas as pd
import sys
import os
# from methyl import set_logger
# from lifelines.datasets import load_waltons
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from decimal import Decimal
# from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from lifelines import CoxPHFitter
import numpy as np
import statistics as st

"""
follow lifeline succession from patients with detected DMR
- needed, metadata:f.e.:
    - TCGA-pipelines_2/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/
            meta_info_druglist_merged_drugs_combined.tsv
    - cols:
        # survivaltime       | years_to_last_follow_up 'bcr_patient_uuid',
- the metilene intersect out table
    - TCGA-pipelines_2/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/
        carboplatin,paclitaxel_cisplatin/female/cutoff_0/
        metilene_complement_intersect.tsv
    # needed infos are vital_status (alive dead),

# use the logrank_test at every position within a DMR to select the best
# resolution between up and down :
# changing to the cox regression test: -> see lifeine hint on https://lifelines.readthedocs.io/en/latest/lifelines.statistics.html
# and further reading on: https://discourse.datamethods.org/t/when-is-log-rank-preferred-over-univariable-cox-regression/2344
If the p-value of the log-rank test is less than a chosen significance level
(such as 0.05), it suggests that there is a statistically significant
difference in survival between the groups, and the Kaplan-Meier estimate can be
considered meaningful.
python /homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/
    tcga_metilene/scripts/create_lifeline_plots.py

lifeline_out_tsv = "{output_path}/{project}/metilene/metilene_output/{drug_combi}/{gender}/{cutoff}/{threshold}/metilene_intersect_lifeline_plot_{DMR}.tsv",

include the mean median difference!:
"""
# # # # # #################################################
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

metilene_intersect = snakemake.input[0]
meta_table = snakemake.input[1]
annot_file = snakemake.input[2]
annot_file_2 = snakemake.input[3]
metilene_out = snakemake.input[4]
DMR = snakemake.wildcards.DMR
lifeline_out_pdf = snakemake.output.lifeline_out_pdf
lifeline_out_tsv = snakemake.output.lifeline_out_tsv
threshold = snakemake.wildcards.threshold
cutoff = snakemake.wildcards.cutoff.split('_')[1]
drug_combi = snakemake.wildcards.drug_combi.replace('_', ';')
project = ', '.join(snakemake.wildcards.project.split('_'))
gender = ', '.join(snakemake.wildcards.gender.split('_'))

print('# snakemake inputs:')
[print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

# # # #######################################################################

###############################################################################
# #                                 test_input                                  #
# ###############################################################################
# # snakemake inputs:
# metilene_intersect = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/metilene_intersect.tsv"
# meta_table = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-LUSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv"
# annot_file = "/scr/palinca/gabor/TCGA-pipeline_3/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# annot_file_2 = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/tcga_metilene/../shared/resources/annot_from_betafile.tsv.gz"
# script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/tcga_metilene/scripts/create_lifeline_plots.py"
# # snakemake output:
# lifeline_out_pdf = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr7_27141798_27144197.pdf"
# lifeline_out_tsv = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr7_27141798_27144197.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_3"
# project = "TCGA-LUSC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "male"
# cutoff = "cutoff_0"
# threshold = "threshold_5"
# DMR = "chr7_27141798_27144197"

# # snakemake inputs:
# metilene_intersect = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/metilene_intersect.tsv"
# meta_table = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-HNSC/metilene/merged_meta_files/cutoff_5/meta_info_druglist_merged_drugs_combined.tsv"
# annot_file = "/scr/palinca/gabor/TCGA-pipeline_3/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# annot_file_2 = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/tcga_metilene/../shared/resources/annot_from_betafile.tsv.gz"
# script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/tcga_metilene/scripts/create_lifeline_plots.py"
# # snakemake output:
# lifeline_out_pdf = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_20/metilene_intersect_lifeline_plot_chr11_122099576_122101477.pdf"
# lifeline_out_tsv = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_20/metilene_intersect_lifeline_plot_chr11_122099576_122101477.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_3"
# project = "TCGA-HNSC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female_male"
# cutoff = "cutoff_5"
# threshold = "threshold_20"
# DMR = "chr11_122099576_122101477"

# # snakemake inputs:
# metilene_intersect = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_intersect.tsv"
# meta_table = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv"
# annot_file = "/scr/palinca/gabor/TCGA-pipeline_3/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# annot_file_2 = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/tcga_metilene/../shared/resources/annot_from_betafile.tsv.gz"
# script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/tcga_metilene/scripts/create_lifeline_plots.py"
# # snakemake output:
# lifeline_out_pdf = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr2_176150267_176153069_temp.pdf"
# lifeline_out_tsv = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr2_176150267_176153069_temp.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_3"
# project = "TCGA-CESC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_0"
# threshold = "threshold_5"
# DMR = "chr2_176150267_176153069"

# snakemake inputs:
# metilene_intersect = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_intersect.tsv"
# meta_table = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv"
# annot_file = "/scr/palinca/gabor/TCGA-pipeline_2/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# annot_file_2 = "../shared/resources/annot_from_betafile.tsv.gz"
# script_file = "../tcga_metilene/scripts/create_lifeline_plots.py"
# # snakemake output:
# lifeline_out_pdf = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr7_83648624_83648822.pdf"
# lifeline_out_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr7_83648624_83648822.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
# project = "TCGA-CESC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_0"
# threshold = "threshold_0"
# DMR = "chr7_83648624_83648822"

# # snakemake inputs:
# metilene_intersect = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_intersect.tsv"
# meta_table = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv"
# annot_file = "/scr/palinca/gabor/TCGA-pipeline_2/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# annot_file_2 = "../shared/resources/annot_from_betafile.tsv.gz"
# metilene_out = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_qval.0.05.out"
# script_file = "../tcga_metilene/scripts/create_lifeline_plots.py"
# # snakemake output:
# lifeline_out_pdf = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr14_85529381_85531409.pdf"
# lifeline_out_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr14_85529381_85531409.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
# project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_0"
# threshold = "threshold_0"
# DMR = "chr14_85529381_85531409"


# # snakemake inputs:
# metilene_intersect = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_intersect.tsv"
# meta_table = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv"
# annot_file = "/scr/palinca/gabor/TCGA-pipeline_2/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# annot_file_2 = "../shared/resources/annot_from_betafile.tsv.gz"
# metilene_out = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_qval.0.05.out"
# script_file = "../tcga_metilene/scripts/create_lifeline_plots.py"
# # snakemake output:
# lifeline_out_pdf = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr5_66828566_66828775.pdf"
# lifeline_out_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr5_66828566_66828775.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
# project = "TCGA-CESC_TCGA-HNSC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_0"
# threshold = "threshold_0"
# DMR = "chr5_66828566_66828775"
###############################################################################
#                                 test_input                                  #
###############################################################################


### triggers:
    # DF['ENSG'] = DF_annot['ENSG'].values[0]
                 # ~~~~~~~~~~~~~~~~~~~~~~~^^^
# IndexError: index 0 is out of bounds for axis 0 with size 0
# -> do not limit the processed gtf on ENSEMBL:
# bedtools intersect -b <(echo "chr20\t58856147\t58856148") -a /scr/palinca/gabor/TCGA-pipeline/metadata/gencode.v36.annotation.gtf.gz
# chr20   HAVANA  gene    58856148        58856148        .       +       .       gene_id "ENSG00000087460.22"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "GNAS";  level 2; tag "ncRNA_host"; havana_gene "OTTHUMG00000033069.19";
# grep ENSG00000087460 /scr/palinca/gabor/TCGA-pipeline/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz                                                                                                         ±[●●][main]
# -> those coordinates can be found in original annotation, but not in the
# ENSEMBL limited processed one

# :
# KeyError: 'chr12_132887020_132888306'

# # snakemake inputs:
# metilene_intersect = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/metilene_intersect.tsv"
# meta_table = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv"
# annot_file = "/scr/palinca/gabor/TCGA-pipeline/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# script_file = "../tcga_metilene/scripts/create_lifeline_plots.py"
# # snakemake output:
# lifeline_out_pdf = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr20_58850101_58856248.pdf"
# lifeline_out_tsv = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr20_58850101_58856248.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline"
# project = "TCGA-CESC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female_male"
# cutoff = "cutoff_0"
# threshold = "threshold_0"
# DMR = "chr20_58850101_58856248"

# (base) [gabor@palinca shared]$ bedtools intersect -b <(echo "chr11\t115505380\t115505381") -a <(awk 'NR>1{print $3"\t"$4"\t"$5"\t"$0}' /scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/data_files/jhu-usc.edu_CESC.HumanMethylation450.1.lvl-3.TCGA-C5-A1BE-01B-11D-A13
# Z-05.gdc_hg38.txt | grep -v "*")
# chr11   115505380       115505381       cg01948062      0.0948425263590688      chr11   115505380       115505381       CADM1;CADM1;CADM1;CADM1;CADM1;CADM1;CADM1;CADM1;CADM1;CADM1;CADM1       protein_coding;protein_coding;protein_coding;protein_coding;protein_coding;
# protein_coding;protein_coding;protein_coding;protein_coding;protein_coding;protein_coding       ENST00000331581.9;ENST00000452722.6;ENST00000536727.4;ENST00000536781.1;ENST00000537058.4;ENST00000537140.4;ENST00000540951.1;ENST00000541434.4;ENST00000542447.5;ENST00000
# 543249.1;ENST00000545380.4      -814;-964;-964;-964;-964;-856;-964;-985;-856;-422;-990  CGI:chr11:115502856-115505163   S_Shore
# # snakemake inputs:
# metilene_intersect = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/metilene_intersect.tsv"
# meta_table = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/merged_meta_files/cutoff_5/meta_info_druglist_merged_drugs_combined.tsv"
# annot_file = "/scr/palinca/gabor/TCGA-pipeline/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# script_file = "../tcga_metilene/scripts/create_lifeline_plots.py"
# # snakemake output:
# lifeline_out_pdf = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_5/metilene_intersect_lifeline_plot_chr11_115503038_115505464.pdf"
# lifeline_out_tsv = "/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_5/metilene_intersect_lifeline_plot_chr11_115503038_115505464.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline"
# project = "TCGA-CESC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female_male"
# cutoff = "cutoff_5"
# threshold = "threshold_5"
# DMR = "chr11_115503038_115505464"
# # ->>>
# > /homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/.snakemake/scripts/tmppd019qrh.create_lifeline_plots.py(192)<module>()
# -> print(e)
# (Pdb) e
# ValueError('Values must be numeric: no strings, datetimes, objects, etc.')
# DONE with
# if DF_plot.value_counts().index.nunique() != 2: instead of --> that covers
# if DF_plot.value_counts().index.nunique() == 1:
# the cases we have just UP or just DOWN or an emtpy Series

#     results = logrank_test(T[UP_bool], T[DOWN_bool], E[UP_bool], E[DOWN_bool])
#                            ~^^^^^^^^^
#     raise ValueError("Cannot index with multidimensional key")
# ValueError: Cannot index with multidimensional key


DF_meta = pd.read_table(meta_table, usecols=['bcr_patient_uuid', 'survivaltime', 'years_to_last_follow_up', 'vital_status'])
DF_metilene = pd.read_table(metilene_intersect, header=[0, 1, 2, 3, 4], index_col=[0, 1, 2, 3], na_values='.').dropna()
# limit the DF_metilene to the recent range:
DF_metilene = DF_metilene.loc[(slice(None), slice(None), slice(None), DMR), :]

# the threshold can be applied here right away since the DF_metilene is already
# filtered on the correct cases:

# threshold % over or under the median
thresh = float(threshold.split('_')[1])
# add the mean of medians directly into the metielen DF:
DF_dead_median = DF_metilene.loc[:, ('dead',slice(None),slice(None),slice(None),slice(None))].median(axis=1)
DF_alive_median = DF_metilene.loc[:, ('alive',slice(None),slice(None),slice(None),slice(None))].median(axis=1)
DF_mean_of_median = pd.concat([DF_alive_median, DF_dead_median], axis=1).mean(axis=1)
DF_mean_of_median.name = 'mean_median'
new_index = DF_metilene.index.to_frame().join(DF_mean_of_median.to_frame()).set_index('mean_median', append=True).index
DF_metilene.index = new_index

def apply_thresh(row):
    mean_median = row.name[-1]
    limit_val = (thresh/100) * mean_median
    upper_limit = mean_median + limit_val
    lower_limit = mean_median - limit_val
    row_thresh = row.apply(lambda x: pd.NA if float(x) < upper_limit and float(x) > lower_limit else x)
    return row_thresh

DF_metilene = DF_metilene.apply(apply_thresh, axis=1)
# # throw away positions where after thresh appliance the start positions just hold
# # NA

DF_metilene.dropna(how='all', inplace=True)

#############################################
# ### IMPORTANT T check!!!
# get the used cases out of the meta_table:
used_cases = [i[1] for i in DF_metilene.columns]
DF_meta = DF_meta.set_index('bcr_patient_uuid').loc[used_cases, :]
DF_meta['T'] = DF_meta['survivaltime'].combine_first(DF_meta['years_to_last_follow_up'])

T = DF_meta['T']
# DF_meta['E'] = DF_meta['vital_status'] == 'dead'
DF_meta['E'] = np.where(DF_meta['vital_status'].values == 'dead', True, False)
# in case not value available for survivaltime AND years_to_last_follow_up, a
# NaN is set in T, :
# '/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_2/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin_paclitaxel/female_male/cutoff_0/metilene_complement_intersect.tsv'
# on case_id:
                                     # vital_status  survivaltime  years_to_last_follow_up   T      E
# bcr_patient_uuid
# 31d8591b-9d80-406c-b61a-7a6544c73466        alive           NaN                      NaN NaN  False
# 7625a27d-a012-4f4d-b637-5f85d7cdae8c         dead           NaN                      NaN NaN   True

### this is not threshold related, its a check, whether we have all T values
### available for the set we are interested in !!!
### therefore the T is created beforehand
# shrink the complete DF_meta table, then extract T and E:
index_to_delete = DF_meta[T.isna()].index
# since we use also DF_metilene we have to reduce it to the cases for which the
# representing whether the “death” was observed or not (alternatively an
# individual can be censored)
# T and E infos are available:
DF_metilene.drop(labels=index_to_delete, axis=1, level=1, inplace=True)
DF_meta = DF_meta[T.notna()]
### this is not threshold related, its a check, whether we have all T values
### available for the set we are interested in !!!
#### IMPORTANT T check!!!
#################################################################

######################### revision, the start position is not estimated with
# cox regression for every start ->

# the groups are representing the cases which are either higher or lower
# methylated than the median of all cases within the DMR, therefore calculate
# the kmf at every postition and take the one with the highest p_value
starts = [i[1] for i in DF_metilene.index]
p_value = 0
pvalue_list = []
DF_starts = []

for start in starts:
    # threshold is already apllied for each start position, the UP and DOWN
    # grouping is performed here
    # limit DF to the start:
    DF_metilene_temp = DF_metilene.loc[ (slice(None), start, slice(None), slice(None)), :]
    # drop beta values which are marked as nans by threshold:
    DF_metilene_temp = DF_metilene_temp.dropna(axis=1)
    # get the mean_median, which was set in threshold def apply_thresh(row):
    mean_median = DF_metilene_temp.index.values[-1][-1]
    cases_kept = [ i[1] for i in DF_metilene_temp.columns]
    # limit the DF_mete to the used cases:
    DF_meta_temp = DF_meta.loc[cases_kept, :]
    # # also T and E must be adjusted here correctly,
    T = DF_meta_temp['T']
    E = DF_meta_temp['vital_status'] == 'dead' #  E is a either boolean or binary array

    # S_temp is supposed to be a Series, since we access one start postitin
    # within the DMR:
    S_temp = DF_metilene_temp.loc[(slice(None), start, slice(None), slice(None)), :].apply(lambda x: 'UP' if float(x) > mean_median else 'DOWN')
    # UP_bool = (S_temp == 'UP').values
    # DOWN_bool = (S_temp == 'DOWN').values
    DF_start = S_temp.to_frame().reset_index()
    DF_start['T'] = T.values
    DF_start['E'] = E.values
    DF_start['UP'] = DF_start[0].apply(lambda x: True if x=='UP' else False)
    try:
        # results = logrank_test(T[UP_bool], T[DOWN_bool], E[UP_bool], E[DOWN_bool])
        # changing to the cox regression test: -> see lifeine hint on https://lifelines.readthedocs.io/en/latest/lifelines.statistics.html
        # and further reading on: https://discourse.datamethods.org/t/when-is-log-rank-preferred-over-univariable-cox-regression/2344
        cph = CoxPHFitter().fit(DF_start.loc[:, ['UP', 'T', 'E']], 'T', 'E')
        p_value = cph.summary['p'].values[0]
        # p_value = results.p_value
        pvalue_list.append(p_value)
        DF_start['start'] = start
        DF_start['p_value'] = p_value
        DF_start['mean_median'] = mean_median
        DF_starts.append(DF_start)
    except Exception as e:
        print(e)
        # if thats not feasable we set an p_value of 1:
        pvalue_list.append(1)
        DF_starts.append(pd.DataFrame())
        # lifelines.exceptions.ConvergenceError: Convergence halted due to
        # matrix inversion problems. Suspicion is high collinearity. Please see
        # the following tips in the lifelines documentation:
        # https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence
        # -in-the-cox-proportional-hazard-modelMatrix is singular.
        # -- in that case, every start position is problematic, because the
        # events (E) holds True or False for every position of the case set,
        # write empty files in that case as well:
        # open(lifeline_out_pdf, 'a').close()
        # pd.DataFrame(columns=['case_id', 'drugs', 'gender', 'projects',
        #                       'UP_or_DOWN', 'beta_value', 'vital_status',
        #                       'survivaltime', 'years_to_last_follow_up', 'T',
        #                       'E', 'mean_median', 'DMR', 'p_value', 'start', 'chr',
        #                       'threshold', 'fst_life_mean', 'scnd_life_mean',
        #                       'plot_type', 'ENSG', 'gene_type', 'gene_status',
        #                       'gene_name']).to_csv(lifeline_out_tsv, sep='\t',
        #                                            index=False)
        # os._exit(0)


# take the lowest p_value, testet through all postitions:
if min(pvalue_list) == 1:
    open(lifeline_out_pdf, 'a').close()
    pd.DataFrame(
        columns=['case_id', 'drugs', 'gender', 'projects', 'UP_or_DOWN',
                  'beta_value', 'vital_status', 'survivaltime',
                  'years_to_last_follow_up', 'T', 'E', 'mean_median', 'DMR',
                  'p_value', 'start', 'chr', 'threshold', 'fst_life_mean',
                  'scnd_life_mean', 'plot_type', 'ENSG', 'gene_type',
                  'gene_status', 'gene_name']).to_csv(lifeline_out_tsv,
                                                      sep='\t', index=False)
    os._exit(0)
plot_index = pvalue_list.index(min(pvalue_list))
# get the best start DF through the minimum pvalue :
DF_start = DF_starts[plot_index]
DF_start.rename({0: 'UP_or_DOWN'}, axis=1, inplace=True)
# p_value = pvalue_list[plot_index]
# DF_beta = DF_metilene.iloc[ plot_index, :].dropna()
# alive_median = DF_beta.loc[('alive', slice(None), slice(None), slice(None), slice(None))].median()
# dead_median = DF_beta.loc[('dead', slice(None), slice(None), slice(None), slice(None))].median()
# mean_median = st.mean([alive_median, dead_median])
# DF_plot = DF_beta.apply(lambda x: 'UP' if float(x) > mean_median else 'DOWN')

p_value = DF_start['p_value'].values[0]
p_value_str = f'p_value = {Decimal(str(min(pvalue_list))):.2e}'
# again, consider the threshold dependend betavalue which are maybe set to NA:
# cases_to_delete = [i[1] for i in DF_metilene.loc[
    # :, DF_metilene.isna().values[0]].columns]


# kmf.fit(T, event_observed=E)  # or, more succinctly, kmf.fit(T, E)

# groups = DF_plot['UP_DOWN']
# in case we cannot create 2 groups (i.e., with the thresh invokening, just UP
# cases are
# present:/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_2/TCGA-HNSC/metilene/metilene_output/carboplatin,paclitaxel_cisplatin_paclitaxel/female_male/cutoff_0/metilene_intersect.tsv
# leads to just up cases with thresh =10: DF_metilene.iloc[plot_index,
# :].median()
# 0.0204743706513952
# (Pdb) limit_val
# 0.0445787454608645
# just write 2 emtpy files:

if DF_start['UP_or_DOWN'].value_counts().index.nunique() != 2:
    open(lifeline_out_pdf, 'a').close()
    pd.DataFrame(
        columns=['case_id', 'drugs', 'gender', 'projects', 'UP_or_DOWN',
                  'beta_value', 'vital_status', 'survivaltime',
                  'years_to_last_follow_up', 'T', 'E', 'mean_median', 'DMR',
                  'p_value', 'start', 'chr', 'threshold', 'fst_life_mean',
                  'scnd_life_mean', 'plot_type', 'ENSG', 'gene_type',
                  'gene_status', 'gene_name']).to_csv(lifeline_out_tsv,
                                                      sep='\t', index=False)
    os._exit(0)
else:
    # cases_to_keep = [i[1] for i in DF_plot.index]
    # DF_meta = DF_meta.loc[cases_to_keep]
    T = DF_start['T']
    E = DF_start['E']
    UP_bool = (DF_start['UP_or_DOWN'] == 'UP').values
    DOWN_bool = (DF_start['UP_or_DOWN'] == 'DOWN').values
    fig, ax = plt.subplots(figsize=(8, 6))
    kmf_UP = KaplanMeierFitter()
    kmf_UP.fit(T[UP_bool], E[UP_bool], label='UP')
    UP_life_expectancy_mean = vars(kmf_UP)['survival_function_'].iloc[:, 0].mean()
    ax = kmf_UP.plot_survival_function(ax=ax)

    kmf_DOWN = KaplanMeierFitter()
    kmf_DOWN.fit(T[DOWN_bool], E[DOWN_bool], label='DOWN')
    DOWN_life_expectancy_mean = vars(kmf_DOWN)['survival_function_'].iloc[:, 0].mean()
    kmf_DOWN.plot_survival_function(ax=ax)

    add_at_risk_counts(kmf_UP, kmf_DOWN, ax=ax)
    ax.set_title(
        f'{p_value_str}, threshold = {threshold.split("_")[1]} \nDMR: {": ".join(DMR.split("_")[:2]) + "-" + DMR.split("_")[2]}, Start: {str(starts[plot_index])}\n{project}, {drug_combi}, {gender}, cutoff={cutoff}')

    plt.tight_layout()
    print(f'saving: {lifeline_out_pdf}')
    plt.savefig(lifeline_out_pdf)
    plt.close()

    # DF = pd.concat([kmf_UP.survival_function_, kmf_DOWN.survival_function_], axis=1)
    # DF_beta = DF_beta.to_frame()
    # DF_beta.columns = ['beta_value']
    # DF_beta = DF_beta.reset_index(level=[0, 2, 3, 4], drop=True)
    # DF_plot = DF_plot.to_frame()
    # DF_plot.columns = ['UP_or_DOWN']
    # DF_plot = DF_plot.reset_index(level=[0, 2, 3, 4]).drop('vital_status', axis=1)
    # DF = pd.concat([DF_plot, DF_beta, DF_meta], axis=1)

    DF_meta.index.name='case_id'
    # add survivaltime and years_to_last_follow_up to the start DF
    DF_start = DF_start.set_index('case_id').join(DF_meta.loc[:, ['survivaltime', 'years_to_last_follow_up']])
    start = starts[plot_index]
    DF_beta = DF_metilene.loc[(slice(None), start, slice(None), slice(None)),:].reset_index(drop=True).T.reset_index().loc[:, ['case_id', 0]].set_index('case_id')
    DF_beta.columns = ['beta_value']
    DF_start = DF_start.join(DF_beta)
    DF_start['DMR'] = DMR
    chrom = DMR.split('_')[0]
    DF_start['chr'] = chrom
    DF_start['threshold'] = thresh
    DF_start['fst_life_mean'] = UP_life_expectancy_mean
    DF_start['scnd_life_mean'] = DOWN_life_expectancy_mean
    DF_start['plot_type'] = 'base_plot'
    DF_start = DF_start.reset_index()
    #####
    ### cannot find :
    # (Pdb) chrom 'chr11' (Pdb) start 115505380 -> DF['ENSG'] = DF_annot['ENSG'].values[0]
    # this can be found in the betavalue tables and its ENST:
     # Chromosome Composite Element REF  Beta_value      Start  ...                                      Transcript_ID                                    Position_to_TSS                 CGI_Coordinate Feature_Type
    # 1614      chr11            cg01948062    0.094843  115505380  ...  ENST00000331581.9;ENST00000452722.6;ENST000005...  -814;-964;-964;-964;-964;-856;-964;-985;-856;-...  CGI:chr11:115502856-115505163      S_Shore
    # temp_data_file = '/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/data_files/jhu-usc.edu_CESC.HumanMethylation450.1.lvl-3.TCGA-C5-A1BE-01B-11D-A13Z-05.gdc_hg38.txt'
    # temp_beta_DF = pd.read_table(temp_data_file)
    # temp_beta_DF = temp_beta_DF.set_index('Chromosome').loc[chrom, :].reset_index()
    # temp_beta_DF[(temp_beta_DF['Start'] <= start) & (temp_beta_DF['End'] > start)]
    # chr11| 115169218| 115504957 | ENSG00000182985 | protein_coding     | KNOWN       | CADM1         ║
    # /scr/palinca/gabor/TCGA-pipeline/metadata_processed/gencode.v36.annotation.gtf_genes.gz
    # this does not lie withing the range of this gene, but when including the
    # ENST field (ENST00000331581) and relate back to the belonging ENSG, it
    # can be found --> ind preprocess, do not filter just on gene, but include
    # on transcript aswell,
    # gene level, with it, also the ENSG00000182985 must be present
    # "transcript" values for:
    # zgrep ENST00000331581 /scr/palinca/gabor/TCGA-pipeline/metadata/gencode.v36.annotation.gtf.gz:
    # exon       CDS        UTR        start_codon stop_codon
    # handle the case if the start lies outside of the ENSG range and report
    # the transcript value
    ###### #
    # # TODO include the ENSGs for the positions found
    # overkill with bedtool,just one position must be found:
    # annot_BT = BedTool.from_dataframe(DF_annot)
    # DF_BT = BedTool.from_dataframe(DF)
    # ####
    DF_annot = pd.read_table(annot_file)
    # # limit the annotatin DF to the right chromosome:
    DF_annot = DF_annot.set_index('chr').loc[chrom,:].reset_index()
    # acces the ENSG:
    DF_annot = DF_annot[(DF_annot['start'] <= start) & (DF_annot['end'] > start)].loc[:, ["ENSG", "gene_type", "gene_status", "gene_name"]]
    if DF_annot.empty:
        DF_annot_2 = pd.read_table(annot_file_2).set_index('chr').loc[chrom,:].reset_index()
        DF_annot_2 = DF_annot_2[(DF_annot_2['start'] <= start) & (DF_annot_2['stop'] > start)].loc[:, ["ENST", "gene_type", "gene_status", "gene_name"]]
        DF_start['ENSG'] = DF_annot_2['ENST'].values[0]
        DF_start['gene_type'] = DF_annot_2['gene_type'].values[0]
        DF_start['gene_status'] = DF_annot_2['gene_status'].values[0]
        DF_start['gene_name'] = DF_annot_2['gene_name'].values[0]
    else:
        DF_start['ENSG'] = DF_annot['ENSG'].values[0]
        DF_start['gene_type'] = DF_annot['gene_type'].values[0]
        DF_start['gene_status'] = DF_annot['gene_status'].values[0]
        DF_start['gene_name'] = DF_annot['gene_name'].values[0]
    print(f'saving: {lifeline_out_tsv}')
    DF_start = DF_start.loc[:,['case_id', 'drugs', 'gender', 'projects', 'UP_or_DOWN', 'beta_value', 'vital_status', 'survivaltime', 'years_to_last_follow_up', 'T', 'E', 'mean_median', 'DMR', 'p_value', 'start', 'chr', 'threshold', 'fst_life_mean', 'scnd_life_mean', 'plot_type', 'ENSG', 'gene_type', 'gene_status', 'gene_name']]
    DF_start.to_csv(lifeline_out_tsv, sep='\t', index=None)

