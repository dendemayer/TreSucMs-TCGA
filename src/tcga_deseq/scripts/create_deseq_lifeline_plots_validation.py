from matplotlib.patches import bbox_artist
import pandas as pd
import numpy as np
import sys
import os
# from methyl import set_logger
# from lifelines.datasets import load_waltons
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from decimal import Decimal
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from lifelines import CoxPHFitter
import statistics as st
# from pandas.core import apply

#######################################################################
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

meta_table = snakemake.input.meta_table
deseq_lifeline_tsv = snakemake.input.deseq_lifeline_tsv
counts = snakemake.input.counts
threshold = snakemake.wildcards.threshold
deseq_lifeline_pdf_UP = snakemake.output.deseq_lifeline_pdf_UP
deseq_lifeline_pdf_DOWN = snakemake.output.deseq_lifeline_pdf_DOWN
deseq_lifeline_tsv_UP = snakemake.output.deseq_lifeline_tsv_UP
deseq_lifeline_tsv_DOWN = snakemake.output.deseq_lifeline_tsv_DOWN
dropped_cases = snakemake.output.dropped_cases
ENSG = snakemake.wildcards.ENSG
count_type = snakemake.wildcards.count_type + '_count'
DRUGS = snakemake.wildcards.drug_combi.split('_')

drug_combi = snakemake.wildcards.drug_combi.replace('_', ';')
cutoff = snakemake.wildcards.cutoff.split('_')[1]
project = ', '.join(snakemake.wildcards.project.split('_'))
gender = ', '.join(snakemake.wildcards.gender.split('_'))

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]
#######################################################################

###############################################################################
#                                 test input                                  #
###############################################################################
# # snakemake inputs:
# meta_table = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv"
# deseq_lifeline_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000249235_lifeline.tsv.gz"
# counts = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/DESeq2_norm_counts_all_cases.gz"
# script_file = "../tcga_deseq/scripts/create_deseq_lifeline_plots_validation.py"
# # snakemake output:
# deseq_lifeline_tsv_UP = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000249235_lifeline_UP_val.tsv.gz"
# deseq_lifeline_tsv_DOWN = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000249235_lifeline_DOWN_val.tsv.gz"
# deseq_lifeline_pdf_UP = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000249235_lifeline_UP_val.pdf"
# deseq_lifeline_pdf_DOWN = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000249235_lifeline_DOWN_val.pdf"
# dropped_cases = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000249235_dropped_cases.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
# project = "TCGA-CESC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female_male"
# cutoff = "cutoff_0"
# threshold = "threshold_0"
# in_de = "DE"
# count_type = "norm"
# ENSG = "ENSG00000249235"
# DRUGS = drug_combi.split('_')
# # snakemake inputs:
# meta_table = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv"
# deseq_lifeline_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000081041_lifeline.tsv.gz"
# counts = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/DESeq2_norm_counts_all_cases.gz"
# script_file = "../tcga_deseq/scripts/create_deseq_lifeline_plots_validation.py"
# # snakemake output:
# deseq_lifeline_tsv_UP = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000081041_lifeline_UP_val.tsv.gz"
# deseq_lifeline_tsv_DOWN = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000081041_lifeline_DOWN_val.tsv.gz"
# deseq_lifeline_pdf_UP = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000081041_lifeline_UP_val.pdf"
# deseq_lifeline_pdf_DOWN = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000081041_lifeline_DOWN_val.pdf"
# dropped_cases = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000081041_dropped_cases.tsv"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
# project = "TCGA-CESC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_0"
# threshold = "threshold_0"
# in_de = "DE"
# count_type = "norm"
# ENSG = "ENSG00000081041"
# DRUGS = drug_combi.split('_')

###############################################################################
#                                 test input                                  #
###############################################################################


def write_empty_files():
    for path in [deseq_lifeline_tsv_UP, deseq_lifeline_tsv_DOWN]:
        pd.DataFrame(columns=['count', 'project_id',
                              'pharmaceutical_therapy_drug_name', 'gender',
                              'vital_status', 'survivaltime',
                              'years_to_last_follow_up', 'T', 'E',
                              'count_type', 'ENSG', 'mean_median', 'threshold',
                              'UP_or_DOWN', 'in_therapy', 'p_value', 'plot_type', 'fst_life_mean', 'scnd_life_mean']).to_csv(
                                  path, sep='\t', index=False)
        print(f'writing empty file {path}')
    for path in [deseq_lifeline_pdf_UP, deseq_lifeline_pdf_DOWN]:
        open(path, 'a').close()
        print(f'writing empty file {path}')
    os._exit(0)

# in case the S_counts is empty, write empty files:
if pd.read_table(counts).empty:
    write_empty_files()
# limit the counts to the ENSG applied:
S_counts = pd.read_table(counts, index_col=0).loc[ENSG,:]
# median = S_counts[S_counts > 0].median()
S_counts_temp = S_counts # to access the count values which are excluded after thresh apllication
S_counts_gt0 = S_counts[S_counts > 0]
# add the vital state infos to the case ids:
dead_alive_info = pd.read_table(meta_table, usecols=['bcr_patient_uuid', 'vital_status']).set_index('bcr_patient_uuid').loc[S_counts.index,:].reset_index()
S_counts_gt0 = S_counts_gt0.to_frame().reset_index().merge(dead_alive_info, how='left').set_index('vital_status')
# estimate the alive and dead medians:
try:
    alive_median = S_counts_gt0.loc['alive', ENSG].median()
except Exception as e:
    print(f'{e}\n to many zeros in counts, setting the median manually to 0')
    alive_median = 0
try:
    dead_median = S_counts_gt0.loc['dead', ENSG].median()
except Exception as e:
    print(f'{e}\n to many zeros in counts, setting the median manually to 0')
    dead_median = 0
mean_median = st.mean([alive_median, dead_median])

thresh = float(threshold.split('_')[1])
limit_val = (thresh/100) * mean_median
upper_limit = mean_median + limit_val
lower_limit = mean_median - limit_val
S_counts = S_counts.apply(lambda x: pd.NA if float(x) < upper_limit and float(x) > lower_limit else x)

cases_dropped = S_counts[S_counts.isna()].index
S_cases_dropped = S_counts_temp.loc[cases_dropped]
if not S_cases_dropped.empty:
    S_cases_dropped.to_csv(dropped_cases, sep='\t')
else:
    pd.DataFrame().to_csv(dropped_cases, sep='\t')

# now drop those cases, do not include them into the survival analyses:
S_counts = S_counts.dropna()
DF_meta = pd.read_table(meta_table, usecols=['bcr_patient_uuid', 'gender', 'project_id', 'survivaltime', 'years_to_last_follow_up', 'vital_status', 'pharmaceutical_therapy_drug_name']).set_index('bcr_patient_uuid').loc[S_counts.index,:]
DF_meta['T'] = DF_meta['survivaltime'].combine_first(DF_meta['years_to_last_follow_up'])
T = DF_meta['T']
# DF_meta['E'] = DF_meta['vital_status'] == 'dead'
DF_meta['E'] = np.where(DF_meta['vital_status'].values == 'dead', True, False)
# take the needed metadata out of the
# although there is metadata already available in the deseq_lifeline_tsv, it
# dous not contain all cases, therefore we need the meta_table again:

DF_to_plot = pd.merge(S_counts.to_frame(), DF_meta, left_index=True, right_index=True, how='left')
DF_to_plot.rename({ENSG: 'count'}, axis=1, inplace=True)
DF_to_plot['count_type'] = count_type
DF_to_plot['ENSG'] = ENSG
DF_to_plot.index.name = 'case_id'
DF_to_plot['mean_median'] = mean_median
DF_to_plot['threshold'] = thresh

DF_to_plot['UP_or_DOWN'] = DF_to_plot['count'].apply(lambda x: 'UP' if x > mean_median else 'DOWN')

def set_if_therapy(drug_invoked):
    temp_val = False
    for drug in DRUGS:
        if drug == drug_invoked:
            temp_val = True
    return temp_val

DF_to_plot['in_therapy'] = DF_to_plot['pharmaceutical_therapy_drug_name'].apply(set_if_therapy)

UP_DF = DF_to_plot[DF_to_plot['UP_or_DOWN'] == 'UP']
DOWN_DF = DF_to_plot[DF_to_plot['UP_or_DOWN'] == 'DOWN']
# both in_therapy cols must contain True and False vaues, otherwise no
# comparison in possible:
if len(UP_DF['in_therapy'].value_counts().index) != 2:
    write_empty_files()
    os._exit(0)

if len(DOWN_DF['in_therapy'].value_counts().index) != 2:
    write_empty_files()
    os._exit(0)

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


### plot up and down
###############################################################################
#                                 DOWN_PLOT                                   #
###############################################################################

T = DOWN_DF['T']
E = DOWN_DF['E']
therapy_bool = DOWN_DF['in_therapy']
fig, ax = plt.subplots(figsize=(8, 6))
kmf_THERAPY = KaplanMeierFitter()
kmf_THERAPY.fit(T[therapy_bool], E[therapy_bool], label='in therapy')
# DOWN_in_therapy_life_median = kmf_THERAPY._median
DOWN_in_therapy_life_mean = vars(kmf_THERAPY)['survival_function_'].iloc[:, 0].mean()
ax = kmf_THERAPY.plot_survival_function(ax=ax)

kmf_NO_THERAPY = KaplanMeierFitter()
kmf_NO_THERAPY.fit(T[~therapy_bool], E[~therapy_bool], label='not in therapy')
# DOWN_not_in_therapy_life_median = kmf_NO_THERAPY._median
DOWN_not_in_therapy_life_mean = vars(kmf_NO_THERAPY)['survival_function_'].iloc[:, 0].mean()
kmf_NO_THERAPY.plot_survival_function(ax=ax)

add_at_risk_counts(kmf_THERAPY, kmf_NO_THERAPY, ax=ax)
p_value_DOWN_str = f'p_value_DOWN = {Decimal(str(p_value_DOWN)):.2e}'

ax.set_title(f'{p_value_DOWN_str}, threshold = {round(thresh)}\nENSG: {ENSG}\n{project}, {drug_combi}\n{project}, {drug_combi}, {gender}, cutoff={cutoff}')

plt.tight_layout()
print(f'saving: {deseq_lifeline_pdf_DOWN}')
plt.savefig(deseq_lifeline_pdf_DOWN)
plt.close()

###############################################################################
#                                      UP_plot                                #
###############################################################################

T = UP_DF['T']
E = UP_DF['E']
therapy_bool = UP_DF['in_therapy']
fig, ax = plt.subplots(figsize=(8, 6))
kmf_THERAPY = KaplanMeierFitter()
kmf_THERAPY.fit(T[therapy_bool], E[therapy_bool], label='in therapy')
# UP_in_therapy_life_median=kmf_THERAPY._median
UP_in_therapy_life_mean = vars(kmf_THERAPY)['survival_function_'].iloc[:, 0].mean()
ax = kmf_THERAPY.plot_survival_function(ax=ax)

kmf_NO_THERAPY = KaplanMeierFitter()
kmf_NO_THERAPY.fit(T[~therapy_bool], E[~therapy_bool], label='not in therapy')
# UP_not_in_therapy_life_median=kmf_NO_THERAPY._median
UP_not_in_therapy_life_mean = vars(kmf_NO_THERAPY)['survival_function_'].iloc[:, 0].mean()
kmf_NO_THERAPY.plot_survival_function(ax=ax)

add_at_risk_counts(kmf_THERAPY, kmf_NO_THERAPY, ax=ax)
p_value_UP_str = f'p_value_UP = {Decimal(str(p_value_UP)):.2e}'
ax.set_title(f'{p_value_UP_str}, threshold = {round(thresh)}\nENSG: {ENSG}\n{project}, {drug_combi}, {gender}, cutoff={cutoff}')

plt.tight_layout()
print(f'saving: {deseq_lifeline_pdf_UP}')
plt.savefig(deseq_lifeline_pdf_UP)
plt.close()

###############################################################################
#                               writing tables                                #
###############################################################################
# 2 kmf means per table are added, the kmf estimation must be therefore ran
# already, so we save the tables after plotting the pdfs


UP_DF = UP_DF.assign(p_value=np.repeat(p_value_UP,len(UP_DF)))
UP_DF['plot_type'] = 'UP_validation'
UP_DF['fst_life_mean'] = UP_in_therapy_life_mean
UP_DF['scnd_life_mean'] = UP_not_in_therapy_life_mean
UP_DF.to_csv(deseq_lifeline_tsv_UP, sep='\t')
print(f'writing {deseq_lifeline_tsv_UP}')

DOWN_DF = DOWN_DF.assign(p_value=np.repeat(p_value_DOWN,len(DOWN_DF)))
DOWN_DF['plot_type'] = 'DOWN_validation'
DOWN_DF['fst_life_mean'] = DOWN_in_therapy_life_mean
DOWN_DF['scnd_life_mean'] = DOWN_not_in_therapy_life_mean
DOWN_DF.to_csv(deseq_lifeline_tsv_DOWN, sep='\t')
print(f'writing {deseq_lifeline_tsv_DOWN}')
