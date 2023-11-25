import pandas as pd
import numpy as np
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
import statistics as st

"""
# including here the deseq result infos, from:
vital_status	gender	project_id	pharmaceutical_therapy_drug_name	bcr_patient_uuid	count	UP_or_DOWN	project_id	pharmaceutical_therapy_drug_name	gender	vital_status	survivaltime	years_to_last_follow_up	T	E	included_in_KM	count_type	p_value	mean_median	threshold	ENSG	plot_type	fst_life_mean	scnd_life_mean
to: (therefore the according deseeq result table is needed as input)
vital_status	gender	project_id	pharmaceutical_therapy_drug_name	bcr_patient_uuid	count	UP_or_DOWN	project_id	pharmaceutical_therapy_drug_name	gender	vital_status	survivaltime	years_to_last_follow_up	T	E	included_in_KM	count_type	p_value	mean_median	threshold	ENSG	plot_type	fst_life_mean	scnd_life_mean
log2FoldChange lfcSE stat pvalue padj

"""

#######################################################################
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

meta_table = snakemake.input.meta_table
summary_table_info = snakemake.input.summary_table_info
deseq_counts = snakemake.input.deseq_counts
deseq_result = snakemake.input.deseq_result
threshold = snakemake.wildcards.threshold
deseq_lifeline_pdf = snakemake.output.deseq_lifeline_pdf
deseq_lifeline_tsv = snakemake.output.deseq_lifeline_tsv
ENSG = snakemake.wildcards.ENSG
count_type = snakemake.wildcards.count_type + '_count'
cutoff = snakemake.wildcards.cutoff.split('_')[1]
drug_combi = snakemake.wildcards.drug_combi.replace('_', ';')
project = ', '.join(snakemake.wildcards.project.split('_'))
gender = ', '.join(snakemake.wildcards.gender.split('_'))

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

#######################################################################

# # #######################################################################

# KeyError: 'alive':
# # snakemake inputs:
# meta_table = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv"
# summary_table_info = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/summary_for_DESeq2_INFO.tsv"
# deseq_counts = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/DESeq2_heatmap_log2fDECREASE_norm_counts.tsv"
# deseq_result = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/DESeq2_log2fDECREASE_result.tsv"
# script_file = "../tcga_deseq/scripts/create_deseq_lifeline_plots.py"
# # snakemake output:
# deseq_lifeline_pdf = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000249235_lifeline.pdf"
# deseq_lifeline_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000249235_lifeline.tsv.gz"
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

# # snakemake inputs:
# meta_table = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/merged_meta_files/cutoff_8/meta_info_druglist_merged_drugs_combined.tsv"
# summary_table_info = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/summary_for_DESeq2_INFO.tsv"
# deseq_counts = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/DESeq2_heatmap_log2fDECREASE_raw_counts.tsv"
# deseq_result = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/DESeq2_log2fDECREASE_result.tsv"
# script_file = "../tcga_deseq/scripts/create_deseq_lifeline_plots.py"
# # snakemake output:
# deseq_lifeline_pdf = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_20/DESeq2_log2f_DECREASE_raw_ENSG00000215533_lifeline.pdf"
# deseq_lifeline_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_20/DESeq2_log2f_DECREASE_raw_ENSG00000215533_lifeline.tsv.gz"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
# project = "TCGA-CESC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_8"
# threshold = "threshold_20"
# in_de = "DE"
# count_type = "raw"
# ENSG = "ENSG00000215533"

# # snakemake inputs:
# meta_table = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv"
# summary_table_info = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/summary_for_DESeq2_INFO.tsv"
# deseq_counts = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/DESeq2_heatmap_log2fDECREASE_norm_counts.tsv"
# script_file = "../tcga_deseq/scripts/create_deseq_lifeline_plots.py"
# deseq_result = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/DESeq2_log2fDECREASE_result.tsv"
# # snakemake output:
# deseq_lifeline_pdf = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000081041_lifeline.pdf"
# deseq_lifeline_tsv = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000081041_lifeline.tsv.gz"
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

# # snakemake inputs:
# meta_table = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_4/TCGA-CESC/DESeq2/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv"
# summary_table_info = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_4/TCGA-CESC/DESeq2/DESeq2_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/summary_for_DESeq2_INFO.tsv"
# deseq_counts = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_4/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/DESeq2_heatmap_log2fDECREASE_norm_counts.tsv"
# script_file = "../tcga_deseq/scripts/create_deseq_lifeline_plots.py"
# # snakemake output:
# deseq_lifeline_pdf = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_4/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000081041_lifeline.pdf"
# deseq_lifeline_tsv = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_4/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000081041_lifeline.tsv.gz"
# # snakemake wildcards:
# output_path = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_4"
# project = "TCGA-CESC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_0"
# threshold = "threshold_0"
# in_de = "DE"
# count_type = "norm"
# ENSG = "ENSG00000081041"
# # #######################################################################

def write_empty_files():
    print(f'writing empty file {deseq_lifeline_pdf}:')
    open(deseq_lifeline_pdf,'a').close()
    print(f'writing empty file {deseq_lifeline_tsv}:')
    pd.DataFrame(columns=['count', 'UP_or_DOWN', 'project_id',
                          'pharmaceutical_therapy_drug_name', 'gender',
                          'vital_status', 'survivaltime',
                          'years_to_last_follow_up', 'T', 'E',
                          'included_in_KM', 'count_type', 'p_value', 'mean_median',
                          'threshold', 'ENSG', 'plot_type', 'fst_life_mean',
                          'scnd_life_mean']).to_csv(deseq_lifeline_tsv,
                                                      sep='\t',
                                               index=False)
DF_counts = pd.read_table(deseq_counts).loc[ENSG,:]
# to be able to filter on 'alive' and 'dead' on median estimation, include
# those infos as MI in the count frame:
DF_counts.index = pd.MultiIndex.from_frame(pd.read_table(summary_table_info, index_col=0).T)

thresh = float(threshold.split('_')[1])

# TODO change here accordingly to alive_median and dead_median
# median = DF_counts[DF_counts > 0].median()
DF_counts_gt0 = DF_counts[DF_counts > 0]
try:
    alive_median = DF_counts_gt0.loc[('alive', slice(None), slice(None), slice(None), slice(None))].median()
except Exception as e:
    print(f'{e}\n to many zeros in counts, setting the median manually to 0')
    alive_median = 0
try:
    dead_median = DF_counts_gt0.loc[('dead', slice(None), slice(None), slice(None), slice(None))].median()
except Exception as e:
    print(f'{e}\n to many zeros in counts, setting the median manually to 0')
    dead_median = 0

mean_median = st.mean([alive_median, dead_median])
limit_val = (thresh/100) * mean_median
upper_limit = mean_median + limit_val
lower_limit = mean_median - limit_val
# bring DF_counts back to the previous form with index is just the ENSG:
DF_counts = DF_counts.reset_index(['vital_status', 'gender', 'project_id', 'pharmaceutical_therapy_drug_name'], drop=True)
DF_counts = DF_counts.apply(lambda x: pd.NA if float(x) < upper_limit and float(x) > lower_limit else x)
# save the cases and their count values which are getting dropped due to the
# threshold parameter:
cases_dropped = DF_counts[DF_counts.isna()].index
# cases_dropped = DF_counts[DF_counts.isna()].reset_index()['bcr_patient_uuid'].to_list()
DF_dropped_cases = pd.read_table(deseq_counts).loc[ENSG, cases_dropped]
# now drop those cases, do not include them into the survival analyses:
DF_counts = DF_counts.dropna()

invoked_cases = pd.read_table(summary_table_info).loc[4][1:].values
# just take case_id included in the deseq run, limit the DF_meta to those cases
# via indexing on invoked_cases:
DF_meta = pd.read_table(meta_table, usecols=['bcr_patient_uuid', 'gender', 'project_id', 'survivaltime', 'years_to_last_follow_up', 'vital_status', 'pharmaceutical_therapy_drug_name']).set_index('bcr_patient_uuid').loc[invoked_cases,:]

DF_meta['T'] = DF_meta['survivaltime'].combine_first(DF_meta['years_to_last_follow_up'])
T = DF_meta['T']
# DF_meta['E'] = DF_meta['vital_status'] == 'dead'
DF_meta['E'] = np.where(DF_meta['vital_status'].values == 'dead', True, False)

# plot the specific ENSG:
DF_to_plot = DF_counts.to_frame()
DF_to_plot['UP_or_DOWN'] = np.where(DF_to_plot[ENSG].values > mean_median, 'UP', 'DOWN')
DF_to_plot = pd.merge(DF_to_plot, DF_meta, left_index=True, right_index=True, how='left')
DF_dropped_cases = pd.merge(DF_dropped_cases, DF_meta, left_index=True, right_index=True, how='left')
T = DF_to_plot['T']
E = DF_to_plot['E']
UP_bool = (DF_to_plot['UP_or_DOWN'] == 'UP').values
DF_to_plot['UP_bool'] = UP_bool
DOWN_bool = (DF_to_plot['UP_or_DOWN'] == 'DOWN').values
DF_to_plot['DOWN_bool'] = DOWN_bool

fig, ax = plt.subplots(figsize=(8, 6))

# before KaplanMeierFitter check whether we have actually dead and alive cases
# to compare against each other:
if len(DF_to_plot['vital_status'].value_counts()) != 2:
    write_empty_files()
    os._exit(0)
# its also mandatory that UP_or_DOWN are binary
if len(DF_to_plot['UP_or_DOWN'].value_counts()) != 2:
    write_empty_files()
    os._exit(0)
kmf_UP = KaplanMeierFitter()
kmf_UP.fit(T[UP_bool], E[UP_bool], label='UP')
# UP_life_expectancy_median = kmf_UP._median
UP_life_expectancy_mean = vars(kmf_UP)['survival_function_'].iloc[:, 0].mean()
ax = kmf_UP.plot_survival_function(ax=ax)

kmf_DOWN = KaplanMeierFitter()
kmf_DOWN.fit(T[DOWN_bool], E[DOWN_bool], label='DOWN')
# DOWN_life_expectancy_median = kmf_DOWN._median
DOWN_life_expectancy_mean = vars(kmf_DOWN)['survival_function_'].iloc[:, 0].mean()

kmf_DOWN.plot_survival_function(ax=ax)

add_at_risk_counts(kmf_UP, kmf_DOWN, ax=ax)

# A very low variance means that the column UP_bool completely determines
# whether a subject dies or not. See
# https://stats.stackexchange.com/questions/11109/how-to-deal-with-perfect-separation-in-logistic-regression.
# ConvergenceWarning: Column UP_bool have very low variance when conditioned on
# death event present or not. This may harm convergence. This could be a form
# of 'complete separation'. For example, try the following code:

# >>> events = DF_to_plot['E'].astype(bool)
# >>> print(DF_to_plot.loc[events, 'UP_bool'].var())
# >>> print(DF_to_plot.loc[~events, 'UP_bool'].var())

# (Pdb) pd.read_table(deseq_counts).loc[ENSG,:].value_counts()
# 0.000000       108
# 3198.009867      1
# 0.722794         1
# 0.945835         1
try:
    cph = CoxPHFitter().fit(DF_to_plot.loc[:, ['UP_bool', 'T', 'E']], 'T', 'E')
except Exception as e:
    print('{e}')
    print(f'writing empty files')
    write_empty_files()
    os._exit(0)

p_value = cph.summary['p'].values[0]
p_value_str = f'p_value = {Decimal(str(p_value)):.2e}'

ax.set_title(f'base plot, {p_value_str}, threshold = {threshold.split("_")[1]}\nENSG: {ENSG}\n{project}, {drug_combi}\n{gender}, cutoff={cutoff}')

plt.tight_layout()
print(f'saving: {deseq_lifeline_pdf}')
plt.savefig(deseq_lifeline_pdf)
plt.close()

# add the info of the lifeexpany of each group:
# prepare to store the table:
print(f'saving {deseq_lifeline_tsv}')
DF_to_plot['included_in_KM'] = True
DF_dropped_cases['included_in_KM'] = False
DF_to_plot = pd.concat([DF_to_plot, DF_dropped_cases])
DF_to_plot.rename({ENSG: 'count'}, axis=1, inplace=True)
DF_to_plot['count_type'] = count_type
DF_to_plot.index.name = 'case_id'
DF_to_plot['p_value'] = p_value
DF_to_plot['mean_median'] = mean_median
DF_to_plot['threshold'] = thresh
DF_to_plot['ENSG'] = ENSG
DF_to_plot['plot_type'] = 'base_plot'
DF_to_plot['fst_life_mean'] = UP_life_expectancy_mean
DF_to_plot['scnd_life_mean'] = DOWN_life_expectancy_mean
result_DF = pd.read_table(deseq_result).loc[ENSG,: ]
# we ommit the baseMean, since it refers to the raw count, but here also
# normalized and transformed counts are used -> could be misleading
DF_to_plot['log2FoldChange'] = result_DF['log2FoldChange']
DF_to_plot['lfcSE'] = result_DF['lfcSE']
DF_to_plot['stat'] = result_DF['stat']
DF_to_plot['pvalue'] = result_DF['pvalue']
DF_to_plot['padj'] = result_DF['padj']
# DF_to_plot.to_csv(deseq_lifeline_tsv, sep='\t')
DF_to_plot.drop(['UP_bool', 'DOWN_bool'], axis=1).to_csv(deseq_lifeline_tsv, sep='\t')
