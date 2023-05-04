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
# from pandas.core import apply

#######################################################################
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

meta_table = snakemake.input.meta_table
summary_table_info = snakemake.input.summary_table_info
deseq_counts = snakemake.input.deseq_counts
threshold = snakemake.wildcards.threshold
deseq_lifeline_pdf = snakemake.output.deseq_lifeline_pdf
deseq_lifeline_tsv = snakemake.output.deseq_lifeline_tsv
ENSG = snakemake.wildcards.ENSG
count_type = snakemake.wildcards.count_type + '_count'

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]
#######################################################################

# # # snakemake inputs:
# meta_table = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC/DESeq2/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv"
# summary_table_info = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC/DESeq2/DESeq2_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/summary_for_DESeq2_INFO.tsv"
# deseq_counts = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/DESeq2_heatmap_log2fDECREASE_norm_counts.tsv"
# script_file = "../tcga_deseq/scripts/create_deseq_lifeline_plots.py"
# # snakemake output:
# deseq_lifeline_pdf = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000249235_lifeline.pdf"
# deseq_lifeline_tsv = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/DESeq2_log2f_DECREASE_norm_ENSG00000249235_lifeline.tsv"
# # snakemake wildcards:
# output_path = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3"
# project = "TCGA-CESC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_0"
# threshold = "threshold_0"
# in_de = "DE"
# count_type = "norm"
# ENSG = "ENSG00000249235"

#######################################################################
def write_empty_files():
    print(f'writing empty file {deseq_lifeline_pdf}:')
    open(deseq_lifeline_pdf,'a').close()
    print(f'writing empty file {deseq_lifeline_tsv}:')
    pd.DataFrame(columns=['case_id', count_type, 'UP_or_DOWN', 'project_id',
                          'pharmaceutical_therapy_drug_name', 'gender',
                          'vital_status', 'survivaltime',
                          'years_to_last_follow_up', 'T', 'E', 'UP_bool',
                          'DOWN_bool', 'included_in_KM', 'p_value', 'median',
                          'threshold', 'ENSG']).to_csv(deseq_lifeline_tsv, sep='\t',
                                               index=False)
DF_counts = pd.read_table(deseq_counts).loc[ENSG,:]

thresh = float(threshold.split('_')[1])

median = DF_counts[DF_counts > 0].median()
limit_val = (thresh/100) * median
upper_limit = median + limit_val
lower_limit = median - limit_val
DF_counts = DF_counts.apply(lambda x: pd.NA if float(x) < upper_limit and float(x) > lower_limit else x)
# save the cases and their count values which are getting dropped due to the
# threshold parameter:
cases_dropped = DF_counts[DF_counts.isna()].index
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
DF_to_plot['UP_or_DOWN'] = np.where(DF_to_plot[ENSG].values > median, 'UP', 'DOWN')
DF_to_plot = pd.merge(DF_to_plot, DF_meta, left_index=True, right_index=True, how='left')
DF_dropped_cases = pd.merge(DF_dropped_cases, DF_meta, left_index=True, right_index=True, how='left')
T = DF_to_plot['T']
E = DF_to_plot['E']
UP_bool = (DF_to_plot['UP_or_DOWN'] == 'UP').values
DF_to_plot['UP_bool'] = UP_bool
DOWN_bool = (DF_to_plot['UP_or_DOWN'] == 'DOWN').values
DF_to_plot['DOWN_bool'] = DOWN_bool

fig, ax = plt.subplots(figsize=(8, 6))

kmf_UP = KaplanMeierFitter()
kmf_UP.fit(T[UP_bool], E[UP_bool], label='UP')
ax = kmf_UP.plot_survival_function(ax=ax)

kmf_DOWN = KaplanMeierFitter()
kmf_DOWN.fit(T[DOWN_bool], E[DOWN_bool], label='DOWN')
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

ax.set_title(f'{p_value_str}, threshold = {threshold.split("_")[1]}\nENSG: {ENSG}')

plt.tight_layout()
print(f'saving: {deseq_lifeline_pdf}')
plt.savefig(deseq_lifeline_pdf)
plt.close()

# prepare to store the table:
print(f'saving {deseq_lifeline_tsv}')
DF_to_plot['included_in_KM'] = True
DF_dropped_cases['included_in_KM'] = False
DF_to_plot = pd.concat([DF_to_plot, DF_dropped_cases])
DF_to_plot.rename({ENSG: count_type}, axis=1, inplace=True)
DF_to_plot.index.name = 'case_id'
DF_to_plot['p_value'] = p_value
DF_to_plot['median'] = median
DF_to_plot['threshold'] = thresh
DF_to_plot['ENSG'] = ENSG
DF_to_plot.to_csv(deseq_lifeline_tsv, sep='\t')