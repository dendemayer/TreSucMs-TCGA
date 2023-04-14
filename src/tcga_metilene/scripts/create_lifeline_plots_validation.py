import sys
import copy
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

"""
validate found DMRs
-> need the start position on whichs betavalues the lifeline plot is based on,
    to be found in start_tsv
-> all betavalues, (pre filtered gender and cutoff )
    to be found in summary_for_metilene, and the complement
thats it, consider also here the threshold when dividing in UP and DOWN in the
merged summary for metilene
"""
#####################
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

meta_table = snakemake.input.meta_table
start_tsv = snakemake.input.start_tsv
summary = snakemake.input.summary
summary_complement = snakemake.input.summary_complement

# # start_tsv is the table which is based on the lifeline_DMR_plot, and holds the
# # extract start value
DMR = snakemake.wildcards.DMR
UP_val_plot = snakemake.output.UP_val_plot
UP_val_tsv =snakemake.output.UP_val_tsv
DOWN_val_plot =snakemake.output.DOWN_val_plot
DOWN_val_tsv =snakemake.output.DOWN_val_tsv
threshold = snakemake.wildcards.threshold
drug_combi = snakemake.wildcards.drug_combi

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

#####################

# # snakemake inputs:
# meta_table = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv"
# start_tsv = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr5_66828566_66828775.tsv"
# summary = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC_TCGA-HNSC/metilene/metilene_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/summary_for_metilene.tsv"
# summary_complement = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC_TCGA-HNSC/metilene/metilene_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/summary_for_metilene_complement.tsv"
# # snakemake output:
# UP_val_plot = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr5_66828566_66828775_UP_val.pdf"
# UP_val_tsv = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr5_66828566_66828775_UP_val.tsv"
# DOWN_val_plot = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr5_66828566_66828775_DOWN_val.pdf"
# DOWN_val_tsv = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr5_66828566_66828775_DOWN_val.tsv"
# # snakemake wildcards:
# output_path = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3"
# project = "TCGA-CESC_TCGA-HNSC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_0"
# threshold = "threshold_0"
# DMR = "chr5_66828566_66828775"


def write_empty_files():
    for path in [UP_val_tsv, DOWN_val_tsv]:
        pd.DataFrame(columns=['case_id', 'group', 'vital_status', 'drugs',
                          'gender', 'projects', 'beta_values', 'T', 'E',
                          'chr_start']).to_csv(path, sep='\t', index=False)
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
DRUGS = drug_combi.split('_')
chr_ = DMR.split('_')[0]
# read in all beta values at the desired start position:
summary_DF = pd.read_table(summary, na_values='.').set_index(['Chromosome', 'Start']).sort_index().loc[(chr_,start),:].dropna()
summary_complement_DF = pd.read_table(summary_complement, na_values='.').set_index(['Chromosome', 'Start']).sort_index().loc[(chr_,start),:].dropna()
DF_summary = pd.concat([summary_DF, summary_complement_DF], axis=1)

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
    median_temp = row.median()
    beta_val_max = row.max()
    limit_val = (thresh/100) * beta_val_max
    upper_limit = median_temp + limit_val
    lower_limit = median_temp - limit_val
    row_thresh = row.apply(lambda x: pd.NA if float(x) < upper_limit and float(x) > lower_limit else x)
    return row_thresh

DF_summary = DF_summary.apply(apply_thresh, axis=1)
DF_summary.dropna(axis=1, inplace=True)

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

DF_temp = DF_summary.T
median = DF_summary.median(axis=1).values[0]
temp_DF = DF_summary.apply(lambda x: 'UP' if float(x) > median else 'DOWN').to_frame()
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

DOWN_DF = final_DF.set_index((chr_, start)).loc['DOWN', :].reset_index().set_index('case_id')
# for later summaries, the col names must be identical, add the chr_start as
# extra col:
DOWN_DF.rename({(chr_, start): 'group'}, axis=1, inplace=True)
UP_DF = final_DF.set_index((chr_, start)).loc['UP', :].reset_index().set_index('case_id')
UP_DF.rename({(chr_, start): 'group'}, axis=1, inplace=True)

def set_if_therapy(drug_invoked):
    temp_val = False
    for drug in DRUGS:
        if drug == drug_invoked:
            temp_val = True
    return temp_val

DOWN_DF['in_therapy'] = DOWN_DF['drugs'].apply(set_if_therapy)
UP_DF['in_therapy'] = UP_DF['drugs'].apply(set_if_therapy)
# both in_therapy cols must contain True and False valse, otherwise no
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


UP_DF['chr_start'] = chr_ + str(start)
UP_DF['p_value'] = p_value_UP
UP_DF.to_csv(UP_val_tsv, sep='\t')
print(f'writing {UP_val_tsv}')

DOWN_DF['chr_start'] = chr_ + str(start)
DOWN_DF['p_value'] = p_value_DOWN
DOWN_DF.to_csv(DOWN_val_tsv, sep='\t')
print(f'writing {DOWN_val_tsv}')

### plot up and down

T = DOWN_DF['T']
E = DOWN_DF['E']
therapy_bool = DOWN_DF['in_therapy']
fig, ax = plt.subplots(figsize=(8, 6))
kmf_THERAPY = KaplanMeierFitter()
kmf_THERAPY.fit(T[therapy_bool], E[therapy_bool], label='in therapy')
ax = kmf_THERAPY.plot_survival_function(ax=ax)

kmf_NO_THERAPY = KaplanMeierFitter()
kmf_NO_THERAPY.fit(T[~therapy_bool], E[~therapy_bool], label='not in therapy')
kmf_NO_THERAPY.plot_survival_function(ax=ax)

add_at_risk_counts(kmf_THERAPY, kmf_NO_THERAPY, ax=ax)
p_value_DOWN_str = f'p_value_DOWN = {Decimal(str(p_value_DOWN)):.2e}'

ax.set_title(f'{p_value_DOWN_str}, threshold = {round(thresh)}\nDMR: {": ".join(DMR.split("_")[:2]) + "-" + DMR.split("_")[2]}\nStart: {start}')

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
ax = kmf_THERAPY.plot_survival_function(ax=ax)

kmf_NO_THERAPY = KaplanMeierFitter()
kmf_NO_THERAPY.fit(T[~therapy_bool], E[~therapy_bool], label='not in therapy')
kmf_NO_THERAPY.plot_survival_function(ax=ax)

add_at_risk_counts(kmf_THERAPY, kmf_NO_THERAPY, ax=ax)
p_value_UP_str = f'p_value_UP = {Decimal(str(p_value_UP)):.2e}'
ax.set_title(f'{p_value_UP_str}, threshold = {round(thresh)}\nDMR: {": ".join(DMR.split("_")[:2]) + "-" + DMR.split("_")[2]}\nStart: {start}')

plt.tight_layout()
print(f'saving: {UP_val_plot}')
plt.savefig(UP_val_plot)
plt.close()
