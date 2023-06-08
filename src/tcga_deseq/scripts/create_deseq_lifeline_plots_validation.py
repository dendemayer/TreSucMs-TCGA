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

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

def write_empty_files():
    for path in [deseq_lifeline_tsv_UP, deseq_lifeline_tsv_DOWN]:
        pd.DataFrame(columns=['count', 'project_id',
                              'pharmaceutical_therapy_drug_name', 'gender',
                              'vital_status', 'survivaltime',
                              'years_to_last_follow_up', 'T', 'E',
                              'count_type', 'ENSG', 'median', 'threshold',
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
median = S_counts[S_counts > 0].median()
S_counts_temp = S_counts # to access the count values which are excluded after thresh apllication
thresh = float(threshold.split('_')[1])
limit_val = (thresh/100) * median
upper_limit = median + limit_val
lower_limit = median - limit_val
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
DF_to_plot['median'] = median
DF_to_plot['threshold'] = thresh

DF_to_plot['UP_or_DOWN'] = DF_to_plot['count'].apply(lambda x: 'UP' if x > median else 'DOWN')

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

ax.set_title(f'{p_value_DOWN_str}, threshold = {round(thresh)}')

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
ax.set_title(f'{p_value_UP_str}, threshold = {round(thresh)}')

plt.tight_layout()
print(f'saving: {deseq_lifeline_pdf_UP}')
plt.savefig(deseq_lifeline_pdf_UP)
plt.close()

###############################################################################
#                               writing tables                                #
###############################################################################
# 2 kmf medians per table are added, the kmf estimation must be therefore ran
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
