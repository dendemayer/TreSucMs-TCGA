import pandas as pd
# from methyl import set_logger
# from lifelines.datasets import load_waltons
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from decimal import Decimal
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts

"""
follow lifeline succession from patients with detected DMR
- needed, metadata:f.e.:
    - TCGA-pipelines_2/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/meta_info_druglist_merged_drugs_combined.tsv
    - cols:
        # survivaltime       | years_to_last_follow_up 'bcr_patient_uuid',
- the metilene intersect out table
    - TCGA-pipelines_2/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_complement_intersect.tsv
    # needed infos are vital_status (alive dead),

# use the logrank_test at every position within a DMR to select the best
# resolution between up and down :
If the p-value of the log-rank test is less than a chosen significance level
(such as 0.05), it suggests that there is a statistically significant
difference in survival between the groups, and the Kaplan-Meier estimate can be
considered meaningful.
"""

meta_table = snakemake.input.meta_table
met_table = snakemake.input.metilene_intersect
DMR = snakemake.wildcards.range
lifeline_out_pdf = snakemake.output.lifeline_out_pdf
lifeline_out_tsv = snakemake.output.lifeline_out_tsv
# '/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_2/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_complement_intersect_lifeline_chr1_78045659_78046029.pdf'
threshold = snakemake.wildcards.threshold

DF_meta = pd.read_table(meta_table, usecols=['bcr_patient_uuid', 'survivaltime', 'years_to_last_follow_up', 'vital_status'])
DF_metilene = pd.read_table(met_table, header=[0,1,2,3,4], index_col=[0,1,2,3], na_values='.').dropna()
# limit the DF_metilene to the recent range:
DF_metilene = DF_metilene.loc[(slice(None), slice(None), slice(None), DMR), :]



#############################################
# get the used cases out of the meta_table:
used_cases = [ i[1] for i in DF_metilene.columns]
DF_meta = DF_meta.set_index('bcr_patient_uuid').loc[used_cases, :]
DF_meta['T'] = DF_meta['survivaltime'].combine_first(DF_meta['years_to_last_follow_up'])

T = DF_meta['T']
DF_meta['E'] = DF_meta['vital_status'] == 'dead'
# in case not value available for survivaltime AND years_to_last_follow_up, a
# NaN is set in T, :
# '/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_2/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin_paclitaxel/female_male/cutoff_0/metilene_complement_intersect.tsv'
# on case_id:
                                     # vital_status  survivaltime  years_to_last_follow_up   T      E
# bcr_patient_uuid
# 31d8591b-9d80-406c-b61a-7a6544c73466        alive           NaN                      NaN NaN  False
# 7625a27d-a012-4f4d-b637-5f85d7cdae8c         dead           NaN                      NaN NaN   True

# shrink the complete DF_meta table, then extract T and E:
index_to_delete = DF_meta[T.isna()].index
# since we use also DF_metilene we have to reduce it to the cases for which the
# T and E infos are available:
DF_metilene.drop(labels=index_to_delete, axis=1, level=1, inplace=True)
DF_meta = DF_meta[T.notna()]
T=DF_meta['T']
E = DF_meta['vital_status'] == 'dead' #  E is a either boolean or binary array
# representing whether the “death” was observed or not (alternatively an
# individual can be censored)
#################################################################

# the groups are representing the cases which are either higher or lower
# methylated than the median of all cases within the DMR, therefore calculate
# the kmf at every postition and take the one with the highest
starts = [i[1] for i in DF_metilene.index]
p_value = 0
pvalue_list = []
for start in starts:
    median = DF_metilene.loc[(slice(None), start, slice(None), slice(None)), :].median(axis=1).values[0]
    try:
        DF_temp = DF_metilene.loc[(slice(None), start, slice(None), slice(None)), :].apply(lambda x: 'UP' if float(x) > median else 'DOWN' )
    except Exception as e:
        breakpoint()
        print(e)

    UP_bool = (DF_temp == 'UP').values
    DOWN_bool = (DF_temp == 'DOWN').values
    results = logrank_test(T[UP_bool], T[DOWN_bool], E[UP_bool], E[DOWN_bool])
    p_value = results.p_value
    pvalue_list.append(p_value)

# take the lowest p_value,

p_value_str = f'p_value = {Decimal(str(min(pvalue_list))):.2e}'
plot_index = pvalue_list.index(min(pvalue_list))
median = DF_metilene.iloc[plot_index, :].median()
# if threshold is not 0, then apply it here, do not use values that are withing
# threshold % over or under the median
thresh = float(threshold.split('_')[1])
try:
    beta_val_max = DF_metilene.iloc[plot_index, :].max()
except Exception as e:
    breakpoint()
    print(e)

limit_val = (thresh/100) * beta_val_max
upper_limit = median + limit_val
lower_limit = median - limit_val

if thresh != 0:
    # DF_metilene.iloc[plot_index, :].apply(lambda: pd.NA if float(x) > median elif
    DF_plot = DF_metilene.iloc[plot_index, :].apply(lambda x: pd.NA if float(x) < upper_limit and float(x) > lower_limit else x)
    DF_plot.dropna(inplace=True)
    # also T and E must be set newly since the length might be shrunken
    try:
        case_ids = [i[1] for i in DF_plot.index]
    except Exception as e:
        breakpoint()
        print(e)
    DF_meta_temp = DF_meta.loc[case_ids, :]
    T = DF_meta_temp['T']
    E = DF_meta_temp['vital_status'] == 'dead' #  E is a either boolean or binary array
    # generate the new p_value:
    DF_plot = DF_plot.apply(lambda x: 'UP' if float(x) > median else 'DOWN' )
    UP_bool = (DF_plot == 'UP').values
    DOWN_bool = (DF_plot == 'DOWN').values
    try:
        results = logrank_test(T[UP_bool], T[DOWN_bool], E[UP_bool], E[DOWN_bool])
    except Exception as e:
        breakpoint()
        print(e)

    p_value = results.p_value
    p_value_str = f'p_value = {Decimal(str(p_value)):.2e}'
else:
    DF_plot = DF_metilene.iloc[plot_index, :].apply(lambda x: 'UP' if float(x) > median else 'DOWN' )

DF_plot.name = 'UP_DOWN'

# kmf.fit(T, event_observed=E)  # or, more succinctly, kmf.fit(T, E)

# groups = DF_plot['UP_DOWN']
# in case we cannot create 2 groups (i.e., with the thresh invokening, just UP
# cases are present:/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_2/TCGA-HNSC/metilene/metilene_output/carboplatin,paclitaxel_cisplatin_paclitaxel/female_male/cutoff_0/metilene_intersect.tsv
# leads to just up cases with thresh =10: DF_metilene.iloc[plot_index, :].median()
# 0.0204743706513952
# (Pdb) limit_val
# 0.0445787454608645
# just write 2 emtpy files:
if DF_plot.value_counts().nunique() == 1:
    open(lifeline_out_pdf, 'a').close()
    pd.DataFrame(columns=['timeline', 'UP', 'DOWN', 'p_value']).to_csv(lifeline_out_tsv, sep='\t')
else:
    UP_bool = (DF_plot == 'UP').values
    DOWN_bool = (DF_plot == 'DOWN').values
    fig, ax = plt.subplots(figsize=(8, 6))

    try:
        kmf_UP = KaplanMeierFitter()
        kmf_UP.fit(T[UP_bool], E[UP_bool], label='UP')
        ax = kmf_UP.plot_survival_function(ax=ax)

        kmf_DOWN = KaplanMeierFitter()
        kmf_DOWN.fit(T[DOWN_bool], E[DOWN_bool], label='DOWN')
        kmf_DOWN.plot_survival_function(ax=ax)
    except Exception as e:
        breakpoint()
        print(e)

    add_at_risk_counts(kmf_UP, kmf_DOWN, ax=ax)
    ax.set_title(f'{p_value_str} \nDMR: {": ".join(DMR.split("_")[:2]) + "-" + DMR.split("_")[2]}\nStart: {str(starts[plot_index])}')

    plt.tight_layout()
    print(f'saving: {lifeline_out_pdf}')
    plt.savefig(lifeline_out_pdf)
    plt.close()

    DF = pd.concat([kmf_UP.survival_function_, kmf_DOWN.survival_function_], axis=1)
    DF['p_value'] = p_value
    DF.to_csv(lifeline_out_tsv , sep='\t')
