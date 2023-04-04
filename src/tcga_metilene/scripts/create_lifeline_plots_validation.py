import pandas as pd
# from methyl import set_logger
# from lifelines.datasets import load_waltons
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from decimal import Decimal
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts


met_table = snakemake.input.metilene_intersect
meta_table = snakemake.input.meta_table
DMR = snakemake.wildcards.range
UP_val_plot = snakemake.output.UP_val_plot
UP_val_tsv =snakemake.output.UP_val_tsv
DOWN_val_plot =snakemake.output.DOWN_val_plot
DOWN_val_tsv =snakemake.output.DOWN_val_tsv
threshold = snakemake.wildcards.threshold
drug_combi = snakemake.wildcards.drug_combi
DRUGS = drug_combi.split('_')

print(f'meta_table = "{meta_table}"')
print(f'met_table = "{met_table}"')
print(f'DMR = "{DMR}"')
print(f'threshold = "{threshold}"')

DF_meta = pd.read_table(meta_table, usecols=['bcr_patient_uuid', 'survivaltime', 'years_to_last_follow_up', 'vital_status', 'pharmaceutical_therapy_drug_name'])
DF_metilene = pd.read_table(met_table, header=[0,1,2,3,4], index_col=[0,1,2,3], na_values='.').dropna()
# limit the DF_metilene to the recent range:
DF_metilene = DF_metilene.loc[(slice(None), slice(None), slice(None), DMR), :]

thresh = float(threshold.split('_')[1])
# set every beta value which is to close to the median to pd.NA:
def apply_thresh(row):
    median = row.median()
    beta_val_max = row.max()
    limit_val = (thresh/100) * beta_val_max
    upper_limit = median + limit_val
    lower_limit = median - limit_val
    row_thresh = row.apply(lambda x: pd.NA if float(x) < upper_limit and float(x) > lower_limit else x)
    return row_thresh

DF_metilene = DF_metilene.apply(apply_thresh, axis=1)

# prepare the DF_meta table, do not limit on case_id but make sure we just
# include cases in the DF_meta and DF_metilene for which we have the T value
# available:
DF_meta['T'] = DF_meta['survivaltime'].combine_first(DF_meta['years_to_last_follow_up'])
T = DF_meta['T']
index_to_delete = DF_meta[T.isna()].index
DF_metilene.drop(labels=index_to_delete, axis=1, level=1, inplace=True)
DF_meta = DF_meta[T.notna()]
breakpoint()

