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
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

meta_table = snakemake.input.meta_table
metilene_intersect = snakemake.input.metilene_intersect
DMR = snakemake.wildcards.DMR
lifeline_out_pdf = snakemake.output.lifeline_out_pdf
lifeline_out_tsv = snakemake.output.lifeline_out_tsv
threshold = snakemake.wildcards.threshold

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]
#######################################################################

def apply_thresh(row):
    median = row.median()
    beta_val_max = row.max()
    limit_val = (thresh/100) * beta_val_max
    upper_limit = median + limit_val
    lower_limit = median - limit_val
    row_thresh = row.apply(lambda x: pd.NA if float(x) < upper_limit and float(x) > lower_limit else x)
    return row_thresh
