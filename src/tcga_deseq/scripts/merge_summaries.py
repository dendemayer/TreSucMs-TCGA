import pandas as pd
import os
import sys

sys.stderr = sys.stdout = open(snakemake.log[0], 'w')

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

summary_table = snakemake.input.summary_table
summary_table_info = snakemake.input.summary_table_info
summary_table_complement = snakemake.input.summary_table_complement
summary_table_complement_info = snakemake.input.summary_table_complement_info

summary_table_both = snakemake.output.summary_table_both
summary_table_both_info = snakemake.output.summary_table_both_info

# write empty files if summary_table is empty:
if pd.read_table(summary_table).empty:
    [pd.DataFrame().to_csv(i) for i in snakemake.output]
    os._exit(0)

DF_summary = pd.read_table(summary_table, sep='\t', index_col=0)
DF_summary_info = pd.read_table(summary_table_info, sep='\t', index_col=0)
DF_summary_complement = pd.read_table(summary_table_complement, sep='\t', index_col=0)
DF_summary_complement_info = pd.read_table(summary_table_complement_info, sep='\t', index_col=0)

DF_summary_both = pd.concat([DF_summary, DF_summary_complement], axis=1)
DF_summary_info_both = pd.concat([DF_summary_info, DF_summary_complement_info], axis=1)

DF_summary_both.to_csv(summary_table_both, sep='\t')
DF_summary_info_both.to_csv(summary_table_both_info, sep='\t')

