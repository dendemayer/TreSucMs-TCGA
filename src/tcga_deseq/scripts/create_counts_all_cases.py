import pandas as pd
import os
import sys
import subprocess

sys.stderr = sys.stdout = open(snakemake.log[0], 'w')

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

summary_table_both = snakemake.input.summary_table_both
summary_table_both_info = snakemake.input.summary_table_both_info
R_path = snakemake.input.R_path
project = snakemake.wildcards.project
deseq_counts = snakemake.output.deseq_counts
deseq_output = os.path.split(deseq_counts[0])[0]


# # if the input table is empty, write empty tables for ever output table which
# # shall be created
if pd.read_table(summary_table_both, sep='\t').empty:
    [open(out, 'w').close() for out in snakemake.output]
    os._exit(0)


# # also if just alive or just dead cases are included, no diff exp analysis is
# # possible:
DF_temp = pd.read_table(summary_table_both_info, sep='\t').T
if DF_temp.iloc[1:,0].nunique() == 1:
    for out in snakemake.output:
        open(out, 'a').close()
    os._exit(0)

sequence = ['Rscript', R_path, summary_table_both, summary_table_both_info, deseq_output, project]
subprocess.check_call(sequence, stdout=sys.stdout, stderr=sys.stderr)
