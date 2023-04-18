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

R_path = snakemake.input.R_path
summary_table = snakemake.input.summary_table
summary_table_info = snakemake.input.summary_table_info
project = snakemake.wildcards.project
summary_out = snakemake.output.summary_out
deseq_output = os.path.split(summary_out)[0]

# if the input table is empty, write empty tables for ever output table which
# shall be created
if pd.read_table(summary_table, sep='\t').empty:
    for out in snakemake.output:
        open(out, 'a').close()
    os._exit(0)


sequence = ['Rscript', R_path, summary_table, summary_table_info, deseq_output, project]
subprocess.check_call(sequence, stdout=sys.stdout, stderr=sys.stderr)
# subprocess.Popen(sequence, stdout=sys.stdout, stderr=sys.stderr)
