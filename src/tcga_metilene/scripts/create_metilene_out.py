import os
import pandas as pd
import subprocess
import sys

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

summary_table = snakemake.input[0]
met_prefix = snakemake.output[1]
gender = snakemake.wildcards.gender
# check whether we have actually alive and dead groups we can compare
# check if alive is available in columns (ommit Start and Chromosome col ->
# iloc[2:,0]), the value count of it must be 2, than
# alive and dead is present, if its 1, then either all cols hold True for alive
# or False
if len(pd.read_table(summary_table, nrows=1).T.reset_index().iloc[2:,0].str.contains('alive').value_counts()) == 1:
    temp_DF = pd.DataFrame()
    for i in snakemake.output:
        temp_DF.to_csv(i)
# if input table is empty, dont start at all:
elif pd.read_table(summary_table).empty:
    temp_DF = pd.DataFrame()
    for i in snakemake.output:
        temp_DF.to_csv(i)
# in case female_male gender, but there is just one available, stop
elif pd.read_table(summary_table, nrows=1).T.reset_index().iloc[2:,0].apply(lambda x: x.split(';')[3]).nunique() == 1 and len(gender.split('_')) == 2:
    temp_DF = pd.DataFrame()
    for i in snakemake.output:
        temp_DF.to_csv(i)
else:
    # m = snakemake.wildcards[5]
    # M = snakemake.wildcards[6]
    # d = snakemake.wildcards[7]
    m = 3  # -m, --mincpgs	Integer	10	The minimum # of CpGs in a DMR
    M = 1000  # -M, --maxdist	Integer	300	The allowed nt distance between two CpGs within a DMR
    d = 0.03  # -d, --minMethDiff	double	0.1	The minimum mean methylation difference for calling DMRs
    resource_path = os.path.join(os.path.split(os.path.split(__file__)[0])[0], 'resources')
    metilene = os.path.join(resource_path, 'metilene_linux64')
    metilene_sorted_out = snakemake.output[0]
    # metilene -m 3 -M 1000 -d 0.03 -a alive -b dead summary_for_metilene.tsv | sort -V -k1,1 -k2,2n > metilene_out_sorted.tsv
    metilene_call = [metilene, f'-m {m}', f'-M {M}', f'-d {d}', f'-t {snakemake.threads}', '-a alive', '-b dead', summary_table , '| sort -V -k1,1 -k2,2n', f'> {metilene_sorted_out}']
    # calling the metilene tool:
    subprocess.check_call(" ".join(metilene_call), stdout=sys.stdout, stderr=sys.stderr, shell=True)
    print('called metilene with the following configuration:')
    print(" ".join(metilene_call))

    metilene_filter = os.path.join(resource_path, 'metilene_output.pl')
    # perl metilene_output.pl -q metilene_out_sorted.tsv -a 'alive' -b 'dead' -c 3 -d 0.03
    met_prefix = met_prefix.replace('_qval.0.05.out', '')
    # metilene_filter_call = [f'perl {metilene_filter}', f'-q {metilene_sorted_out}' ,'-a alive', '-b dead', f'-c {m}', f'-d {d}']
    metilene_filter_call = [f'perl {metilene_filter}', f'-q {metilene_sorted_out}' ,'-a alive', '-b dead', f'-c {m}', f'-d {d}', f'-o {met_prefix}']

    # calling the filter step:
    subprocess.check_call(" ".join(metilene_filter_call), stdout=sys.stdout, stderr=sys.stderr, shell=True)
    print('called metilene-filter with the following configuration:')
    print(" ".join(metilene_filter_call))
