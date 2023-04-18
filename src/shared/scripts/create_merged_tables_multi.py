import sys
import pandas as pd
sys.stderr = sys.stdout = open(snakemake.log[0], 'w')
pd.concat([pd.read_table(i) for i in snakemake.input]).to_csv(
    snakemake.output[0], sep='\t', index=False)
