import pandas as pd
pd.concat([pd.read_table(i) for i in snakemake.input]).to_csv(
    snakemake.output[0], sep='\t', index=False)
