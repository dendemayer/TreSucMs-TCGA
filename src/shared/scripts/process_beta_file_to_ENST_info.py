import pandas as pd
import re

# prepare the header sorted file:
# head -1 jhu-usc.edu_CESC.HumanMethylation450.1.lvl-3.TCGA-C5-A1BE-01B-11D-A13Z-05.gdc_hg38 > header.txt
# sort -V -k3,3 -k4,4n <(sed '1d' jhu-usc.edu_CESC.HumanMethylation450.1.lvl-3.TCGA-C5-A1BE-01B-11D-A13Z-05.gdc_hg38.txt) >> header_sorted.txt
DF = pd.read_table('/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/resources/header_sorted.txt').loc[:, ['Chromosome', 'Start', 'End', 'Gene_Type', 'Gene_Symbol', 'Transcript_ID', 'Feature_Type']]
DF = DF[DF['Chromosome'] != '*']
# exclude multiple values for Gene_Type etc in one row:
DF['Gene_Type'] = DF['Gene_Type'].apply(lambda x: ';'.join(sorted(list(set(x.split(';'))))))
DF['Gene_Symbol'] = DF['Gene_Symbol'].apply(lambda x: ';'.join(sorted(list(set(x.split(';'))))))
# addititonally shorten the ENST00000331581.xx notation before the digit point
DF['Transcript_ID'] = DF['Transcript_ID'].apply(lambda x: ';'.join([re.search('ENST\d+', i).group(0) for i in  sorted(list(set(x.split(';')))) if i != '.']))
# duplicated cooedinates can be dropped:, before:
# [480457 rows x 7 columns]
# after:
# [480436 rows x 7 columns]
DF.drop_duplicates(['Chromosome', 'Start', 'End'], inplace=True)
DF.rename({'Chromosome': 'chr', 'Start': 'start', 'End': 'stop','Gene_Type': 'gene_type', 'Gene_Symbol': 'gene_name', 'Transcript_ID': 'ENST'}, axis=1, inplace=True)
DF['gene_status'] = '-'
DF.to_csv('/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/resources/annot_from_betafile.tsv.gz', sep='\t', index=False)
