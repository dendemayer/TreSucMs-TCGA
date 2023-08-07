import sys
import re
import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], 'w')

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

annot_input = snakemake.input.annot_input[0]
annot_output = snakemake.output.annot_output[0]


# ###############################################################################
# #                                 test_input
# # ###############################################################################
# # snakemake inputs:
# annot_input = "/scr/palinca/gabor/TCGA-pipeline/metadata/gencode.v36.annotation.gtf.gz"
# # snakemake output:
# annot_output = "/scr/palinca/gabor/TCGA-pipeline/metadata_processed/gencode.v36.annotation.gtf_genes_transcripts.gz"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline"
# ###############################################################################
# #                                 test_input
# ###############################################################################

# col 8 has different col numbers in gene and exon rows, first get the gene and
# parse out ENSG, gene_type, gene_status and gene_name, than add exon infos:
DF_annot = pd.read_table(annot_input, skiprows=5, header=None).set_index([2]).sort_index().loc['gene', :].reset_index(drop=True).set_index([0, 3, 4]).loc[:, 8].to_frame()
DF_annot['ENSG'] = DF_annot[8].apply(lambda x: re.search('ENSG\d*', x.split(';')[0]).group(0))
DF_annot['gene_type'] = DF_annot[8].apply(lambda x: re.search('".*"', x.split(';')[1]).group(0).replace('"', ''))
DF_annot['gene_status'] = DF_annot[8].apply(lambda x: re.search('".*"', x.split(';')[2]).group(0).replace('"', ''))
DF_annot['gene_name'] = DF_annot[8].apply(lambda x: re.search('".*"', x.split(';')[3]).group(0).replace('"', ''))
# add the ENSG annotations:
DF_annot = DF_annot.reset_index().drop(8,axis=1).rename({0: 'chr', 3: 'start', 4: 'stop'}, axis=1)
DF_annot['gtf_row_type'] = 'gene'

# get rows with transcr informatin:
#### exon and transcr hold same ENSG information:
# zgrep -w -E 'transcript' gencode.v36.annotation.gtf.gz | grep ENSEMBL | awk '{match($0, "ENSG[0-9]+", a); print a[0]}' | sort -u > all_transcript_ENSG.txt
# zgrep -w -E 'exon' gencode.v36.annotation.gtf.gz | grep ENSEMBL | awk '{match($0, "ENSG[0-9]+", a); print a[0]}' | sort -u > all_exon_ENSG.txt
# cat all_exon_ENSG.txt all_transcript_ENSG.txt | sort | uniq -u
DF_annot_transcr = pd.read_table(annot_input, skiprows=5, header=None).set_index(2).sort_index().loc['transcript', :].reset_index(drop=True).set_index([0, 3, 4]).loc[:, 8].to_frame()

DF_annot_transcr['ENSG'] = DF_annot_transcr[8].apply(lambda x: re.search('ENSG\d*', x.split(';')[0]).group(0))
# in contrast to the gene DF, the gene_type is in the 2cnd col, not the first:
DF_annot_transcr['gene_type'] = DF_annot_transcr[8].apply(lambda x: re.search('".*"', x.split(';')[2]).group(0).replace('"', ''))
DF_annot_transcr['gene_status'] = DF_annot_transcr[8].apply(lambda x: re.search('".*"', x.split(';')[3]).group(0).replace('"', ''))
DF_annot_transcr['gene_name'] = DF_annot_transcr[8].apply(lambda x: re.search('".*"', x.split(';')[4]).group(0).replace('"', ''))
DF_annot_transcr['gtf_row_type'] = 'transcript'
# add the ENSG annotations:
DF_annot_transcr = DF_annot_transcr.reset_index().drop(8,axis=1).rename({0: 'chr', 3: 'start', 4: 'stop'}, axis=1)
DF_final = pd.concat([DF_annot, DF_annot_transcr]).drop_duplicates()
print(f'saving {annot_output}')
DF_final.to_csv(annot_output, index=False, sep='\t')