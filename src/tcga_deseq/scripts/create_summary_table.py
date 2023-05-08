import pandas as pd
import sys
import os

# from shared.scripts.create_merged_tables import PROJECT
"""
creating the summary tables (input for DESeq2)
data table:
/scr/dings/PEVO/NEW_downloads_3/DEseq_31_8/TCGA-CESC/carboplatin_carboplatin,paclitaxel_cisplatin/DRUG_combi_complement_summary_dead_alive.tsv.
    0   1   2
ENSG00000000003 2541    3382
ENSG00000000005 0   0   8   3
INFO table.
/scr/dings/PEVO/NEW_downloads_3/DEseq_31_8/TCGA-CESC/carboplatin_carboplatin,paclitaxel_cisplatin/DRUG_combi_complement_summary_dead_alive_INFO.tsv
        0       1       2       3       4
vital_status    alive   alive   alive   alive
gender  female  female  female  female  femal
PROJECT TCGA-CESC       TCGA-CESC       TCGA-
drugnames       cisplatinBgemcitabin    erbit
case_id e711b88e-fa0e-4781-8549-7a53ed3a13a1
"""

sys.stderr = sys.stdout = open(snakemake.log[0], 'w')

meta_table = snakemake.input.meta_table
gender = snakemake.wildcards.gender
genders = gender.split('_')
drug_combi = snakemake.wildcards.drug_combi
drugs = drug_combi.split('_')
OUTPUT_PATH = snakemake.wildcards.output_path
PROJECT = snakemake.wildcards.project
cutoff = snakemake.wildcards.cutoff
summary_table = snakemake.output.summary_table
summary_table_info = snakemake.output.summary_table_info
summary_table_complement = snakemake.output.summary_table_complement
summary_table_complement_info = snakemake.output.summary_table_complement_info

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

DF_meta = pd.read_table(meta_table).set_index('bcr_patient_uuid')
# filter on gender, don't use direct loc indexing, if one gender not available,
# exception occurs:
index_list = []
for g in genders:
    index_list.extend(DF_meta[DF_meta['gender'] == g].index.tolist())

DF_meta = DF_meta.loc[index_list, :]

# filter on drug selection, both normal and complement:
index_list = []
for d in drugs:
    index_list.extend(DF_meta[DF_meta['pharmaceutical_therapy_drug_name'] == d].index.tolist())

DF_meta_complement = DF_meta[~DF_meta.index.isin(index_list)]
DF_meta_complement = DF_meta_complement.sort_values(by=['vital_status', 'gender', 'project_id', 'pharmaceutical_therapy_drug_name', 'bcr_patient_uuid'])

DF_meta = DF_meta.loc[index_list, :]
DF_meta = DF_meta.sort_values(by=['vital_status', 'gender', 'project_id', 'pharmaceutical_therapy_drug_name', 'bcr_patient_uuid'])

# for both, add the abs path to the refered filename
DF_meta_complement['full_path'] = OUTPUT_PATH + '/' + DF_meta_complement['project_id'] + '/DESeq2/data_files/' + DF_meta_complement['filename']
DF_meta['full_path'] = OUTPUT_PATH + '/' + DF_meta['project_id'] + '/DESeq2/data_files/' + DF_meta['filename']

full_path_complement_list = DF_meta_complement['full_path'].values.tolist()
full_path_list = DF_meta['full_path'].values.tolist()

def create_data_summary(full_path_list):
    DF_summary = pd.DataFrame()
    for file in full_path_list:
        DF_temp = pd.read_table(file, header=None)[pd.read_table(file, header=None)[0].str.contains('ENSG')]
        DF_temp[0] = DF_temp[0].apply(lambda x: x.split('.')[0])
        DF_temp = DF_temp.set_index(0)
        DF_summary = pd.concat([DF_summary, DF_temp], axis=1)
    DF_summary.columns = range(0, DF_summary.shape[1])
    DF_summary.index.name = None
    return DF_summary

# write the datatables
DF_summary = create_data_summary(full_path_list)
DF_summary.to_csv(summary_table, sep='\t')
DF_summary_complement = create_data_summary(full_path_complement_list)
DF_summary_complement.to_csv(summary_table_complement, sep='\t')

# creating the belongin info tables:
DF_meta.reset_index().loc[:, ['vital_status', 'gender', 'project_id', 'pharmaceutical_therapy_drug_name', 'bcr_patient_uuid']].T.to_csv(summary_table_info, sep='\t')
DF_meta_complement_T = DF_meta_complement.reset_index().loc[:, ['vital_status', 'gender', 'project_id', 'pharmaceutical_therapy_drug_name', 'bcr_patient_uuid']].T
# why is the renaming step here?
# DF_meta_complement_T = DF_meta_complement_T.rename({'project_id': 'PROJECT', 'pharmaceutical_therapy_drug_name': 'drugnames', 'bcr_patient_uuid': 'case_id'})
DF_meta_complement_T.to_csv(summary_table_complement_info, sep='\t')
