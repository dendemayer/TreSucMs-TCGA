import pandas as pd
from natsort import natsort_keygen
"""
this script is suitable for single and multiprojects
the invoked data files are always
# filtering the summary tables on projects(project_id),
-> a filtering on the projects is not necessary, the are already limited to
them
drugs(pharmaceutical_therapy_drug_name), gender(gender)
"""
out_file = snakemake.output[0]
out_file_complement = snakemake.output[1]
meta_file = snakemake.input[0]
DF_meta = pd.read_table(meta_file)
OUTPUT_PATH = snakemake.wildcards[0]
projects = snakemake.wildcards[1].split('_')
drugs = snakemake.wildcards[2].split('_')
genders = snakemake.wildcards[3].split('_')
cutoff = float(snakemake.wildcards[4].split('_')[1])

# DF_meta.set_index(['project_id', 'pharmaceutical_therapy_drug_name',
# 'gender'])
DF_meta_filtered = pd.DataFrame()
for gender in genders:
    DF_meta_filtered = pd.concat(
        [DF_meta_filtered, DF_meta[DF_meta['gender'] == gender]])
# place here the cutoff part, s.t. it is included within the complement and the
# normal summary tables:
"""
the metadata we need for the final summary table are:
vital_state, case_id, drugs,  gender, project_id
invoke the cutoff parameter right here:
    that means:
        if the cutoff parameter is set f.e. to 2 (years), then a patient
        deceased after 2.2 years is set to alive status, to stay consistent
        with the kaplan meier plots, the survivaltime is set to null and the
        years_to_last_follow_up is set to the survivaltime
example:
from:
 id                                   | vital_status | survivaltime       | years_to_last_follow_up ║
 4a2ac5b2-3d7c-4d02-9356-6a70596e5316 | dead         | 2.271232876712329  |                         ║
to:
 id                                   | vital_status | survivaltime       | years_to_last_follow_up ║
 4a2ac5b2-3d7c-4d02-9356-6a70596e5316 | alive        |                    | 2.271232876712329       ║

"""
if cutoff > 0:
    # first get the id's of the cases in which the change shall be performed:
    # dead cases in which survivaltime is greater then the cutoff
    try:
        DF_temp = DF_meta_filtered[DF_meta_filtered['vital_status'] == 'dead']
        id_list = DF_temp[DF_temp['survivaltime'] > cutoff]['bcr_patient_uuid'].to_list()
        # then change vital state to alive
        DF_meta_filtered.set_index('bcr_patient_uuid', inplace=True)
        DF_meta_filtered.loc[id_list, 'vital_status'] = 'alive'
        DF_meta_filtered.loc[
            id_list, 'years_to_last_follow_up'] = DF_meta_filtered.loc[
                id_list, 'survivaltime']
        DF_meta_filtered.loc[id_list, 'survivaltime'] = pd.NA
        DF_meta_filtered.reset_index(inplace=True)
    except Exception as e:
        print(f'while cutoff {cutoff} invokening, no cases found suitable for it')
        print('continuing')


# DF_meta_filtered filtered on the right gender, use this DF to proceed with
# further filtering (drugs):
DF_meta_filtered_2 = pd.DataFrame()
for drug in drugs:
    DF_meta_filtered_2 = pd.concat(
        [DF_meta_filtered_2,
         DF_meta_filtered[
             DF_meta_filtered['pharmaceutical_therapy_drug_name'] == drug]])

# DF_meta_filtered_2 -> contains the drugs apllied, now drop every item out of
# DF_meta_filtered which apears in DF_meta_filtered_2, thats the complement
# table than (DF_meta_filtered)

DF_meta_filtered.set_index('bcr_patient_uuid', inplace=True)
DF_meta_filtered_2.set_index('bcr_patient_uuid', inplace=True)
DF_complement = DF_meta_filtered.drop(DF_meta_filtered_2.index)


# # Chromosome  Start
# # dead;262311a0-7a10-4934-8fba-11ee581ad738;cisplatin;female;TCGA-CESC
# # chr1    15865   0.937932741729432
# # the data tables look like this:
# # Composite Element REF   Beta_value  Chromosome  Start   End  annotations...
# # cg00000029  0.438094185668971   chr16   53434200    53434201
# # # where to find the datatables: this summary table is metilene specific, so
# # that can be set fix
# # OUTPUT_PATH/project_id/metilene/data_files/jhu-usc.edu_CESC.HumanMethylation450.9.lvl-3.TCGA-DS-A0VM-01A-11D-A10W-05.gdc_hg38.txt
DF_meta_filtered_2['file_path'] = OUTPUT_PATH + '/' + DF_meta_filtered_2['project_id'] + '/metilene/data_files/' + DF_meta_filtered_2['filename']
DF_complement['file_path'] = OUTPUT_PATH + '/' + DF_complement['project_id'] + '/metilene/data_files/' + DF_complement['filename']
def return_data_DF(Series):
    # create the col_name:
    name_list = []
    name_list.append(Series.loc['vital_status'])
    name_list.append(Series.name)
    name_list.append(Series.loc['pharmaceutical_therapy_drug_name'])
    name_list.append(Series.loc['gender'])
    name_list.append(Series.loc['project_id'])
    col_name = ';'.join(name_list)
    return pd.read_table(
        Series.loc['file_path']).set_index(
            ['Chromosome', 'Start']).rename(
                {'Beta_value': col_name}, axis=1)[col_name]

def write_final_DFs(DF, out_file):
    DFs = []
    for index in DF.index:
        DFs.append(return_data_DF(DF.loc[index, :]))
    if len(DFs) > 0:
        DF_final = pd.concat(DFs, axis=1)
        DF_final = DF_final.sort_values(by=['Chromosome', 'Start'], key=natsort_keygen())
        DF_final.drop(('*', -1), inplace=True)
    else:
        DF_final = pd.DataFrame()
    # print(f'writing {out_file}')
    DF_final.to_csv(out_file, sep='\t')

write_final_DFs(DF_meta_filtered_2, out_file)
write_final_DFs(DF_complement, out_file_complement)

# # metilene input:
# # The input consists of a single __SORTED__ (for genomic positions)
# # tab-separated file. It must contain a header line of the format:
# # | chr | pos | g1_xxx | g1_xxx | [...] | g2_xxx | g2_xxx | [...] |
# # create a list of DF with index chr, pos , multi index col consisting of
