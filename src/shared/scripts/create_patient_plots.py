import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import subprocess
from PyPDF2 import PdfMerger
import sys

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

meta_table = snakemake.input.meta_table

plot_file_age = snakemake.output.plot_file_age
plot_file_age_in_therapy = snakemake.output.plot_file_age_in_therapy
plot_file_age_not_in_therapy = snakemake.output.plot_file_age_not_in_therapy

plot_file_survival = snakemake.output.plot_file_survival
plot_file_survival_in_therapy = snakemake.output.plot_file_survival_in_therapy
plot_file_survival_not_in_therapy = snakemake.output.plot_file_survival_not_in_therapy

out_md_vital = snakemake.output.out_md_vital
out_pdf_vital = snakemake.output.out_pdf_vital
final_pdf = snakemake.output.final_pdf

drug_str = snakemake.params.drug_str
project = snakemake.wildcards.project
cutoff = snakemake.wildcards.cutoff

###############################################################################
#                               test input set                                #
###############################################################################
# meta_table = '/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined.tsv'

# plot_file_age = '/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined_age.pdf'
# plot_file_age_in_therapy = '/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined_age_in_therapy.pdf'
# plot_file_age_not_in_therapy = '/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined_age_not_in_therapy.pdf'

# plot_file_survival = '/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined_survival.pdf'
# plot_file_survival_in_therapy = '/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined_survival_in_therapy.pdf'
# plot_file_survival_not_in_therapy = '/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined_survival_not_in_therapy.pdf'

# out_md_vital = '/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined_vital.md'
# out_pdf_vital = '/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined_vital.pdf'
# final_pdf = '/scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined_final.pdf'
###############################################################################
#                               test input set                                #
###############################################################################
project_str = ' + '.join(project.split('_'))
drugs = drug_str.split('_')
drug_header = '; '.join(drug_str.split('_'))
cutoff = cutoff.split('_')[1]

meta_DF = pd.read_table(meta_table)
meta_DF = meta_DF.loc[:, ['vital_status', 'gender', 'survivaltime', 'years_to_last_follow_up', 'age_at_diagnosis', 'project_id', 'pharmaceutical_therapy_drug_name']]
meta_DF['Survivaltime'] = meta_DF.loc[:, 'survivaltime'].fillna(meta_DF.loc[:, 'years_to_last_follow_up'])
meta_DF.sort_values(['vital_status', 'gender'], ascending=False, inplace=True)

# add a hue col for in therapy or not in therapy
meta_DF.set_index('pharmaceutical_therapy_drug_name').loc[drugs, :]
meta_DF['in_therapy'] = [True if (i in drugs) else False for i in meta_DF['pharmaceutical_therapy_drug_name'].tolist()]

meta_add_proj = meta_DF.copy(deep=True)
meta_add_proj['project_id'] = project_str
meta_add_proj = pd.concat([meta_DF, meta_add_proj])
meta_add_proj_gen = meta_add_proj.copy(deep=True)
meta_add_proj_gen['gender'] = 'female + male'
meta_add_proj_gen = pd.concat([meta_add_proj, meta_add_proj_gen])

sns.set_style("whitegrid")
g = sns.displot(meta_add_proj_gen[meta_add_proj_gen['in_therapy']],
                x="age_at_diagnosis", col="project_id", row="gender",
                binwidth=len(meta_DF['gender'].value_counts()),
                height=len(meta_DF['project_id'].value_counts()) + 0.7,
                facet_kws=dict(margin_titles=True), hue='vital_status',
                kde=True)
plt.subplots_adjust(top=0.9)
g.fig.suptitle("Age at Therapy start,\npatients in therapy", y=0.99)
plt.savefig(plot_file_age_in_therapy)
plt.close('all')

sns.set_style("whitegrid")
g = sns.displot(meta_add_proj_gen[~meta_add_proj_gen['in_therapy']],
                x="age_at_diagnosis", col="project_id", row="gender",
                binwidth=len(meta_DF['gender'].value_counts()),
                height=len(meta_DF['project_id'].value_counts()) + 0.7,
                facet_kws=dict(margin_titles=True), hue='vital_status',
                kde=True)
plt.subplots_adjust(top=0.9)
g.fig.suptitle("Age at Therapy start,\npatients not in therapy", y=0.99)
plt.savefig(plot_file_age_not_in_therapy)
plt.close('all')

sns.set_style("whitegrid")
g = sns.displot(meta_add_proj_gen,
                x="age_at_diagnosis", col="project_id", row="gender",
                binwidth=len(meta_DF['gender'].value_counts()),
                height=len(meta_DF['project_id'].value_counts()) + 0.7,
                facet_kws=dict(margin_titles=True), hue='vital_status',
                kde=True)
plt.subplots_adjust(top=0.9)
g.fig.suptitle("Age at Therapy start,\nall patients", y=0.99)
plt.savefig(plot_file_age)
plt.close('all')

sns.set_style("whitegrid")
g = sns.displot(meta_add_proj_gen[meta_add_proj_gen['in_therapy']],
                x="Survivaltime", col="project_id", row="gender",
                binwidth=len(meta_DF['gender'].value_counts()),
                height=len(meta_DF['project_id'].value_counts()) + 0.7,
                facet_kws=dict(margin_titles=True), hue='vital_status',
                kde=True)
plt.subplots_adjust(top=0.9)
g.fig.suptitle("Survivaltime at Therapy start,\npatients in therapy", y=0.99)
plt.savefig(plot_file_survival_in_therapy)
plt.close('all')

sns.set_style("whitegrid")
g = sns.displot(meta_add_proj_gen[~meta_add_proj_gen['in_therapy']],
                x="Survivaltime", col="project_id", row="gender",
                binwidth=len(meta_DF['gender'].value_counts()),
                height=len(meta_DF['project_id'].value_counts()) + 0.7,
                facet_kws=dict(margin_titles=True), hue='vital_status',
                kde=True)
plt.subplots_adjust(top=0.9)
g.fig.suptitle("Survivaltime at Therapy start,\npatients not in therapy", y=0.99)
plt.savefig(plot_file_survival_not_in_therapy)
plt.close('all')

sns.set_style("whitegrid")
g = sns.displot(meta_add_proj_gen,
                x="Survivaltime", col="project_id", row="gender",
                binwidth=len(meta_DF['gender'].value_counts()),
                height=len(meta_DF['project_id'].value_counts()) + 0.7,
                facet_kws=dict(margin_titles=True), hue='vital_status',
                kde=True)
plt.subplots_adjust(top=0.9)
g.fig.suptitle("Survivaltime at Therapy start,\nall patients", y=0.99)
plt.savefig(plot_file_survival)
plt.close('all')

meta_DF['in_therapy'] = ['in therapy' if (i in drugs) else 'not in therapy' for i in meta_DF['pharmaceutical_therapy_drug_name'].tolist()]

with open(out_md_vital, 'w') as f:
    f.write(f'# Patients included of {project_str}, with therapy {drug_header} and Cutoff {cutoff}: \n- {len(meta_DF)} patients in total\n\n# vital grouped:\n\n')
meta_DF.loc[:, ['vital_status', 'gender', 'project_id', 'in_therapy']].value_counts(subset=['vital_status']).to_frame().reset_index().rename({0: '#'}, axis=1).set_index(['vital_status']).to_markdown(out_md_vital, mode='a')

with open(out_md_vital, 'a') as f:
    f.write('\n\n# therapy grouped:\n\n')
meta_DF.loc[:, ['vital_status', 'gender', 'project_id', 'in_therapy']].value_counts(subset=['in_therapy']).to_frame().reset_index().rename({0: '#'}, axis=1).set_index(['in_therapy']).to_markdown(out_md_vital, mode='a')

with open(out_md_vital, 'a') as f:
    f.write('\n\n# vital and therapy grouped:\n\n')
meta_DF.loc[:, ['vital_status', 'gender', 'project_id', 'in_therapy']].value_counts(subset=['vital_status', 'in_therapy']).to_frame().reset_index().rename({0: '#'}, axis=1).set_index(['vital_status']).to_markdown(out_md_vital, mode='a')

with open(out_md_vital, 'a') as f:
    f.write('\n\n# vital and gender grouped:\n\n')
meta_DF.loc[:, ['vital_status', 'gender', 'project_id', 'in_therapy']].value_counts(subset=['vital_status', 'gender']).to_frame().reset_index().rename({0: '#'}, axis=1).set_index(['vital_status']).to_markdown(out_md_vital, mode='a')

with open(out_md_vital, 'a') as f:
    f.write('\n\n# vital, gender and project grouped:\n\n')
meta_DF.loc[:, ['vital_status', 'gender', 'project_id', 'in_therapy']].value_counts(subset=['vital_status', 'gender', 'project_id']).to_frame().reset_index().rename({0: '#'}, axis=1).set_index(['vital_status']).to_markdown(out_md_vital, mode='a')

with open(out_md_vital, 'a') as f:
    f.write('\n\n# vital, gender, project and therapy grouped:\n\n')
meta_DF.loc[:, ['vital_status', 'gender', 'project_id', 'in_therapy']].value_counts(subset=['vital_status', 'gender', 'project_id', 'in_therapy']).to_frame().reset_index().rename({0: '#'}, axis=1).set_index(['vital_status']).to_markdown(out_md_vital, mode='a')


sequence = ['pandoc', out_md_vital, '-t', 'latex', '-o', out_pdf_vital ]
# such that pandoc finds the iftex.sty, chdir to this sourcefile, than go one
# up and into resources/iftex
os.chdir(os.path.join(os.path.dirname(__file__), os.pardir, 'resources', 'iftex'))
call  = subprocess.check_call(sequence, stdout=sys.stdout, stderr=sys.stderr)
# call  = subprocess.check_call(sequence)

merger = PdfMerger()

pdfs_to_merge = [out_pdf_vital, plot_file_age_in_therapy, plot_file_age_not_in_therapy, plot_file_age, plot_file_survival_in_therapy, plot_file_survival_not_in_therapy, plot_file_survival]

for pdf in pdfs_to_merge:
    merger.append(pdf)

merger.write(final_pdf)
merger.close()
