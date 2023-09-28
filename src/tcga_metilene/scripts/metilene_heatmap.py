import pandas as pd
import sys
import os
import subprocess

"""
the actual input for R is build dynamically, after the requested heatmap output
is created, the temporary files (beta_val_inname and info_inname) are removed
-> temp input for are must be named with the specific DMR, parallel instances
of other DMR requests would interfere otherwise
"""

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print('# snakemake inputs:')
[print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

metilene_intersect = snakemake.input.metilene_intersect
metilene_out = snakemake.input.metilene_out
R_heatmap_script = snakemake.input.R_heatmap_script

# # start_tsv is the table which is based on the lifeline_DMR_plot, and holds the
# # extract start value
pdf_heatmap_out = snakemake.output.pdf_heatmap_out
output_path = snakemake.wildcards.output_path
project = snakemake.wildcards.project
gender = snakemake.wildcards.gender
cutoff = snakemake.wildcards.cutoff
DMR = snakemake.wildcards.DMR
drug_combi = snakemake.wildcards.drug_combi

###############################################################################
#                                  test set                                   #
###############################################################################

# # snakemake inputs:
# metilene_intersect = "/scr/palinca/gabor/TCGA-pipeline_4/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_intersect.tsv"
# metilene_out = "/scr/palinca/gabor/TCGA-pipeline_4/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_qval.0.05.out"
# R_heatmap_script = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/tcga_metilene/scripts/metilene_heatmap.R"
# # snakemake output:
# pdf_heatmap_out = "/scr/palinca/gabor/TCGA-pipeline_4/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_intersect_heatmaps_beta_value_chr7_27163810_27167288.pdf"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_4"
# project = "TCGA-CESC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_0"
# DMR = "chr7_27163810_27167288"

# import countmatrix as:
#                 │ 0       │ 1     │ 2
# ENSG00000000003 │ 4134    │ 4765  │ 3825
# ENSG00000000005 │ 1       │ 0     │ 1
# ENSG00000000419 │ 2418    │ 2831  │ 1854
# and the belonging header:
#
#                    │ 0                  │ 1                  │ 2
# vital_status       │ alive              │ alive              │ alive
# gender             │ female             │ female             │ female
# project_id         │ TCGA-CESC          │ TCGA-CESC          │ TCGA-CESC
# pharmaceutical_the…│ carboplatin,paclit…│ carboplatin,paclit…│ carboplatin,paclit…
# bcr_patient_uuid   │ 63dcfa32-c69e-45ad…│ afb5442e-0249-475b…│ e2ff7cfd-a101-45b9…

# found for every positions of the DMRs

# therefore needed is the metilene_intersect.tsv table:
# vital_status │           │           │                    │ dead               │ dead
# case_id      │           │           │                    │ b726f18b-b996-4978…│ d3ff2aeb-46be-4f5
# drugs        │           │           │                    │ carboplatin        │ carboplatin
# gender       │           │           │                    │ female             │ female
# projects     │           │           │                    │ TCGA-CESC          │ TCGA-CESC
# Chromosome   │ Start     │ End       │ region             │                    │
# chr6         │ 166530703 │ 166530704 │ chr6_166530702_166…│ 0.499550999056649  │ 0.493254585682986
# chr6         │ 166530795 │ 166530796 │ chr6_166530702_166…│ 0.502649598280444  │ 0.466584973276633
# chr7         │ 27112964  │ 27112965  │ chr7_27112963_2711…│ 0.289101666827651  │ 0.125578700648676
# chr7         │ 27113280  │ 27113281  │ chr7_27112963_2711…│ 0.357887952854563  │ 0.186961366011444

cutoff = cutoff.split('_')[1]

met_int_DF = pd.read_table(metilene_intersect, index_col=[0, 1, 2, 3], header=[0, 1, 2, 3, 4], na_values="-").dropna(how='all')
# limit the beta_val to plot right away to the recent DMR
met_int_DF = met_int_DF.loc[(slice(None), slice(None), slice(None), DMR ), :]
# selecting increasing and decreasing methylations
# sorting on: vital_status, project, gender:
met_int_DF = met_int_DF.sort_index(axis=1, level=[0, 4, 3])
# the row notation is going to be chr_start
met_int_DF = met_int_DF.reset_index()
met_int_DF['chr_start'] = met_int_DF['Chromosome'] + '_' + met_int_DF['Start'].astype('str')
met_int_DF = met_int_DF.set_index('chr_start').drop(['Chromosome', 'Start', 'End', 'region'], axis=1)
cols = met_int_DF.columns
cols = cols.to_frame().reset_index(drop=True).loc[:, ['vital_status', 'gender', 'projects', 'drugs']]
# just keep gender, projects, drugs, if they are multivariat, otherwise their
# annotation is not neccessary in the heatmap:
for col in cols.columns:
    if len(cols[col].value_counts()) == 1:
        cols.drop(col,axis=1,inplace=True)

cols = cols.T
# writing temporary input files for R, after the subprocess R instance ran
# succesfully, those tempfiles can be deleted
header_temp = pdf_heatmap_out.replace('.pdf', '_header_temp.tsv')
print(f'saving temp file for metilene_heatmap.R {header_temp}')
cols.to_csv(header_temp, sep='\t')
met_int_DF.columns = range(0, met_int_DF.shape[1])
beta_val_temp = pdf_heatmap_out.replace('.pdf', '_input_temp.tsv')
print(f'saving temp file for metilene_heatmap.R {beta_val_temp}')
# breakpoint() # metilene needs empty values to be written as '-', but thats
# problematic in R pheatmap, change them to "" or delete the entyre row of all
# are NaNs
# this frame is handed over to R, in case there are values missing for some
# patients, the still would be NAs, write 0, s.t. they can be plottet correctly
met_int_DF.to_csv(beta_val_temp, sep='\t', na_rep="0")


###
# handing over  both files needed for the heatmap, and infos for the main plot:
# project, gender, cutoff, DMR, drug_combi -> wildcards
# q_value, mmd, CpGs  -> out of metilene out
###
# getting the met_out vars:
met_out = pd.read_table(metilene_out, index_col=None, header=None)
met_out.columns = ['chr', 'start', 'stop', 'q-value',  'mean methylation difference', '#CpGs', 'mean alive', 'mean dead']
met_out = met_out.set_index(['chr', 'start', 'stop'])
pos_indexer = DMR.split('_')
# make int out of the postitions, for correct indexing the met out frame
pos_indexer[1:3] = [int(i) for i in pos_indexer[1:3]]
# index the met_out frame for the recent DMR:
met_out = met_out.loc[tuple(pos_indexer), : ]
q_value, mean_methylation_difference, CpGs, mean_alive, mean_dead = tuple([str(i) for i in met_out])
CpGs =str(int(float(CpGs)))
sequence = ['Rscript', R_heatmap_script, beta_val_temp, header_temp, pdf_heatmap_out, project, gender, drug_combi, cutoff,  DMR, q_value, mean_methylation_difference, CpGs]
call  = subprocess.check_call(sequence, stdout=sys.stdout, stderr=sys.stderr)

os.remove(beta_val_temp)
os.remove(header_temp)
