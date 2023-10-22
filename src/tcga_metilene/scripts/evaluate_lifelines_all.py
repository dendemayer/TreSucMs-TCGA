import os
import pandas as pd
from PyPDF2 import PdfMerger
import re
import sys

###############################################################################
#                              snakemake inputs                               #
###############################################################################

sys.stderr = sys.stdout = open(snakemake.log[0], 'w')

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

metilene_lifeline_aggregated = snakemake.input[0]
metilene_lifeline_eval = snakemake.output[0] # metilene_lifeline_eval
metilene_lifeline_eval_pdfs = snakemake.output[1] # metilene_lifeline_eval_pdfs
###############################################################################

###############################################################################
#                                  test set                                   #
###############################################################################

## snakemake inputs:
#metilene_lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_aggregated.tsv.gz"
#evaluate_script = "../tcga_metilene/scripts/evaluate_lifelines_all.py"
## snakemake output:
#metilene_lifeline_eval = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated.tsv.gz"
#metilene_lifeline_eval_pdfs = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated.pdf"
## snakemake wildcards:
#output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
#project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
#drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
#gender = "male"
#cutoff = "cutoff_5"

### snakemake inputs: metilene_lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_aggregated.tsv.gz"
#evaluate_script = "../tcga_metilene/scripts/evaluate_lifelines_all.py"
## snakemake output:
#metilene_lifeline_eval = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated.tsv.gz"
#metilene_lifeline_eval_pdfs = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated.pdf"
## snakemake wildcards:
#output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
#project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
#drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
#gender = "female_male"
#cutoff = "cutoff_0"

## snakemake inputs:
##KeyError: 'chr6_28336674'
#metilene_lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_aggregated.tsv.gz"
#evaluate_script = "../tcga_metilene/scripts/evaluate_lifelines_all.py"
## snakemake output:
#metilene_lifeline_eval = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated.tsv.gz"
#metilene_lifeline_eval_pdfs = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated.pdf"
## snakemake wildcards:
#output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
#project = "TCGA-HNSC"
#drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
#gender = "male"
#cutoff = "cutoff_5"


# # snakemake inputs:
#metilene_lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_aggregated.tsv.gz"
## snakemake output:
#metilene_lifeline_eval = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated.tsv.gz"
#metilene_lifeline_eval_pdfs = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated.pdf"
## snakemake wildcards:
#output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
#project = "TCGA-HNSC"
#drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
#gender = "female"
#cutoff = "cutoff_0"

# # snakemake inputs:
# metilene_lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_aggregated.tsv.gz"
# # snakemake output:
# metilene_lifeline_eval = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated.tsv.gz"
# metilene_lifeline_eval_pdfs = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated.pdf"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
# project = "TCGA-CESC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_5"

# # snakemake inputs:
# metilene_lifeline_aggregated = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_4/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene2/metilene2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5/metilene2_lifelines_aggregated.tsv.gz"
# # snakemake output:
# metilene_lifeline_eval = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_4/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene2/metilene2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5/metilene2_lifelines_evaluated.tsv.gz"
# # snakemake wildcards:
# output_path = "/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_4"
# project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_0"
# end = "tsv.gz"
###############################################################################
#                                  test set                                   #
###############################################################################

DF_eval = pd.read_table(metilene_lifeline_aggregated)
if DF_eval.empty:
    DF_eval.to_csv(metilene_lifeline_eval, sep='\t')
    open(metilene_lifeline_eval_pdfs, 'a').close()
    os._exit(0)

## starting point are the scored base_plots, sort them p_sum wise ->
## sorted after the combined p_value of the scored plots (2 out of 3)
#DF_base_sort_temp = DF_eval.set_index(['plot_type', 'scored']).sort_index().loc[('base_plot', True), :].sort_values('p_sum').reset_index()
#DF_base_sort_temp.set_index(['chr_start', 'threshold'], inplace=True)
#DF_eval = DF_eval.set_index(DF_base_sort_temp.index.names).loc[DF_base_sort_temp.index,:]
## for each chr_start take the best result:
#DF_eval_final = pd.concat([pd.concat([DF_eval.loc[start,:].reset_index().iloc[:3,:], pd.DataFrame({'chr_start':[start]})], axis=1).fillna(start) for start in starts]).set_index('chr_start')

#DONE: take the best combination out of the mixture of all thresholds applied:
# limit on the relevant, scored plots:
DF_eval_cross = DF_eval[DF_eval['scored']]
starts = DF_eval_cross.reset_index()['chr_start'].value_counts().index.to_list()
# make a kartesian cross product between the base plots and the validation
# plots, -> inspecting all combinations of available p-valies in dependence of
# the threshold invoked
DF_eval_cross = DF_eval_cross.loc[:, ['chr_start', 'ENSG', 'plot_type', 'threshold', 'p_value', 'scored']].merge(DF_eval.loc[:, ['chr_start', 'ENSG', 'plot_type', 'threshold', 'p_value', 'scored']], how='cross').set_index(['scored_x', 'scored_y']).loc[(True, True), :]
# sum up all possible p_value sums
DF_eval_cross['p_sum_cross'] = DF_eval_cross['p_value_x'] + DF_eval_cross['p_value_y']
# limit on the base_plots
DF_base_temp = DF_eval_cross.set_index('plot_type_x').loc['base_plot'].sort_values('p_sum_cross')
# just take the minimum for base and val plot of same starts: indexing
# chr_start_x and chr_start_y, concatenating the final, best p_sum_cross:
DF_base_temp_cross = pd.concat([ DF_base_temp[DF_base_temp['plot_type_y'] != 'base_plot'].reset_index().set_index(['chr_start_x', 'chr_start_y']).loc[(start, start), :].sort_values('p_sum_cross').iloc[0, :].to_frame() for start in starts ], axis=1).T
# index the start DF (DF_eval_cross) directly on the gathered combinations:
DF_base_plot = DF_base_temp_cross.reset_index().loc[:, ['level_1', 'plot_type_x', 'threshold_x', 'p_sum_cross']]
DF_base_plot.columns =  ['chr_start', 'plot_type', 'threshold', 'p_sum_cross']
DF_val_plot = DF_base_temp_cross.reset_index().loc[:, ['level_1', 'plot_type_y', 'threshold_y', 'p_sum_cross']]
DF_val_plot.columns =   ['chr_start', 'plot_type', 'threshold', 'p_sum_cross']

DF_base_val = pd.concat([DF_base_plot, DF_val_plot]).sort_values(['chr_start', 'plot_type'], ascending=False).set_index(['chr_start', 'plot_type', 'threshold'])
DF_eval = DF_eval.set_index(['chr_start', 'plot_type', 'threshold'])

DF_eval_final = DF_base_val.join(DF_eval)
# ordering p_sum_cross wise:
DF_eval_final = DF_eval_final.sort_values('p_sum_cross').reset_index()
# ordering each 2er postition pair to base_plot val_plot:
chr_start_list = DF_eval_final['chr_start'].drop_duplicates().to_list()
DF_eval_final = pd.concat([ DF_eval_final.set_index('chr_start').loc[chr_start,:].sort_values('plot_type', ascending=False).reset_index() for chr_start in chr_start_list])

# the not scored plot are included aswell:
# indexing with the complement of the CMP type set in the DF_eval_final, out of
# that, take the best p_value:
complement_dict = {'DOWN': 'UP', 'UP': 'DOWN'}
# for every chr_start found, concat the according validation plot:
DF_eval_final = pd.concat([pd.concat([DF_eval_final.set_index('chr_start').loc[chr_start, :].reset_index(), DF_eval.loc[(chr_start, complement_dict[DF_eval_final.set_index('chr_start').loc[chr_start, 'CMP'].drop_duplicates().values[0]] + '_validation', slice(None)),:].reset_index().sort_values('p_value').iloc[0,:].to_frame().T]) for chr_start in chr_start_list])


### the succession shall be base_plot, scored plot, not_scored_plot
##DF_eval_final = DF_eval_final.reset_index()
### get the base plot index:
##base_index = DF_eval_final.set_index('plot_type', append=True).loc[(slice(None),'base_plot'), :].reset_index('plot_type').index
### get everything not base plot (val plot) but scored true:
##scored_val_index = DF_eval_final[ DF_eval_final['plot_type'] != 'base_plot' ].set_index('scored', append=True).loc[(slice(None), True), :].reset_index('scored').index
### get the rest
##not_scored_val_index = DF_eval_final[ DF_eval_final['plot_type'] != 'base_plot' ].set_index('scored', append=True).loc[(slice(None), False), :].reset_index('scored').index
##final_index = []
##for i,j,k in zip(base_index, scored_val_index, not_scored_val_index):
##    final_index.extend([i,j,k])


print(f'saving table {metilene_lifeline_eval}')
DF_eval_final.to_csv(metilene_lifeline_eval, sep='\t', index=None)

pdfs_to_merge = [re.sub('tsv.*', 'pdf', i) for i in  DF_eval_final['file_path'].values.tolist()]

print(f'saving {metilene_lifeline_eval_pdfs}')

# do not open to many pdfs at once:
# OSError: [Errno 24] Too many open files: '/scr/palinca/gabor/TCGA-pipeline/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_5/metilene_intersect_lifeline_plot_chr22_44172456_44172932.pdf'
# (Pdb) len(pdfs_to_merge)
# 1095
# after the slicing the sum must be 1095 again!
slice_length = 100
mod = len(pdfs_to_merge) % slice_length
pdfs_to_merge_lists = []
start_ind = 0
slices = int(len(pdfs_to_merge) / slice_length)
if len(pdfs_to_merge) > slice_length:
    # divide the pdfs to merge in slice_lengther steps, get the rest:
    for i in range(slices):
        pdfs_to_merge_lists.append(pdfs_to_merge[start_ind: start_ind+slice_length])
        start_ind += slice_length
    # add the odd rest if something is left:
    if mod != 0:
        pdfs_to_merge_lists.append(pdfs_to_merge[slice_length*slices:])


# with pdf merger just the pdfs_to_merge list is needed:
# loop through the slices, prepend already merged pdfs before the single pdfs
if len(pdfs_to_merge) > slice_length:
    # create all the temporary pdfs slice merge filenames:
    temp_merge_name = [metilene_lifeline_eval_pdfs.replace('.pdf', f'_temp_{int(i)}.pdf') for i in range(slices+1)]
    first_pdf = True
    for i in range(len(pdfs_to_merge_lists)):
        merger = PdfMerger()
        first_slice = True
        for pdf in pdfs_to_merge_lists[i]:
            # before any slice prepend already merged pdfs, unless the first merge
            if first_slice and not first_pdf:
                first_slice = False
                first_pdf = False
                #
                merger.append(temp_merge_name[i-1])
                merger.append(pdf)
            else:
                first_pdf = False
                first_slice = False
                merger.append(pdf)
        if i != max(range(len(pdfs_to_merge_lists))):
            merger.write(temp_merge_name[i])
            merger.close()
        else:
            merger.write(metilene_lifeline_eval_pdfs)
            merger.close()
else:
    merger = PdfMerger()
    for pdf in pdfs_to_merge:
        merger.append(pdf)
    merger.write(metilene_lifeline_eval_pdfs)
    merger.close()
