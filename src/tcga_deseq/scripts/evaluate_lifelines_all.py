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

deseq_lifeline_aggregated = snakemake.input[0]
deseq_lifeline_eval = snakemake.output[0] # deseq_lifeline_eval
deseq_lifeline_eval_pdfs = snakemake.output[1] # deseq_lifeline_eval_pdfs
count_type = snakemake.wildcards.count_type
###############################################################################
###############################################################################
#                                  test set                                   #
###############################################################################

# snakemake inputs:
# deseq_lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_aggregated.tsv.gz"
# script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/tcga_deseq/scripts/evaluate_lifelines_all.py"
# # snakemake output:
# deseq_lifeline_eval = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated_norm_count.tsv.gz"
# deseq_lifeline_eval_pdfs = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated_norm_count.pdf"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_3"
# project = "TCGA-HNSC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_0"
# count_type = "norm_count"

DF_eval = pd.read_table(deseq_lifeline_aggregated)
if DF_eval.empty:
    DF_eval.to_csv(deseq_lifeline_eval)
    open(deseq_lifeline_eval_pdfs, 'a').close()
    os._exit(0)

DF_eval = DF_eval[DF_eval['count_type'] == count_type]
DF_eval_cross = DF_eval[DF_eval['scored']]
# make a kartesian cross product between the base plots and the validation
# plots, -> inspecting all combinations of available p-values in dependence of
# the threshold invoked
# one evaluation step for each count type -> the DF can be limited here to the
# respective counttype requested through the wildcard
DF_eval_cross = DF_eval_cross.set_index('count_type').loc[count_type, :].reset_index()
ENSGs = DF_eval_cross.reset_index()['ENSG'].value_counts().index.to_list()
DF_eval_cross = DF_eval_cross.loc[:, ['ENSG', 'plot_type', 'threshold', 'p_value', 'scored']].merge(DF_eval.loc[:, ['ENSG', 'plot_type', 'threshold', 'p_value', 'scored']], how='cross').set_index(['scored_x', 'scored_y']).loc[(True, True), :]
# sum up all possible p_value sums
DF_eval_cross['p_sum_cross'] = DF_eval_cross['p_value_x'] + DF_eval_cross['p_value_y']
# limit on the base_plots
DF_base_temp = DF_eval_cross.set_index('plot_type_x').loc['base_plot'].sort_values('p_sum_cross')
# just take the minimum for base and val plot of same ENSGs: indexing
# ENSG_x and ENSG_y, concatenating the final, best p_sum_cross:
DF_base_temp_cross = pd.concat([ DF_base_temp[DF_base_temp['plot_type_y'] != 'base_plot'].reset_index().set_index(['ENSG_x', 'ENSG_y']).loc[(ensg, ensg), :].sort_values('p_sum_cross').iloc[0, :].to_frame() for ensg in ENSGs ], axis=1).T
# index the start DF (DF_eval_cross) directly on the gathered combinations:
DF_base_plot = DF_base_temp_cross.reset_index().loc[:, ['level_1', 'plot_type_x', 'threshold_x', 'p_sum_cross']]
DF_base_plot.columns =  ['ENSG', 'plot_type', 'threshold', 'p_sum_cross']
DF_val_plot = DF_base_temp_cross.reset_index().loc[:, ['level_1', 'plot_type_y', 'threshold_y', 'p_sum_cross']]
DF_val_plot.columns =   ['ENSG', 'plot_type', 'threshold', 'p_sum_cross']

#DF_base_val = pd.concat([DF_base_plot, DF_val_plot]).sort_values(['ENSG', 'plot_type'], ascending=False).set_index(['ENSG', 'plot_type', 'threshold']).reset_index()
DF_base_val = pd.concat([DF_base_plot, DF_val_plot]).sort_values(['ENSG', 'plot_type'], ascending=False).set_index(['ENSG', 'plot_type', 'threshold'])
DF_eval = DF_eval.drop_duplicates().set_index(['ENSG', 'plot_type', 'threshold'])

DF_eval_final = DF_base_val.join(DF_eval)
DF_eval_final = DF_eval_final.sort_values('p_sum_cross').reset_index()
# ordering each 2er postition pair to base_plot val_plot:
ensg_list = DF_eval_final['ENSG'].drop_duplicates().to_list()
DF_eval_final = pd.concat([ DF_eval_final.set_index('ENSG').loc[ensg,:].sort_values('plot_type', ascending=False).reset_index() for ensg in ensg_list])

# the not scored plot are included aswell:
# indexing with the complement of the CMP type set in the DF_eval_final, out of
# that, take the best p_value:
complement_dict = {'DOWN': 'UP', 'UP': 'DOWN'}
# for every chr_start found, concat the according validation plot:
DF_eval_final = DF_eval_final.drop_duplicates()
DF_eval_final = pd.concat([pd.concat([DF_eval_final.set_index('ENSG').loc[ensg, :].reset_index(), DF_eval.loc[(ensg, complement_dict[DF_eval_final.set_index('ENSG').loc[ensg, 'CMP'].drop_duplicates().values[0]] + '_validation', slice(None)),:].reset_index().sort_values('p_value').iloc[0,:].to_frame().T]) for ensg in ensg_list])

print(f'saving table {deseq_lifeline_eval}')
DF_eval_final.to_csv(deseq_lifeline_eval, sep='\t', index=None)

pdfs_to_merge = [re.sub('tsv.*', 'pdf', i) for i in  DF_eval_final['file_path'].values.tolist()]

print(f'saving {deseq_lifeline_eval_pdfs}')

# do not open to many pdfs at once:
# OSError: [Errno 24] Too many open files: '/scr/palinca/gabor/TCGA-pipeline/TCGA-HNSC/deseq/deseq_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_5/deseq_intersect_lifeline_plot_chr22_44172456_44172932.pdf'
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
    temp_merge_name = [deseq_lifeline_eval_pdfs.replace('.pdf', f'_temp_{int(i)}.pdf') for i in range(slices+1)]
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
            merger.write(deseq_lifeline_eval_pdfs)
            merger.close()
else:
    merger = PdfMerger()
    for pdf in pdfs_to_merge:
        merger.append(pdf)
    merger.write(deseq_lifeline_eval_pdfs)
    merger.close()

