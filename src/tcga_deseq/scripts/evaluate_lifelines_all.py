import os
import pandas as pd
from PyPDF2 import PdfMerger
import re
import sys

###############################################################################
#                              snakemake inputs                               #
###############################################################################

if 'snakemake' in dir():
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
    deseq_lifeline_eval_combined_pvals = snakemake.output[2]
    count_type = snakemake.wildcards.count_type
###############################################################################
###############################################################################
#                                  test set                                   #
###############################################################################
else:
    # snakemake inputs:
    deseq_lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin_cisplatin,gemcitabine/male/cutoff_0/threshold_0/DESeq2_lifelines_aggregated.tsv.gz"
    script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../tcga_deseq/scripts/evaluate_lifelines_all.py"
    # snakemake output:
    deseq_lifeline_eval = "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin_cisplatin,gemcitabine/male/cutoff_0/threshold_0/DESeq2_lifelines_evaluated-norm_count.tsv.gz"
    deseq_lifeline_eval_pdfs = "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin_cisplatin,gemcitabine/male/cutoff_0/threshold_0/DESeq2_lifelines_evaluated-norm_count.pdf"
    deseq_lifeline_eval_combined_pvals = "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin_cisplatin,gemcitabine/male/cutoff_0/threshold_0/DESeq2_lifelines_evaluated-norm_count_pval_combined.tsv.gz"
    # snakemake wildcards:
    output_path = "/scr/palinca/gabor/TCGA-pipeline_8"
    project = "TCGA-HNSC"
    drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin_cisplatin,gemcitabine"
    gender = "male"
    cutoff = "cutoff_0"
    count_type = "norm_count"
    # # snakemake inputs:
    # deseq_lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_aggregated.tsv.gz"
    # script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../tcga_deseq/scripts/evaluate_lifelines_all.py"
    # # snakemake output:
    # deseq_lifeline_eval = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count_test.tsv.gz"
    # deseq_lifeline_eval_pdfs = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count_test.pdf"
    # deseq_lifeline_eval_combined_pvals = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count_pval_combined_test.tsv.gz"
    # # snakemake wildcards:
    # output_path = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod"
    # project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
    # drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
    # gender = "female_male"
    # cutoff = "cutoff_0"
    # count_type = "norm_count"

DF_eval = pd.read_table(deseq_lifeline_aggregated)
if DF_eval.empty:
    DF_eval.to_csv(deseq_lifeline_eval, sep='\t')
    DF_eval.to_csv(deseq_lifeline_eval_combined_pvals, sep='\t')
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

################################################################################
#                           calculation all pval combinations                  #
################################################################################
DF_eval_base = DF_eval_cross[DF_eval_cross['plot_type'] == 'base_plot']
DF_eval_val = DF_eval_cross[DF_eval_cross['plot_type'] != 'base_plot']
# ensg = ENSGs[0]
DF_eval_temp = pd.DataFrame()
for ensg in ENSGs:
    try:
        # renaming also the p_val respectively to p_value_life_base and
        # p_value_life_val:
        temp_DF = DF_eval_base.set_index('ENSG').loc[ensg, ['threshold', 'p_value_life']].merge(DF_eval_val.set_index('ENSG').loc[ensg, ['threshold', 'p_value_life']], how='cross', suffixes=('_base', '_val'))
    except AttributeError:
        # it occurs that just one ENSG is available, then take the two p_vals
        # available, one for base and one for val:
        # print(f'ensg: {ensg}')
        base_temp = DF_eval_base.set_index('ENSG').loc[ensg, ['threshold', 'p_value_life']]
        val_temp = DF_eval_val.set_index('ENSG').loc[ensg, ['threshold', 'p_value_life']]
        # the renaming of the p_value must be done for the Series as well:
        val_temp.rename({'p_value_life': 'p_value_life_val', 'threshold': 'threshold_val'}, inplace=True)
        base_temp.rename({'p_value_life': 'p_value_life_base', 'threshold': 'threshold_base'}, inplace=True)
        temp_DF = pd.concat([base_temp.to_frame().T, val_temp.to_frame().T], axis=1)
        temp_DF = temp_DF.reset_index(drop=True)
    temp_DF['ENSG']= ensg
    temp_DF['p_value_life_sum'] = temp_DF['p_value_life_base'] + temp_DF['p_value_life_val']
    temp_DF['p_value_life_prod'] = temp_DF['p_value_life_base'] * temp_DF['p_value_life_val']
    tmp_cmp = DF_eval_val.set_index('ENSG').loc[ensg,:]['CMP']
    if isinstance(tmp_cmp, pd.Series):
        tmp_cmp = tmp_cmp.drop_duplicates().values[0]
    DF_eval_temp = pd.concat([DF_eval_temp, temp_DF])

print(f'saving table {deseq_lifeline_eval_combined_pvals}')
DF_eval_temp.to_csv(deseq_lifeline_eval_combined_pvals, sep='\t', index=None)
################################################################################
#                     calculating all pvalue combinations                      #
################################################################################


################################################################################
#                              old pval selection                              #
################################################################################
# # DF_eval_cross = DF_eval_cross.loc[:, ['ENSG', 'plot_type', 'threshold', 'p_value_life', 'scored']].merge(DF_eval.loc[:, ['ENSG', 'plot_type', 'threshold', 'p_value_life', 'scored']], how='cross').set_index(['scored_x', 'scored_y']).loc[(True, True), :]
# # sum up all possible p_value sums
# # DF_eval_cross['p_sum_cross'] = DF_eval_cross['p_value_life_x'] + DF_eval_cross['p_value_life_y']
# # DF_eval_cross['p_prod_cross'] = DF_eval_cross['p_value_life_x'] * DF_eval_cross['p_value_life_y']
# # limit on the base_plots
# DF_base_temp = DF_eval_cross.set_index('plot_type_x').loc['base_plot'].sort_values('p_sum_cross')
# # just take the minimum for base and val plot of same ENSGs: indexing
# # ENSG_x and ENSG_y, concatenating the final, best p_sum_cross:
# DF_base_temp_cross = pd.concat([ DF_base_temp[DF_base_temp['plot_type_y'] != 'base_plot'].reset_index().set_index(['ENSG_x', 'ENSG_y']).loc[(ensg, ensg), :].sort_values('p_sum_cross').iloc[0, :].to_frame() for ensg in ENSGs ], axis=1).T
# # index the start DF (DF_eval_cross) directly on the gathered combinations:
# DF_base_plot = DF_base_temp_cross.reset_index().loc[:, ['level_1', 'plot_type_x', 'threshold_x', 'p_sum_cross']]
# DF_base_plot.columns =  ['ENSG', 'plot_type', 'threshold', 'p_sum_cross']
# DF_val_plot = DF_base_temp_cross.reset_index().loc[:, ['level_1', 'plot_type_y', 'threshold_y', 'p_sum_cross']]
# DF_val_plot.columns =   ['ENSG', 'plot_type', 'threshold', 'p_sum_cross']

# #DF_base_val = pd.concat([DF_base_plot, DF_val_plot]).sort_values(['ENSG', 'plot_type'], ascending=False).set_index(['ENSG', 'plot_type', 'threshold']).reset_index()
# DF_base_val = pd.concat([DF_base_plot, DF_val_plot]).sort_values(['ENSG', 'plot_type'], ascending=False).set_index(['ENSG', 'plot_type', 'threshold'])
# DF_eval = DF_eval.drop_duplicates().set_index(['ENSG', 'plot_type', 'threshold'])

# DF_eval_final = DF_base_val.join(DF_eval)
# DF_eval_final = DF_eval_final.sort_values('p_sum_cross').reset_index()
# # ordering each 2er postition pair to base_plot val_plot:
# ensg_list = DF_eval_final['ENSG'].drop_duplicates().to_list()
# DF_eval_final = pd.concat([ DF_eval_final.set_index('ENSG').loc[ensg,:].sort_values('plot_type', ascending=False).reset_index() for ensg in ensg_list])

# # the not scored plots are included aswell:
# # indexing with the complement of the CMP type set in the DF_eval_final, out of
# # that, take the best p_value:
# complement_dict = {'DOWN': 'UP', 'UP': 'DOWN'}
# # for every chr_start found, concat the according validation plot:
# DF_eval_final = DF_eval_final.drop_duplicates()
# DF_eval_final = pd.concat([pd.concat([DF_eval_final.set_index('ENSG').loc[ensg, :].reset_index(), DF_eval.loc[(ensg, complement_dict[DF_eval_final.set_index('ENSG').loc[ensg, 'CMP'].drop_duplicates().values[0]] + '_validation', slice(None)),:].reset_index().sort_values('p_value_life').iloc[0,:].to_frame().T]) for ensg in ensg_list])
################################################################################
#                              old pval selection                              #
################################################################################

complement_dict = {'DOWN': 'UP', 'UP': 'DOWN'}
DF_final_final = pd.DataFrame()
for ensg in ENSGs:
    cmp = DF_eval_cross.set_index('ENSG').loc[ensg,:]['CMP'].drop_duplicates().values[0]
    if isinstance(DF_eval_cross.set_index(['ENSG', 'plot_type']).loc[ensg, 'base_plot'], pd.DataFrame):
        try:
            base_DF = DF_eval_cross.set_index(['ENSG', 'plot_type']).loc[ensg, 'base_plot'].sort_values('p_value_life').iloc[0,:].to_frame().T
            val_DF = DF_eval_cross.set_index(['ENSG', 'plot_type']).loc[ensg, f'{cmp}_validation'].sort_values('p_value_life').iloc[0,:].to_frame().T
            # DF_eval_cross.set_index(['ENSG', 'plot_type']).loc[ensg, :]
        except Exception as e:
            print(f'{ensg}, {e}')
        try:
            compl_DF = DF_eval.set_index(['ENSG', 'plot_type']).loc[ensg, f'{complement_dict[cmp]}_validation'].sort_values('p_value_life').iloc[0,:].to_frame().T
        except KeyError:
            compl_DF = pd.DataFrame()
        temp_DF = pd.concat([base_DF, val_DF, compl_DF])
        temp_DF['p_val_life_prod_scored'] = temp_DF[temp_DF['scored'].values].loc[:, 'p_value_life'].prod()
        temp_DF['p_val_life_sum_scored'] = temp_DF[temp_DF['scored'].values].loc[:, 'p_value_life'].sum()
        DF_final_final = pd.concat([DF_final_final, temp_DF])
    else:
        try:
            base_DF = DF_eval_cross.set_index(['ENSG', 'plot_type']).loc[ensg, 'base_plot']
            val_DF = DF_eval_cross.set_index(['ENSG', 'plot_type']).loc[ensg, f'{cmp}_validation']
        except Exception as e:
            print(f'{ensg}, {e}')
        try:
            compl_DF = DF_eval.set_index(['ENSG', 'plot_type']).loc[ensg, f'{complement_dict[cmp]}_validation']
        except KeyError:
            compl_DF = pd.DataFrame()
        temp_DF = pd.concat([base_DF.to_frame().T, val_DF.to_frame().T, compl_DF.to_frame().T])
        # print(temp_DF)
        temp_DF['p_val_life_prod_scored'] = temp_DF[temp_DF['scored'].values].loc[:, 'p_value_life'].prod()
        temp_DF['p_val_life_sum_scored'] = temp_DF[temp_DF['scored'].values].loc[:, 'p_value_life'].sum()
        DF_final_final = pd.concat([DF_final_final, temp_DF])
# old index:
# Index(['ENSG', 'plot_type', 'threshold', 'p_sum_cross', 'count_type', 'CMP', 'p_value_life', 'fst_life_mean', 'scnd_life_mean', 'life_mean_diff', 'file_path', 'p_sum', 'p_prod', 'scored', 'log2FoldChange', 'lfcSE', 'stat', 'p_value_deseq', 'padj']
DF_eval_final = DF_final_final
DF_eval_final = DF_eval_final.sort_values('p_val_life_prod_scored')


# this sorting can rearrange the index of plot types, loc that correctly again:
ENSGs = DF_eval_final.reset_index(1).index.drop_duplicates()

DF_final_final = pd.DataFrame()
for ensg in ENSGs:
    cmp = DF_eval_cross.set_index('ENSG').loc[ensg,:]['CMP'].drop_duplicates().values[0]
    if isinstance(DF_eval_cross.set_index(['ENSG', 'plot_type']).loc[ensg, 'base_plot'], pd.DataFrame):
        try:
            base_DF = DF_eval_cross.set_index(['ENSG', 'plot_type']).loc[ensg, 'base_plot'].sort_values('p_value_life').iloc[0,:].to_frame().T
            val_DF = DF_eval_cross.set_index(['ENSG', 'plot_type']).loc[ensg, f'{cmp}_validation'].sort_values('p_value_life').iloc[0,:].to_frame().T
            # DF_eval_cross.set_index(['ENSG', 'plot_type']).loc[ensg, :]
        except Exception as e:
            print(f'{ensg}, {e}')
        try:
            compl_DF = DF_eval.set_index(['ENSG', 'plot_type']).loc[ensg, f'{complement_dict[cmp]}_validation'].sort_values('p_value_life').iloc[0,:].to_frame().T
        except KeyError:
            compl_DF = pd.DataFrame()
        temp_DF = pd.concat([base_DF, val_DF, compl_DF])
        temp_DF['p_val_life_prod_scored'] = temp_DF[temp_DF['scored'].values].loc[:, 'p_value_life'].prod()
        temp_DF['p_val_life_sum_scored'] = temp_DF[temp_DF['scored'].values].loc[:, 'p_value_life'].sum()
        DF_final_final = pd.concat([DF_final_final, temp_DF])
    else:
        try:
            base_DF = DF_eval_cross.set_index(['ENSG', 'plot_type']).loc[ensg, 'base_plot']
            val_DF = DF_eval_cross.set_index(['ENSG', 'plot_type']).loc[ensg, f'{cmp}_validation']
            # DF_eval_cross.set_index(['ENSG', 'plot_type']).loc[ensg, :]
        except Exception as e:
            print(f'{ensg}, {e}')
        try:
            compl_DF = DF_eval.set_index(['ENSG', 'plot_type']).loc[ensg, f'{complement_dict[cmp]}_validation']
        except KeyError:
            compl_DF = pd.DataFrame()
        temp_DF = pd.concat([base_DF.to_frame().T, val_DF.to_frame().T, compl_DF.to_frame().T])
        temp_DF['p_val_life_prod_scored'] = temp_DF[temp_DF['scored'].values].loc[:, 'p_value_life'].prod()
        temp_DF['p_val_life_sum_scored'] = temp_DF[temp_DF['scored'].values].loc[:, 'p_value_life'].sum()
        DF_final_final = pd.concat([DF_final_final, temp_DF])

DF_eval_final = DF_final_final
# # base_plots:
# base_temp = DF_eval_final.set_index('scored', append=True).loc[(index_ensg, 'base_plot', True), :]
# base_temp = base_temp.reset_index([1,2]).rename({'level_1': 'plot_type'}, axis=1)
# # belonging val plot, eiterh up or down validated:
# val_temp = DF_eval_final[~(DF_eval_final.reset_index(1)['level_1'] == 'base_plot').values & (DF_eval_final['scored'] == True).values]
# val_temp = val_temp.reset_index(1).rename({'level_1': 'plot_type'}, axis=1)
# # the not scored plots:
# rest_temp = DF_eval_final[DF_eval_final['scored'] == False]
# rest_temp = rest_temp.reset_index(1).rename({'level_1': 'plot_type'}, axis=1)


# DF_temp = []
# for ensg in index_ensg:
#     try:
#         DF_temp.append(pd.concat([base_temp.loc[ensg,:], val_temp.loc[ensg, :], rest_temp.loc[ensg, :]], axis=1).T)
#     except Exception as e:
#         continue

# DF_eval_final = pd.concat(DF_temp)




# to keep the sorting of the base plot and val plot (scored) and val plot (not
# scored) get the pdfs_to_merge right here, before new ordering of DF:
pdfs_to_merge = [re.sub('tsv.*', 'pdf', i) for i in  DF_eval_final['file_path'].values.tolist()]
DF_eval_final = DF_eval_final.reset_index().rename({'level_0': 'ENSG', 'level_1': 'plot_type'}, axis=1)
DF_eval_final = DF_eval_final.sort_values('p_val_life_prod_scored')
print(f'saving table {deseq_lifeline_eval}')
DF_eval_final.to_csv(deseq_lifeline_eval, sep='\t', index=None)


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

