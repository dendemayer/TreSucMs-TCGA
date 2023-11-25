import pandas as pd
import glob
import PyPDF2
import os
import sys
import re

"""
merging here the diffplots in the combinde treshold str dirs
additionally adding here for metilene the chr-start based diffplots
(not in snakefile input to keep compatibility between metilene and deseq)
sorting like:
DESeq2:
aggregated:
    - base
    - UP
    - DOWN
evaluated:
    - base
    - UP
    - DOWN
metilene
aggregated:
    - base_DMR
    - base_chr-start
    - UP_DMR
    - UP_chr-start
    - DOWN_DMR
    - DOWN_chr-start
evaluated:
    - base_DMR
    - base_chr-start
    - UP_DMR
    - UP_chr-start
    - DOWN_DMR
    - DOWN_chr-start
"""

sys.stderr = sys.stdout = open(snakemake.log[0], 'w')

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

print('# snakemake params:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.params.items()]

plot_diffs_aggr = snakemake.input.plot_diffs_aggr
plot_diffs_eval = snakemake.input.plot_diffs_eval
aggr_plots = snakemake.output.aggr_plots
pipeline = snakemake.wildcards.pipeline

###############################################################################
#                             test_set                                        #
###############################################################################
# # snakemake inputs:
# plot_diffs_aggr = ["/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_diffs_UP_validation-beta_vals.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_diffs_DOWN_validation-beta_vals.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_diffs_base_plot-beta_vals.pdf"]
# plot_diffs_eval = ["/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs_UP_validation-beta_vals.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs_DOWN_validation-beta_vals.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs_base_plot-beta_vals.pdf"]
# # snakemake output:
# aggr_plots = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_aggr+eval_diffs_merged.pdf"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_7"
# project = "TCGA-CESC"
# pipeline = "metilene"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female_male"
# cutoff = "cutoff_5"
# threshold_str = "threshold_0_threshold_5_threshold_10_threshold_20"
# snakemake params:

# # snakemake inputs:
# plot_diffs_aggr = ["/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_diffs_UP_validation-beta_vals.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_diffs_DOWN_validation-beta_vals.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_diffs_base_plot-beta_vals.pdf"]
# plot_diffs_eval = ["/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs_UP_validation-beta_vals.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs_DOWN_validation-beta_vals.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs_base_plot-beta_vals.pdf"]
# # snakemake output:
# aggr_plots = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_aggr+eval_diffs_merged.pdf"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_7"
# project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
# pipeline = "metilene"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female_male"
# cutoff = "cutoff_8"
# threshold_str = "threshold_0_threshold_5_threshold_10_threshold_20"
###############################################################################
#                             test_set                                        #
###############################################################################

plot_diffs_aggr = plot_diffs_aggr + plot_diffs_eval
# add the metilene chr-start diff plots if they are present:
if pipeline == 'metilene':
    # adding the metilene specific chr-start diff plots
    chr_start_diffs = [i for i in glob.glob(f'{os.path.split(aggr_plots)[0]}/*chr-start.pdf') if os.path.getsize(i) > 1]
    plot_diffs_aggr = plot_diffs_aggr + chr_start_diffs

plot_diffs_aggr = [i for i in plot_diffs_aggr if os.path.getsize(i) > 1]
if len(plot_diffs_aggr) == 0:
    print(f'no plots to aggregate, writing empty {aggr_plots}')
    open(aggr_plots, 'a').close()
    os._exit(0)

# parsing out UP_validation, DOWN_validation, or base_plot:
DF = pd.DataFrame({'full_path': plot_diffs_aggr, 'type': ['_'.join(os.path.split(plot)[1].split('-')[0].split('_')[-2:]) for plot in plot_diffs_aggr]})

# adding the aggregated or evaluated info:
DF['aggregation'] = DF['full_path'].apply(lambda x: 'evaluated' if re.search('eval', os.path.split(x)[1]) else 'aggregated')

# metilene gets an extra col for chr-start or DMR:
if pipeline == 'metline':
    DF['region'] = DF['full_path'].apply(lambda x: 'chr-start' if re.search('chr-start', os.path.split(x)[1]) else 'DMR')
    # order it directly also based on the region type
    pdfs_to_merge = pd.concat([DF.set_index(['aggregation', 'type', 'region']).loc[('aggregated', ['base_plot', 'UP_validation', 'DOWN_validation']), :], DF.set_index(['aggregation', 'type', 'region']).loc[('evaluated', ['base_plot', 'UP_validation', 'DOWN_validation']), :]])['full_path'].values.tolist()
else:
    #ordering without region type for deseq:
   pdfs_to_merge = DF.set_index(['aggregation', 'type']).loc[(['aggregated', 'evaluated'], ['base_plot', 'UP_validation', 'DOWN_validation'], ), :]['full_path'].tolist()


# ordering, metilene:
# aggregation type            region
# aggregated  base_plot       DMR        /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...
#                             chr-start  /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...
#             UP_validation   DMR        /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...
#                             chr-start  /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...
#             DOWN_validation DMR        /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...
#                             chr-start  /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T... evaluated   base_plot       DMR        /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...
#                             chr-start  /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...
#             UP_validation   DMR        /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...
#                             chr-start  /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...
#             DOWN_validation DMR        /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...
#                             chr-start  /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...

# ordering deseq:
# aggregation type
# aggregated  base_plot        /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...
#             UP_validation    /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...
#             DOWN_validation  /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...
# evaluated   base_plot        /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...
#             UP_validation    /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...
#             DOWN_validation  /scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_T...

print(f'writing pdfs to merge:\n{pdfs_to_merge}\nto:\n{aggr_plots}')
merger = PyPDF2.PdfMerger()
[merger.append(pdf) for pdf in pdfs_to_merge]
merger.write(aggr_plots)
merger.close()
