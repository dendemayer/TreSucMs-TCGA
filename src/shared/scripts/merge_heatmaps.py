import pandas as pd
import PyPDF2
import os
import sys
from natsort import natsort_keygen
import re

sys.stderr = sys.stdout = open(snakemake.log[0], 'w')

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

print('# snakemake params:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.params.items()]

heatmaps = snakemake.input.heatmaps
heatmaps_merged = snakemake.output.heatmaps_merged
pipeline = snakemake.wildcards.pipeline
###############################################################################
#                             test_set                                        #
###############################################################################

# heatmaps = ["/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr12_16605757_16606859.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr13_36431314_36432474.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr13_36431314_36432474.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr12_16605757_16606859.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr19_57769432_57770082.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr13_36431314_36432474.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr6_3848399_3850498.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chrX_16711802_16712814.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr8_140238482_140238568.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr17_4969192_4969308.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr7_95396094_95397625.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr14_85529381_85531409.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr11_115503038_115505464.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr7_121328978_121329025.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chrX_72305564_72308502.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr11_115503038_115505464.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr14_85529381_85531409.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr11_115503038_115505464.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr12_16605757_16606859.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr11_130121549_130122473.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chrX_16711802_16712814.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr16_51150443_51153014.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr7_121328978_121329025.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr8_140238482_140238568.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr19_57769432_57770082.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr8_140238482_140238568.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr17_4969192_4969308.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr11_130121549_130122473.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chrX_72305564_72308502.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr16_51150443_51153014.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr7_121328978_121329025.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr7_83648624_83648822.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr11_115503038_115505464.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr11_130121549_130122473.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr16_51150443_51153014.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr11_6926309_6926793.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr19_57769432_57770082.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr7_95396094_95397625.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr7_95396094_95397625.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr14_85529381_85531409.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr6_26183279_26184017.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr19_57769432_57770082.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr16_54932558_54939135.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr20_58850459_58852587.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr11_6926309_6926793.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr6_26183279_26184017.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr6_26183279_26184017.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr17_4969192_4969308.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr8_140238482_140238568.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr11_130121549_130122473.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr11_6926309_6926793.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr14_85529381_85531409.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chrX_16711802_16712814.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr11_6926309_6926793.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chrX_72305564_72308502.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr13_102399793_102403137.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr13_36431314_36432474.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr6_26183279_26184017.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr16_51150443_51153014.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr13_102399793_102403137.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr13_102399793_102403137.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr20_58850459_58852587.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr16_54932558_54939135.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr7_83648624_83648822.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr20_58850459_58852587.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr12_16605757_16606859.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr6_3848399_3850498.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chrX_72305564_72308502.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr7_95396094_95397625.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr6_3848399_3850498.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr7_83648624_83648822.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr7_83648624_83648822.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr6_3848399_3850498.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chrX_16711802_16712814.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr20_58850459_58852587.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr17_4969192_4969308.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_lineplot_median_beta_value_chr16_54932558_54939135.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_violinplot_beta_value_chr16_54932558_54939135.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr7_121328978_121329025.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_heatmaps_beta_value_chr13_102399793_102403137.pdf"]
# heatmaps_merged = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_beta_value_merged.pdf"
# output_path = "/scr/palinca/gabor/TCGA-pipeline_7"
# project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
# pipeline = "metilene"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_8"
# # snakemake inputs:
# heatmaps = ["/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/DESeq2_heatmap_log2fINCREASE_nt_counts.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DEn,paclitaxel_cisplatin/male/cutoff_5/DESeq2_heatmap_log2fDECREASE_nt_counts.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/DESeq2_heatmap_log2fINCREASE_TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/DESeq2_heatmap_log2fDECREASE_raw_counts.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboptoff_5/DESeq2_heatmap_log2fINCREASE_norm_counts.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/DESeq2_heatmap_log2fDECREASE_norm_counts.pdf", "/scr/palincaESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/DESeq2_heatmap_log2fINCREASE_vsd_counts.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/DECREASE_vsd_counts.pdf"]
# # snakemake output:
# heatmaps_merged = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/DESeq2_heatmap_merged.pdf"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_7"
# project = "TCGA-HNSC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "male"
# cutoff = "cutoff_5"
###############################################################################
#                             test_set                                        #
###############################################################################

heatmaps = [i for i in heatmaps if os.path.getsize(i) > 1]
if len(heatmaps) == 0:
    open(heatmaps_merged, 'a').close()
    os._exit(0)
# in case of metilene, wie want to sort the 4 different plots types
# (boxplot|heatmaps|lineplot|violinplot)  DMR wise:

if pipeline == 'metilene':
    # box= '/scr/palinca/gabor/TCGA-pipeline_7/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/metilene_intersect_boxplot_beta_value_chr6_3848399_3850498.pdf'
    # parsing out DMR of full filename:
    # DMR = re.search(r'chr.*\.pdf', os.path.split(box)[1]).group(0).strip('.pdf')
    # parsing out the plot_type
    # plot_type = re.search(r'int.*plot', os.path.split(box)[1]).group(0).split('_')[1]
    DF = pd.concat([pd.Series([heatmap, re.search(r'chr.*\.pdf$', os.path.split(heatmap)[1]).group(0).strip('.pdf'), re.search(r'intersect.*_beta', os.path.split(heatmap)[1]).group(0).replace('intersect_', '').replace('_beta', '')], index=['full_path', 'DMR', 'type']) for heatmap in heatmaps], axis=1).T.sort_values(['DMR', 'type']).head(50)
    heatmaps = DF.sort_values(['DMR', 'type'], key=natsort_keygen()).loc[:, 'full_path'].values.tolist()


merger = PyPDF2.PdfMerger()
[merger.append(pdf) for pdf in heatmaps]
merger.write(heatmaps_merged)
merger.close()