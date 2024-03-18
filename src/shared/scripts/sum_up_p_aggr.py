import os
# from matplotlib.patches import bbox_artist
import matplotlib.transforms as transforms
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import sys
import re
from scipy.stats import gmean

"""
visualize the impact of the threshold parameter, based on the variances of
p_value_life and the life_mean_diff
"""

if 'snakemake' in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], 'w')

    print('# snakemake inputs:')
    [ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

    print('# snakemake output:')
    [ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

    print('# snakemake wildcards:')
    [ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

    plot_diffs_table = snakemake.input.plot_diffs_table
    p_val_prod_sum = snakemake.output.p_val_prod_sum[0]
    p_val_prod_sum_table = snakemake.output.p_val_prod_sum_table[0]
    # project = snakemake.wildcards.project
    # projects = ', '.join(snakemake.wildcards.project.split('_'))
    # gender = ', '.join(snakemake.wildcards.gender.split('_'))
    # gender = snakemake.wildcards.gender
    # drugs = '; '.join(snakemake.wildcards.drug_combi.split('_'))
    drug_combi = snakemake.wildcards.drug_combi
    pipeline = snakemake.wildcards.pipeline
    threshold_str = snakemake.wildcards.threshold_str
    count_type = snakemake.wildcards.count_type
    # plot_type = snakemake.wildcards.plot_type
else:
    ###########################################################################
    #                                test data                                #
    ###########################################################################
    # snakemake inputs:
    plot_diffs_table = "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz"
    # snakemake output:
    p_val_prod_sum = "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0-cutoff_5-cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene-beta_vals_p_prod_sum.pdf"
    p_val_prod_sum_table = "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0-cutoff_5-cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene-beta_vals_p_prod_sum.tsv.gz"
    # snakemake wildcards:
    output_path = "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod"
    pipeline = "metilene"
    drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
    threshold_str = "threshold_0_threshold_5_threshold_10_threshold_20"
    count_type = "beta_vals"
    # # snakemake inputs:
    # plot_diffs_table = "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz", "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs-norm_count_p_aggr.tsv.gz"
    # # snakemake output:
    # p_val_prod_sum = "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0-cutoff_5-cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2-norm_count_p_prod_sum.pdf"
    # p_val_prod_sum_table = "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0-cutoff_5-cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2-norm_count_p_prod_sum.tsv.gz"
    # # snakemake wildcards:
    # output_path = "/scr/palinca/gabor/TCGA-pipeline_10_pval_prod"
    # pipeline = "DESeq2"
    # drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
    # threshold_str = "threshold_0_threshold_5_threshold_10_threshold_20"
    # count_type = "norm_count"
    ###########################################################################
    #                                test data                                #
    ###########################################################################

"""
same script for different pipelines, parse out the pipeline name
"""
# those are the tables which are invoked into the aggregated DF
file_list = [i for i in plot_diffs_table if not pd.read_table(i).empty]
if len(file_list) == 0:
    open(p_val_prod_sum,'w').close()
    pd.DataFrame().to_csv(p_val_prod_sum_table, sep='\t')
    os._exit(0)

DF_list = []

# additionally calculate the geometric mean
for table in file_list:
    """
    same p_val_prod_sum for all 3 val plots, take the base_plot rows, throw away cols which are plot type specific, kepp all which are shared over the 3 plot_types:
    # the upper bound shall be  0.0025 -> just include those p_prod in the
    # ranking an sum those up
    (Pdb) 0.05*0.05
    0.0025000000000000005
    """
    trait = ''
    if pipeline == 'DESeq2':
        trait = 'ENSG'
        temp_DF = pd.read_table(table).set_index('plot_type').loc['base_plot', [trait, 'count_type', 'CMP', 'file_path', 'log2FoldChange', 'lfcSE', 'stat', 'p_value_deseq', 'padj', 'p_val_life_prod_scored', 'p_val_life_sum_scored']].reset_index()
    else:
        trait = 'chromosome - start'
        temp_DF = pd.read_table(table).set_index('plot_type', drop=True).loc['base_plot', ['ENSG', trait, 'DMR', 'CMP', 'file_path', 'mean_methylation_difference', 'q-value', '#CpGs', 'mean_alive', 'mean_dead', 'p_val_life_prod_scored', 'p_val_life_sum_scored']].reset_index()
    # boolean of values which shall be included in ranking
    temp_DF['included_in_prod_sum'] =  temp_DF['p_val_life_prod_scored'].apply(lambda x: True if x<0.0025 else False)
    temp_DF['dist_to_0.0025'] = 0.0025 - temp_DF.loc[temp_DF['included_in_prod_sum'],'p_val_life_prod_scored']
    temp_DF['prod_sum'] = temp_DF.loc[temp_DF['included_in_prod_sum'], 'dist_to_0.0025'].sum()
    temp_DF['geometric_mean'] = gmean(temp_DF.loc[temp_DF['included_in_prod_sum'], 'dist_to_0.0025'])
    temp_DF['nr_included_in_prod_sum'] = len(temp_DF[temp_DF['included_in_prod_sum']])
    temp_cutoff = table.split(os.path.sep)[-3]
    temp_sex = table.split(os.path.sep)[-4]
    temp_pipeline = table.split(os.path.sep)[-7]
    temp_project = table.split(os.path.sep)[-8]
    temp_DF['pipeline'] = temp_pipeline
    temp_DF['project'] = temp_project
    temp_DF['cutoff'] = temp_cutoff
    temp_DF['sex'] = temp_sex
    temp_DF['input_file'] = table
    DF_list.append(temp_DF)

DF = pd.concat(DF_list)
# DF.sort_values('prod_sum', ascending=False, inplace=True)
DF.sort_values('nr_included_in_prod_sum', ascending=False, inplace=True)
DF.to_csv(p_val_prod_sum_table, sep='\t', index=None)
DF= DF[DF['included_in_prod_sum']]
DF = DF.reset_index(drop=True)
DF['project + sex + cutoff'] = DF['project'] + ' + ' + DF['sex'] + ' + ' + DF['cutoff']
DF['project + sex'] = DF['project'] + ' + ' + DF['sex']
# sort str len wise:
sort_len_index = DF['project + sex'].str.len().sort_values(ascending=False).index
DF_temp = DF.loc[sort_len_index,:]
# now sort cutoff block wise, (project + sex)
# block_index = DF.loc[:, 'project + sex'].drop_duplicates().values
# DF_temp = pd.concat([DF.set_index('project + sex', drop=False).loc[i,:].sort_values('cutoff').reset_index(drop=True) for i in block_index])
# get the hue sorting:
# hue_order = DF_temp['project + sex + cutoff'].drop_duplicates().values
hue_order = DF_temp['project + sex'].drop_duplicates().values

fig, ax = plt.subplots(1,2,sharey=True, figsize=(11,5))
# palette=reversed(sns.cubehelix_palette()),

sns.barplot(DF, x='nr_included_in_prod_sum', y='project + sex + cutoff', ax=ax[0], hue='project + sex', hue_order=hue_order, legend=False, palette=reversed(sns.cubehelix_palette()))
ax[0].set_xlabel('number of\ngenes')
ax[0].xaxis.grid()
for i in ax[0].containers:
    ax[0].bar_label(i,)
g = sns.barplot(DF, x='geometric_mean', y='project + sex + cutoff', ax=ax[1], hue='project + sex', legend=True, hue_order=hue_order, palette=reversed(sns.cubehelix_palette()))
# sns.move_legend(g, 'lower center', bbox_to_anchor=(.55, .45))
sns.move_legend(g, loc=1, bbox_to_anchor=(5.14, 0.45))
ax[1].set_xlabel('geometric mean of\ndistances to 0.0025')
ax[1].xaxis.grid()
# for i in ax[1].containers:
    # ax[1].bar_label(i,)
ax[0].set_ylabel('Set combinations')
ax[1].set_ylabel('')
plt.suptitle(f'{pipeline}: Summary of the p-value validation products for all sets\n')
plt.subplots_adjust(left=0.43, right=0.63, top=0.9, bottom=0.15)
print(f'saving plot:{p_val_prod_sum}')
plt.savefig(p_val_prod_sum)
DF_temp = DF.drop_duplicates('project + sex + cutoff').loc[:, ['project + sex + cutoff', 'nr_included_in_prod_sum', 'prod_sum', 'cutoff', 'project', 'sex']]
DF_temp['str_len'] = DF_temp['project + sex + cutoff'].str.len()
DF_temp = DF_temp.sort_values(['str_len', 'nr_included_in_prod_sum'], ascending=False)
DF_temp = DF_temp.loc[:, ['project', 'sex', 'cutoff', 'nr_included_in_prod_sum', 'prod_sum']]
DF_temp['prod_sum'] = DF_temp['prod_sum'].round(3)
DF_temp.to_latex(p_val_prod_sum.replace('pdf', 'tex'), header=['project', 'sex', 'cutoff', 'number of invoked genes', 'sum of validation products'], index=None)
plt.cla()
plt.clf()

