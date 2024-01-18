import os
# from matplotlib.patches import bbox_artist
import matplotlib.transforms as transforms
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import sys
import re

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

    lifelines_evaluated = snakemake.input.lifelines_evaluated
    plot_diffs = snakemake.output.plot_diffs
    plot_diffs_table = snakemake.output.plot_diffs_table
    project = snakemake.wildcards.project
    # projects = ', '.join(snakemake.wildcards.project.split('_'))
    # gender = ', '.join(snakemake.wildcards.gender.split('_'))
    gender = snakemake.wildcards.gender
    # drugs = '; '.join(snakemake.wildcards.drug_combi.split('_'))
    drug_combi = snakemake.wildcards.drug_combi
    pipeline = snakemake.wildcards.pipeline
    # plot_type = snakemake.wildcards.plot_type
    # cutoff = snakemake.wildcards.cutoff.split('_')[1]
    cutoff = snakemake.wildcards.cutoff
    count_type = snakemake.wildcards.count_type
else:
    # snakemake inputs:
    lifelines_evaluated = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz"
    script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../shared/scripts/plot_thresh_diff_p_aggr.py"
    # snakemake output:
    plot_diffs = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.pdf"
    plot_diffs_table = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs-beta_vals_p_aggr.tsv.gz"
    # snakemake wildcards:
    output_path = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod"
    project = "TCGA-HNSC"
    pipeline = "metilene"
    drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
    gender = "female"
    cutoff = "cutoff_5"
    threshold_str = "threshold_0_threshold_5_threshold_10_threshold_20"
    count_type = "beta_vals"
    # # snakemake inputs:
    # lifelines_evaluated = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz"
    # script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../shared/scripts/plot_thresh_diff_p_aggr.py"
    # # snakemake output:
    # plot_diffs = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs_base_plot-beta_vals_p_aggr.pdf"
    # plot_diffs_table = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs_base_plot-beta_vals_p_aggr.tsv.gz"
    # # snakemake wildcards:
    # output_path = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod"
    # project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
    # pipeline = "metilene"
    # drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
    # gender = "female_male"
    # cutoff = "cutoff_0"
    # threshold_str = "threshold_0_threshold_5_threshold_10_threshold_20"
    # count_type = "beta_vals"
    # # snakemake inputs:
    # lifelines_evaluated = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_lifelines_evaluated-norm_count.tsv.gz"
    # script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../shared/scripts/plot_thresh_diff.py"
    # # snakemake output:
    # plot_diffs = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs_base_plot-norm_count_p_aggr.pdf"
    # plot_diffs_table = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/DESeq2_plot_eval_diffs_base_plot-norm_count_p_aggr.tsv.gz"
    # # snakemake wildcards:
    # output_path = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod"
    # project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
    # pipeline = "DESeq2"
    # drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
    # gender = "female"
    # cutoff = "cutoff_0"
    # threshold_str = "threshold_0_threshold_5_threshold_10_threshold_20"
    # count_type = "norm_count"

projects = ', '.join(project.split('_'))
gender = ', '.join(gender.split('_'))
drugs = '; '.join(drug_combi.split('_'))
cutoff = cutoff.split('_')[1]
###############################################################################
#                                 test input                                  #
###############################################################################

###############################################################################
#                    input to estimate the plot dimansions                    #
###############################################################################
# # (Pdb) len(base_DF_pval_sort['chr_start'].drop_duplicates())
# # 443
# # snakemake inputs:
# lifelines_evaluated = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_aggregated.tsv.gz"
# script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/tcga_metilene/../shared/scripts/plot_thresh_diff.py"
# # snakemake output:
# plot_diffs = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_diffs_base_plot-beta_vals.pdf"
# plot_diffs_table = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_diffs_base_plot-beta_vals.tsv.gz"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_3"
# project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
# pipeline = "metilene"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female_male"
# cutoff = "cutoff_5"
# plot_type = "base_plot"
# count_type = "beta_vals"
# plot_type = "base_plot"
# projects = project
# drugs = drug_combi

# plot with few 11  datapoints:
# good fig size:

# (Pdb) len(base_DF_pval_sort['chr_start'].drop_duplicates())
# 21
#
# # snakemake inputs:
# lifelines_evaluated = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_aggregated.tsv.gz"
# script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/tcga_metilene/../shared/scripts/plot_thresh_diff.py"
# # snakemake output:
# plot_diffs = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_diffs_base_plot.pdf"
# plot_diffs_table = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_diffs_base_plot.tsv.gz"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_3"
# project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
# pipeline = "metilene"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female"
# cutoff = "cutoff_0"
# plot_type = "base_plot"
# projects = project
# drugs = drug_combi

# (Pdb) num_cols
# 442
# # figsize = ( num_cols * 0.05, 8)
# lifelines_evaluated = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_aggregated.tsv.gz"
# script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/tcga_metilene/../shared/scripts/plot_thresh_diff.py"
# # snakemake output:
# plot_diffs = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_diffs_base_plot.pdf"
# plot_diffs_table = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_diffs_base_plot.tsv.gz"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_3"
# project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
# pipeline = "metilene"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female_male"
# cutoff = "cutoff_5"
# plot_type = "base_plot"
# projects = project
# drugs = drug_combi

#bei 21 -> 10.5,  bei 442 -> 22.1
# 435 -> 50
###############################################################################
#                    input to estimate the plot dimensions                    #
###############################################################################



DF = pd.read_table(lifelines_evaluated)
if DF.empty:
    open(plot_diffs, 'w').close()
    DF.to_csv(plot_diffs_table, sep='\t')
    # open(metilene_lifelines_evaluated, 'a').close()
    os._exit(0)
DF = DF.set_index('plot_type')
DF.drop_duplicates('file_path', inplace=True)

def plot_dynamical_metilene(DF, out_file, out_file_table):
    """
    using here the chromosome - start as feature for dissemination
    """
    base_DF = DF
    base_DF['chromosome - start'] = base_DF['chr_start'].apply(lambda x: x.split('_')[0] + ' - ' + x.split('_')[1])
    # p_val_sorted_DMR = base_DF.set_index('threshold').loc[0, :].sort_values(what_to_plot)['DMR'].to_list()
    # try:
    #     if re.search('lifelines_evaluated',os.path.split(lifelines_evaluated)[1]):
    #         # in the evaluated case, the best p_value_life is already given
    #         p_val_sorted_DMR = base_DF.sort_values('p_value_life')['DMR'].to_list()
    #     else:
    #         p_val_sorted_DMR = base_DF.set_index('threshold').loc[0, :].sort_values('p_value_life')['DMR'].to_list()
    # except Exception as e:
    #     print(f'exception {e}\nnot enough values for plot, writing empty files {out_file} and {out_file_table}')
    #     pd.DataFrame().to_csv(out_file_table, sep='\t')
    #     open(out_file, 'w').close()
    #     os._exit(0)
    # base_DF_pval_sort = base_DF.reset_index().set_index('DMR').loc[p_val_sorted_DMR, :]
    base_DF_pval_sort = base_DF.reset_index().sort_values('p_val_life_prod_scored')
    num_cols = len(base_DF_pval_sort.reset_index()['chr_start'].drop_duplicates())
    # figsize = (num_cols * 0.8, num_rows * 0.02)
    # figsize in width , height!
    # figsize = ( num_cols * 0.05, 8)
    # taking f(x) = 0.028*x+9.921 for the width (which is in that plot the
    # uniqe rows of feature
    # figsize = (0.028*num_cols+9.921, 8)
    figsize = (0.09541062801932366*num_cols + 8.496376811594203, 8)
    # figsize = ((-0.001 * num_cols +0.522) , 8)
    # figsize = (num_cols, num_rows)
    fig, ax = plt.subplots(2,1,figsize=figsize)
    fig.suptitle(f'aggregation over base and validation plot, sorted over product of p values, chr-start positions ({count_type})\n{projects}, {gender},\n {drugs}, cutoff={cutoff}\n')
    # just plot the UP plots if there are also UP plots available:
    if 'UP' in base_DF_pval_sort.set_index('plot_type').loc['base_plot']['CMP'].values:
        up_chr_start = base_DF_pval_sort.set_index(['plot_type', 'CMP']).sort_index().loc[('base_plot', 'UP'), :]['chr_start'].values
        temp_UP = base_DF_pval_sort.set_index('chr_start').loc[up_chr_start,:].reset_index()
        temp_UP = temp_UP.set_index('plot_type', drop=False)
        temp_UP.loc['DOWN_validation','plot_type'] = 'validation, NOT scored'
        temp_UP = temp_UP.reset_index(drop=True)
        temp_UP['p_prod'] = 'p_val_prod'
        sns.lineplot(ax=ax[0], data=temp_UP, x='chr_start', y='p_value_life', marker='.', hue='plot_type', hue_order=['validation, NOT scored', 'UP_validation', 'base_plot'], palette=sns.color_palette()[:3])
        sns.lineplot(ax=ax[0], data=temp_UP, x='chr_start', y='p_val_life_prod_scored', marker='.', legend='auto', hue='p_prod', palette=[sns.cubehelix_palette()[4]])
        ax[0].grid()
        ax[0].set_yscale('log')
        ax[0].set_xticks(ax[0].get_xticks(),ax[0].get_xticklabels(), rotation=90)
        ax[0].title.set_text('UP plots')
        ax[0].axhline(y=0.05, color='red')
        trans = transforms.blended_transform_factory(ax[0].get_yticklabels()[0].get_transform(), ax[0].transData)
        ax[0].text(0.01, 0.05, "0.05", color="black", transform=trans, ha="left", va="bottom")
    if 'DOWN' in base_DF_pval_sort.set_index('plot_type').loc['base_plot']['CMP'].values:
        down_chr_start = base_DF_pval_sort.set_index(['plot_type', 'CMP']).sort_index().loc[('base_plot', 'DOWN'), :]['chr_start'].values
        temp_DOWN = base_DF_pval_sort.set_index('chr_start').loc[down_chr_start,:].reset_index()
        temp_DOWN = temp_DOWN.set_index('plot_type', drop=False)
        temp_DOWN.loc['UP_validation','plot_type'] = 'validation, NOT scored'
        temp_DOWN = temp_DOWN.reset_index(drop=True)
        temp_DOWN['p_prod'] = 'p_val_prod'
        sns.lineplot(ax=ax[1], data=temp_DOWN, x='chr_start', y='p_value_life', marker='.', hue='plot_type', hue_order=['validation, NOT scored', 'DOWN_validation', 'base_plot'], palette=sns.color_palette()[:3])
        sns.lineplot(ax=ax[1], data=temp_DOWN, x='chr_start', y='p_val_life_prod_scored', marker='.', legend='auto', hue='p_prod', palette=[sns.cubehelix_palette()[4]])
        ax[1].grid()
        ax[1].title.set_text('DOWN plots')
        ax[1].axhline(y=0.05, color='red')
        ax[1].set_yscale('log')
        trans = transforms.blended_transform_factory(ax[1].get_yticklabels()[0].get_transform(), ax[1].transData)
        ax[1].text(0.01, 0.05, "0.05", color="black", transform=trans, ha="left", va="bottom")
    plt.xticks(rotation=90)
    plt.xlim(left=-1)
    plt.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.3, hspace=1.2)
    print(f'saving plot {out_file}')
    plt.savefig(out_file, bbox_inches='tight')
    print(f'saving table {out_file_table}')
    base_DF_pval_sort.to_csv(out_file_table, sep='\t')

def plot_dynamical_deseq(DF, out_file, out_file_table, count_type):
    """
    using here the ENSG as feature for dissemination
    """
    base_DF = DF
    if re.search('lifelines_evaluated',os.path.split(lifelines_evaluated)[1]):
        if len(base_DF.value_counts('ENSG').value_counts() > 1):
            # ENSGs might occur doubled, since it is possible that in one thresh no
            # scoring happended, but in another it was, drop those doubled ENSGs which
            # are not scored
            base_DF = base_DF.reset_index()
            ensg_to_delete = base_DF.value_counts('ENSG')[base_DF.value_counts('ENSG') > 1 ].index
            num_index = base_DF.set_index(['ENSG', 'scored'], append=True).reset_index(0).loc[(ensg_to_delete, True), 'level_0'].values
            # base_DF = base_DF.drop(num_index).set_index('plot_type')
    # in deseq, the count_type wildcard is available:
    base_DF = base_DF[base_DF['count_type'] == count_type]
    try:
        base_DF.sort_values('p_val_life_prod_scored', inplace=True)
        # p_val_sorted_ENSG = base_DF['ENSG'].drop_duplicates().values.tolist()
    except Exception as e:
        print(f'exception {e}\nnot enough values for plot, writing empty files {out_file} and {out_file_table}')
        pd.DataFrame().to_csv(out_file_table, sep='\t')
        open(out_file, 'w').close()
        os._exit(0)
    # base_DF_pval_sort = base_DF.set_index('ENSG').loc[p_val_sorted_ENSG, :]
    base_DF_pval_sort = base_DF
    # num_rows, num_cols = base_DF_pval_sort.shape
    num_cols = len(base_DF_pval_sort.reset_index()['ENSG'].drop_duplicates())
    # figsize = ( 0.028*num_cols+9.921, 8)
    figsize = (0.09541062801932366*num_cols + 8.496376811594203, 10)
    # figsize = (num_cols * 0.8, num_rows * 0.02)
    fig, ax = plt.subplots(2,1,figsize=figsize)
    fig.suptitle(f'aggregation over base and validation plot, sorted over product of p values, ENSGs ({count_type})\n{projects}, {gender},\n {drugs}, cutoff={cutoff}\n')
    # ## plotting the UP values and the not Scored validation values
    # do not take all CMP == UP like this:
    # temp_UP = base_DF_pval_sort[base_DF_pval_sort['CMP'] == 'UP'].copy(deep=True)
    # there are a few validation plots which have the oposite CMP, but they are not scored anyways, select through the base UP and belonging ENSG:
    up_ensg = base_DF_pval_sort.set_index(['plot_type', 'CMP']).sort_index().loc[('base_plot', 'UP'), :]['ENSG'].values
    temp_UP = base_DF_pval_sort.set_index('ENSG').loc[up_ensg,:].reset_index()
    temp_UP = temp_UP.set_index('plot_type', drop=False)
    temp_UP.loc['DOWN_validation','plot_type'] = 'validation, NOT scored'
    temp_UP = temp_UP.reset_index(drop=True)
    temp_UP['p_prod'] = 'p_val_prod'
    sns.lineplot(ax=ax[0], data=temp_UP, x='ENSG', y='p_value_life', marker='.', hue='plot_type', hue_order=['validation, NOT scored', 'UP_validation', 'base_plot'], palette=sns.color_palette()[:3])
    sns.lineplot(ax=ax[0], data=temp_UP, x='ENSG', y='p_val_life_prod_scored', marker='.', legend='auto', hue='p_prod', palette=[sns.cubehelix_palette()[4]])
    ax[0].grid()
    ax[0].set_yscale('log')
    ax[0].set_xticks(ax[0].get_xticks(),ax[0].get_xticklabels(), rotation=90)
    ax[0].title.set_text('UP plots')
    ax[0].axhline(y=0.05, color='red')
    trans = transforms.blended_transform_factory(ax[0].get_yticklabels()[0].get_transform(), ax[0].transData)
    ax[0].text(0.01, 0.05, "0.05", color="black", transform=trans, ha="left", va="bottom")

    down_ensg = base_DF_pval_sort.set_index(['plot_type', 'CMP']).sort_index().loc[('base_plot', 'DOWN'), :]['ENSG'].values
    temp_DOWN = base_DF_pval_sort.set_index('ENSG').loc[down_ensg,:].reset_index()
    temp_DOWN = temp_DOWN.set_index('plot_type', drop=False)
    temp_DOWN.loc['UP_validation','plot_type'] = 'validation, NOT scored'
    temp_DOWN = temp_DOWN.reset_index(drop=True)
    temp_DOWN['p_prod'] = 'p_val_prod'
    sns.lineplot(ax=ax[1], data=temp_DOWN, x='ENSG', y='p_value_life', marker='.', hue='plot_type', hue_order=['validation, NOT scored', 'DOWN_validation', 'base_plot'], palette=sns.color_palette()[:3])
    sns.lineplot(ax=ax[1], data=temp_DOWN, x='ENSG', y='p_val_life_prod_scored', marker='.', legend='auto', hue='p_prod', palette=[sns.cubehelix_palette()[4]])
    ax[1].grid()
    ax[1].set_yscale('log')
    ax[1].title.set_text('DOWN plots')
    ax[1].axhline(y=0.05, color='red')
    trans = transforms.blended_transform_factory(ax[1].get_yticklabels()[0].get_transform(), ax[1].transData)
    ax[1].text(0.01, 0.05, "0.05", color="black", transform=trans, ha="left", va="bottom")
    plt.xticks(rotation=90)
    plt.xlim(left=-1)
    plt.subplots_adjust(left=0.06, right=0.98, top=0.88, bottom=0.3, hspace=1.3)
    print(f'saving plot {out_file}')
    plt.savefig(out_file, bbox_inches='tight')
    print(f'saving table {out_file_table}')
    base_DF_pval_sort.to_csv(out_file_table, sep='\t')


if pipeline == 'metilene':
    plot_dynamical_metilene(DF, plot_diffs, plot_diffs_table)
elif pipeline == 'DESeq2':
    plot_dynamical_deseq(DF, plot_diffs, plot_diffs_table, count_type)
# plot_dynamical(base_DF, 'p_value_life', plot_diffs_p_val , plot_diffs_p_val_table)
