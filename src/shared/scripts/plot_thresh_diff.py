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

    lifeline_aggregated = snakemake.input.lifeline_aggregated
    plot_diffs = snakemake.output.plot_diffs
    plot_diffs_table = snakemake.output.plot_diffs_table
    project = snakemake.wildcards.project
    # projects = ', '.join(snakemake.wildcards.project.split('_'))
    # gender = ', '.join(snakemake.wildcards.gender.split('_'))
    gender = snakemake.wildcards.gender
    # drugs = '; '.join(snakemake.wildcards.drug_combi.split('_'))
    drug_combi = snakemake.wildcards.drug_combi
    pipeline = snakemake.wildcards.pipeline
    plot_type = snakemake.wildcards.plot_type
    # cutoff = snakemake.wildcards.cutoff.split('_')[1]
    cutoff = snakemake.wildcards.cutoff
    count_type = snakemake.wildcards.count_type
else:
    # snakemake inputs:
    lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz"
    script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../shared/scripts/plot_thresh_diff.py"
    # snakemake output:
    plot_diffs = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs_base_plot-beta_vals.pdf"
    plot_diffs_table = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-HNSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_8/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs_base_plot-beta_vals.tsv.gz"
    # snakemake wildcards:
    output_path = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod"
    project = "TCGA-HNSC"
    pipeline = "metilene"
    drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
    gender = "female"
    cutoff = "cutoff_8"
    threshold_str = "threshold_0_threshold_5_threshold_10_threshold_20"
    plot_type = "base_plot"
    count_type = "beta_vals"
    # # snakemake inputs:
    # lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_evaluated-beta_vals.tsv.gz"
    # script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../shared/scripts/plot_thresh_diff.py"
    # # snakemake output:
    # plot_diffs = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs_DOWN_validation-beta_vals.pdf"
    # plot_diffs_table = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_eval_diffs_DOWN_validation-beta_vals.tsv.gz"
    # # snakemake wildcards:
    # output_path = "/scr/palinca/gabor/TCGA-pipeline_7_pval_prod"
    # project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
    # pipeline = "metilene"
    # drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
    # gender = "female_male"
    # cutoff = "cutoff_0"
    # threshold_str = "threshold_0_threshold_5_threshold_10_threshold_20"
    # plot_type = "DOWN_validation"
    # count_type = "beta_vals"
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
# lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_aggregated.tsv.gz"
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
# lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_aggregated.tsv.gz"
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
# lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_5/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_aggregated.tsv.gz"
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



DF = pd.read_table(lifeline_aggregated)
if DF.empty:
    open(plot_diffs, 'w').close()
    DF.to_csv(plot_diffs_table, sep='\t')
    # open(metilene_lifeline_aggregated, 'a').close()
    os._exit(0)
DF = DF.set_index('plot_type')
DF.drop_duplicates('file_path', inplace=True)

def plot_dynamical_metilene(DF, out_file, out_file_table, plot_type):
    """
    using here the chromosome - start as feature for dissemination
    """
    thresh_len = len(DF['threshold'].value_counts())
    base_DF = DF.loc[plot_type, :]
    # base_DF['chromosome - start'] = base_DF['DMR'].apply(lambda x: x.split('_')[0] + ' - ' + x.split('_')[1])
    base_DF['chromosome - start'] = base_DF['chr_start'].apply(lambda x: x.split('_')[0] + ' - ' + x.split('_')[1])
    # p_val_sorted_DMR = base_DF.set_index('threshold').loc[0, :].sort_values(what_to_plot)['DMR'].to_list()
    try:
        if re.search('lifelines_evaluated',os.path.split(lifeline_aggregated)[1]):
            # in the evaluated case, the best p_value_life is already given
            p_val_sorted_DMR = base_DF.sort_values('p_value_life')['DMR'].to_list()
        else:
            p_val_sorted_DMR = base_DF.set_index('threshold').loc[0, :].sort_values('p_value_life')['DMR'].to_list()
    except Exception as e:
        print(f'exception {e}\nnot enough values for plot, writing empty files {out_file} and {out_file_table}')
        pd.DataFrame().to_csv(out_file_table, sep='\t')
        open(out_file, 'w').close()
        os._exit(0)
    base_DF_pval_sort = base_DF.set_index('DMR').loc[p_val_sorted_DMR, :]
    num_cols = len(base_DF_pval_sort.reset_index()['DMR'].drop_duplicates())
    # figsize = (num_cols * 0.8, num_rows * 0.02)
    # figsize in width , height!
    # figsize = ( num_cols * 0.05, 8)
    # taking f(x) = 0.028*x+9.921 for the width (which is in that plot the
    # uniqe rows of feature
    # figsize = (0.028*num_cols+9.921, 8)
    figsize = (0.09541062801932366*num_cols + 8.496376811594203, 8)
    # figsize = ((-0.001 * num_cols +0.522) , 8)
    # figsize = (num_cols, num_rows)
    fig, ax = plt.subplots(2,1,figsize=figsize, sharex=True)
    if re.search('lifelines_evaluated',os.path.split(lifeline_aggregated)[1]):
        fig.suptitle(f'evaluated life mean diffs and p-values for all DMRs\np-value sorted\n{projects}, {gender}, {drugs}, cutoff={cutoff}\nplot type={" ".join(plot_type.split("_"))}')
    else:
        fig.suptitle(f'aggregated life mean diffs and p-values for all DMRs\np-value sorted for threshold=0\n{projects}, {gender}, {drugs}, cutoff={cutoff}\nplot type={" ".join(plot_type.split("_"))}')
    sns.lineplot(ax=ax[0], data=base_DF_pval_sort, x='DMR', y='p_value_life', marker='.', hue='threshold', palette=sns.cubehelix_palette()[:thresh_len][::-1])
    sns.lineplot(ax=ax[1], data=base_DF_pval_sort, x='DMR', y='life_mean_diff', marker='.', hue='threshold', legend=False, palette=sns.cubehelix_palette()[:thresh_len][::-1])
    ax[0].grid()
    ax[0].set_yscale('log')
    ax[0].axhline(y=0.05)
    trans = transforms.blended_transform_factory(ax[0].get_yticklabels()[0].get_transform(), ax[0].transData)
    ax[0].text(0.01, 0.05, "0.05", color="black", transform=trans, ha="left", va="bottom")

    ax[1].axhline(y=0)
    ax[1].grid()
    plt.xticks(rotation=90)
    plt.xlim(left=-1)
    plt.subplots_adjust(left=0.06, right=0.98, top=0.86, bottom=0.3)
    print(f'saving plot {out_file}')
    plt.savefig(out_file, bbox_inches='tight')
    print(f'saving table {out_file_table}')
    base_DF_pval_sort.to_csv(out_file_table, sep='\t')
    plt.clf()
    plt.cla()
    ###########################################################################
    #                      make the plot chr-start wise:                      #
    ###########################################################################
    num_cols = len(base_DF_pval_sort.reset_index()['chr_start'].drop_duplicates())
    figsize = (0.09541062801932366*num_cols + 8.496376811594203, 8)
    fig, ax = plt.subplots(2,1,figsize=figsize, sharex=True)
    if re.search('lifelines_evaluated',os.path.split(lifeline_aggregated)[1]):
        fig.suptitle(f'evaluated life mean diffs and p-values for all chromosome - start positions\np-value sorted for threshold=0\n{projects}, {gender}, {drugs}, cutoff={cutoff}\nplot type={" ".join(plot_type.split("_"))}')
        sns.lineplot(ax=ax[0], data=base_DF_pval_sort, x='chromosome - start', y='p_value_life', marker='.', palette=sns.cubehelix_palette()[:thresh_len][::-1])
        sns.lineplot(ax=ax[1], data=base_DF_pval_sort, x='chromosome - start', y='life_mean_diff', marker='.', legend=False, palette=sns.cubehelix_palette()[:thresh_len][::-1])
    else:
        fig.suptitle(f'aggregated life mean diffs and p-values for all chromosome - start positions\np-value sorted for threshold=0\n{projects}, {gender}, {drugs}, cutoff={cutoff}\nplot type={" ".join(plot_type.split("_"))}')
        sns.lineplot(ax=ax[0], data=base_DF_pval_sort, x='chromosome - start', y='p_value_life', marker='.', hue='threshold', palette=sns.cubehelix_palette()[:thresh_len][::-1])
        sns.lineplot(ax=ax[1], data=base_DF_pval_sort, x='chromosome - start', y='life_mean_diff', marker='.', hue='threshold', legend=False, palette=sns.cubehelix_palette()[:thresh_len][::-1])
    ax[0].grid()
    ax[0].set_yscale('log')
    ax[0].axhline(y=0.05)
    trans = transforms.blended_transform_factory(ax[0].get_yticklabels()[0].get_transform(), ax[0].transData)
    ax[0].text(0.01, 0.05, "0.05", color="black", transform=trans, ha="left", va="bottom")
    ax[1].axhline(y=0)
    ax[1].grid()
    chr_start_list = base_DF_pval_sort.loc[:, 'chr_start'].values.tolist()
    base_DF_pval_sort.reset_index().set_index('chr_start', inplace=True)
    base_DF_pval_sort = base_DF_pval_sort.reset_index().drop_duplicates(subset=['DMR', 'chr_start', 'threshold']).set_index('chr_start')
    # include the range index in the frame and position the threshold x
    # coordinate based on that index value rather than the str chr-start value
    # to make an offset for centering the label annotation
    # including index:
    base_DF_pval_sort = base_DF_pval_sort.reset_index().reset_index().set_index('chr_start')
    if re.search('lifelines_evaluated',os.path.split(lifeline_aggregated)[1]):
        [ax[1].text(base_DF_pval_sort.loc[chr_start, 'index'] - 0.1 , base_DF_pval_sort.loc[chr_start, 'life_mean_diff'] + 0.01, base_DF_pval_sort.loc[chr_start, 'threshold']) for chr_start in chr_start_list]
    plt.xticks(rotation=90)
    plt.xlim(left=-1)
    plt.subplots_adjust(left=0.06, right=0.98, top=0.86, bottom=0.3)
    out_file = out_file.replace('.pdf', '_chr-start.pdf')
    print(f'saving plot {out_file}')
    plt.savefig(out_file, bbox_inches='tight')

def plot_dynamical_deseq(DF, out_file, out_file_table, plot_type, count_type):
    """
    using here the ENSG as feature for dissemination
    """
    thresh_len = len(DF['threshold'].value_counts())
    # limit the plot to one count type
    base_DF = DF.loc[plot_type, :]
    if re.search('lifelines_evaluated',os.path.split(lifeline_aggregated)[1]):
        if len(base_DF.value_counts('ENSG').value_counts() > 1):
            # ENSGs might occur doubled, since it is possible that in one thresh no
            # scoring happended, but in another it was, drop those doubled ENSGs which
            # are not scored
            base_DF = base_DF.reset_index()
            ensg_to_delete = base_DF.value_counts('ENSG')[base_DF.value_counts('ENSG') > 1 ].index
            num_index = base_DF.set_index(['ENSG', 'scored'], append=True).reset_index(0).loc[(ensg_to_delete, True), 'level_0'].values
            base_DF = base_DF.drop(num_index).set_index('plot_type')
    # in deseq, the count_type wildcard is available:
    base_DF = base_DF[base_DF['count_type'] == count_type]
    try:
        if re.search('lifelines_evaluated',os.path.split(lifeline_aggregated)[1]):
            # in the evaluated case, the best p_value_life is already given
            p_val_sorted_ENSG = base_DF.sort_values('p_value_life')['ENSG'].to_list()
        else:
            p_val_sorted_ENSG = base_DF.set_index('threshold').loc[0, :].sort_values('p_value_life')['ENSG'].to_list()
    except Exception as e:
        print(f'exception {e}\nnot enough values for plot, writing empty files {out_file} and {out_file_table}')
        pd.DataFrame().to_csv(out_file_table, sep='\t')
        open(out_file, 'w').close()
        os._exit(0)
    base_DF_pval_sort = base_DF.set_index('ENSG').loc[p_val_sorted_ENSG, :]
    # num_rows, num_cols = base_DF_pval_sort.shape
    num_cols = len(base_DF_pval_sort.reset_index()['ENSG'].drop_duplicates())
    # figsize = ( 0.028*num_cols+9.921, 8)
    figsize = (0.09541062801932366*num_cols + 8.496376811594203, 8)
    # figsize = (num_cols * 0.8, num_rows * 0.02)
    fig, ax = plt.subplots(2,1,figsize=figsize, sharex=True)
    if re.search('lifelines_evaluated',os.path.split(lifeline_aggregated)[1]):
        fig.suptitle(f'evaluated life mean diffs and p-adjusted for all ENSGs ({count_type})\np-value sorted\n{projects}, {gender}, {drugs}, cutoff={cutoff}\nplot type={" ".join(plot_type.split("_"))}')
        # in evaluated, no hue, and therefore no palette is needed
        sns.lineplot(ax=ax[0], data=base_DF_pval_sort, x='ENSG', y='p_value_life', marker='.')
        sns.lineplot(ax=ax[1], data=base_DF_pval_sort, x='ENSG', y='life_mean_diff', marker='.', legend=False)
    else:
        fig.suptitle(f'aggregated life mean diffs and p-adjusted for all ENSGs ({count_type})\np-value sorted for threshold=0\n{projects}, {gender}, {drugs}, cutoff={cutoff}\nplot type={" ".join(plot_type.split("_"))}')
        sns.lineplot(ax=ax[0], data=base_DF_pval_sort, x='ENSG', y='p_value_life', marker='.', hue='threshold', palette=sns.cubehelix_palette()[:thresh_len][::-1])
        sns.lineplot(ax=ax[1], data=base_DF_pval_sort, x='ENSG', y='life_mean_diff', marker='.', hue='threshold', legend=False, palette=sns.cubehelix_palette()[:thresh_len][::-1])
    ax[1].axhline(y=0)
    ax[0].axhline(y=0.05)
    ax[0].grid()
    ax[0].set_yscale('log')
    trans = transforms.blended_transform_factory(ax[0].get_yticklabels()[0].get_transform(), ax[0].transData)
    ax[0].text(0.01, 0.05, "0.05", color="black", transform=trans, ha="left", va="bottom")
    ax[1].grid()
    # ENSG:
    # base_DF_pval_sort.iloc[0, :].name
    # p_val
    # base_DF_pval_sort.iloc[0, 4]
    # thresh:
    # base_DF_pval_sort.iloc[0, 0]
    # life_mean_diff:
    # base_DF_pval_sort.iloc[0, 7]
    # ax[1].text(base_DF_pval_sort.iloc[0, :].name, base_DF_pval_sort.iloc[0, 7], base_DF_pval_sort.iloc[0, 0])
    base_DF_pval_sort = base_DF_pval_sort.reset_index().reset_index().set_index('ENSG')
    # [ax[0].text(base_DF_pval_sort.loc[ensg, 'index'] - 0.3, base_DF_pval_sort.loc[ensg, 'p_value_life'] + 0.05, base_DF_pval_sort.loc[ensg, 'threshold']) for ensg in base_DF_pval_sort.index.tolist()]
    # adding the threshold value to the marker points:
    if re.search('lifelines_evaluated',os.path.split(lifeline_aggregated)[1]):
        [ax[1].text(base_DF_pval_sort.loc[ensg, 'index'] - 0.25, base_DF_pval_sort.loc[ensg, 'life_mean_diff'] + 0.02, base_DF_pval_sort.loc[ensg, 'threshold']) for ensg in base_DF_pval_sort.index.tolist()]
        # for ensg in base_DF_pval_sort.index.tolist():
        #     try:
        #         ax[1].text(base_DF_pval_sort.loc[ensg, 'index'] - 0.25, base_DF_pval_sort.loc[ensg, 'life_mean_diff'] + 0.02, base_DF_pval_sort.loc[ensg, 'threshold'])
        #     except Exception as e:
        #         print(e)
        #         breakpoint()
    plt.xticks(rotation=90)
    plt.xlim(left=-1)
    plt.subplots_adjust(left=0.06, right=0.98, top=0.86, bottom=0.3)
    print(f'saving plot {out_file}')
    plt.savefig(out_file, bbox_inches='tight')
    print(f'saving table {out_file_table}')
    base_DF_pval_sort.to_csv(out_file_table, sep='\t')


if pipeline == 'metilene':
    plot_dynamical_metilene(DF, plot_diffs, plot_diffs_table, plot_type)
elif pipeline == 'DESeq2':
    plot_dynamical_deseq(DF, plot_diffs, plot_diffs_table, plot_type, count_type)
# plot_dynamical(base_DF, 'p_value_life', plot_diffs_p_val , plot_diffs_p_val_table)
