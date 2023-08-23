import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import sys

sys.stdout = sys.stderr = open(snakemake.log[0], 'w')

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

metilene_lifeline_aggregated = snakemake.input[0]
plot_diffs = snakemake.output.plot_diffs[0]
plot_diffs_table = snakemake.output.plot_diffs_table[0]

projects = ', '.join(snakemake.wildcards.project.split('_'))
gender = ', '.join(snakemake.wildcards.gender.split('_'))
drugs = '; '.join(snakemake.wildcards.drug_combi.split('_'))

"""
visualize the differences of the threshold parameter, based on p_value and the
life_mean_diff
"""
###############################################################################
#                                 test input                                  #
###############################################################################
# # snakemake inputs:
# metilene_lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_aggregated.tsv.gz"
# script_file = "../tcga_metilene/scripts/plot_thresh_diff.py"
# # snakemake output:
# plot_diffs = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_diffs.pdf"
# plot_diffs_table = "/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_plot_diffs.tsv.gz"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_2"
# project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "female_male"
# cutoff = "cutoff_0"

# projects = project
# gender = gender
# drugs = drug_combi
###############################################################################
#                                 test input                                  #
###############################################################################


DF = pd.read_table(metilene_lifeline_aggregated)
if DF.empty:
    open(plot_diffs, 'w').close()
    DF.to_csv(plot_diffs_table)
    # open(metilene_lifeline_aggregated, 'a').close()
    os._exit(0)
DF = DF.set_index('plot_type')

base_DF = DF.loc['base_plot', :]


# base_DF = base_DF.dropna(subset='ENSG')
# [284 rows x 16 columns]
base_DF['chromosome - start'] = base_DF['DMR'].apply(lambda x: x.split('_')[0] + ' - ' + x.split('_')[1])
# plt.show()
# get the thresh length:
thresh_len = len(base_DF['threshold'].value_counts())
def plot_dynamical(base_DF, out_file, out_file_table):
    # p_val_sorted_DMR = base_DF.set_index('threshold').loc[0, :].sort_values(what_to_plot)['DMR'].to_list()
    p_val_sorted_DMR = base_DF.set_index('threshold').loc[0, :].sort_values('p_value')['DMR'].to_list()
    base_DF_pval_sort = base_DF.set_index('DMR').loc[p_val_sorted_DMR, :]
    # [233 rows x 15 columns]
    num_rows, num_cols = base_DF_pval_sort.shape
    figsize = (num_cols * 0.8, num_rows * 0.02)
    fig, ax = plt.subplots(2,1,figsize=figsize, sharex=True)
    fig.suptitle(f'life mean diffs and p-values for all DMRs\np-value sorted for threshold=0\n{projects}, {gender}, {drugs}')
    # sns.barplot(data=base_DF, x='ENSG', y=what_to_plot, hue='threshold')
    # sns.lineplot(ax=ax[0], data=base_DF_pval_sort, x='chromosome - start', y='p_value', hue='threshold', palette=sns.color_palette('magma'))
    sns.lineplot(ax=ax[0], data=base_DF_pval_sort, x='chromosome - start', y='p_value', hue='threshold', palette=sns.cubehelix_palette()[:thresh_len][::-1])
    sns.lineplot(ax=ax[1], data=base_DF_pval_sort, x='chromosome - start', y='life_mean_diff', hue='threshold', legend=False, palette=sns.cubehelix_palette()[:thresh_len][::-1])
    ax[1].axhline(y=0)
    # plt.title('life mean diffs for all DMRs')
    plt.xticks(rotation=90)
    plt.subplots_adjust(left=0.06, right=0.98, top=0.86, bottom=0.3)
    # new_position = ax.get_position().translated(0, 0.07)  # Adjust the second value for vertical movement
    # ax.set_position(new_position)
    print(f'saving plot {out_file}')
    plt.savefig(out_file)
    print(f'saving table {out_file_table}')
    base_DF_pval_sort.loc[:, ['ENSG', 'chromosome - start', 'life_mean_diff', 'p_value']].to_csv(out_file_table, sep='\t')
    # in addition make a single q-value plot to compare meth diff significance
    # with kaplan meier fct estimation -> not part of snakemake output,
    # additional info...
    fig, ax = plt.subplots(figsize=figsize)
    fig.suptitle(f'q-values from metilene for the found DMRs')
    sns.lineplot(data=base_DF_pval_sort, x='chromosome - start', y='q-value')
    plt.xticks(rotation=90)
    plt.subplots_adjust(left=0.06, right=0.98, top=0.86, bottom=0.3)
    plt.savefig(out_file.replace('.pdf', '_q_values.pdf'))


plot_dynamical(base_DF, plot_diffs, plot_diffs_table)
# plot_dynamical(base_DF, 'p_value', plot_diffs_p_val , plot_diffs_p_val_table)
