import pandas as pd
import os
import seaborn as sns
from matplotlib import pyplot as plt
import sys

sys.stdout = sys.stderr = open(snakemake.log[0], 'w')

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

metilene_intersect = snakemake.input.metilene_intersect
pdf_boxplot_out = snakemake.output.pdf_boxplot_out
pdf_lineplot_out = snakemake.output.pdf_lineplot_out
DMR = snakemake.wildcards.DMR

# make use of MI:
# # (Pdb) DF_DMR.index.names
# # FrozenList(['Chromosome', 'Start', 'End', 'region'])
# # (Pdb) DF_DMR.columns.names
# # FrozenList(['vital_status', 'case_id', 'drugs', 'gender', 'projects'])
DF_DMR = pd.read_table(metilene_intersect, header=[0, 1, 2, 3, 4], index_col=[0, 1, 2, 3], na_values='.')
# out of the MI index parse the regions:
# # (Pdb) DF_DMR.index.names
# # FrozenList(['Chromosome', 'Start', 'End', 'region'])

# range_list = pd.Index(
#     [i[3] for i in DF_DMR.index.tolist()]).value_counts().index.to_list()
project_list = pd.Index(
    [i[4] for i in DF_DMR.columns]).value_counts().index.to_list()
palette_len = len(project_list) * 2


# for plotting access each DMR seperately:
# do not loop over every range out of the intersect table since the DMR are
# requested already by the input function:
# def return_plot_DMR_regions_plot(metilene_intersect_tables): in modules/bed_intersect_metilene.py
# for range_ in range_list:
# limit the DF to the recent region:
try:
    DF_to_plot = DF_DMR.loc[(slice(None), slice(None), slice(None), DMR), :]
except Exception as e:
    # KeyError('chr19_58228367_58228578') in meta table TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin_paclitaxel/female_male/cutoff_0/metilene_complement_intersect.tsv
 # vital_status |          |          |                         | alive              | dead               | dead
 # chr19        | 58228368 | 58228369 | chr19_58228367_58228578 | .                  | .                  | .
 # chr19        | 58228412 | 58228413 | chr19_58228367_58228578 | .                  | .                  | .
 # chr19        | 58228437 | 58228438 | chr19_58228367_58228578 | .                  | .                  | .
 # chr19        | 58228439 | 58228440 | chr19_58228367_58228578 | .                  | .                  | .
 # --> found ranges for cases for which we dont have the
    print(e)
    print(DF_DMR)
    # breakpoint()
    open(pdf_boxplot_out, 'a').close()
    open(pdf_lineplot_out, 'a').close()
    os._exit(0)
DF_to_plot_median = DF_to_plot.groupby( by=['vital_status', 'projects'], axis=1).median().reset_index('Start')
projects_list = [ i[4] for i in DF_to_plot.columns]
vital_array = [ i[0] for i in DF_to_plot.columns]
new_co_MI = pd.MultiIndex.from_arrays([vital_array, projects_list], names=('vital_status', 'projects'))
DF_to_plot.columns = new_co_MI
DF_to_plot.reset_index(level=[0,2,3], drop=True, inplace=True)
DF_to_plot = DF_to_plot.T.reset_index()
DF_to_plot['hue'] = DF_to_plot['vital_status'] + '_' + DF_to_plot['projects']
hue = DF_to_plot['hue'].to_list()
DF_to_plot = DF_to_plot.iloc[:, 2:-1]
hue = hue * DF_to_plot.shape[1]
DF_to_plot = DF_to_plot.melt()
DF_to_plot.rename({'value': 'beta_value'}, axis=1, inplace=True)
DF_to_plot['hue'] = hue
DF_to_plot.sort_values(by=['Start', 'hue'], inplace=True)
range_title = DMR.split('_')
range_title = f'{range_title[0]}: {range_title[1]}-{range_title[2]}'
sns.set(style="ticks", palette=sns.color_palette('coolwarm_r', palette_len))
plot = sns.boxplot(x='Start', y='beta_value', hue='hue', data=DF_to_plot)
plt.title(f'Range: {range_title}')
plot.set_xticklabels(plot.get_xticklabels(), rotation=45)
plt.legend(title='vital state and projects', bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.tight_layout()
print(f'saving {pdf_boxplot_out}')
plot.figure.savefig(pdf_boxplot_out)
plt.clf()
plt.close()
DF_to_plot = DF_DMR.loc[(slice(None), slice(None), slice(None), DMR), :]
DF_to_plot.columns = new_co_MI
DF_to_plot.reset_index(level=[0, 2, 3], drop=True, inplace=True)
DF_to_plot = DF_to_plot.T.reset_index()
DF_to_plot['hue'] = DF_to_plot['vital_status'] + '_' + DF_to_plot['projects']
DF_to_plot = DF_to_plot.drop(['vital_status', 'projects'], axis=1)
DF_to_plot = DF_to_plot.groupby('hue').median()
hue = DF_to_plot.index.to_list()
hue = hue * DF_to_plot.shape[1]
DF_to_plot = DF_to_plot.melt()
DF_to_plot['hue'] = hue
DF_to_plot.sort_values(by=['Start', 'hue'], inplace=True)
DF_to_plot.rename({'value': 'beta_value_median'}, axis=1, inplace=True)
#     # ### plot just the median as linegraph:
#     DF_line_plot = DF_to_plot.groupby(
#         ['vital_state', 'project', 'Start']).median(
#     ).reset_index(
#         ['Start', 'vital_state', 'project']).rename(
#         {'mean_beta_value': 'median of means beta value'},
#         axis=1)
#     DF_line_plot['vital_proj'] = DF_line_plot['vital_state']\
#         + '_' + DF_line_plot['project']

plt.figure(figsize=(15, 5))
#     # plt.rcParams["figure.autolayout"] = True
#     # DF_line_plot['Start'] =
#     # pd.Categorical(DF_line_plot['Start'])
plt.figure()
plt.ticklabel_format(
    style='plain', axis='x', useOffset=False)
# # style='vital_proj', markers=True
plot = sns.lineplot(
    data=DF_to_plot,
    x="Start",
    marker='o',
    y="beta_value_median", #  linewidth=3,
    hue="hue")
plt.title(f'Range: {range_title}')
plt.setp(plot.get_xticklabels(), rotation=45)

legend = plt.legend(
    title='vital state of all projects',
    bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
for legobj in legend.legend_handles:
    legobj.set_linewidth(3.0)
#     # plt.legend(
#         # title='vital state of all projects',
#         # bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#     # plt.legend(
#     # title='vital state of all projects', bbox_to_anchor=(1,
#     # 1))
plt.tight_layout()
print(f'saving {pdf_lineplot_out}')
plot.figure.savefig(pdf_lineplot_out)

plt.clf()
plt.close()
