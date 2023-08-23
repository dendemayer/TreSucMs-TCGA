import os
import pandas as pd


def return_bed_interesect_metilene_files(OUTPUT_PATH, PROJECTS, DRUGS,
                                         cutoffs):
    file_list = []
    drugs = '_'.join(DRUGS)
    # if len(PROJECT) > 1:
    #     PROJECT.append('_'.join(sorted([x.upper() for x in PROJECT])))
    for project in PROJECTS:
        for gender in ['female', 'male', 'female_male']:
            for cutoff in cutoffs:
                cutoff = f'cutoff_{str(cutoff)}'
                file_list.append(
                    os.path.join(
                        OUTPUT_PATH, project,
                        'metilene/metilene_output', drugs, gender,
                        cutoff, 'metilene_intersect.tsv'))
    return file_list


def return_plot_DMR_regions_plot(metilene_intersect_tables):
    """
    we need to check whether there are actually DMR in the intersected files,
    if the intersect table is empty, do not request a plot based on that
    table
    else, hand over file names with pattern like:
    metilene_intersect_boxplot_beta_value_chr12_132887020_132888306.pdf
    """

    DMR_dict = {}
    for x in metilene_intersect_tables:
        if not pd.read_table(x).empty:
            DMR_dict.update(
                {x: pd.read_table(
                    x, usecols=[3],
                    skiprows=5)['region'].value_counts().index.to_list()})
    # to each intersect table, we hold the DMR belonging to it, now populate
    # the pdf out filenames:
    metilene_plots = []
    for filename in DMR_dict.keys():
        for DMR in DMR_dict[filename]:
            metilene_plots.append(
                filename.replace('.tsv', f'_boxplot_beta_value_{DMR}.pdf'))
            metilene_plots.append(
                filename.replace(
                    '.tsv', f'_lineplot_median_beta_value_{DMR}.pdf'))
    return metilene_plots


def return_DMR_merge_plots(OUTPUT_PATH, PROJECTS, DRUGS, cutoffs, threshold):
    file_list = []
    drugs = '_'.join(DRUGS)
    # if len(PROJECT) > 1:
    #     PROJECT.append('_'.join(sorted([x.upper() for x in PROJECT])))
    thresholds = [f'threshold_{str(i)}' for i in threshold]
    cutoffs = [f'cutoff_{str(i)}' for i in cutoffs]
    for project in PROJECTS:
        for gender in ['female', 'male', 'female_male']:
            for cutoff in cutoffs:
                for threshold in thresholds:
                    file_list.append(
                        os.path.join(
                            OUTPUT_PATH, project,
                            'metilene/metilene_output', drugs, gender,
                            cutoff, threshold, 'All_plots_merged.pdf'))
    return file_list
