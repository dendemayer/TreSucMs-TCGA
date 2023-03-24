import os
import pandas as pd


def return_bed_interesect_metilene_files(OUTPUT_PATH, PROJECT, DRUGS, cutoffs):
    file_list = []
    drugs = '_'.join(DRUGS)
    if len(PROJECT) > 1:
        PROJECT.append('_'.join(sorted([x.upper() for x in PROJECT])))
    for project in PROJECT:
        for gender in ['female', 'male', 'female_male']:
            for complement in ['', '_complement']:
                for cutoff in cutoffs:
                    cutoff = f'cutoff_{str(cutoff)}'
                    file_list.append(
                        os.path.join(
                            OUTPUT_PATH, project,
                            'metilene/metilene_output', drugs, gender,
                            cutoff, f'metilene{complement}_intersect.tsv'))
    return file_list


def return_plot_DMR_regions_plot(metilene_intersect_tables):
    """
    we need to check whether there are actually DMR in the intersected files,
    if the intersect table is empty, do not request a plot based on that
    table
    else, hand over file names with pattern like:
    metilene_intersect_boxplot_beta_value_chr12_132887020_132888306.pdf
    metilene_complement_intersect{}.tsv
    """

    DMR_dict = {}
    for x in metilene_intersect_tables:
        if not pd.read_table(x).empty:
            DMR_dict.update(
                {x: pd.read_table(
                    x, usecols=[3], skiprows=5)['region'].value_counts(
                    ).index.to_list()})
    # to each intersect table, we hold the DMR belonging to it, now populate
    # the pdf out filenames:
    pdf_file_list = []
    for filename in DMR_dict.keys():
        for DMR in DMR_dict[filename]:
            pdf_file_list.append(
                filename.replace('.tsv', f'_boxplot_beta_value_{DMR}.pdf'))
            pdf_file_list.append(
                filename.replace(
                    '.tsv', f'_lineplot_median_beta_value_{DMR}.pdf'))

    return pdf_file_list
