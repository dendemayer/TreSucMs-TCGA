import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import os
import seaborn as sns
# import glob
from create_matrix_new import set_logger


# fct 12 in main
def drug_frequency(OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID,
                   DRUGS_title):
    """
    :param: OUTPUT_PATH: path for DESeq2 pipeline outputs
    :type: OUTPUT_PATH: str
    :param: DRUGS: applied drug(s)
    :type: DRUGS: list of str
    :param: SCRIPT_PATH: path to the DESeq2_pipeline repo
    :type: SCRIPT_PATH: str
    :param: PROJECT_DRUG_UUID: hash holding a project to the UUID of the\
    belonging drugtable
    :type: PROJECT_DRUG_UUID: dict

    create an overview of all available drugs, of the applied projects in
    DESeq2 output dir

    .. _function_12:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 12 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 12

    """
    # ##### create the drug frequencies in the PROJECT_title, aggregate all
    # applied projects
    # DRUGS_title = '_'.join(sorted(map(str.lower, DRUGS)))
    project_list = []
    for project in PROJECT_DRUG_UUID:
        project_list.append(project)
    PROJECT_title = '_'.join(sorted(map(str.upper, project_list)))
    # print("\nOUTPUT_PATH:\n", OUTPUT_PATH)
    DF_list = []
    for project in PROJECT_DRUG_UUID:
        DF_list.append(
            pd.read_csv(
                os.path.join(
                    OUTPUT_PATH, project,
                    'DF_3t_both_with_DRUG_combi_frequency.tsv'),
                sep='\t', index_col=0))
    DF_complete = pd.concat(DF_list, axis=1)
    file_name = os.path.join(
        OUTPUT_PATH, PROJECT_title, 'DRUG_combi_frequency.tsv')
    # to make is visually usefull, drop every row, where 1 is max:
    DF_complete = DF_complete[DF_complete.max(axis=1) > 1]
    DF_complete.to_csv(file_name, sep='\t')
    logger = set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title)
    log_file = file_name.replace(OUTPUT_PATH + os.path.sep, '')
    logger.info('projects_drug_frequency_12:\t{}'.format(log_file))

    # #### horizontal barplot of drugs
    ax = DF_complete.plot.barh(figsize=(10, 8), legend=True,
                               title='frequency counts drugnames')
    plt.tight_layout()
    fig = ax.get_figure()
    file_name = os.path.join(OUTPUT_PATH, PROJECT_title,
                             'DRUG_combi_frequency.pdf')
    fig.savefig(os.path.join(file_name))
    log_file = file_name.replace(OUTPUT_PATH + os.path.sep, '')
    logger.info('projects_drug_frequency_12:\t{}'.format(log_file))
    plt.close()

    # ##### heatmap of drugs
    # grid_kws = {"height_ratios": (.05, .9), "hspace": .3}
    # f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws)
    # # ax = sns.heatmap(flights, ax=ax,
    #                  # cbar_ax=cbar_ax,
    #                  # cbar_kws={"orientation": "horizontal"})
    # grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
    # fig, (ax, cbar_ax) = plt.subplots(1,2, gridspec_kw=grid_kws, figsize=(12,
    # 14))
    fig, ax = plt.subplots(figsize=(14, 14))
    # fig, ax = plt.subplots()
    label_list = []
    for i in DF_complete.columns:
        label_list.append(i.split('-')[1])

    sns.heatmap(
        DF_complete,
        xticklabels=label_list,
        ax=ax,
        yticklabels=1,
        annot=True,
        fmt='g').set_title(
        'DRUG_combi_frequency for\n{}\n'.format(PROJECT_title),
        fontsize=18)
    # cbar_kws={'label': 'number of Drugs'}
    # .set_xticklabels(ax.get_xticklabels(), rotation=90)

    plt.xlabel('Projects', fontsize=18)
    plt.ylabel('Drugs', fontsize=18)
    plt.tick_params(which='major', labelsize=12, labelbottom=False,
                    bottom=False, top=False, labeltop=True)
    # plt.xticks(rotation=90)
    plt.tight_layout()
    file_name = os.path.join(
        OUTPUT_PATH,
        PROJECT_title,
        'DRUG_combi_frequency_heatmap.pdf')
    plt.savefig(file_name)
    log_file = file_name.replace(OUTPUT_PATH + os.path.sep, '')
    logger.info('projects_drug_frequency_12:\t{}'.format(log_file))
    plt.clf()
    plt.close()
    # #######################################################################
    # create the drug frequencies of all single projects in the outputdir for
    # the overall table for the interactive choice of the drugs


# fct 14 in main
def drug_frequency_all_single_projects(
        OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID, DRUGS_title):
    """
    :param: OUTPUT_PATH: path for DESeq2 pipeline outputs
    :type: OUTPUT_PATH: str
    :param: DRUGS: applied drug(s)
    :type: DRUGS: list of str
    :param: SCRIPT_PATH: path to the DESeq2_pipeline repo
    :type: SCRIPT_PATH: str
    :param: PROJECT_DRUG_UUID: hash holding a project to the UUID of the\
    belonging drugtable
    :type: PROJECT_DRUG_UUID: dict

    take every drug frequency out of the single project folders, therefore walk
    in the OUTPUT_PATH/TCGA-[2..4] folders...

    .. _function_14:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 14 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 14

    """
    # DRUGS_title = '_'.join(sorted(map(str.lower, DRUGS)))
    project_list = []
    for project in PROJECT_DRUG_UUID:
        project_list.append(project)
    PROJECT_title = '_'.join(sorted(map(str.upper, project_list)))

    pattern = ['TCGA-????/', 'TCGA-???/', 'TCGA-??/']
    path_list = []

    for pattern in pattern:
        for path in Path(OUTPUT_PATH).glob(pattern):
            path_list.append(path)

    DF_drug_frequency = []
    for path in path_list:
        path = os.path.join(path, 'DF_3t_both_with_DRUG_combi_frequency.tsv')
        DF_drug_frequency.append(pd.read_csv(path, sep='\t', index_col=0))
    DF_drug_frequency = pd.concat(DF_drug_frequency, axis=1)
    sorted_cols = sorted(DF_drug_frequency.columns)
    DF_drug_frequency_sorted = DF_drug_frequency.loc[:, sorted_cols]

    file_name = os.path.join(
        OUTPUT_PATH, 'DRUG_combi_frequency_all_sorted.tsv')
    DF_drug_frequency_sorted.to_csv(file_name, sep='\t')
    logger = set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title)
    log_file = file_name.replace(OUTPUT_PATH + os.path.sep, '')
    logger.info('projects_drug_frequency_all_14:\t{}'.format(log_file))

    file_name = os.path.join(
        OUTPUT_PATH, 'DRUG_combi_frequency_all.tsv')
    DF_drug_frequency.to_csv(file_name, sep='\t')
    log_file = file_name.replace(OUTPUT_PATH + os.path.sep, '')
    logger.info('projects_drug_frequency_all_14:\t{}'.format(log_file))

    fig, ax = plt.subplots(figsize=(15, 200))
    # fig, ax = plt.subplots(figsize=(DF_complete.shape[0]*0.3,  # (b,h)
    # DF_complete.shape[1]*10))
    # fig, ax = plt.subplots()
    sns.heatmap(
        DF_drug_frequency, ax=ax, yticklabels=1,
        annot=True).set_title(
            'DRUG_combi_frequency for every single Project in {}'.format(
                OUTPUT_PATH))

    plt.tight_layout()
    file_name = os.path.join(OUTPUT_PATH, 'DRUG_combi_frequency_heatmap.pdf')
    plt.savefig(file_name)
    log_file = file_name.replace(OUTPUT_PATH + os.path.sep, '')
    logger.info('projects_drug_frequency_all_14:\t{}'.format(log_file))
    plt.clf()
    plt.close()

    # ax = DF_drug_frequency.plot.barh(figsize=(10, 8), legend=True,
    # title='frequency counts drugnames')
    # plt.tight_layout()
    # fig = ax.get_figure()
    # fig.savefig(os.path.join(OUTPUT_PATH, 'DF_3t_both_with_DRUG_combi.pdf'))
    # plt.close()

    # DF_drug_frequency.name = PROJECT
    # DF_drug_frequency.to_csv(
    # os.path.join( OUTPUT_PATH, 'DF_3t_both_with_DRUG_combi_frequency.tsv'),
    # sep='\t', na_rep='NaN')


if __name__ == '__main__':
    drug_frequency()
