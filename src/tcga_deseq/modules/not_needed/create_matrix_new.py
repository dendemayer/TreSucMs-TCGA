# import glob
# import gzip
# import logging
# import pdb
# import sys
# from pathlib import Path
import glob
# import json
import logging
# import mygene
import os
import pandas as pd
import re
import shutil
import subprocess       # to call R
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
# import numpy as np
# import shutil  # needed to copy meta_info.dat in reproducible run
# ############### heres the automated api stuff
# whats needed: 3 tables,
# * one with metainfo (gender, vital_status, case_id) it is highly desired that
# here the therapeutics are also included, but yet not possible, therefore:
# * one with the therapeutic_agents, linked to case_id
# * MANIFEST infos from the MANIFEST are actually not needed, but could be
# important for later comparison of completeness of the dataset used


# fct 5 in main
def create_summary_table(OUTPUT_PATH, PROJECT, DRUGS_title, drugs):
    """
    :param: OUTPUT_PATH: path for metilene pipeline outputs
    :type: OUTPUT_PATH: str
    :param: PROJECT: list of projects chosen
    :type: PROJECT: list of str
    :param: DRUGS_title: merged drug title out of multiple drugs
    :type: DRUGS_title: str

    creating a summary beta value table with every single methyl file info
    aggregated, cols are multinamed to keep info like:

        * vital_status
        * case_id
        * gender
        * therapeutics
        * project

    summary_table filename:
        complete_summary.tsv

    .. _function_5:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 5 to your call,
        # example:
        $ python main_metilene.py -p TCGA-CESC -d cisplatin -f 5
    ### new api attempt:
    Composite Element …| Chromosome | Start     | alive;f5d90...;
                                                  cisplatin;female;TCGA-CESC
    cg00000029         | chr16      | 53434200  | 0.228704153418621
    cg00000108         | chr3       | 37417715  | NA
    # not filtered on DRUGS and gender:
    -> to keep the compatibility with new and old api style, create a complete
    summary
        - summary_DF.tsv in the form of:
            genes           | 6ff12a54-10da-4941…| 748aada0-448f-4523…
            ENSG00000000003 | 3423               | 2162
            ENSG00000000005 | 0                  | 0

    -> filtered on gender AND drugs
    -> to keep the compatibility with new and old api style, create filtered
    summary in the kind of :
        the multiindex names:
            - DRUG_combi_female_summary_dead_alive_INFO.tsv
            - DRUG_combi_male_summary_dead_alive_INFO.tsv
            - DRUG_combi_summary_dead_alive_INFO.tsv
        content:
                     | 0                  | 1                  | ...
        vital_status | alive              | alive              | ...
        gender       | female             | female             | ...
        PROJECT      | TCGA-HNSC          | TCGA-HNSC          | ...
        drugnames    | carboplatin,paclit…| carboplatin,paclit…| ...
        case_id      | b74519ff-b338-45e9…| aa10d6da-ba20-43e8…| ...

        the count data summary names:
            - DRUG_combi_female_summary_dead_alive.tsv
            - DRUG_combi_male_summary_dead_alive.tsv
            - DRUG_combi_summary_dead_alive.tsv

        content:
                        | 0    | 1    | 2
        ENSG00000000003 | 1459 | 1094 | 1449
        ENSG00000000005 | 0    | 0    | 0
        ENSG00000000419 | 1239 | 2058 | 1698

    """
    # just put those meth file in the summary table which are in the drugtable,
    # use the "bcr_patient_barcode" str to grep the files
    # print('in fct 5')
    data_dir = PROJECT + '_data_files'
    data_dir = os.path.join(OUTPUT_PATH, PROJECT, data_dir)
    dir_list = os.listdir(data_dir)

    try:
        merged_df = pd.read_csv(os.path.join(
            OUTPUT_PATH, PROJECT,
            'meta_info_druglist_merged_drugs_combined.tsv'), sep='\t')
    except FileNotFoundError as e:
        print(e)
        print("run fct create_merged_metatable(OUTPUT_PATH, PROJECT)\n")
        print("or fct sep_down_data_files(OUTPUT_PATH, PROJECT) before\n")
        print("exiting now")
        os._exit(1)
    merged_df.rename({'bcr_patient_barcode_x': 'bcr_patient_barcode'}, axis=1,
                     inplace=True)

    # with the complement of drugs, also nan drugs are included, that
    # could raise an error, set those values to 'Not Available'
    # sequence item 2: expected str instance, float found
    index_to_str = merged_df[
        merged_df['pharmaceutical_therapy_drug_name'].isnull()].index
    merged_df.loc[
        index_to_str, ['pharmaceutical_therapy_drug_name']] = 'Not Available'
    # limit the files to read in to the files we actually have meta data
    # provided through meta_info_druglist_merged_drugs_combined.tsv
    # this is best accomplished through creating a set
    dir_set = set(dir_list).intersection(set(merged_df['filename'].to_list()))
    dir_list = list(dir_set)
    dir_path_list = [os.path.join(data_dir, dir_) for dir_ in dir_list]
    first = True

    def transform_ENSG(count_table, merged_df, PROJECT):
        '''
        - selecting just rows from the htseq countfiles which include ENSG
          pattern
        - discarding the .[:digit:] notation from ENSG notation
        - example: ENSG00000000003.13 -> ENSG00000000003
        '''
        temp_DF = pd.read_csv(count_table, sep='\t', header=None,
                              names=['ENSG', 'count'])
        # change the index s.t. just ENSG rows are included and .digit
        # patterns are stripped
        filtered_ENSG = temp_DF[temp_DF[
            'ENSG'].str.contains(
                '^ENSG')]['ENSG'].apply(lambda x: x.split('.')[0])
        temp_DF = temp_DF.loc[filtered_ENSG.index, :]
        temp_DF['ENSG'] = filtered_ENSG
        temp_DF.set_index('ENSG', inplace=True)
        # get the metainfo to the datafile.
        # we need:
        # vital_status gender      PROJECT     drugnames   case_id
        compl_meta = merged_df[merged_df['filename'] == os.path.basename(
            count_table)]
        vital_status = compl_meta['vital_status'].values[0]
        gender = compl_meta['gender'].values[0]
        drugs = compl_meta['pharmaceutical_therapy_drug_name'].values[0]
        case_id = compl_meta['case_ids'].values[0]
        compl_meta = merged_df[merged_df[
            'filename'] == os.path.basename(count_table)]
        multiindex = pd.MultiIndex.from_arrays(
            [[vital_status], [gender], [PROJECT], [drugs], [case_id]], names=(
                'vital_status', 'gender', 'project', 'drugs', 'case_id'))
        temp_DF.columns = multiindex
        return temp_DF
    first = True
    summary_DF = pd.DataFrame()
    for count_table in dir_path_list:
        if first:
            summary_DF = transform_ENSG(count_table, merged_df, PROJECT)
            first = False
        else:
            add_DF = transform_ENSG(count_table, merged_df, PROJECT)
            # join to the existing DF:
            summary_DF = pd.concat([summary_DF, add_DF], axis=1)

    # check whether we have both, or just female or male:
    temp_set = set(summary_DF.columns.to_frame()['gender'].values)
    gender_list = list(temp_set)
    # save summary_tables and their INFO (MultiIndex header) according to what
    # genders are available:
    if len(gender_list) == 2:
        gender_list = ['female', 'male', 'both']

    logger = set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
    for gender in gender_list:
        path_to_save = os.path.join(OUTPUT_PATH, PROJECT, DRUGS_title,
                                    'DRUGS_combi_' + gender +
                                    '_summary_dead_alive.tsv')
        path_to_INFO = os.path.join(OUTPUT_PATH, PROJECT, DRUGS_title,
                                    'DRUGS_combi_' + gender +
                                    '_summary_dead_alive_INFO.tsv')
        sum_temp_DF = summary_DF
        sum_temp_DF.index.names = ['']
        # TODO : save here the complete table with no filtering in the
        # PROJECT dir
        if gender == 'both':
            # constraint to the drug selection,
            temp_DF = pd.DataFrame()
            for drug in drugs:
                temp_DF = pd.concat(
                    [temp_DF, sum_temp_DF.loc[:, (
                        slice(None), slice(None), slice(None), drug)]], axis=1)
            sum_temp_DF = temp_DF
            # first save the multiindex.
            multii_DF = sum_temp_DF.columns.to_frame().reset_index(drop=True).T
            try:
                multii_DF.to_csv(path_to_INFO, sep='\t')
                logger.info('create_summary_table_5:\t{}'.format(
                    path_to_INFO.replace(OUTPUT_PATH + os.path.sep, '')))
            except Exception as e:
                print(f'Exception occured: {e}, could not save {path_to_INFO}')
                print('exiting')
                os._exit(0)
            # if gender == 'both', no gender constraining is needed:
            sum_temp_DF.columns = pd.RangeIndex(len(sum_temp_DF.T.index))
            try:
                sum_temp_DF.to_csv(path_to_save, sep='\t')
                os.path.sep
                logger.info('create_summary_table_5:\t{}'.format(
                    path_to_save.replace(OUTPUT_PATH + os.path.sep, '')))
            except Exception as e:
                print(f'Exception occured: {e}, could not save {path_to_save}')
                print('exiting')
                os._exit(0)
        else:
            # constrain DF selection to gender:
            sum_temp_DF = sum_temp_DF.loc[:, (slice(None), gender)]
            # constrain to the right drug selection:
            temp_DF = pd.DataFrame()
            for drug in drugs:
                try:
                    temp_DF = pd.concat(
                        [temp_DF, sum_temp_DF.loc[:, (
                            slice(None), slice(None), slice(None), drug)]],
                        axis=1)
                except KeyError:
                    continue
            sum_temp_DF = temp_DF
            # save the filtered MI
            multii_DF = sum_temp_DF.columns.to_frame().reset_index(drop=True).T
            try:
                multii_DF.to_csv(path_to_INFO, sep='\t')
                logger.info('create_summary_table_5:\t{}'.format(
                    path_to_INFO.replace(OUTPUT_PATH + os.path.sep, '')))
            except Exception as e:
                print(f'Exception occured: {e}, could not save {path_to_INFO}')
                print('exiting')
                os._exit(0)
            sum_temp_DF.columns = pd.RangeIndex(len(sum_temp_DF.T.index))
            try:
                sum_temp_DF.to_csv(path_to_save, sep='\t')
                logger.info('create_summary_table_5:\t{}'.format(
                    path_to_save.replace(OUTPUT_PATH + os.path.sep, '')))
            except Exception as e:
                print(f'Exception occured: {e}, could not save {path_to_save}')
                print('exiting')
                os._exit(0)

    # the complete summary table is needed in PROJECT dir without any
    # filtering on gender or drugs
    summary_DF.index.name = 'genes'
    summary_DF.columns.name = ''
    summary_DF = summary_DF.T.reset_index([0, 1, 2, 3], drop=True).T
    log_path = os.path.join(PROJECT, 'summary_DF.tsv')
    try:
        summary_DF.to_csv(os.path.join(OUTPUT_PATH, log_path), sep='\t')
        logger.info('create_summary_table_5:\t{}'.format(log_path))
    except Exception as e:
        print(f'Exception occured: {e}, could not save {log_path}')
        print('exiting')
        os._exit(0)


# fct 7 in main and 9
def provide_DESeq2_table(PROJECT, OUTPUT_PATH, DRUGS, SCRIPT_PATH, logger,
                         cutoff, DRUGS_title):
    """
    :param: PROJECT: list of projects chosen
    :type: PROJECT: list of str
    :param: OUTPUT_PATH: path for DESeq2 pipeline outputs
    :type: OUTPUT_PATH: str
    :param: DRUGS: applied drug(s)
    :type: DRUGS: list of str
    :param: SCRIPT_PATH: path to the DESeq2_pipeline repo
    :type: SCRIPT_PATH: str
    :param: logger: the adjustet logger with the right filehandler
    :type: logger: logging instance

    filtering case_id according to DRUG query
    providing new table for DESeq2 analysis
    dependend on distinguishable factors of the tables provided, a
    singlefactorial (at least differences in vital state) or a mutlifactorial
    run in DESeq2 is performed (gender, therapy or project )

    .. _function_7:
    .. _function_9:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 7 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 7

        # when choosing multiple projects, call:
        $ python main_deseq.py -p TCGA-CESC  -p TCGA-HNSC -d cisplatin -f 9
    """
    # since the drugs are not filtered out until here, write here new drugdirs
    # and save all subsequent files in them..
    # inserting the drugname within the filename is not needed anymore

    # adjust Drugquery to an single str (ordered, commaseperated)
    # temp save the list, we need to loop through it...
    #
    # DRUGS_title = '_'.join(sorted(map(str.lower, DRUGS)))

    # #### for single PROJECTS:
    if isinstance(PROJECT, str):
        PROJECT_title = PROJECT
        # create a new deseq dir:
        DF_DRUG_combi = pd.read_csv(
            os.path.join(OUTPUT_PATH, PROJECT,
                         'DF_3t_both_with_DRUG_combi.tsv'),
            sep='\t')
    else:  # if PROJECT is no str object, a hash longer 1 is handed over,
        # the keys in it are the projectnames...
        project_list = []
        first = True
        for project in PROJECT:
            if first:
                DF_DRUG_combi = pd.read_csv(
                    os.path.join(
                        OUTPUT_PATH,
                        project, 'DF_3t_both_with_DRUG_combi.tsv'), sep='\t')
                project_list.append(project)
                first = False
            else:
                DF_DRUG_combi = pd.concat(
                    [DF_DRUG_combi, pd.read_csv(
                        os.path.join(OUTPUT_PATH, project,
                                     'DF_3t_both_with_DRUG_combi.tsv'),
                        sep='\t')])
                project_list.append(project)
            PROJECT_title = '_'.join(sorted(map(str.upper, project_list)))
        try:
            log_path = os.path.join(
                PROJECT_title, 'DF_3t_both_with_DRUG_combi_all_projects.tsv')
            DF_DRUG_combi.iloc[:, 1:].to_csv(os.path.join(
                OUTPUT_PATH, log_path), sep='\t', index=None)
            logger.info(f'deseq_out_multi_9:\t{log_path}')
        except Exception as e:
            print(f'error thrown:{e}, could not save:{log_path}')
            print('exiting now')
            os._exit(0)

    try:
        log_path = os.path.join(
            PROJECT_title, 'DF_3t_both_with_DRUG_combi_all_projects.tsv')
        DF_DRUG_combi.iloc[:, 1:].to_csv(os.path.join(
            OUTPUT_PATH, log_path), sep='\t', index=None)
        logger.info(f'deseq_out_multi_7:\t{log_path}')
    except Exception as e:
        print(f'error thrown:{e}, could not save:{log_path}')
        print('exiting now')
        os._exit(0)

    print("creating drug directory: ", os.path.join(OUTPUT_PATH, PROJECT_title,
                                                    DRUGS_title))
    os.makedirs(os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title),
                exist_ok=True)
    # create also the complement of the DRUG selection:
    # a selection could be performed like that:
    # DF_DRUG_combi.sort_values(by=['vital_status', 'gender'], inplace=True)
    # but we take the values we need in the actual order and concatenate them,
    # that way we make sure that no other expressions like 'not available' are
    # handed over:

    DF_DRUG_combi = pd.concat([
        DF_DRUG_combi.loc[DF_DRUG_combi['gender'] == 'female', :],
        DF_DRUG_combi.loc[DF_DRUG_combi['gender'] == 'male', :]
    ])
    # here we include the cutoff parameter, include dead cases which have a
    # survivaltime greater the the cutoff parameter(alive cases are here
    # automatically not invoked, since they have NaN in survivaltime):
    # before, they must be redeclared as alive:
    if cutoff > 0:
        temp_DF = (DF_DRUG_combi[DF_DRUG_combi['survivaltime'] > cutoff])
        DF_save = temp_DF
        DF_save.insert(len(temp_DF.columns), 'cut_off', cutoff)
        # those cases are additionally added, due to the chosen cutoff
        # parameter log them:
        file_name = 'cutoff_cases_add_' + str(cutoff) + '.tsv'
        if not temp_DF.empty:
            DF_save.to_csv(
                os.path.join(
                    OUTPUT_PATH,
                    PROJECT_title,
                    file_name), sep='\t', index=False)
            if isinstance(PROJECT, str):
                logger.info(
                    'deseq_out_single_7:\t{}'.format(
                        os.path.join(PROJECT_title, file_name)))
            else:
                logger.info(
                    'deseq_out_multi_9:\t{}'.format(
                        os.path.join(PROJECT_title, file_name)))

        DF_DRUG_combi.loc[temp_DF.index, 'vital_status'] = 'alive'
    DF_DRUG_combi = pd.concat([
        DF_DRUG_combi.loc[DF_DRUG_combi['vital_status'] == 'alive', :],
        DF_DRUG_combi.loc[DF_DRUG_combi['vital_status'] == 'dead', :],
    ])

    # #######################################################################
    # here the filtering according to DRUGS is performed
    # for every drug in DRUGS, save every case_id
    # ########################################################################
    # we include several therapies, concat the filter:
    filter_list = []
    for drug in DRUGS:
        filt = DF_DRUG_combi['drugnames'] == drug
        filter_list.append(filt)

    filt_DF = pd.concat(filter_list, axis=1)
    filt = filt_DF.any(axis=1)
    DF_DRUG_combi_temp = DF_DRUG_combi[filt]
    DF_DRUG_combi_temp_COMPL = DF_DRUG_combi[~filt]

    # add the cases to alive which are included due to the cutoff parameter
    # consider the drug selection

    if cutoff > 0:
        filter_list = []
        for drug in DRUGS:
            filt = DF_save['drugnames'] == drug
            filter_list.append(filt)

        filt_DF = pd.concat(filter_list, axis=1)
        filt = filt_DF.any(axis=1)
        DF_save = DF_save[filt]
        DF_save_COMPL = DF_save[~filt]

        # save the drugfiltered cases, which are added to alive due to the
        # chosen cutoff parameter
        file_name = 'cutoff_cases_add_' + str(cutoff) + '.tsv'
        if not DF_save.empty:
            DF_save.to_csv(
                os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                             file_name), sep='\t', index=False)
            DF_save_COMPL.to_csv(
                os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                             'COMPLEMENT', file_name), sep='\t', index=False)
            if isinstance(PROJECT, str):
                logger.info(
                    'deseq_out_single_7:\t{}'.format(
                        os.path.join(PROJECT_title, DRUGS_title, file_name)))
                logger.info(
                    'deseq_out_single_7:\t{}'.format(
                        os.path.join(PROJECT_title, DRUGS_title, 'COMPLEMENT',
                                     file_name)))
            else:
                logger.info(
                    'deseq_out_multi_9:\t{}'.format(
                        os.path.join(PROJECT_title, DRUGS_title, file_name)))
                logger.info(
                    'deseq_out_multi_9:\t{}'.format(
                        os.path.join(PROJECT_title, DRUGS_title, 'COMPLEMENT',
                                     file_name)))
    DF_DRUG_combi = DF_DRUG_combi_temp
    DF_complement = DF_DRUG_combi_temp_COMPL

    # # reduce to the cols of interest and set index to case_id for later join
    # # with the raw counttable
    # ####################
    DF_DRUG_combi = DF_DRUG_combi.loc[:, ["vital_status", "gender", "PROJECT",
                                          "drugnames",
                                          "case_id"]].set_index('case_id',
                                                                drop=False)
    DF_complement = DF_complement.loc[:, ["vital_status", "gender", "PROJECT",
                                          "drugnames",
                                          "case_id"]].set_index('case_id',
                                                                drop=False)

    # #### at this point, we have the drug filtered summarytable and the
    # complement of it, also create of both the male an female variants.
    # be aware that not necessarily every gender is present in the tables...
    DF_DRUG_combi_male = DF_DRUG_combi[DF_DRUG_combi['gender'] == 'male']
    DF_DRUG_combi_female = DF_DRUG_combi[DF_DRUG_combi['gender'] == 'female']

    DF_complement_male = DF_complement[DF_complement['gender'] == 'male']
    DF_complement_female = DF_complement[DF_complement['gender'] == 'female']

    # ########################################################################

    # load the summary_DF.tsv there are the raw counts and get the case_id
    # needed via a merge
    def create_count_tables(DF_input):
        '''
        join every applied table with the raw count table, summary_DF.tsv ->
        reduce to the cases which are filtered in the combi tables
        '''
        if isinstance(PROJECT, str):
            summary_DF = pd.read_csv(os.path.join(
                OUTPUT_PATH, PROJECT_title, 'summary_DF.tsv'), sep='\t',
                index_col='genes')
        else:
            first = True
            for project in PROJECT:
                if first:
                    summary_DF = pd.read_csv(os.path.join(
                        OUTPUT_PATH, project, 'summary_DF.tsv'), sep='\t',
                        index_col='genes')
                    first = False
                else:
                    summary_DF = pd.concat(
                        [summary_DF, pd.read_csv(
                            os.path.join(
                                OUTPUT_PATH, project, 'summary_DF.tsv'),
                            sep='\t', index_col='genes')], axis=1)

        # for joining with case_id field transpose summary_DF temporarily:
        summary_DF = summary_DF.transpose()
        summary_DF.index.set_names('case_id', inplace=True)
        DF_joined = pd.merge(DF_input, summary_DF, how='inner',
                             left_index=True, right_index=True)
        DF_joined.set_index('vital_status', inplace=True)
        DF_joined = DF_joined.transpose()
        # print("\nDF_joined:\n", DF_joined)
        return DF_joined

    def call_DESEQ(key, DF_joined):
        '''
        key: provide a suffix for deseq output creation
        DF_joined: DF with ENSG and count data to be analysed by DESeq2
        '''
        # ###############
        DF_joined = DF_joined.T.reset_index().T
        # in case of a multiproject, we need to sort newly:
        if not isinstance(PROJECT, str):
            DF_joined = DF_joined.T.sort_values(
                by=['vital_status', 'drugnames', 'PROJECT', 'gender']).T
        # divide info rows from actual counts:
        DF_info = DF_joined.iloc[0:5, ]
        DF_joined = DF_joined.iloc[5:-1, ]
        # in case the key is 'DRUG_combi' or 'complement', the gender shall be
        # unique = 2, otherwise we produce the same results in f.e. DRUG_combi
        # and DRUG_combi_female in TCGA-OV
        if key == 'DRUG_combi' or key == 'complement':
            if len(DF_info.loc['gender'].unique()) == 1:
                return
        # there have to be at least 2 entries for every vital_status in DF_info
        # and also 2 factors in col vital_status
        for i in DF_info.T.value_counts('vital_status'):
            if i < 2:
                print('not enough entries for DESeq2 run: {}'.format(
                    DF_info.T.value_counts('vital_status')))
                return
        if len(DF_info.T.value_counts('vital_status')) == 1:
            print('just one factor present: {}'.format(
                DF_info.T.value_counts('vital_status')))
            print('no DESeq2 run can be performed')
            return

        count_table_log = key + '_summary_dead_alive.tsv'
        count_table = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                   count_table_log)
        info_table_log = key + '_summary_dead_alive_INFO.tsv'
        info_table = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                  info_table_log)
        deseq_output = os.path.join(OUTPUT_PATH, PROJECT_title,
                                    DRUGS_title, 'DESeq2_out_' + key)
        deseq_output_log = 'DESeq2_out_' + key
        if isinstance(PROJECT, str):
            logger.info(
                'deseq_out_single_7:\t{}'.format(
                    os.path.join(PROJECT_title, DRUGS_title, count_table_log)))
            logger.info(
                'deseq_out_single_7:\t{}'.format(
                    os.path.join(PROJECT_title, DRUGS_title, info_table_log)))
        else:
            logger.info(
                'deseq_out_multi_9:\t{}'.format(
                    os.path.join(PROJECT_title, DRUGS_title, count_table_log)))
            logger.info(
                'deseq_out_multi_9:\t{}'.format(
                    os.path.join(PROJECT_title, DRUGS_title, info_table_log)))

        # now change the ENSG description with help of table
        # ENSG_symbol.tsv to gene symbol, if present
        # DF_ENSG_symbol = pd.read_csv(
        #     os.path.join(SCRIPT_PATH, 'ENSG_symbol.tsv'),
        #     sep='\t', header=None)
        # DF_ENSG_symbol.iloc[:, 1].fillna(DF_ENSG_symbol.iloc[:, 0],
        #                                  # inplace=True)
        # DF_ENSG_symbol = DF_ENSG_symbol.iloc[:-1, :]
        # DF_joined.set_index(DF_ENSG_symbol.iloc[:, 1], inplace=True)
        # # call the script
        # check length in counttable, when DRUG is not found, this table is
        # empty
        # and would give an error in rscript:

        if DF_joined.shape[1] != 0:
            # with introducing complement also spaces and all sorts of chars
            # are introduced:
            # Note: levels of factors in the design contain characters other
            # than letters, numbers, '_' and '.'. It is recommended (but not
            # required) to use only letters, numbers, and delimiters '_' or
            # '.', as these are safe characters for column names in R. [This is
            # a message, not a warning or an error]
            # this leads to problems in the deseq design->
            # # invalid class “DESeqDataSet” object: levels of factors in the
            # design have non-unique level names after make.names() is
            # applied. best to only use letters and numbers for levels of
            # factors in the design
            # -> replace every critical char with with '_'
            # -> also cut str at specific length:
            # -> do that just for complement selection:
            # it's better to set every complement drug to the same name for
            # DESeq2, s.t. not for every single name a design is described:
            if re.search('.*complement', key):
                DF_info.loc['drugnames', :] = 'complement_drug'
                # DF_info.loc['drugnames', :] = DF_info.loc[
                #     'drugnames', :].str.replace(' ', '_', regex=False)
                # DF_info.loc['drugnames', :] = DF_info.loc[
                #     'drugnames', :].str.replace('  ', '_', regex=False)
                # DF_info.loc['drugnames', :] = DF_info.loc[
                #     'drugnames', :].str.replace('/', '_', regex=False)
                # DF_info.loc['drugnames', :] = DF_info.loc[
                #     'drugnames', :].str.replace('[', '_', regex=False)
                # DF_info.loc['drugnames', :] = DF_info.loc[
                #     'drugnames', :].str.replace(']', '_', regex=False)
                # DF_info.loc['drugnames', :] = DF_info.loc[
                #     'drugnames', :].str.replace('(', '_', regex=False)
                # DF_info.loc['drugnames', :] = DF_info.loc[
                #     'drugnames', :].str.replace(')', '_', regex=False)
                # DF_info.loc['drugnames', :] = DF_info.loc[
                #     'drugnames', :].str.replace('-', '_', regex=False)
                # DF_info.loc['drugnames', :] = DF_info.loc[
                #     'drugnames', :].str.replace('+', '_', regex=False)
                # DF_info.loc['drugnames', :] = DF_info.loc[
                #     'drugnames', :].str.replace(',', '_', regex=False)
                # DF_info.loc['drugnames', :] = DF_info.loc[
                #     'drugnames', :].str.replace('.', '_', regex=False)
                # DF_info.loc['drugnames', :] = DF_info.loc[
                #     'drugnames', :].str.replace('..', '_', regex=False)
                # -> also cut str at specific length:
            # max_length = int(8)
            # length_filter = DF_info.loc[
            #     'drugnames', :].str.len() > max_length
            # DF_info.loc[
            #     'drugnames', length_filter] = DF_info.loc[
            #         'drugnames', length_filter].apply(
            #             lambda x: x[:max_length])
            DF_joined.to_csv(count_table, sep='\t')
            DF_info.to_csv(info_table, sep='\t')
            os.makedirs(deseq_output, exist_ok=True)
            print("\ncount_table:\n", count_table)
            print("\ninfo_table:\n", info_table)
            R_path = os.path.join(SCRIPT_PATH, 'DESeq2_diffexp_multifactor.R')
            # R_path = os.path.join(OUTPUT_PATH,
            # 'DESeq2_diffexp_single_cancer.R')
            print("now calling: \n", R_path, "\nand\n", count_table,
                  "\t(arg[1])\nand\n")
            # to invoke the r version given by the conda env, use Rscript:
            Rscript = 'Rscript'
            sequence = [Rscript, R_path, count_table, info_table, deseq_output,
                        PROJECT_title]
            # print("\nsequence:\n", sequence)
            # print("\ndeseq_output:\n", deseq_output)
            subprocess.check_call(sequence)

            # ### log at this point every R output:
            if isinstance(PROJECT, str):
                for res_file in \
                        glob.glob(os.path.join(deseq_output, 'DESeq2_*')):
                    logger.info('deseq_out_single_7:\t{}'.format(
                        os.path.join(PROJECT_title, DRUGS_title,
                                     deseq_output_log,
                                     os.path.basename(res_file))))
            else:
                for res_file in \
                        glob.glob(os.path.join(deseq_output, 'DESeq2_*')):
                    logger.info('deseq_out_multi_9:\t{}'.format(
                        os.path.join(PROJECT_title, DRUGS_title,
                                     deseq_output_log,
                                     os.path.basename(res_file))))

    # ############################################## end of def call_DESEQ

    DF_DRUG_combi = DF_DRUG_combi.sort_values(by=['vital_status', 'gender'])
    DF_DRUG_combi_female = DF_DRUG_combi_female.sort_values(by=['vital_status',
                                                                'gender'])
    DF_DRUG_combi_male = DF_DRUG_combi_male.sort_values(by=['vital_status',
                                                            'gender'])
    DF_DRUG_complement = DF_complement.sort_values(by=['vital_status',
                                                       'gender'])
    DF_DRUG_complement_male = DF_complement_male.sort_values(
        by=['vital_status', 'gender'])
    DF_DRUG_complement_female = DF_complement_female.sort_values(
        by=['vital_status', 'gender'])

    DF_hash = {'DRUG_combi_complement': DF_DRUG_complement,
               'DRUG_combi_male_complement': DF_DRUG_complement_male,
               'DRUG_combi_female_complement': DF_DRUG_complement_female,
               'DRUG_combi': DF_DRUG_combi,
               'DRUG_combi_female': DF_DRUG_combi_female,
               'DRUG_combi_male': DF_DRUG_combi_male}

    for key in DF_hash:
        DF_hash[key] = create_count_tables(DF_hash[key])

    for key in DF_hash:
        call_DESEQ(key, DF_hash[key])


# fct 8 (single proj) and 10 (multi proj)in main
def create_statistics_from_DESeq2_tables(
        OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT, logger, DRUGS_title):
    '''
    :param: OUTPUT_PATH: path for DESeq2 pipeline outputs
    :type: OUTPUT_PATH: str
    :param: DRUGS: applied drug(s)
    :type: DRUGS: list of str
    :param: SCRIPT_PATH: path to the DESeq2_pipeline repo
    :type: SCRIPT_PATH: str
    :param: PROJECT: list of projects chosen
    :type: PROJECT: list of str
    :param: logger: the adjustet logger with the right filehandler
    :type: logger: logging instance

    adding some statistics to the result output tables from DESeq2
    table is saved in
    OUTPUT_PATH/PROJECT_title/DRUGS_title/DESeq2_out***/results_statistics.tsv

    .. _function_8:
    .. _function_10:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 8 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 8

        # when choosing multiple projects, call:
        $ python main_deseq.py -p TCGA-CESC  -p TCGA-HNSC -d cisplatin -f 10

    '''
    # depending on the datatype of project the tablepath is build:

    # DRUGS_title = '_'.join(sorted(map(str.lower, DRUGS)))
    file_name_list = []
    table_name_list = []
    if isinstance(PROJECT, dict):  # a multi project is handed over as dict
        PROJECT_title = []
        for project in PROJECT:
            PROJECT_title.append(project)
        PROJECT_title = '_'.join(sorted(PROJECT))
    else:
        PROJECT_title = PROJECT

    # ########## START of DRUG frequencies
    # #### create here the the prject specific
    # DF_3t_both_with_DRUG_combi_frequency.pdf/tsv
    # and the project aggregated
    # DRUG_combi_frequency_all.tsv
    # DRUG_combi_frequency_all_sorted.tsv
    # DRUG_combi_frequency_heatmap.pdf
    if isinstance(PROJECT, str):
        path_to_meta = os.path.join(OUTPUT_PATH, PROJECT,
                                    'meta_info_druglist_merged_drugs_combined')
        path_to_freq = os.path.join(OUTPUT_PATH, PROJECT,
                                    'DF_3t_both_with_DRUG_combi_frequency')
        log_file = os.path.join(PROJECT,
                                'DF_3t_both_with_DRUG_combi_frequency')
        meta_DF = pd.read_csv(path_to_meta + '.tsv', sep='\t')
        freq_table = meta_DF['pharmaceutical_therapy_drug_name'].value_counts()
        freq_table.name = PROJECT
        try:
            freq_table.to_csv(path_to_freq + '.tsv', sep='\t')
            logger.info('deseq_stats_single_8:\t{}'.format(log_file + '.tsv'))
        except Exception as e:
            print(f'Exception occured: {e}, could not save {path_to_freq}')
            print('exiting')
            os._exit(0)
        # dataset consists of 3 input datas, drugnames (S.index), frequencies
        # (s.values) and y_pos -> ticks as long as Series len(S)
        freq_table = freq_table.sort_values()
        drug_freq = freq_table.values
        drug_names = freq_table.index
        y_pos = np.arange(len(freq_table))
        plt.cla()
        plt.clf()
        plt.figure(figsize=(7, 5))
        plt.yticks(y_pos, drug_names)
        plt.barh(y_pos, drug_freq, label='Drugfrequencies')
        plt.grid(axis='x')
        plt.title(f'Drug frequencies for {PROJECT}')
        plt.legend(['drug_names'], loc=4)
        plt.tight_layout()
        try:
            plt.savefig(path_to_freq + '.pdf')
            logger.info('deseq_stats_single_8:\t{}'.format(log_file + '.pdf'))
        except Exception as e:
            print(f'Exception occured: {e}, could not save {path_to_freq}')
            print('exiting')
            os._exit(0)
    else:
        temp_freq_DF = pd.DataFrame()
        for project in PROJECT:
            temp_freq_DF = pd.concat([temp_freq_DF, pd.read_csv(os.path.join(
                OUTPUT_PATH, project,
                'DF_3t_both_with_DRUG_combi_frequency.tsv'),
                sep='\t').set_index('Unnamed: 0')], axis=1)
        # save the drug frequency table for the aggregated prjects:
        table_path = os.path.join(OUTPUT_PATH, PROJECT_title,
                                  'DRUG_combi_frequency.tsv')
        log_file = os.path.join(PROJECT_title, 'DRUG_combi_frequency.tsv')
        try:
            temp_freq_DF.index.name = ''
            temp_freq_DF.to_csv(table_path, sep='\t')
            logger.info('deseq_stats_multi_10:\t{}'.format(log_file))
        except Exception as e:
            print(f'Exception occured: {e}, could not save {table_path}')
            print('exiting')
            os._exit(0)
        # save the heatmap of the aggregated projects.
        heatmap_path = os.path.join(OUTPUT_PATH, PROJECT_title,
                                    'DRUG_combi_frequency_heatmap.pdf')
        log_file = os.path.join(PROJECT_title,
                                'DRUG_combi_frequency_heatmap.pdf')
        y = temp_freq_DF.index.tolist()
        x = temp_freq_DF.columns.tolist()
        fig = px.imshow(temp_freq_DF, x=x, y=y,
                        color_continuous_scale='Viridis', aspect="auto")
        fig.update_xaxes(side="top")
        projects = [project for project in PROJECT]
        projects = ', '.join(projects)
        fig.update_layout(
            title_text=f'available drug combinations for: {projects}')
        try:
            fig.write_image(heatmap_path)
            logger.info('deseq_stats_multi_10:\t{}'.format(log_file))
        except Exception as e:
            print(f'Exception occured: {e}, could not save {heatmap_path}')
            print('exiting')
            os._exit(0)
    # ########## END of DRUG frequencies

    glob_path = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                             'DESeq2_out*')

    for path in glob.glob(glob_path):
        table_name_list.append(os.path.join(path, 'DESeq2_results.tsv'))
        file_name_list.append(
            os.path.join(path, 'results_statistics.tsv'))
    if len(table_name_list) == 0:
        return

    # print("\ntable_name_list:\n", table_name_list)
    # print("\nfile_name_list:\n", file_name_list)
    for table, file_name in zip(table_name_list, file_name_list):
        file_log = file_name.replace(OUTPUT_PATH + os.path.sep, '')
        DESeq_res_DF = pd.read_csv(table, sep='\t')
        # print("\nDESeq_res_DF:\n", DESeq_res_DF)
        # abs log2folgchange > 1; pvalue < 0.05; padj < 0.05
        filt_log = abs(DESeq_res_DF['log2FoldChange']) > 1
        filt_pval = abs(DESeq_res_DF['pvalue']) < 0.05
        filt_adj_pval = abs(DESeq_res_DF['padj']) < 0.05
        S_log = filt_log.value_counts()
        S_pval = filt_pval.value_counts()
        S_adj_pval = filt_adj_pval.value_counts()
        S_log.name = 'abs log2FoldChange > 1'
        S_pval.name = 'pvalue < 0.05'
        S_adj_pval.name = 'padj < 0.05'
        DF_final = pd.concat([S_log, S_pval, S_adj_pval], axis=1)
        # adding project name to the single project table s.t. they are
        # distinguishable in the report file:
        if isinstance(PROJECT, str):
            DF_final.index.name = PROJECT

        if isinstance(PROJECT, dict):
            MF_name = []
            for i in PROJECT:
                MF_name.append(i.split('-')[1])
            MF_name = '_'.join(MF_name)
            DF_final.index.name = MF_name
        print("\nDF_final:\n", DF_final)
        DF_final.to_csv(file_name, sep='\t')
        if isinstance(PROJECT, str):
            logger.info('deseq_stats_single_8:\t{}'.format(file_log))
        else:
            logger.info('deseq_stats_multi_10:\t{}'.format(file_log))

        DF_res = DESeq_res_DF.rename({'Unnamed: 0': 'ENSG_name'}, axis=1)
        DF_res = DF_res.set_index('ENSG_name')
        DF_res_sorted_IN = DF_res.sort_values(
            by='log2FoldChange').dropna()
        DF_res_sorted_IN = DF_res_sorted_IN[DF_res_sorted_IN['pvalue'] < 0.05]
        # print("\nDF_res_sorted_IN, filtered pvalue:\n", DF_res_sorted_IN)
        DF_res_sorted_IN = DF_res_sorted_IN[DF_res_sorted_IN['padj'] < 0.05]
        # print("\nDF_res_sorted_IN, filtered padj:\n", DF_res_sorted_IN)
        index_INCREASE = DF_res_sorted_IN.nsmallest(
            n=10, columns='log2FoldChange').index
        index_DECREASE = DF_res_sorted_IN.nlargest(
            n=10, columns='log2FoldChange').index
        # print('ENSG list for PROJECT_title: ', self.PROJECT_title)
        ENSG_list = {'INCREASE': list(index_INCREASE),
                     'DECREASE': list(index_DECREASE)}
        print("\nENSG_list:\n", ENSG_list)
        # save those 10 highest and 10 smallest ENSG, together with all info
        # out of the deseq result table plus mygene stuff
        # the ENSG name is available throug the countfiles, get the genenames
        # from the annotation table rather with the mygene lib:
        DF_annot = pd.read_csv(
            os.path.join(
                OUTPUT_PATH, 'gencode.v36.annotation.gtf.gz'),
            sep='\t', skiprows=5, header=None)
        DF_annot = DF_annot[(DF_annot[2] == 'gene')]
        # filtering out the ENSG and create a new col for it:
        DF_annot['ENSG'] = DF_annot[8].apply(
            lambda x: re.search(r'ENSG\d+', x).group() if re.search(r'ENSG\d+',
                                                                    x) else x)
        DF_annot.set_index('ENSG', inplace=True)
        DF_annot[8].apply(lambda x: re.search('".*"', x.split(';')[1]).group())
        DF_annot[8].apply(
            lambda x: re.search('".*"', x.split(';')[1]).group().strip('"'))
        # filtering out the gene_type (called 'type_of_gene')
        DF_annot['type_of_gene'] = DF_annot[8].apply(
            lambda x: re.search('".*"', x.split(';')[1]).group().strip('"'))
        # for a better representation in the REPORT.pdf, replace _ with -,
        # s.t. Carriage returns are introduced within the tables
        DF_annot['type_of_gene'] = DF_annot[
            'type_of_gene'].str.replace('_', '-')
        # filtering out the 'symbol' -> gene symbol in annot file called
        # gene_name must be named 'symbol'
        DF_annot['symbol'] = DF_annot[8].apply(
            lambda x: re.search('".*"', x.split(';')[3]).group().strip('"'))
        ginfo_DF_merged = pd.merge(
            DF_res, DF_annot, left_index=True, right_index=True)
        for ENSG_key in ENSG_list:  # ENSG_key denotes INCREASE or DECREASE
            # print("\nginfo_DF_merged:\n", ginfo_DF_merged)
            file_name_ENSG_genes = file_name.replace(
                'statistics.tsv',
                ENSG_key + '_' + '10_lgfch_ENSG_and_Gene_info.tsv')
            print('saving the ENSG with Gene Info at: {}'.format(
                file_name_ENSG_genes))
            ginfo_DF_merged.loc[ENSG_list[ENSG_key], :].to_csv(
                file_name_ENSG_genes, sep='\t')
            file_log = file_name_ENSG_genes.replace(
                OUTPUT_PATH + os.path.sep, '')
            if isinstance(PROJECT, str):
                logger.info('deseq_stats_single_8:\t{}'.format(file_log))
            else:
                logger.info('deseq_stats_multi_10:\t{}'.format(file_log))
        # additionally save the complete results with the annotations:
        file_log = os.path.join(os.path.split(file_name)[0],
                                'DESeq2_results_with_annotation.tsv')
        file_log = file_log.replace(OUTPUT_PATH + os.path.sep, '')
        ginfo_DF_merged.index.name = 'ENSG'
        ginfo_DF_merged.to_csv(os.path.join(OUTPUT_PATH, file_log), sep='\t')
        if isinstance(PROJECT, str):
            logger.info('deseq_stats_single_8:\t{}'.format(file_log))
        else:
            logger.info('deseq_stats_multi_10:\t{}'.format(file_log))


def set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title):
    '''
    :param: OUTPUT_PATH: path for DESeq2 pipeline outputs
    :type: OUTPUT_PATH: str
    :param: PROJECT_title: merged project title out of multiple projects
    :type: PROJECT_title: str
    :param: DRUGS_title: merged drug title out of multiple drugs
    :type: DRUGS_title: str

    every paths and options are set, configure here the logfiles, with which
    the snakemake config files are going to be created we create 2 loggers, in
    case just one project is applied the logs are written in
    PROJECT/DRUGS_title/test_log.log in case multi project is applied, the logs
    are written in PROJECT_title/DRUGS_title/test_log.log with that, it is
    clear which config file shall be createt out of the logfiles present in one
    outputpath (the drugs path must therefore be created from the first fct, to
    write the log file also, the dir of the logfile must be logged, s.t.
    snakemake knows where the input file for the final snakemake configuratioin
    file is located
    '''

    # in that dir)
    #   setting a level
    #   adding a filehandler
    #   adding a formatter to the filehandler
    # logging.basicConfig(
    # filename=os.path.join(
    # OUTPUT_PATH,
    # PROJECT_title,
    # DRUGS_title,
    # 'logfile.tsv'),
    # level=logging.INFO,
    # format='%(asctime)s\t%(levelname)s\t%(message)s')

    logger = logging.getLogger('logger')
    if logger.hasHandlers():
        logger.removeHandler(logger.handlers[0])
    formatter = logging.Formatter(
        '%(asctime)s\t%(levelname)s\t%(name)s\t%(message)s\t')
    handler = logging.FileHandler(
        os.path.join(
            OUTPUT_PATH, PROJECT_title, DRUGS_title, 'test_log.log'))
    handler.setFormatter(formatter)
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    return logger
    # if(logger_single.hasHandlers()):
    # logger_single.removeHandler(logger_single.handlers[0])
    # handler = logging.FileHandler(os.path.join(
    # OUTPUT_PATH, PROJECT, DRUGS_title,'test_log.log'))
    # handler.setFormatter(formatter)
    # logger_single.addHandler(handler)

    # # logger_multi = logging.getLogger('multi_logger')
    # # logger_multi.setLevel(logging.INFO)


# fct 15 in main
def create_snake_config(OUTPUT_PATH, PROJECT_title, DRUGS_title, PROJECT_list,
                        DRUGS, SCRIPT_PATH, cutoff):
    '''
    :param: OUTPUT_PATH: path for DESeq2 pipeline outputs
    :type: OUTPUT_PATH: str
    :param: PROJECT_title: concatenated str of multiple projects
    :type: PROJECT_title: str
    :param: DRUGS_title: concatenated str of multiple drugs
    :type: DRUGS_title: str
    :param: PROJECT_list: list of applied projects
    :type: PROJECT_list: list of str
    :param: DRUGS: applied drug(s)
    :type: DRUGS: list of str
    :param: SCRIPT_PATH: path to the DESeq2_pipeline repo
    :type: SCRIPT_PATH: str
    :param: cutoff: path to the DESeq2_pipeline repo
    :type: cutoff: float

    out of the log files in PROJECT_title/DRUGS_title/test_log.log
    parse out all outputfiles of the applied run and create
    PROJECT_title/DRUGS_title/snakemake_config.yaml


    .. _function_15:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 15 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 15
    '''

    # create a list to all log files:
    file_name_list = []
    # if just one project is applied, the PROJECT_title is the same as
    # PROJECT_list[0] and is the only project to add:
    if len(PROJECT_list) == 1:
        conf_file_logger = os.path.join(PROJECT_title, DRUGS_title,
                                        'snakemake_config.yaml')
        conf_file = os.path.join(
            OUTPUT_PATH, PROJECT_title, DRUGS_title, 'snakemake_config.yaml')
        file_name_list.append(
            os.path.join(
                OUTPUT_PATH, PROJECT_title, DRUGS_title, 'test_log.log'))
    # # #else walk through single project list and add the combination of them
    # as # # last DF to the DF_list
    else:
        conf_file_logger = os.path.join(PROJECT_title, DRUGS_title,
                                        'snakemake_config.yaml')
        conf_file = os.path.join(
            OUTPUT_PATH, PROJECT_title, DRUGS_title, 'snakemake_config.yaml')
        for project in PROJECT_list:
            file_name_list.append(os.path.join(OUTPUT_PATH, project,
                                               DRUGS_title,
                                               'test_log.log'))
        file_name_list.append(os.path.join(OUTPUT_PATH, PROJECT_title,
                                           DRUGS_title, 'test_log.log'))
    # print("\nfile_name_list:\n", file_name_list)
    # the location for the conffile has to be added, before the config_file is
    # actually written:
    logger = set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title)
    logger.info('config_file_15:\t{}'.format(conf_file_logger))

    temp_dict = {}
    for log_file in file_name_list:
        with open(log_file, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    break
                else:
                    line_list = line.strip().split(sep='\t')
                    key = line_list[3]
                    value = line_list[4]
                    if key not in temp_dict:
                        temp_dict.update({key: [value]})
                    else:
                        temp_dict[key].append(value)
                        # repetitions from multiple runs can be cleared here:
                        temp_dict[key] = list(set(temp_dict[key]))
    # print("\ntemp_dict:\n", temp_dict)
    # also add DRUGS_title and PROJECTS_title to conf:
    temp_dict.update({'DRUGS:': DRUGS})
    temp_dict.update({'PROJECTS:': PROJECT_list})
    temp_dict.update({'cutoff:': [str(cutoff)]})

    # now create the config.yaml:
    # key:
    #   - value
    # if one project is applied, config_file is placed in project/drugs, else
    # in project_title/drugs
    with open(conf_file, 'w') as f:
        for key in temp_dict:
            f.write(key + '\n')
            for value in temp_dict[key]:
                f.write('  - ' + value + '\n')

    print('\na snakemake configuration file has been created ')
    print('with the following projects:\n')
    print(PROJECT_list)
    # additionally cp that file with the description in the name in the script
    # dir:
    name_copied = os.path.join(SCRIPT_PATH, os.path.os.pardir, 'Snakes',
                               'snakemake_config_' + PROJECT_title +
                               '_' + DRUGS_title + '.yaml')

    os.system('cp ' + conf_file + ' ' + name_copied)
    print('\na copy of the Snakemake configuration is created in your ' +
          'SCRIPT_PATH:\n\n{}'.format(name_copied))


def create_log_for_A(OUTPUT_PATH, PROJECT_list, PROJECT_title, DRUGS_title):
    '''
    :param: OUTPUT_PATH: path for DESeq2 pipeline outputs
    :type: OUTPUT_PATH: str
    :param: PROJECT_list: list of applied projects
    :type: PROJECT_list: list of str
    :param: PROJECT_title: concatenated str of multiple projects
    :type: PROJECT_title: str
    :param: DRUGS_title: concatenated str of multiple drugs
    :type: DRUGS_title: str

    in case the -A option is set without the -D, the log file has to be createt
    in advance. With that, it is simultaneously checked if the files needed
    (downloaded datafiles, etc) are actually already present
    '''
    os.makedirs(os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title),
                exist_ok=True)
    logger = set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title)

    for PROJECT in PROJECT_list:
        # check for the meta_info.dat in PROJECT/
        if os.path.exists(os.path.join(OUTPUT_PATH, PROJECT, 'meta_info.dat')):
            logger.info(
                'meta_info_1:\t{}'.format(
                    os.path.join(PROJECT, 'meta_info.dat')))
        else:
            print('could not find {}, have you performed the download steps\
                  with -D option in advance?'.format(os.path.join(
                      OUTPUT_PATH, PROJECT, 'meta_info.dat')))
            os._exit(0)
        #
        # checking the log_file not possible, it's not created yet, it's ok to
        # just create it here, since the only info it will hold here is that it
        # is created in the new DRUGS_title dir, which is true...
        os.makedirs(os.path.join(
            OUTPUT_PATH, PROJECT, DRUGS_title), exist_ok=True)
        logger.info(
            'log_file:\t{}'.format(
                os.path.join(PROJECT, DRUGS_title, 'test_log.log')))

        # check for 'nationwidechildrens.org_clinical_drug_cesc.txt'
        # and the 'nationwidechildrens.org_clinical_patient_cesc.txt'
        # in PROJECT/ dir:
        infix = PROJECT.replace('TCGA-', '').lower()
        drug_table_name = os.path.join(
            OUTPUT_PATH, PROJECT,
            'nationwidechildrens.org_clinical_drug_' + infix + '.txt')
        if os.path.exists(drug_table_name):
            logger.info(
                'meta_tables_2:\t{}'.format(
                    drug_table_name.replace(OUTPUT_PATH + os.path.sep, '')))
        else:
            print('could not find {}, have you performed the download steps\
                  with -D option in advance?'.format(drug_table_name))
            os._exit(0)
        patient_table_name = os.path.join(
            OUTPUT_PATH, PROJECT,
            'nationwidechildrens.org_clinical_patient_' + infix + '.txt')
        if os.path.exists(patient_table_name):
            logger.info(
                'meta_tables_2:\t{}'.format(
                    patient_table_name.replace(OUTPUT_PATH + os.path.sep, '')))
        else:
            print('could not find {}, have you performed the download steps\
                  with -D option in advance?'.format(patient_table_name))
            os._exit(0)

        # check for sep_down files, first check existence of dir (OUTPUT_PATH,
        # PROJECT, 'TCGA-****_data_files'), then glob every file out of it:
        data_dir = os.path.join(OUTPUT_PATH, PROJECT, PROJECT + '_data_files')
        if os.path.exists(data_dir):
            for data_file in glob.glob(os.path.join(data_dir,
                                                    '*' + os.path.sep + '*')):
                logger.info('manifest_3:\t{}'.format(
                    data_file.replace(OUTPUT_PATH + os.path.sep, '')))
        else:
            print('could not find {}, have you performed the download steps\
                  with -D option in advance?'.format(data_dir))
            os._exit(0)
        # check for MANIFEST.txt:
        if os.path.exists(os.path.join(OUTPUT_PATH, PROJECT, 'MANIFEST.txt')):
            logger.info('manifest_3:\t{}'.format(os.path.join(
                PROJECT, 'MANIFEST.txt')))
        else:
            print('could not find {}, have you performed the download steps\
                  with -D option in advance?'.format(os.path.join(
                PROJECT, 'MANIFEST.txt')))
            os._exit(0)

        # check for 'DF_3t_left.tsv' and 'DF_3t_right.tsv' in PROJECT, they are
        # not crucial for the consecutive analyses, but if present, are logged:
        for table_iter in ['DF_3t_left.tsv', 'DF_3t_right.tsv']:
            if os.path.exists(os.path.join(OUTPUT_PATH, PROJECT, table_iter)):
                logger.info('merged_tables_4:\t{}'.format(os.path.join(
                    PROJECT, table_iter)))

        # check 'DF_3t_both.tsv' for in PROJECT
        if os.path.exists(os.path.join(OUTPUT_PATH, PROJECT,
                                       'DF_3t_both.tsv')):
            logger.info('merged_tables_4:\t{}'.format(os.path.join(
                PROJECT, 'DF_3t_both.tsv')))
        else:
            print('could not find {}, have you performed the download steps\
                  with -D option in advance?'.format(os.path.join(
                PROJECT, 'DF_3t_both.tsv')))
            os._exit(0)

        # check for summary_DF.tsv in PROJECT/
        if os.path.exists(os.path.join(OUTPUT_PATH, PROJECT,
                                       'summary_DF.tsv')):
            logger.info('summary_tables_5:\t{}'.format(os.path.join(
                PROJECT, 'summary_DF.tsv')))
        else:
            print('could not find {}, have you performed the download steps\
                  with -D option in advance?'. format(os.path.join(
                PROJECT, 'summary_DF.tsv')))
            os._exit(0)

        # check for ['DF_3t_both_with_DRUG_combi.tsv',
        # 'DF_3t_both_with_DRUG_combi_frequency.pdf',
        # 'DF_3t_both_with_DRUG_combi_frequency.tsv', ]
        for table_iter in ['DF_3t_both_with_DRUG_combi.tsv',
                           'DF_3t_both_with_DRUG_combi_frequency.pdf',
                           'DF_3t_both_with_DRUG_combi_frequency.tsv']:
            if os.path.exists(os.path.join(OUTPUT_PATH, PROJECT, table_iter)):
                logger.info('DRUG_combi_6:\t{}'.format(os.path.join(
                    PROJECT, table_iter)))
            else:
                print('could not find {}, have you performed the download\
                      steps with -D option in advance?'. format(os.path.join(
                    PROJECT, table_iter)))
                os._exit(0)


def snake_meta(PROJECT, OUTPUT_PATH, DRUGS_title, SCRIPT_PATH,
               PROJECT_title):
    '''
    :param: PROJECT: list of projects chosen
    :type: PROJECT: list of str
    :param: FILE_TYPE: type of raw data to download from TCGA
    :type: FILE_TYPE: str
    :param: OUTPUT_PATH: path for DESeq2 pipeline outputs
    :type: OUTPUT_PATH: str
    :param: DRUGS_title: concatenated str of multiple drugs
    :type: DRUGS_title: str
    :param: SCRIPT_PATH: path to the DESeq2_pipeline repo
    :type: SCRIPT_PATH: str
    :param: PROJECT_title: concatenated str of multiple projects
    :type: PROJECT_title: str

    copy the respective
    SCRIPT_PATH/Snakes/meta_infos/PROJECT_title/DRUGS_title*/meta_info.dat
    into the actual OUTPUT_PATH/PROJECT/ path
    '''
    # its ok to glob the DRUGS_title_cutoff value, the meta_info.dat ist not
    # affected by cutoff
    temp_dir = PROJECT_title + '_' + DRUGS_title + '*'
    meta_source = glob.glob(os.path.join(SCRIPT_PATH, 'Snakes', 'meta_infos',
                                         temp_dir,
                                         PROJECT + '_meta_info.dat'))[0]
    dest = os.path.join(OUTPUT_PATH, PROJECT, 'meta_info.dat')
    shutil.copy2(meta_source, dest)


def apply_list(value):
    # print(drug)
    if re.search("palixtaxel", value):
        return "paclitaxel"
    if re.search("premetrexed", value):
        return "pemetrexed"
    if re.search("ironotecan", value):
        return "irinotecan"
    if (re.search("5-flurouracil", value) or
        re.search("5-fu", value) or re.search("5fu", value) or
        re.search("5 fu", value) or
        re.search("fluorouracil", value) or
        re.search("fluorouracil", value) or
            re.search("fluoruracil", value)):
        return "5-fluorouracil"
    if re.search("cisplatinum", value):
        return "cisplatin"
    if re.search("vinorelbin", value):
        return "vinorelbine"

    # ##harmonize synonyms:
    if re.search("cisplatin-xrt", value):
        return "cisplatin"
    if re.search("carboplatinum", value):
        return "carboplatin"
    if re.search("taxotere", value):
        return "docetaxel"
    if re.search("taxol", value) or re.search("abraxane", value):
        return "paclitaxel"
    if re.search("navelbine", value):
        return "vinorelbine"
    if re.search("cetuximab", value):
        return "erbitux"
    if re.search("erlotinib", value):
        return "tarceva"
    if re.search("panitumumab", value):
        return "vectibix"
    if re.search("paraplatin", value):
        return "carboplatin"
    if re.search("vepesid", value):
        return "etoposide"
    if re.search("gemzar", value):
        return "gemcitabine"
    if re.search("ironotecan", value):
        return "irinotecan"
    if re.search("alimta", value):
        return "pemetrexed"
    # for BRCA merge arimidex and anastrozole:
    if re.search("anastrozole", value):
        return "arimidex"
    else:
        return value
# ###############end of new automated api section
