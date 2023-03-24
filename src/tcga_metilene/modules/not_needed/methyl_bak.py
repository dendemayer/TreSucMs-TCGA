#!/usr/bin/env python3.7

# import json
# import subprocess       # to call R
# import tarfile
from matplotlib import pyplot as plt
from natsort import natsorted
import glob
import logging
import numpy as np  # for aranging REF at x
import os
import pandas as pd
import re
import requests
import seaborn as sns
import shutil
import json
import tarfile


# fct 1 in main
# this will be held available in the repo:
def download_GDC_manifest(PROJECT, FILE_TYPE, OUTPUT_PATH, DRUGS_title,
                          SCRIPT_PATH, PROJECT_title, snakerun, api_manifest):
    """
    new api change:
        # this fct is not needed anymore, the api manifest will be provided
        # within the pipeline repos

    :param: PROJECT: list of projects chosen
    :type: PROJECT: list of str
    :param: FILE_TYPE: type of raw data to download from TCGA
    :type: FILE_TYPE: str
    :param: OUTPUT_PATH: path for metilene pipeline outputs
    :type: OUTPUT_PATH: str
    :param: DRUGS_title: merged drug title out of multiple drugs
    :type: DRUGS_title: str
    :param: SCRIPT_PATH: path to the metilene_pipeline repo
    :type: SCRIPT_PATH: str
    :param: PROJECT_title: merged project title out of multiple projects
    :type: PROJECT_title: str
    :param: snakerun: the adjustet logger with the right filehandler
    :type: logger: logging instance

    creating the meta_info.dat file in your OUTPUT_PATH/PROJECT directory

    .. _function_1:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 1 to your call,
        # example:
        $ python main_metilene.py -p TCGA-CESC -d cisplatin -f 1

    """
    if snakerun:
        set_logger.snake_meta(PROJECT, FILE_TYPE, OUTPUT_PATH, DRUGS_title,
                              SCRIPT_PATH, PROJECT_title)
        return
    # download Data Release 31.0
    # https://docs.gdc.cancer.gov/Data/Release_Notes/gdc_manifest_20211029_data_release_31.0_active.tsv.gz
    url_gdc_manifest = ('https://docs.gdc.cancer.gov/Data/Release_' +
                        'Notes/' + api_manifest)
    print(f"\ngdc_manifest from:\n{url_gdc_manifest}\nis downloaded" +
          f" to\t{OUTPUT_PATH}")
    filename_manifest = api_manifest
    dest_manifest = os.path.join(OUTPUT_PATH, filename_manifest)
    if os.path.exists(dest_manifest):
        print(f'{dest_manifest} already loaded')
        return
    # download the data:
    data = requests.get(url_gdc_manifest)
    logger = set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
    with open(dest_manifest, 'wb') as f:
        try:
            f.write(data.content)
            logger.info(f'download_GDC_manifest_1:\t{filename_manifest}')
        except BaseException:
            print('problems while downloading meta_info.dat')
            print('is https://portal.gdc.cancer.gov/ online?')
            print('exiting the program')
            os._exit(0)


# fct 2 in main
# those aux files are downloaded with help of
# gdc_manifest_20211029_data_release_31.0_active.tsv.gz
def download_clinical_tables(PROJECT, OUTPUT_PATH, SCRIPT_PATH, DRUGS_title,
                             api_manifest):
    """
    api change:
        clinical tables according to the projects are loaded in
        OUTPUT_PATH/PROJECT/aux_files/*:
        nationwidechildrens.org_auxiliary_cesc.txt
        nationwidechildrens.org_biospecimen_aliquot_cesc.txt
        nationwidechildrens.org_biospecimen_analyte_cesc.txt
        nationwidechildrens.org_biospecimen_diagnostic_slides_cesc.txt
        nationwidechildrens.org_biospecimen_portion_cesc.txt
        nationwidechildrens.org_biospecimen_protocol_cesc.txt
        nationwidechildrens.org_biospecimen_sample_cesc.txt
        nationwidechildrens.org_biospecimen_shipment_portion_cesc.txt
        nationwidechildrens.org_biospecimen_slide_cesc.txt
        nationwidechildrens.org_clinical_drug_cesc.txt
        nationwidechildrens.org_clinical_follow_up_v2.0_cesc.txt
        nationwidechildrens.org_clinical_follow_up_v4.0_cesc.txt
        nationwidechildrens.org_clinical_follow_up_v4.0_nte_cesc.txt
        nationwidechildrens.org_clinical_nte_cesc.txt
        nationwidechildrens.org_clinical_omf_v4.0_cesc.txt
        nationwidechildrens.org_clinical_patient_cesc.txt
        nationwidechildrens.org_clinical_radiation_cesc.txt
        nationwidechildrens.org_ssf_normal_controls_cesc.txt
        nationwidechildrens.org_ssf_tumor_samples_cesc.txt
    :param: UUID: unique file identifier of the meta table
    :type: UUID: str
    :param: PROJECT: list of projects chosen
    :type: PROJECT: list of str
    :param: OUTPUT_PATH: path for metilene pipeline outputs
    :type: OUTPUT_PATH: str
    :param: DRUGS_title: merged drug title out of multiple drugs
    :type: DRUGS_title: str

    with the UUID the clinical tables will be downloaded in the
    OUTPUT_PATH/PROJECT:

        * nationwidechildrens.org_clinical_patient_****.txt
        * nationwidechildrens.org_clinical_drug_****.txt

    .. _function_2:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 2 to your call,
        # example:
        $ python main_metilene.py -p TCGA-CESC -d cisplatin -f 2
    """
    # if not os.path.exists(os.path.join(OUTPUT_PATH, PROJECT)):
    #     print("{} dir not existent, call meta_filter fct before".format(
    #         os.path.join(OUTPUT_PATH, PROJECT)))
    #     os._exit(0)

    # get the UUID of all patient related meta tables out of the gdc_manifest:
    # read in the table:

    filename_manifest = api_manifest
    gdc_manifest_DF = pd.read_csv(
        os.path.join(
            SCRIPT_PATH, filename_manifest),
        sep='\t').set_index('project_id').loc[PROJECT, :].reset_index()
    # the additional metadata tables begin with childre, filter them and list
    # the belonging UUID's:

    # the list holds UUID to filename:
    # array([['07998a19-ef5d-4377-8121-6527fc736ec6',
    #         'nationwidechildrens.org_biospecimen_diagnostic_slides_cesc.txt'],
    #     ['120a9a8c-adad-4b13-8773-d28b7b841276',
    #         'nationwidechildrens.org_auxiliary_cesc.txt'], ....
    UUID_list = gdc_manifest_DF[
        gdc_manifest_DF[
            'filename'].str.contains(
                'children.*.txt', regex=True)]['id'].to_numpy().tolist()

    data_endpt = "https://api.gdc.cancer.gov/data/"
    # data_endpt = "https://api.gdc.cancer.gov/data/{}".format(UUID_list)
    params = {"ids": UUID_list}

    response = requests.post(data_endpt, data=json.dumps(params),
                             headers={"Content-Type": "application/json"})

    # The file name can be found in the header within the
    # Content-Disposition key.
    try:
        response_head_cd = response.headers["Content-Disposition"]
    except KeyError:
        print('TCGA GDC dataportal not reachable at the moment\
            check the website: https://portal.gdc.cancer.gov/\
            for potentially network maintenance news and try again later')
        os._exit(0)

    logger = set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
    file_logger = re.findall("filename=(.+)", response_head_cd)[0]
    file_name = os.path.join(OUTPUT_PATH, PROJECT, file_logger)

    print("\nadditional meta info table:\n", file_name)
    logger.info('download_clinical_tables_2:\t{}'.format(
        os.path.join(PROJECT, file_logger)))

    with open(file_name, "wb") as output_file:
        output_file.write(response.content)

    # unpack all auxiliary files and put them into PROJECT/aux_files/ dir:
    # file_name holds the just downloaded tar.gz, composed out of date and
    # other numerical flags, ex:gdc_download_20220524_114836.301313.tar.gz
    # this filename is generic and changes with every download, do not log it!
    # instead log the content from this archive.
    # after extracting, the generic .gz file can be deleted
    os.makedirs(os.path.join(OUTPUT_PATH, PROJECT, 'aux_files'), exist_ok=True)
    tar = tarfile.open(file_name, 'r:gz')
    tar.extractall(os.path.join(OUTPUT_PATH, PROJECT, 'aux_files'))
    tar.close()
    list_aux_tables = glob.glob(os.path.join(OUTPUT_PATH, PROJECT,
                                             'aux_files', '*', '*'))
    # now mv every table one dir above and delete the remaining empty dir:
    for aux_table in list_aux_tables:
        src = aux_table
        dest = os.path.join(
            os.path.split(
                aux_table)[0], os.path.pardir, os.path.split(aux_table)[1])
        os.replace(src, dest)
    compl_ls = os.listdir(os.path.join(OUTPUT_PATH, PROJECT, 'aux_files'))
    for element in compl_ls:
        dir_full_path = os.path.join(OUTPUT_PATH, PROJECT, 'aux_files',
                                     element)
        if os.path.isdir(dir_full_path):
            os.removedirs(dir_full_path)
    # log every file within aux_files:
    logger = set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
    for aux_file in os.listdir(os.path.join(
            OUTPUT_PATH, PROJECT, 'aux_files')):
        logger.info('download_clinical_tables_2\t{}'.format(
            os.path.join(PROJECT, 'aux_files', aux_file)))
    os.remove(file_name)


# #####################################
# trying to download the methyl files separately depending on the UUID given
# in the druglist tabel:


# fct 3 in main
def create_merged_metatable(OUTPUT_PATH, PROJECT, DRUGS_title, workflow_type,
                            SCRIPT_PATH):
    """
    :param: OUTPUT_PATH: path for pipeline outputs
    :type: OUTPUT_PATH: str
    :param: PROJECT: list of projects chosen
    :type: PROJECT: list of str
    :param: DRUGS_title: merged drug title out of multiple drugs
    :type: DRUGS_title: str
    :workflow_type: either methylation data or rna-count data
    :type: workflow_type: str
    :param: SCRIPT_PATH: path to the pipeline script
    :type: SCRIPT_PATH: str

    meta_info.dat and nationwidechildren are linked through:
    meta_info.dat -> associated_entities.0.case_id
    nationwide... -> bcr_patient_uuid
    whats needed is col "id" out of meta_info.dat for sep
    downloadeta_info.dat
    and nationwidechildren are linked through:

    .. _function_3:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 3 to your call,
        # example:
        $ python main_metilene.py -p TCGA-CESC -d cisplatin -f 3
    """
    path_aux = os.path.join(OUTPUT_PATH, PROJECT, 'aux_files')

    # aliquot table contains info about which file is used,
    # to be combined with manifest TCGA-CESC.tsv (mani_DF)
    aliq_table_path = glob.glob(path_aux + '/*aliquot*')[0]
    aliq_DF = pd.read_csv(aliq_table_path, sep='\t').loc[2:, :]
    aliq_DF['bcr_patient_uuid'] = aliq_DF['bcr_patient_uuid'].str.lower()
    aliq_DF['bcr_aliquot_uuid'] = aliq_DF['bcr_aliquot_uuid'].str.lower()

    # ## merge aliquot ids from TCGA-CESC.tsv (mani_DF) with
    # nationwidechildrens.org_biospecimen_aliquot_cesc
    manifest_file = os.path.join(SCRIPT_PATH, 'GCv36_Manifests',
                                 f'{PROJECT}.tsv')
    mani_DF = pd.read_csv(manifest_file, sep='\t')  # join on 'aliquot_ids'
    # this mani_DF can be filtered right away on the workflow_type:
    workflow_filter = mani_DF['workflow_type'] == workflow_type
    mani_DF = mani_DF[workflow_filter]
    mani_DF.rename({'aliquot_ids': 'bcr_aliquot_uuid'}, axis=1, inplace=True)
    mani_DF['bcr_aliquot_uuid'] = mani_DF['bcr_aliquot_uuid'].str.lower()

    # drug_table contains info which therapy is used
    # nationwidechildrens.org_clinical_drug_cesc.txt
    drug_table_path = glob.glob(path_aux + '/*drug*')[0]
    drug_DF = pd.read_csv(drug_table_path, sep='\t').loc[2:, :]
    drug_DF['bcr_patient_uuid'] = drug_DF['bcr_patient_uuid'].str.lower()

    # sample table contains info about tumor site:
    sample_table_path = glob.glob(path_aux + '/*biospecimen_sample*')[0]
    sample_DF = pd.read_csv(sample_table_path, sep='\t').loc[2:, :]
    sample_DF['bcr_sample_uuid'] = sample_DF['bcr_sample_uuid'].str.lower()
    sample_DF['bcr_patient_uuid'] = sample_DF['bcr_patient_uuid'].str.lower()
    # constrain samples to use on primary tumor samples:
    primary_tumor_filter = sample_DF['sample_type'] == 'Primary Tumor'
    sample_DF = sample_DF[primary_tumor_filter]

    # follow_up table (vital_DF) contains info like vital status and survival
    # times:
    # nationwidechildrens.org_clinical_follow_up_v4.0_cesc
    vital_table_path = glob.glob(path_aux + '/*clinical_follow_up_v4*')[0]
    vital_DF = pd.read_csv(vital_table_path, sep='\t').loc[2:, :]
    vital_DF['bcr_patient_uuid'] = vital_DF['bcr_patient_uuid'].str.lower()

    # nationwidechildrens.org_clinical_patient_cesc.txt'
    # patient_DF contains the gender info
    patient_table_path = glob.glob(path_aux + '/*patient*')[0]
    patient_DF = pd.read_csv(patient_table_path, sep='\t').loc[2:, :]
    patient_DF['bcr_patient_uuid'] = patient_DF['bcr_patient_uuid'].str.lower()
    # it occurs that the col naming for age_at_diagnosis varies, check that:
    # nationwidechildrens.org_clinical_patient_cesc_columns -> age_at_diagnosis
    # nationwidechildrens.org_clinical_patient_hnsc_columns ->
    # age_at_initial_pathologic_diagnosis
    # -> harmonize the col to age_at_diagnosis
    col_str = ' '.join(patient_DF.columns.to_list())
    age_at_dia = re.findall('age_at.*?_diagnosis', col_str)[0]
    patient_DF.rename({age_at_dia: 'age_at_diagnosis'}, axis=1, inplace=True)

    # [312 rows x 28 columns]:
    aliq_mani_DF = pd.merge(aliq_DF, mani_DF, on='bcr_aliquot_uuid')
    # ## END merge aliquot ids from TCGA-CESC.tsv with nation...aliquot..tsv ##

    # ######################## concat the drugnames in drug_DF:
    # this can be done before any merging steps with other tables:
    # the drugs are in col: pharmaceutical_therapy_drug_name
    # -> nupis must be converted to a single appearance, the drugs are
    # concantenatet via comma, correct misspelled drugnames:
    drug_DF['pharmaceutical_therapy_drug_name'] = drug_DF[
        'pharmaceutical_therapy_drug_name'].str.lower()
    # those nupis are from the drug table!
    not_uniq_filt = drug_DF['bcr_patient_barcode'].value_counts() > 1
    not_uniq_patient_index = not_uniq_filt[not_uniq_filt].index
    # out of the multiple drugs, drop duplicates, make a sorted list:
    drug_DF = drug_DF.set_index('bcr_patient_barcode', drop=False)
    for nupi in not_uniq_patient_index:
        sorted_drug_list = sorted(drug_DF.loc[nupi, :][
                'pharmaceutical_therapy_drug_name'].drop_duplicates(
                ).to_list())
        # correct the drugs within the list:
        print(f'sorted_drug_list before correction:\n{sorted_drug_list}')
        for index, value in enumerate(sorted_drug_list):
            sorted_drug_list[index] = apply_list(value)
        print(f'sorted_drug_list after correction:\n{sorted_drug_list}\n')
        # now the duplicates for the nupi can be droped and the remaining
        # value for the drug can be the concatenated sorted, corrected drug
        # string:
        # from: https://stackoverflow.com/questions/
        # 13035764/remove-pandas-rows-with-duplicate-indices :
        # df3 = df3[~df3.index.duplicated(keep='first')]
        # drug_DF.loc[nupi, :] = drug_DF.loc[nupi, :][~drug_DF.loc[nupi,
        # :].index.duplicated(keep='first')]
        # -> the "request" on the index asks 'which one is duplicated?', since
        # the first occurence is not 'yet' duplicated, this is False, every
        # followed repetition is a duplication and gives True, the complement
        # of this list accesses the first entry:
        # (Pdb) drug_DF.loc[nupi, :].index.duplicated()
        # array([False,  True,  True,  True,  True,  True,  True])
        # (Pdb) ~drug_DF.loc[nupi, :].index.duplicated()
        # array([ True, False, False, False, False, False, False])
        temp_nupi_row = drug_DF.loc[
            nupi, :][~drug_DF.loc[nupi, :].index.duplicated()]
        temp_nupi_row[
            'pharmaceutical_therapy_drug_name'] = ','.join(sorted_drug_list)
        # now this multiple index can be deleted and then concatenated with the
        # temp row.
        drug_DF.drop(nupi, inplace=True)
        drug_DF = pd.concat([drug_DF, temp_nupi_row])

    # combined drug therapies are now joined in one row per patient, temp save
    # that table:
    drug_out_path = os.path.join(OUTPUT_PATH, PROJECT,
                                 'meta_info_druglist_merged.tsv')
    drug_DF.to_csv(drug_out_path, sep='\t', index=False)
    logger = set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
    logger.info('create_merged_metatable_3:\t{}'.format(
        os.path.join(PROJECT, 'meta_info_druglist_merged.tsv')))
    print(f'saved meta_info_druglist_merged.tsv to {drug_out_path}')

    # # ## merge  aliq_mani_DF with drugs_DF on bcr_patient_uuid
    aliq_mani_drug_DF = pd.merge(aliq_mani_DF, drug_DF, on='bcr_patient_uuid')

    aliq_mani_drug_sample_DF = pd.merge(aliq_mani_drug_DF, sample_DF,
                                        on='bcr_patient_uuid')
    # join everything else on col 'bcr_patient_uuid
    # [141 rows x 122 columns]:
    complete_DF = aliq_mani_drug_sample_DF.merge(vital_DF,
                                                 on='bcr_patient_uuid')
    # [144 rows x 224 columns]
    complete_DF = aliq_mani_drug_sample_DF.merge(patient_DF,
                                                 on='bcr_patient_uuid')
    # lowercase the vital_status and gender infos:
    complete_DF['vital_status'] = complete_DF['vital_status'].str.lower()
    complete_DF['gender'] = complete_DF['gender'].str.lower()
    complete_path = os.path.join(
        OUTPUT_PATH, PROJECT, 'meta_info_druglist_merged_drugs_combined.tsv')
    # in rare cases duplicates arose, they lead to the same filename and can
    # be dropped here:
    complete_DF.drop_duplicates(subset='bcr_patient_uuid', inplace=True)
    complete_DF.to_csv(complete_path, sep='\t', index=False)
    logger.info('create_merged_metatable_3:\t{}'.format(
        os.path.join(PROJECT, 'meta_info_druglist_merged_drugs_combined.tsv')))
    print('saved meta_info_druglist_merged_drugs_combined.tsv ')
    print(f'in {complete_path}')


# fct 4 in main
def sep_down_data_files(OUTPUT_PATH, PROJECT, DRUGS_title):
    """
    :param: OUTPUT_PATH: path for metilene pipeline outputs
    :type: OUTPUT_PATH: str
    :param: PROJECT: list of projects chosen
    :type: PROJECT: list of str
    :param: DRUGS_title: merged drug title out of multiple drugs
    :type: DRUGS_title: str

    with help of 'id' col of 'meta_info_druglist_merged_drugs_combined.tsv',

    .. _function_4:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 4 to your call,
        # example:
        $ python main_metilene.py -p TCGA-CESC -d cisplatin -f 4
    """
    # # if not os.path.exists(os.path.join(OUTPUT_PATH, PROJECT)):
    #     # print(
    #         # "{} dir not existent, call meta_filter() -f 1,".format(
    #             # os.path.join(OUTPUT_PATH, PROJECT)),
    #         # end='')
    #     # print("download_clinical_drug_table and create_merged_metatable fct
    #     ")
    #     # os._exit(0)
    # # open the allready collapsed meta_info_merged:

    logger = set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
    OUTPUT_PATH_temp = os.path.join(OUTPUT_PATH, PROJECT)
    try:
        merged_df = pd.read_csv(os.path.join(
            OUTPUT_PATH_temp, 'meta_info_druglist_merged_drugs_combined.tsv'),
            sep='\t')
    except FileNotFoundError as e:
        print(e)
        print("run fct create_merged_metatable(OUTPUT_PATH, PROJECT)")
        os._exit(1)
    else:
        ids_to_download = list(merged_df['id'])
        # print("\nids_to_download:\n", ids_to_download)
        temp_dir = PROJECT + '_data_files'
        temp_dir = os.path.join(OUTPUT_PATH_temp, temp_dir)
        print("\ndownload in:\n", temp_dir)
        # TODO if dir allready exist request the large download!
        # # if os.path.exists(temp_dir):
        #     # # this step is expensive in time and storage, if the data dir
        #     # # already
        #     # # exist, return
        #     # print("subdir already exists, NO new datafiles are downloaded")
        #     # return
        # else:
        os.makedirs(temp_dir, exist_ok=True)
        for id_ in ids_to_download:
            data_endpt = "https://api.gdc.cancer.gov/data/{}".format(id_)
            try:
                response = requests.get(
                    data_endpt, headers={"Content-Type": "application/json"})
            except ConnectionError:
                print('Remote end closed connection without response')
                print(' check the website: https://portal.gdc.cancer.gov/')
                print(' for potentially network maintenance news and try ')
                print(' again later')
                os._exit(0)
            # The file name can be found in the header within the
            # Content-Disposition key.
            try:
                response_head_cd = response.headers["Content-Disposition"]
            except KeyError:
                print('TCGA GDC dataportal not reachable at the moment')
                print(' check the website: https://portal.gdc.cancer.gov/')
                print(' for potentially network maintenance news and try ')
                print(' again later')
                os._exit(0)
            file_name = re.findall("filename=(.+)", response_head_cd)[0]
            file_name = os.path.join(temp_dir, file_name)
            print("\nfile_name:\n", os.path.basename(file_name))
            try:
                with open(file_name, "wb") as output_file:
                    output_file.write(response.content)
            except ConnectionError:
                print('Remote end closed connection without response')
                print(' check the website: https://portal.gdc.cancer.gov/')
                print(' for potentially network maintenance news and try ')
                print(' again later')
                os._exit(0)
        # now log every downloaded methyl data file:
        for data_file in glob.glob(os.path.join(temp_dir, '*')):
            data_file = data_file.replace(OUTPUT_PATH + os.path.sep, '')
            logger.info('sep_down_data_files_4:\t{}'.format(
                os.path.join(data_file)))
        # #######################################check for consistancy

# #########################


# fct 5 in main
def create_summary_table(OUTPUT_PATH, PROJECT, DRUGS_title):
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

    since the MANIFEST is not available anymore (api request through separate
    download reveals no MANIFEST) the information which methylation file
    belongs to which patient comes from "bcr_patient_barcode" col out of
    meta_info_druglist_merged.tsv (before merged out of meta_info.dat +
    nationwidechildrens.org_clinical_drug_cesc.txt) the case_id can also be
    taken from meta_info_druglist_merged (col 'bcr_patient_uuid' is the
    case_id)

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
    open_files_ids = list(merged_df['bcr_patient_barcode'])

    dir_list_2 = {}
    # to relate back filename to bcr_patient_barcode (without extra parsing)
    # keep them as hash:
    # # now grep the filenames out of the dir_list, ignore hidden files:
    # and check for consistency in table dimensions!
    for file_name in dir_list:
        for of in open_files_ids:
            if re.search(of, file_name) and not file_name.startswith('.'):
                # print('found {}'.format(file_name))
                dir_list_2.update({file_name: of})
    # print("\ndir_list_2:\n", dir_list_2)
    # build up the whole DF starting with first el of dir_list_2
    # print("\ndir_list_2:\n", dir_list_2)
    summary_DF = pd.DataFrame()
    count = 0
    for i in dir_list_2:  # i is the id of the file
        if count == 0:
            # print('first iteration', i)
            summary_DF = pd.read_csv(os.path.join(
                data_dir, i), sep='\t')
            # important: check shape: some of the tables are inconsistent!
            if summary_DF.shape != (485577, 11):
                print('shape is {}, skipping table {}'.format(summary_DF.shape,
                                                              i))
                continue
            summary_DF = summary_DF.loc[:, [
                'Composite Element REF', 'Chromosome', 'Start', 'Beta_value']]
            print(
                'starting with table:\n{}\nto complete_summary.tsv'.format(i))
            # rename Beta_value with vital_status;case_id;drug_combi;PROJECT
            # goto the belonging patient id out of dir_list_2 hash
            filt = merged_df['bcr_patient_barcode'] == dir_list_2[i]
            # DONE include here the gender info!
            new_col_list = merged_df[filt].loc[:, [
                'vital_status',
                'case_ids',
                'pharmaceutical_therapy_drug_name',
                'gender']].values[0]
            new_col_list = np.append(new_col_list, PROJECT)
            # print(new_col_list)
            new_col_list = ';'.join(new_col_list)
            # print("\nnew_col_list:\n", new_col_list)
            summary_DF.rename(columns={'Beta_value': new_col_list},
                              inplace=True)
            # print("\nsummary_DF:\n", summary_DF)
            count = count + 1
        else:
            print('adding table:\n{}\nto complete_summary.tsv'.format(i))
            # merge here at col Composite Element REF
            temp_DF = pd.read_csv(
                os.path.join(data_dir, i), sep='\t')
            # important: check shape: some of the tables are inconsistent!
            if temp_DF.shape != (485577, 11):
                print('shape is {}, skipping table {}'.format(
                    summary_DF.shape, i))
                continue
            temp_DF = temp_DF.loc[:, [
                'Composite Element REF', 'Beta_value']]
            filt = merged_df['bcr_patient_barcode'] == dir_list_2[i]
            # DONE include here the gender info!
            new_col_list = merged_df[filt].loc[:, [
                'vital_status',
                'case_ids',
                'pharmaceutical_therapy_drug_name',
                'gender']].values[0]
            new_col_list = np.append(new_col_list, PROJECT)
            # print(new_col_list)
            new_col_list = ';'.join(new_col_list)
            temp_DF.rename(
                columns={'Beta_value': new_col_list}, inplace=True)
            # merge is by default a inner join
            summary_DF = pd.merge(
                summary_DF, temp_DF, on='Composite Element REF')
            # print("\nsummary_DF:\n", summary_DF)
            count = count + 1

    # print("\nsummary_DF:\n", summary_DF)
    logger = set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
    summary_DF.to_csv(
        os.path.join(
            OUTPUT_PATH, PROJECT, 'complete_summary.tsv'),
        sep='\t', index=False, na_rep='NA')
    logger.info('create_summary_table_5:\t{}'.format(
        os.path.join(PROJECT, 'complete_summary.tsv')))


# fct. 6 in main
def create_table_for_metilene(OUTPUT_PATH, PROJECT, DRUGS, DRUGS_title,
                              met_dir, cutoff):
    """
    :param: OUTPUT_PATH: path for metilene pipeline outputs
    :type: OUTPUT_PATH: str
    :param: PROJECT: list of projects chosen
    :type: PROJECT: list of str
    :param: DRUGS: applied drug(s)
    :type: DRUGS: list of str
    :param: DRUGS_title: merged drug title out of multiple drugs
    :type: DRUGS_title: str
    :param: met_dir: results from the metilene run are saved there
    :type: met_dir: str
    :param: cutoff: convert dead cases to alive cases if they outlive the\
        cutoff parameter (in years)
    :type: cutoff: float

    creates a filtered summary table depending on the drugs
    header structure like:

    vital_status;case_id;druglist;gender;PROJECT

    file for metilene input:

    summary_dead_alive_dropped_NA_dropped_REF_sorted.tsv

    for each gender possibility an additional subdir is created:
        OUTPUT_PATH/PROJECT/DRUGS_title/gender

    whereby gender is either one of male, female or both

    help_file:

    dropped_REF_sorted.tsv

    changes made by cutoff parameter are invoked here:

    .. _function_6:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 6 to your call,
        # example:
        $ python main_metilene.py -p TCGA-CESC -d cisplatin -f 6
    """
    logger = set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
    # DONE: provide tables if gender and their combinations can be applied
    # read in table, right away with multiindex, REF, Chromosome and Start
    summary_DF = pd.read_csv(
        os.path.join(OUTPUT_PATH, PROJECT, 'complete_summary.tsv'), sep='\t',
        index_col=[0, 1, 2])
    # drop every row with chromosome name *
    summary_DF.drop(summary_DF.loc[(slice(None), '*'), :].index, inplace=True)
    # sort first chromosome wise, then start wise:
    # summary_DF.sort_index(level=[1, 2], inplace=True) # ## not sufficient,
    # would sort like: chr19, chr2, natsort is needed:
    # # WICHTIG, sort first Chromosome, then start: -> therefore make col
    # # 'Chromosome' categorial first, sort it then
    temp_chromosome_start_DF = summary_DF.reset_index(
        level=['Chromosome', 'Start'])
    temp_chromosome_start_DF = summary_DF.reset_index(
        level=['Chromosome', 'Start'])
    temp_chromosome_col = temp_chromosome_start_DF.loc[:, 'Chromosome']
    # prepare the categorial col, tell how the nonunique factors shall be
    # sorted:
    temp_chromosome_start_DF.loc[:, 'Chromosome'] = pd.Categorical(
        temp_chromosome_col, ordered=True, categories=natsorted(
            temp_chromosome_col.unique()))
    # do the actual sorting:
    temp_chromosome_start_DF = temp_chromosome_start_DF.sort_values(
        by=["Chromosome", "Start"], ascending=[True, True])
    # DF is correctly sorted now, set chr and start as index again:)
    summary_DF = temp_chromosome_start_DF.set_index(['Chromosome', 'Start'],
                                                    append=True)

    # drop rows which contain only NA's
    summary_DF.dropna(how='all', inplace=True)
    # build the Multiindex out of the col names:
    MI_dict = {'vital_state': [], 'UUID': [], 'DRUGS': [], 'gender': [],
               'PROJECT': []}
    col_name_list = list(summary_DF.columns)

    # ####################################################
    # DONE include changes made by cutoff parameter here:
    # #################

    # #################### start include cutoff parameter:
    # open the 'meta_info_druglist_merged_drugs_combined.tsv' table, it holds
    # the information needed to invoke the cutoff parameter:
    if cutoff > 0:
        meta_drugs_DF = pd.read_csv(
            os.path.join(
                OUTPUT_PATH, PROJECT,
                'meta_info_druglist_merged_drugs_combined.tsv'), sep='\t',
            index_col=0)
        # all alives and not numeric values can be deleted from the caselist:
        filt = meta_drugs_DF.loc[:, 'death_days_to'].str.isnumeric()
        # # filt = (meta_drugs_DF.loc[
        #     # :, 'death_days_to'] != '[Not Applicable]' or meta_drugs_DF.loc[
        #         # :, 'death_days_to'] != '[Discrepancy]')
        temp_meta = meta_drugs_DF[filt]
        # temp_meta.loc[:, [
        # 'cases.0.demographic.vital_status', 'death_days_to']]
        filt = temp_meta['vital_status'] == 'dead'
        temp_meta = temp_meta[filt]
        # make years out of last_concact_days_to col:
        temp_meta.loc[:, 'death_days_to'] = temp_meta.loc[
            :, 'death_days_to'].apply(lambda x: float(x)/365)
        meta_cases = temp_meta.index.to_numpy()
        cutoff_dict = {}
        for meta_case in meta_cases:
            for col_name in col_name_list:
                match = re.search(meta_case, col_name)
                if match:
                    # check if cases lived longer than cutoff:
                    # if temp_meta.loc[meta_case, 'death_days_to'] is a Series:
                    # ValueError: The truth value of a Series is ambiguous. Use
                    # a.empty, a.bool(), a.item(), a.any() or a .all().
                    if temp_meta.loc[meta_case, 'death_days_to'] > cutoff:
                        # include them to the alive cases->
                        col_name.split(';', 1)
                        index = col_name_list.index(col_name)
                        col_new = ';'.join(
                            ['alive', col_name.split(';', 1)[1]])
                        col_name_list[index] = col_new
                        days_to_death = temp_meta.loc[
                            meta_case, 'death_days_to']
                        cutoff_dict.update({col_new: [days_to_death, cutoff]})
        # #####################################################################
        # #### the header is changed and dead cases which outlive the cutoff
        # are converted to alive
        # #### with this, the invokation of cutoff is completed,
        # #####################################################################
        # survivaltime are changed to alive cases, log those cases:
        if len(cutoff_dict.values()) != 0:
            cutoff_add_DF = pd.DataFrame().from_dict(cutoff_dict)
            cutoff_add_DF.index = ['survivaltime', 'cutoff']
            cut_col_list = cutoff_add_DF.columns
            cut_case_list = [x.split(';')[1] for x in cut_col_list]
            temp = cutoff_add_DF.T
            temp['case_id'] = cut_case_list
            temp = temp.set_index('case_id')
            cutoff_add_DF = temp
            cutoff_add_DF = cutoff_add_DF.join(meta_drugs_DF)
            cutoff_add_DF = cutoff_add_DF.loc[
                :, [
                    'id',
                    'gender',
                    'vital_status',
                    'pharmaceutical_therapy_drug_name',
                    'survivaltime', 'last_contact_days_to',
                    # 'cases.0.demographic.year_of_birth',
                    # 'cases.0.demographic.year_of_death',
                    'age_at_diagnosis', 'cutoff']]
            cutoff_add_DF.insert(
                len(cutoff_add_DF.columns)-1, 'PROJECT', PROJECT)
            cutoff_add_DF.rename(
                columns={'id': 'UUID', 'cases.0.demographic.gender': 'gender',
                         'cases.0.demographic.vital_status': 'vital_status',
                         'pharmaceutical_therapy_drug_name': 'drugnames'},
                inplace=True)
            cutoff_add_DF['age_at_diagnosis'] = cutoff_add_DF[
                'age_at_diagnosis'].apply(lambda x: x/365)

            file_name = os.path.join(
                OUTPUT_PATH, PROJECT,
                'cutoff_cases_add_' + str(cutoff) + '.tsv')
            log_name = os.path.join(
                PROJECT, 'cutoff_cases_add_' + str(cutoff) + '.tsv')
            cutoff_add_DF.to_csv(file_name, sep='\t')
            logger.info(
                'create_table_for_metilene_6:\t{}'.format(log_name))
            # also save the DRUG filtered table:
            drug_filt = cutoff_add_DF['drugnames'].isin(DRUGS)
            cutoff_add_DF_drug = cutoff_add_DF[drug_filt]
            file_name = os.path.join(
                OUTPUT_PATH,
                PROJECT,
                DRUGS_title, 'cutoff_cases_add_' + str(cutoff) + '.tsv')
            log_name = os.path.join(
                PROJECT,
                DRUGS_title, 'cutoff_cases_add_' + str(cutoff) + '.tsv')
            cutoff_add_DF_drug.to_csv(file_name, sep='\t')
            logger.info(
                'create_table_for_metilene_6:\t{}'.format(log_name))

    # ################# end include cutoff parameter:

    for col in col_name_list:
        MI_dict['vital_state'].append(col.split(';')[0])
        MI_dict['UUID'].append(col.split(';')[1])
        MI_dict['DRUGS'].append(col.split(';')[2])
        MI_dict['gender'].append(col.split(';')[3])
        MI_dict['PROJECT'].append(col.split(';')[4])
    # the multiindex can't be build directly from the dict, first create DF out
    # of dict, than use from_frame method:
    MI_index = pd.MultiIndex.from_frame(pd.DataFrame(MI_dict))
    summary_DF = summary_DF.T.set_index(MI_index).T
    # search for col_names matching the DRUG pattern exactly:
    # NEW: also create the complement of it:
    summary_DF_complement = summary_DF.loc[
        :, ~summary_DF.columns.isin(DRUGS, level='DRUGS')]
    summary_DF = summary_DF.loc[:, summary_DF.columns.isin(DRUGS,
                                                           level='DRUGS')]
    # check for the genders
    # values, counts = np.unique(
    #     summary_DF.columns.isin(
    #         ['female'], level='gender'), return_counts=True)
    gender_list = []
    for gender in ['female', 'male']:
        values = np.unique(
            summary_DF.columns.isin([gender], level='gender'))
        if True in values:
            gender_list.append(gender)
    if len(gender_list) == 2:
        gender_list.append('both')
    # here just single projects are processed, name PROJECT_title to single
    PROJECT_title = PROJECT
    # for export, bring the header in that form:
    # Chromosome  Start   alive;aa10d6da-...;carboplatin...;female;TCGA-HNSC
    # chr1    15865   0.844660062689115
    # and also save the belonging REF in:
    # Composite Element REF
    # cg13869341
    # cg14008030
    for gender in gender_list:
        path_gen = os.path.join(
            OUTPUT_PATH, PROJECT_title, DRUGS_title, gender)
        os.makedirs(os.path.join(path_gen, met_dir), exist_ok=True)
        new_col_list = ['Chromosome', 'Start']
        # for both, write the whole table
        if gender == 'both':
            temp_REF_DF = summary_DF.reset_index(
                'Composite Element REF')['Composite Element REF']
            temp_REF_DF.to_csv(os.path.join(
                path_gen, 'dropped_REF_sorted.tsv'), sep='\t', index=False)
            # the REF are now saved to the DF_dropped_REF_sorted.tsv
            # this index can be dropped now:
            temp_DF = summary_DF.reset_index(
                ['Chromosome', 'Start']).reset_index('Composite Element REF',
                                                     drop=True)
            for col_i in range(2, len(temp_DF.columns)):
                new_col_list.append(';'.join(temp_DF.columns[col_i]))
            temp_DF.columns = new_col_list
            summary_name = \
                'summary_dead_alive_dropped_NA_dropped_REF_sorted.tsv'
            temp_DF.to_csv(os.path.join(path_gen, summary_name),
                           sep='\t', index=False)
            logger.info(
                'create_table_for_metilene_6:\t{}'.format(
                    os.path.join(
                        PROJECT_title, DRUGS_title, gender, summary_name)))
            logger.info('create_table_for_metilene_6:\t{}'.format(
                os.path.join(PROJECT_title, DRUGS_title, gender,
                             'dropped_REF_sorted.tsv')))
            # to_csv(os.path.join(path_gen, 'summary_DF'), sep='\t')
        # or write the gender limited parts:
        else:
            # summary_DF.loc[:, (slice(None), slice(None), slice(None),
            # gender)].to_csv(os.path.join(path_gen, 'summary_DF'), sep='\t')
            temp_REF_DF = summary_DF.loc[
                :, (slice(None), slice(None), slice(None), gender)
            ].reset_index(
                'Composite Element REF')['Composite Element REF']
            temp_REF_DF.to_csv(os.path.join(
                path_gen, 'dropped_REF_sorted.tsv'), sep='\t', index=False)
            # the REF are now saved to the DF_dropped_REF_sorted.tsv
            # this index can be dropped now:
            temp_DF = summary_DF.loc[
                :, (slice(None), slice(None), slice(None), gender)
            ].reset_index(
                ['Chromosome', 'Start']).reset_index(
                    'Composite Element REF', drop=True)
            for col_i in range(2, len(temp_DF.columns)):
                new_col_list.append(';'.join(temp_DF.columns[col_i]))
            temp_DF.columns = new_col_list
            temp_DF.to_csv(os.path.join(
                path_gen,
                'summary_dead_alive_dropped_NA_dropped_REF_sorted.tsv'
            ), sep='\t', index=False)
            summary_name = \
                'summary_dead_alive_dropped_NA_dropped_REF_sorted.tsv'
            logger.info('create_table_for_metilene_6:\t{}'.format(
                        os.path.join(
                            PROJECT_title, DRUGS_title, gender, summary_name)))
            logger.info('create_table_for_metilene_6:\t{}'.format(
                        os.path.join(
                            PROJECT_title,
                            DRUGS_title, gender, 'dropped_REF_sorted.tsv')))
        # TODO: how to handle NA values?


# fct 7 and 12 in main
def provide_metilene_table(OUTPUT_PATH, DRUGS, SCRIPT_PATH,
                           met_opt_list, met_dir, PROJECT, DRUGS_title):
    r'''
    taking: 'summary_dead_alive_dropped_NA_dropped_REF_sorted.tsv'
    Chromosome	Start	alive;adb7c5c8-4afd-40dc-89f1-...;cisplatin;female;
    TCGA-CESC
    chr1	15865	0.83049028997308
    in project dir, project tag is in header, not in aggregated dir

    * depending on options -m, -M, and -d create new dir like:
        * "m\_3\_-M\_1000\_-d\_0.03" and write the perl script filtered
            metilene output in that dir
    * adjust the -d option in metilene also in the filter perl script,
        otherwise those postitions are filtered out of the metilene output
    * create in the metilene dir:
        * from metilene_linux64:
            * 'met_output_sorted.txt'
        * from metilene_output.pl:
            * metilene_qval.0.05.bedgraph
            * metilene_qval.0.05.out
            * metilene_qval.0.05.pdf

    .. _function_7:
    .. _function_12:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 7 to your call,
        # example:
        $ python main_metilene.py -p TCGA-CESC -d cisplatin -f 7

        # when choosing multiple projects, call:
        $ python main_metilene.py -p TCGA-CESC  -p TCGA-HNSC -d cisplatin -f 12

    '''
    # plots for all projects:
    if isinstance(PROJECT, str):
        PROJECT_title = PROJECT
    else:
        project_list = []
        for project in PROJECT:
            project_list.append(project)
        PROJECT_title = '_'.join(sorted(map(str.upper, project_list)))

    gender_list = ['female', 'male', 'both']
    for gender in gender_list:
        # check whether the gender dir exists, else continue
        if not os.path.exists(
            os.path.join(
                OUTPUT_PATH, PROJECT_title, DRUGS_title, gender)):
            continue
        else:
            met_dir_2 = os.path.join(gender, met_dir)
            os.chdir(os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                  met_dir_2))
            # Rscript_path = os.path.join(SCRIPT_PATH, 'metilene_output.R')
            Perl_path = os.path.join(SCRIPT_PATH, 'metilene_output.pl')
            metilene_path = os.path.join(SCRIPT_PATH, 'metilene_linux64')
            # -d, --minMethDiff <n>
            # minimum mean methylation difference (default:0.100000)
            met_opt = ' ' + ' '.join(met_opt_list)
            pre_opt_metilene = met_opt + ' -a alive -b dead '
            input_path = os.path.join(
                os.pardir,
                'summary_dead_alive_dropped_NA_dropped_REF_sorted.tsv')
            # important sort flags:
            # -V, --version-sort
            #          natural sort of (version) numbers within text
            # -n, --numeric-sort
            #          compare according to string numerical value
            post_opt_metilene = ' | sort -k1,1V -k2,2n > ' + \
                'met_output_sorted.txt'
            # the R script just takes 10 cols, scip col 7 and 8: not
            # neccessary, perl calls the R script...
            # call_awk = r"awk '{print
            # $1,'\t',$2,'\t'$3,'\t',$4,'\t',$5,'\t',$6, '\t',$9,'\t',$10}'
            # met_output_sorted.txt > met_output_for_R.txt" call_R_output =
            # Rscript_path + ' met_output_for_R.txt R.out.pdf'

            call_metilene = metilene_path + pre_opt_metilene + input_path \
                + post_opt_metilene
            print("\ncall_metilene:\n", call_metilene)
            try:
                os.system(call_metilene)
            except BaseException:
                print('problems with the metilene binary')
                os._exit(0)

            logger = set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title)
            if isinstance(PROJECT, str):
                logger.info(
                    'provide_metilene_table_7:\t{}'.format(
                        os.path.join(
                            PROJECT_title, DRUGS_title, met_dir_2,
                            'met_output_sorted.txt')))
            else:
                logger.info(
                    'provide_metilene_table_12:\t{}'.format(
                        os.path.join(
                            PROJECT_title, DRUGS_title, met_dir_2,
                            'met_output_sorted.txt')))

            call_perl_output = Perl_path + ' -q ' + 'met_output_sorted.txt' \
                + ' -c 3 -a alive -b dead'

            # re.search returns a Match object, group(0) gives the entire match
            # searching for the -d option and getting then the values:
            match = re.search(r'(-d\s*\w.*\w|-d\s*\w)', met_opt)
            if match:
                call_perl_output = call_perl_output + ' ' + match.group(0)
                # print(match.group(0))
            print("\ncall_perl_output:\n", call_perl_output)
            try:
                os.system(call_perl_output)
            except BaseException:
                print('problems with calling metilene_output.pl')
                os._exit(0)
            # here we have to separate the pdf output:
            pdf_separate = 'pdfseparate metilene_qval.0.05.pdf\
                metilene_qval.0.05_%d.pdf '
            os.system(pdf_separate)
            for perl_met_out in glob.glob(
                os.path.join(
                    OUTPUT_PATH,
                    PROJECT_title, DRUGS_title, met_dir_2, 'metilene_qval.*')):
                perl_met_out = perl_met_out.replace(
                    OUTPUT_PATH + os.path.sep, '')
                if isinstance(PROJECT, str):
                    logger.info(
                        'provide_metilene_table_7:\t{}'.format(perl_met_out))
                else:
                    logger.info(
                        'provide_metilene_table_12:\t{}'.format(perl_met_out))


# fct 8 and  13 in main
def call_bed_intersect(OUTPUT_PATH, DRUGS, met_opt_list, met_dir, PROJECT,
                       DRUGS_title):
    """
    cut a new table out of the files summary_dead_alive_dropped_NA_dropped_REF
    _sorted.tsv and dropped_REF_sorted.tsv
    insert new col with 'End' position of the REF
    create intersect_reads.tsv
    save the headings in intersect_header.tsv

    .. _function_8:
    .. _function_13:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 8 to your call,
        # example:
        $ python main_metilene.py -p TCGA-CESC -d cisplatin -f 8

        # when choosing multiple projects, call:
        $ python main_metilene.py -p TCGA-CESC  -p TCGA-HNSC -d cisplatin -f 13
    """
    # change the location where the right ranges are located from metilene,
    # this is in the options dir

    # change wd s.t. metilene saves the results in OUTPUT:
    # determine whether a single or multiproject is applied by the PROJECT
    # data-type:

    if isinstance(PROJECT, str):
        PROJECT_title = PROJECT
    else:
        project_list = []
        for project in PROJECT:
            project_list.append(project)
        PROJECT_title = '_'.join(sorted(map(str.upper, project_list)))

    ######
    gender_list = ['female', 'male', 'both']
    for gender in gender_list:
        if not os.path.exists(os.path.join(OUTPUT_PATH, PROJECT_title,
                                           DRUGS_title, gender)):
            continue
        else:
            os.chdir(os.path.join(OUTPUT_PATH, PROJECT_title,
                                  DRUGS_title, gender))
            summary_DF_dead_alive_sorted = pd.read_csv(
                'summary_dead_alive_dropped_NA_dropped_REF_sorted.tsv',
                sep='\t')
            DF_dropped_REF_sorted = pd.read_csv('dropped_REF_sorted.tsv',
                                                sep='\t')
            DF_dropped_REF_sorted.columns = ['REF']
            # print("\n summary_DF_dead_alive_sorted:\n",
            # summary_DF_dead_alive_sorted)
            # print("\nDF_dropped_REF_sorted:\n", DF_dropped_REF_sorted)
            # TODO: take original End values!!!!
            # the original end values are in the raw data files, some have
            # equal start ending?!
            # summary_DF_dead_alive_sorted.insert(
            #     loc=2, column='End',
            #     value=summary_DF_dead_alive_sorted['Start'].apply(
            #         lambda x: x + 1))
            summary_DF_dead_alive_sorted.insert(
                loc=2, column='End',
                value=summary_DF_dead_alive_sorted['Start'])

            summary_DF_dead_alive_sorted.insert(
                loc=3, column='REF', value=DF_dropped_REF_sorted['REF'])
            # print("\n summary_DF_dead_alive_sorted.iloc[:, 0:4]:\n",
            #       summary_DF_dead_alive_sorted.iloc[:, 0:4])
            summary_DF_dead_alive_sorted.to_csv(
                'intersect_reads.tsv', sep='\t', index=False, header=False)
            # # for some reason on some projects a .1 is added, like
            # TCGA-HNSC.1,
            # # make sure this gets lost:
            # #
            # header_temp = []
            # for col in summary_DF_dead_alive_sorted.columns:
            #     header_temp.append(col.replace('.1', ''))
            # print("\nheader_temp:\n", header_temp)
            header_temp = summary_DF_dead_alive_sorted.columns
            with open('intersect_header.tsv', 'w') as f:
                for item in header_temp:
                    f.write("%s\t" % item)
            logger = set_logger(
                OUTPUT_PATH,
                PROJECT_title,
                DRUGS_title)
            if isinstance(PROJECT, str):
                logger.info('call_bed_intersect_8:\t{}'.format(
                    os.path.join(PROJECT_title, DRUGS_title, gender,
                                 'intersect_header.tsv')))
                logger.info('call_bed_intersect_8:\t{}'.format(
                    os.path.join(PROJECT_title, DRUGS_title, gender,
                                 'intersect_reads.tsv')))
            else:
                logger.info('call_bed_intersect_13:\t{}'.format(
                    os.path.join(PROJECT_title, DRUGS_title, gender,
                                 'intersect_header.tsv')))
                logger.info('call_bed_intersect_13:\t{}'.format(
                    os.path.join(PROJECT_title, DRUGS_title, gender,
                                 'intersect_reads.tsv')))

            # bedtools intersect -a intersect_reads.tsv -b
            # metilene_qval.0.05.out call_bed_intersect = 'bedtools intersect
            # -a intersect_reads.tsv -b metilene_qval.0.05.out >
            # intersect_out.tsv'
            call_bed_intersect = \
                'bedtools intersect -a intersect_reads.tsv -b '\
                + os.path.join(met_dir,
                               'metilene_qval.0.05.out') + ' > ' \
                + os.path.join(met_dir, 'intersect_out.tsv')

            print("\ncall_bed_intersect:\n", call_bed_intersect)
            os.system(call_bed_intersect)
            if isinstance(PROJECT, str):
                logger.info('call_bed_intersect_8:\t{}'.format(
                    os.path.join(PROJECT_title, DRUGS_title, gender, met_dir,
                                 'intersect_out.tsv')))
            else:
                logger.info('call_bed_intersect_13:\t{}'.format(
                    os.path.join(PROJECT_title, DRUGS_title, gender, met_dir,
                                 'intersect_out.tsv')))
            # ### take those metilene results of the single runs in the project
            # dirs and put them here together:


# fct 9 and 14 in main
def create_plots(OUTPUT_PATH,
                 DRUGS, met_opt_list, met_dir, PROJECT, DRUGS_title):
    """
    Creating plots for the metilene results

    * aggregated plots for every applied project
    * seperated plots for each project
        * lineplots for every value of every position
        * boxplots for each position
            * filtering out those boxplot where all alives among all projects
                are higher as deads or vice versa
            * this leads to the ranking of good or best candidates

    .. _function_9:
    .. _function_14:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 9 to your call,
        # example:
        $ python main_metilene.py -p TCGA-CESC -d cisplatin -f 9

        # when choosing multiple projects, call:
        $ python main_metilene.py -p TCGA-CESC  -p TCGA-HNSC -d cisplatin -f 14
    """
    # plots for all projects:
    if isinstance(PROJECT, str):
        PROJECT_title = PROJECT
    else:
        project_list = []
        for project in PROJECT:
            project_list.append(project)
        PROJECT_title = '_'.join(sorted(map(str.upper, project_list)))

# ################# START multiproject part
    gender_list = ['female', 'male', 'both']
    for gender in gender_list:
        if not os.path.exists(os.path.join(OUTPUT_PATH, PROJECT_title,
                                           DRUGS_title, gender)):
            continue
        else:
            os.chdir(os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                  gender))
            # aggregating for every project:
            # check which projectnames are in header. create a set:
            intersect_DF_header = list(pd.read_csv('intersect_header.tsv',
                                                   sep='\t').columns)[:-1]
            # intersect_DF = pd.read_csv(os.path.join(
            #     OUTPUT_PATH,'intersect_out.tsv'), header=None, sep='\t')
            intersect_DF = pd.read_csv(os.path.join(
                met_dir, 'intersect_out.tsv'), header=None, sep='\t')
            intersect_DF.columns = intersect_DF_header
            # print("\nintersect_DF:\n", intersect_DF)
            # print("\nintersect_DF.columns:\n", intersect_DF.columns)
            intersect_DF.set_index(
                ['Chromosome', 'Start', 'End', 'REF'], inplace=True)
            # print("\nintersect_DF:\n", intersect_DF)
            # transforming the DF for re index with colnames, first split cols
            # at ; vital_status;uuid;drugs;project
            intersect_DF = intersect_DF.T
            vital_status_list = []
            uuid_list = []
            drug_list = []
            project_list = []
            for i in intersect_DF_header:
                # print("\ni:\n", i)
                if len(i.split(';')) == 5:
                    temp_list = i.split(';')
                    vital_status_list.append(temp_list[0])
                    uuid_list.append(temp_list[1])
                    drug_list.append(temp_list[2])
                    project_list.append(temp_list[4])
            # PROBLEM here with set index!
            # print("\nintersect_DF_header:\n", intersect_DF_header)
            intersect_DF.set_index(
                [vital_status_list, uuid_list, drug_list, project_list],
                inplace=True)
            # print("\nintersect_DF:\n", intersect_DF)
            intersect_DF = intersect_DF.T
            # print("\nintersect_DF.loc[:, ('alive', 'TCGA-CESC')]:\n",
            #       intersect_DF.loc[:, ('alive', slice(None), slice(None),
            #                            ['TCGA-CESC', 'TCGA-HNSC'])])
            project_set = list(set(project_list))
            # print("\nproject_set:\n", project_set)
            # the plot shall include a dead and a alive graph for every
            # project, check if both vital states are available for all
            # projects:
            project_to_delete = []

            for vital_state in ['dead', 'alive']:
                # print('check for vital_state:', vital_state)
                for project in project_set:
                    # print("\nproject:\n", project)
                    temp_DF = intersect_DF.loc[:, (vital_state, slice(None),
                                                   slice(None), project)]
                    if temp_DF.empty:
                        print('delete', project, 'no vital_state', vital_state)
                        project_to_delete.append(project)
            if len(project_to_delete) != 0:
                project_to_delete = list(set(project_to_delete))
                # print("\nproject_to_delete:\n", project_to_delete)
                # Remove an item by value: remove()
                for element in project_to_delete:
                    project_set.remove(element)
                # print("\nproject_set:\n", project_set)
                # print("\nproject_to_delete:\n", project_to_delete)

            # plot table projec specific:
            # REF  beta_median_dead  beta_median_alive  abs_diff
            # cg12848345          0.228446           0.351432  0.122986
            # cg23202979          0.120083           0.370412  0.250329

            # first the medians:
            # print("\nintersect_DF:\n", intersect_DF)
            # stocks.groupby('Symbol').Close.mean())
            intersect_DF_median = pd.DataFrame()
            for project in project_set:
                for vital_state in ['dead', 'alive']:
                    # print("\nproject:\n", project)
                    temp_DF = intersect_DF.loc[:, (vital_state, slice(None),
                                                   slice(None), project)]
                    # print("\ntemp_DF:\n", vital_state, project, temp_DF)
                    # temp_DF is here a multiindexed series, cast to DataFrame
                    # to be able to change the colname
                    temp_DF = pd.DataFrame(temp_DF.median(axis=1))
                    new_col_name = project + '_beta-median_' + vital_state
                    temp_DF.columns = [[vital_state], [new_col_name]]
                    # print("\ntemp_DF:\n", temp_DF)
                    # concatenate median series to new DF
                    if intersect_DF_median.empty:
                        intersect_DF_median = temp_DF
                    else:
                        intersect_DF_median = pd.concat([intersect_DF_median,
                                                        temp_DF], axis=1)
            # in case all series are empty, nothing is to plot, return:
            if intersect_DF_median.empty:
                print('nothing to plot, returning')
                return

            # print("\nintersect_DF_median:\n", intersect_DF_median)
            # to refer in plot to 'REF', drop first 3 levels,
            # then make REF to a col
            intersect_DF_median.reset_index(
                level=['Chromosome', 'Start', 'End'], drop=True, inplace=True)
            intersect_DF_median.reset_index(inplace=True)
            # print("\nintersect_DF_median:\n", intersect_DF_median)
            # add up the medians for alive resp dead across the projects and
            # create a abs_diff_mean col
            intersect_DF_median[('abs_diff_median', 'all_projects')] = abs(
                intersect_DF_median.loc[:, 'dead'].sum(axis=1) -
                intersect_DF_median.loc[:, 'alive'].sum(axis=1))
            # print("\nintersect_DF_median:\n", intersect_DF_median)
            # ## plot direct with pandas:
            ax = plt.gca()
            for colnr in range(1, len(intersect_DF_median.columns)):
                # print("\nintersect_DF_median.columns[colnr]:\n",
                # intersect_DF_median.columns[colnr])
                intersect_DF_median.plot(
                    kind='line',
                    x='REF', y=intersect_DF_median.columns[colnr], ax=ax)
            file_name = os.path.join(
                os.getcwd(), met_dir, 'plot_beta_median_line_graph.pdf')
            plt.savefig(file_name)
            print('plot saved in:', file_name)
            logger = set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title)
            if isinstance(PROJECT, str):
                logger.info(
                    'create_plots_9:\t{}'.format(
                        os.path.join(
                            PROJECT_title,
                            DRUGS_title,
                            gender,
                            met_dir,
                            'plot_beta_median_line_graph.pdf')))
            else:
                logger.info(
                    'create_plots_14:\t{}'.format(
                        os.path.join(
                            PROJECT_title,
                            DRUGS_title,
                            gender,
                            met_dir,
                            'plot_beta_median_line_graph.pdf')))
            plt.clf()

            # ######## plot with matplotlib:
            # drop first level of cols, then rename remaining one level cols
            temp_DF_median = intersect_DF_median.copy()
            temp_DF_median.columns = temp_DF_median.columns.droplevel()
            temp_DF_median.rename(
                columns={
                    '': 'REF',
                    'all_projects': 'abs_diff_de_al_all_projects'},
                inplace=True)
            # print("\ntemp_DF_median:\n", temp_DF_median)
            sorted_DF = temp_DF_median.sort_values(by='REF')
            # print("\nsorted_DF:\n", sorted_DF)
            counter = 0
            DF_list = [temp_DF_median, sorted_DF]
            for temp_DF_median in DF_list:
                # create the range, where ref shall be placed at x:
                x_indexes = np.arange(len(temp_DF_median))
                # print("\nx_indexes:\n", x_indexes)
                # print("\ntemp_DF_median.iloc[:, 1:-1]:\n",
                #       temp_DF_median.iloc[:, 0].values)
                plt.figure(figsize=(14, 5))
                plt.plot(x_indexes, temp_DF_median.iloc[:, 1:])

                plt.legend(
                    list(temp_DF_median.iloc[:, 1:].columns),
                    loc='upper right', bbox_to_anchor=(1.3, 1))

                plt.xticks(
                    ticks=x_indexes, labels=temp_DF_median.iloc[:, 0].values,
                    rotation='vertical')
                plt.title("All projects beta median line graph")
                plt.xlabel("REF")
                plt.ylabel("beta value")
                # Pad margins so that markers don't get clipped by the axes
                # plt.margins(0.2, 0.2)
                # Tweak spacing to prevent clipping of tick-labels
                plt.subplots_adjust(bottom=0.25, right=0.8)
                # plt.tight_layout()
                # plt.grid(True)
                file_name = ''
                if counter == 0:
                    file_name = os.path.join(
                        os.getcwd(),
                        met_dir, 'plot_beta_median_line_graph_matplot.pdf')
                else:
                    file_name = os.path.join(
                        os.getcwd(),
                        met_dir, 'plot_beta_median_' +
                        'line_graph_matplot_REF_sort.pdf')
                counter = counter + 1
                plt.savefig(file_name)
                file_logger = file_name.replace(OUTPUT_PATH + os.path.sep, '')
                if isinstance(PROJECT, str):
                    logger.info('create_plots_9:\t{}'.format(file_logger))
                else:
                    logger.info('create_plots_14:\t{}'.format(file_logger))
                # plt.show()
                plt.clf()
                print('plot saved in:', file_name)

            # ####TODO deprecation error,
            # ###### boxplot = DF_box_plot.boxplot(by='vital_state') not ok!!!
            # # BOXPLOTS:
            # # take first 3 REF medians with highest diff:
            # # print("\nintersect_DF_median:\n", intersect_DF_median)
            # intersect_DF_median.sort_values(
            #     by=('abs_diff_median', 'all_projects'),
            #     ascending=False, inplace=True)
            # max_diff_REF = list(intersect_DF_median.iloc[:3, 0])
            # # with that, access the intersect_DF rows:
            # DF_box_plot = intersect_DF.loc[(
            #     slice(None), slice(None), slice(None), max_diff_REF), :]
            # DF_box_plot.reset_index(
            #     level=['Chromosome', 'Start', 'End'], drop=True,
            #     inplace=True)
            # DF_box_plot = DF_box_plot.T
            # # print("\nDF_box_plot:\n", DF_box_plot)
            # DF_box_plot = DF_box_plot.reset_index(level=0)
            # DF_box_plot.rename(columns={'level_0': 'vital_state'},
            # inplace=True)
            # plt.figure()
            # boxplot = DF_box_plot.boxplot(by='vital_state')
            # file_name = os.path.join(OUTPUT_PATH, met_dir,
            #                          'boxplot_beta_median_highest_diff.pdf')
            # plt.savefig(file_name)
            # print('plot saved in:', file_name)
            # plt.clf()

            # ##########BOXPLOTS separated by project
            # filter out the candidates

            # print("\nintersect_DF:\n", intersect_DF)
            # here, the Start col is still present
            print("\nproject_to_delete:\n", project_to_delete)
            print("\nproject_set:\n", project_set)
            # access the columns with projects that have both, alive and dead
            # cases:

            intersect_DF_complete = intersect_DF.loc[:, (
                ['alive', 'dead'], slice(None), slice(None), project_set)]
            intersect_DF_complete.columns.names = [
                'vital_state', 'UUID', 'drug', 'project']
            intersect_DF_complete_Position = intersect_DF_complete.copy()

            # ######################### START check the positions:
            intersect_DF_complete_Position.reset_index(
                level=['Chromosome', 'End'], drop=True, inplace=True)
            intersect_DF_complete_Position = intersect_DF_complete_Position.T
            # print("\nintersect_DF_complete_Position:\n",
            # intersect_DF_complete_Position)
            intersect_DF_complete_Position.reset_index(
                level=['UUID', 'drug'], drop=True, inplace=True)
            intersect_DF_complete_Position.reset_index(level=['vital_state',
                                                              'project'])
            # intersect_DF_complete_Position =
            # intersect_DF_complete_Position.iloc[:, :4]
            # print("\nintersect_DF_complete_Position:\n",
            # intersect_DF_complete_Position)
            count = 0
            # create new folders for candidates and not candidate boxplots:
            os.makedirs(os.path.join(os.getcwd(), met_dir, 'boxplot_rest'),
                        exist_ok=True)
            os.makedirs(os.path.join(os.getcwd(), met_dir, 'candidates'),
                        exist_ok=True)
            for column in intersect_DF_complete_Position.columns:
                count = count + 1
                # ### check if medians of all dead are higher or medians of all
                # dead are lower than alive medians
                check_median_DF = \
                    intersect_DF_complete_Position.loc[
                        :, column].groupby(
                            level=['vital_state', 'project']).median()
                alive_median_list = list(check_median_DF.loc['alive', :])
                dead_median_list = list(check_median_DF.loc['dead', :])
                all_higher = True
                all_lower = True
                for i in dead_median_list:
                    if i > min(alive_median_list):
                        all_lower = False
                    if i < max(alive_median_list):
                        all_higher = False
                if all_higher or all_lower:
                    file_name = 'boxplot_project_vital_state_POS_' + str(
                        column[0]) + '_' + str(column[1]) + '_candidate.pdf'
                    file_name = os.path.join(
                        os.getcwd(), met_dir, 'candidates', file_name)
                else:
                    file_name = 'boxplot_project_vital_state_POS_' + str(
                        column[0]) + '_' + str(column[1]) + '.pdf'
                    file_name = os.path.join(
                        os.getcwd(), met_dir, 'boxplot_rest', file_name)
                plt.figure()
                temp_DF = intersect_DF_complete_Position.reset_index()
                temp_DF_bench = temp_DF.groupby(
                    ['project', 'vital_state']).median()
                temp_DF_bench = temp_DF_bench.loc[:, (column[0], column[1])]
                median_diffs = []
                for project in project_set:
                    alive_med = temp_DF_bench.loc[(project, 'alive')]
                    dead_med = temp_DF_bench.loc[(project, 'dead')]
                    median_diffs.append(abs(alive_med - dead_med))
                projects_diff = round(float(sum(median_diffs)), ndigits=4)
                file_name = file_name.replace(
                    '.pdf', '_' + str(projects_diff) + '.pdf')
                # temp_DF.grou
                # print("\nintersect_DF_complete_Position:\n", temp_DF)
                sns.set(style="ticks", palette="bright")
                # sns.despine(offset=10, trim=True)
                plot = sns.boxplot(
                    x="project", y=column, hue="vital_state", palette=[
                        "m", "g"], data=temp_DF)

                # Seaborn box plot returns a matplotlib axes instance.
                # Unlike pyplot itself, which has a method plt.title(), the
                # corresponding argument for an axes is ax.set_title().
                # Therefore you need to call
                # sns.boxplot('Day', 'Count', data= gg).set_title('lalala')
                # Of course you could also use the returned axes instance to
                # make it more readable:
                # ax = sns.boxplot('Day', 'Count', data= gg)
                # ax.set_title('lalala')
                # ax.set_ylabel('lololo')
                print("\nsaving positional boxplots in:\n", file_name)
                plt_title = 'REF:    ' + \
                    column[1] + '\nposition: ' + str(column[0])
                plot.set_title(plt_title)
                plot.set_ylabel('beta_values')
                plot.figure.savefig(file_name)
                # also save the data this plot is based on:
                filename_table = file_name.replace('.pdf', '.tsv')
                intersect_DF_complete_Position.loc[:, (
                    column[0], column[1])].to_csv(
                    filename_table, sep='\t')
                file_logger = file_name.replace(OUTPUT_PATH + os.path.sep, '')
                file_logger_table = file_logger.replace('.pdf', '.tsv')

                if isinstance(PROJECT, str):
                    logger.info('create_plots_9:\t{}'.format(file_logger))
                    logger.info(
                        'create_plots_9:\t{}'.format(file_logger_table))
                else:
                    logger.info('create_plots_14:\t{}'.format(file_logger))
                    logger.info(
                        'create_plots_14:\t{}'.format(file_logger_table))
                plt.clf()
                # ####################### check positions END
        # #################################END of multiproject part

    # START    # plots for the single projects:
    # else:
    # OUTPUT_PATH = os.path.join(OUTPUT_PATH, PROJECT, DRUGS_title, met_dir)
    for gender in gender_list:
        if not os.path.exists(os.path.join(OUTPUT_PATH, PROJECT_title,
                                           DRUGS_title, gender)):
            continue
        else:
            os.chdir(os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                  gender, met_dir))
            print("\nos.getcwd():\n", os.getcwd())
            np.warnings.filterwarnings(
                'error', category=np.VisibleDeprecationWarning)
            # change wd s.t. metilene saves the results in OUTPUT:
            # when creating a list out of .columns, last, extra element in
            # named 'Unnamed: 106', ommitt that when asigning the colnames:
            intersect_DF_header = list(
                pd.read_csv(
                    os.path.join(
                        os.getcwd(),
                        os.pardir,
                        'intersect_header.tsv'),
                    sep='\t').columns)[
                :-1]
            # check if metilene output is empty:
            if os.stat('intersect_out.tsv').st_size == 0:
                print(
                    'intersect_out.tsv is empty, no plots for ',
                    PROJECT_title)
                return
            intersect_DF = pd.read_csv(
                'intersect_out.tsv', header=None, sep='\t')
            intersect_DF.columns = intersect_DF_header
            print("\nintersect_DF:\n", intersect_DF)

            # changing to multiindex:
            intersect_DF.set_index(
                ['Chromosome', 'Start', 'End', 'REF'],
                inplace=True)
            # print("\nintersect_DF:\n", intersect_DF)
            # transforming the DF for re index with colnames, first split cols
            # at ; vital_status;uuid;DRUGS_title
            intersect_DF = intersect_DF.T
            vital_status_list = []
            uuid_list = []
            drug_list = []
            project_list = []
            for i in intersect_DF_header:
                print("\ni:\n", i)
                if len(i.split(';')) == 5:
                    temp_list = i.split(';')
                    vital_status_list.append(temp_list[0])
                    uuid_list.append(temp_list[1])
                    drug_list.append(temp_list[2])
                    project_list.append(temp_list[4])
            print("\nvital_status_list:\n", vital_status_list)
            intersect_DF.set_index([vital_status_list, uuid_list, drug_list,
                                    project_list], inplace=True)
            intersect_DF = intersect_DF.T
            print("\nintersect_DF:\n", intersect_DF)

            # the plot shall include a dead and a alive median graph for every
            # project, check if both vital states are available:

            # ############# MEDIAN:

            intersect_DF_median = pd.DataFrame()
            for vital_state in ['dead', 'alive']:
                temp_DF = intersect_DF.loc[:, vital_state]
                if temp_DF.empty:
                    print(
                        'no vital_state',
                        vital_state,
                        'return to next Project')
                    return
                temp_DF = pd.DataFrame(temp_DF.median(axis=1))
                new_col_name = PROJECT_title + '_beta-median_' + vital_state
                temp_DF.columns = [new_col_name]
                # print("\ntemp_DF:\n", temp_DF)
                # concatenate median series to new DF
                if intersect_DF_median.empty:
                    intersect_DF_median = temp_DF
                else:
                    intersect_DF_median = pd.concat([intersect_DF_median,
                                                    temp_DF], axis=1)
            # ############# plot graph:
            # to refer in plot to 'REF', drop first 3 levels, then make REF to
            # col
            intersect_DF_median['abs_diff_median'] = abs(
                intersect_DF_median.iloc[:, 0
                                         ] - intersect_DF_median.iloc[:, 1])
            # print("\nintersect_DF_median:\n", intersect_DF_median)
            intersect_DF_median.reset_index(
                level=['Chromosome', 'Start', 'End'], drop=True, inplace=True)
            intersect_DF_median.reset_index(inplace=True)
            # print("\nintersect_DF_median:\n", intersect_DF_median)

            ax = plt.gca()
            for colnr in range(1, len(intersect_DF_median.columns)):
                # print("\nintersect_DF_median.columns[colnr]:\n",
                #       intersect_DF_median.columns[colnr])
                intersect_DF_median.plot(
                    kind='line',
                    x='REF', y=intersect_DF_median.columns[colnr], ax=ax)
            file_name = 'plot_median_diff.pdf'
            plt.savefig(file_name)
            if isinstance(PROJECT, str):
                logger.info('create_plots_9:\t{}'.format(
                    os.path.join(PROJECT_title, DRUGS_title, gender, met_dir,
                                 file_name)))
            else:
                logger.info('create_plots_14:\t{}'.format(
                    os.path.join(PROJECT_title, DRUGS_title, gender,
                                 met_dir, file_name)))
            print('plot saved in: {}: {}'.format(os.getcwd(), file_name))
            # plt.clf() clears the entire current figure with all its axes,
            # but leaves the window opened, such that it may be reused for
            # other plots.
            plt.clf()

            # BOXPLOTS:
            # take first 3 REF medians with highest diff:
            intersect_DF_median.sort_values(
                by='abs_diff_median', ascending=False, inplace=True)
            max_diff_REF = list(intersect_DF_median.iloc[:3, 0])
            # with that, access the intersect_DF rows:
            DF_box_plot = intersect_DF.loc[
                (slice(None), slice(None), slice(None), max_diff_REF), :]
            DF_box_plot.reset_index(
                level=['Chromosome', 'Start', 'End'], drop=True, inplace=True)
            DF_box_plot = DF_box_plot.T
            print("\nDF_box_plot:\n", DF_box_plot)
            DF_box_plot = DF_box_plot.reset_index(level=0)
            DF_box_plot.rename(
                columns={
                    'level_0': 'vital_state'},
                inplace=True)
            print("\nDF_box_plot:\n", DF_box_plot)
            DF_box_plot.reset_index(drop=True, inplace=True)
            print("\nDF_box_plot:\n", DF_box_plot)
            print("\nDF_box_plot.columns:\n", list(DF_box_plot.columns)[1:])
            DF_box_plot = pd.melt(
                DF_box_plot, id_vars=[
                    'vital_state'], value_vars=list(DF_box_plot.columns)[1:])
            print("\nDF_box_plot.columns:\n", list(DF_box_plot.columns))
            plt.figure()
            # boxplot = DF_box_plot.boxplot(by='vital_state', layout=(2, 3))
            plot = sns.boxplot(
                x='REF', y='value', hue='vital_state', palette=[
                    'm', 'g'], data=DF_box_plot)
            file_name = 'boxplot_beta_median_2.pdf'
            # plt.savefig(file_name)

            plot.set_title('highest median diffs of ' + DRUGS_title)
            plot.set_ylabel('beta_value')
            plot.figure.savefig(file_name)
            if isinstance(PROJECT, str):
                logger.info(
                    'create_plots_9:\t{}'.format(
                        os.path.join(
                            PROJECT_title,
                            DRUGS_title, gender,
                            met_dir,
                            file_name)))
            else:
                logger.info(
                    'create_plots_14:\t{}'.format(
                        os.path.join(
                            PROJECT_title,
                            DRUGS_title, gender,
                            met_dir,
                            file_name)))

            print('plot saved in: {}: {}'.format(os.getcwd(), file_name))
            plt.clf()

            # ############MEAN:

            intersect_DF_mean = pd.DataFrame()
            for vital_state in ['dead', 'alive']:
                temp_DF = intersect_DF.loc[:, vital_state]
                if temp_DF.empty:
                    print(
                        'no vital_state',
                        vital_state,
                        'return to next Project')
                    return
                temp_DF = pd.DataFrame(temp_DF.mean(axis=1))
                new_col_name = PROJECT_title + '_beta-mean_' + vital_state
                temp_DF.columns = [new_col_name]
                # print("\ntemp_DF:\n", temp_DF)
                # concatenate median series to new DF
                if intersect_DF_mean.empty:
                    intersect_DF_mean = temp_DF
                else:
                    intersect_DF_mean = pd.concat(
                        [intersect_DF_mean, temp_DF], axis=1)
            # ############# plot graph:
            # to refer in plot to 'REF', drop first 3 levels, then make REF to
            # col
            intersect_DF_mean['abs_diff_mean'] = abs(
                intersect_DF_mean.iloc[:, 0] - intersect_DF_mean.iloc[:, 1])
            # print("\nintersect_DF_mean:\n", intersect_DF_mean)
            intersect_DF_mean.reset_index(
                level=['Chromosome', 'Start', 'End'], drop=True, inplace=True)
            intersect_DF_mean.reset_index(inplace=True)
            # print("\nintersect_DF_mean:\n", intersect_DF_mean)
            ax = plt.gca()
            for colnr in range(1, len(intersect_DF_mean.columns)):
                # print("\nintersect_DF_mean.columns[colnr]:\n",
                #       intersect_DF_mean.columns[colnr])
                intersect_DF_mean.plot(
                    kind='line',
                    x='REF', y=intersect_DF_mean.columns[colnr], ax=ax)
            file_name = 'plot_mean_diff.pdf'
            plt.savefig(file_name)
            if isinstance(PROJECT, str):
                logger.info(
                    'create_plots_9:\t{}'.format(
                        os.path.join(
                            PROJECT_title,
                            DRUGS_title, gender,
                            met_dir,
                            file_name)))
            else:
                logger.info(
                    'create_plots_14:\t{}'.format(
                        os.path.join(
                            PROJECT_title,
                            DRUGS_title, gender,
                            met_dir,
                            file_name)))
            print('plot saved in:', file_name)
            plt.clf()
    # ###############

    # # print("\nintersect_DF:\n", intersect_DF)
    # # we want to plot each betavalue for a REF, f.e. cg02584756
    # # cg14207654 cg10511249, for each REF we compare dead vs. alive

    # print(intersect_DF.loc[(slice(None), slice(None), slice(None),
    #                         'cg25999578'), :])
    # print("\nintersect_DF.index:\n", intersect_DF.index)
    # intersect_DF.sort_index(inplace=True)
    # print("\nintersect_DF.index:\n", intersect_DF.index)
    # print(intersect_DF.loc[(
    #     slice(None), slice(None), slice(None), 'cg12848345'), :])

    # alive_cols = []
    # # create the alive and dead graph_DF:
    # for i in intersect_DF_header:
    #     if re.search('alive',i):
    #         alive_cols.append(i)
    # alive_graph_DF = intersect_DF.loc[:,
    # intersect_DF.columns.isin(alive_cols)]
    # alive_graph_DF = alive_graph_DF.median(axis=1)

    # dead_cols = []
    # for i in intersect_DF_header:
    #     if re.search('dead',i):
    #         dead_cols.append(i)
    # dead_graph_DF = intersect_DF.loc[:, intersect_DF.columns.isin(dead_cols)]
    # dead_graph_DF = dead_graph_DF.median(axis=1)
    # graph_DF = pd.concat([intersect_DF['REF'], dead_graph_DF,
    # alive_graph_DF],
    #                     axis=1)
    # graph_DF.columns = ['REF', 'beta_median_dead', 'beta_median_alive']
    # graph_DF['abs_diff'] = abs(graph_DF['beta_median_dead'] -
    # graph_DF['beta_median_alive'])

    # print("\ngraph_DF:\n", graph_DF)
    # # plt.xkcd()

    # # gca stands for 'get current axis'
    # ax = plt.gca()
    # graph_DF.plot(kind='line',x='REF',y='beta_median_dead', ax=ax)
    # graph_DF.plot(kind='line',x='REF',y='beta_median_alive', color='red',
    # ax=ax)
    # graph_DF.plot(kind='line',x='REF',y='abs_diff', color='green', ax=ax)

    # plt.savefig(os.path.join(OUTPUT_PATH,'plot_test_diff.pdf'))


# # ###########fct 10 and 15 in main
def create_mean_boxplot(OUTPUT_PATH, DRUGS, met_opt_list, met_dir, PROJECT,
                        DRUGS_title):
    """
    - Load the REF positions out of the intersect.out table and the DMR ranges\
    out of metilene_qval.0.05.out

    - build the means for every case and make boxplot

    cols of metilene_qval.0.05.bedgraph:
        chr | start | stop | mean methylation difference

    cols of metilene_qval.0.05.out:
        chr | start | stop | q-value | mean methylation difference | #CpGs\
            | p(MWU) | p (2D KS) | mean g1 | mean g2 |

    - we use the metilene_qval.0.05.out since it holds more info for further\
    analysis

    content of boxplot:
        - mean meth differences plotted for every project, per project\
        comparison of dead and alive

        - now sum up the difference of median between dead and alive per\
        project an add up those differences for all projects available, the\
        higher the difference, the better the result

        - -> save those pairs, sort, and flag them accordingly in the\
            filename, s.t. they can be included in the report, if they\
            are good candidates

    - additionally, for every range (box)plottet, make a linegraph for the\
    medians of the betavalues of every case, at every position and for every\
    project


    .. _function_10:
    .. _function_15:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 10 to your call,
        # example:
        $ python main_metilene.py -p TCGA-CESC -d cisplatin -f 10

        # when choosing multiple projects, call:
        $ python main_metilene.py -p TCGA-CESC  -p TCGA-HNSC -d cisplatin -f 15

    """
    # plots for all projects:
    if isinstance(PROJECT, str):
        PROJECT_title = PROJECT
    else:
        project_list = []
        for project in PROJECT:
            project_list.append(project)
        PROJECT_title = '_'.join(sorted(map(str.upper, project_list)))
    print('calling create_mean_boxplot with PROJECT=', PROJECT_title)

    gender_list = ['female', 'male', 'both']
    for gender in gender_list:
        if not os.path.exists(os.path.join(OUTPUT_PATH, PROJECT_title,
                                           DRUGS_title, gender)):
            continue
        else:
            os.chdir(os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                  gender, met_dir))

            DF_range = pd.read_csv(
                'metilene_qval.0.05.out', sep='\t', header=None)
            # DF_range = pd.read_csv(os.path.join(OUTPUT_PATH,
            #                                     'own_filtered_met_output.txt'),
            #                        sep='\t', header=None)
            # print("\nDF_range:\n", DF_range)

            # list made of min, max position list pairs:
            list_min_max = []
            for i in range(0, len(DF_range)):
                list_min_max.append(list(DF_range.iloc[i, [1, 2]]))

            with open(os.path.join(os.pardir, 'intersect_header.tsv'
                                   ), 'r') as f:
                intersect_DF_header = f.readline().strip().split('\t')
            # # intersect_DF_header = list(pd.read_csv(
            #     # os.path.join(
            #         # os.pardir, 'intersect_header.tsv'),
            #         # sep='\t').columns)[:-1]
            intersect_DF = pd.read_csv(
                'intersect_out.tsv', header=None, sep='\t')
            intersect_DF.columns = intersect_DF_header
            # print("\nintersect_DF:\n", intersect_DF)
            # print("\nintersect_DF.columns:\n", intersect_DF.columns)
            intersect_DF.set_index(
                ['Chromosome', 'Start', 'End', 'REF'],
                inplace=True)
            # print("\nintersect_DF:\n", intersect_DF)
            # transforming the DF for re index with colnames, first split cols
            # at ; vital_status;uuid;drugs;project
            intersect_DF = intersect_DF.T
            vital_status_list = []
            uuid_list = []
            drug_list = []
            project_list = []
            for i in intersect_DF_header:
                # print("\ni:\n", i)
                if len(i.split(';')) == 5:  # ignore first fields, chromosome,
                    # start, ref...
                    temp_list = i.split(';')
                    vital_status_list.append(temp_list[0])
                    uuid_list.append(temp_list[1])
                    drug_list.append(temp_list[2])
                    project_list.append(temp_list[4])
            intersect_DF.set_index(
                [vital_status_list, uuid_list, drug_list, project_list],
                inplace=True)
            intersect_DF = intersect_DF.T
            # print("\nintersect_DF.loc[:, ('alive', 'TCGA-CESC')]:\n",
            #       intersect_DF.loc[:, ('alive', slice(None), slice(None),
            #                            ['TCGA-CESC', 'TCGA-HNSC'])])
            project_set = list(set(project_list))
            # print("\nproject_set:\n", project_set)
            # the plot shall include a dead and a alive graph for every
            # project, check if both vital states are available for all
            # projects:
            project_to_delete = []

            for vital_state in ['dead', 'alive']:
                # print('check for vital_state:', vital_state)
                for project in project_set:
                    # print("\nproject:\n", project)
                    temp_DF = intersect_DF.loc[
                        :, (vital_state, slice(None), slice(None), project)]
                    if temp_DF.empty:
                        print('delete', project, 'no vital_state', vital_state)
                        project_to_delete.append(project)
            if len(project_to_delete) != 0:
                project_to_delete = list(set(project_to_delete))
                # Remove an item by value: remove()
                for element in project_to_delete:
                    project_set.remove(element)
            print("\nproject_set:\n", project_set)
            print("\nproject_to_delete:\n", project_to_delete)
            intersect_DF.reset_index(inplace=True)
            # print("\nintersect_DF:\n", intersect_DF)

            os.makedirs('boxplot_range', exist_ok=True)
            # create a new dataframe for every min max pair fitting in range:
            # also load the metilene_qval.0.05.out to access the pvalues and
            # numbers of CpG's found, needed on REPORT creation
            met_out_DF = pd.read_csv(
                os.path.join(
                    OUTPUT_PATH,
                    PROJECT_title,
                    DRUGS_title, gender,
                    met_dir, 'metilene_qval.0.05.out'), sep='\t', header=None)
            # | chr | start | stop | q-value | mean methylation difference |
            # CpGs | mean g1 | mean g2 |
            met_out_DF.columns = ['chr', 'start', 'stop', 'q-value',
                                  'mean_methylation_difference', 'CpGs',
                                  'mean_alive', 'mean_dead']
            met_out_DF = met_out_DF.set_index(['chr', 'start', 'stop'])
            DF_table_index = 0
            for i in list_min_max:
                # we count the index to access the chromosome belonging to that
                # range, located in the DF_range table
                filt = (
                    intersect_DF.loc[:, 'Start'] >= i[0]) & (
                        intersect_DF.loc[:, 'End'] <= i[1])
                if intersect_DF[filt].empty:
                    continue
                else:
                    print('checking range {} - {}'.format(i[0], i[1]))
                    # print(intersect_DF[filt].iloc[:, range(0, 6)])
                    # make means of every col, make dataframe of multiindex
                    # mean series with resetting all index levels to cols, then
                    # drop the uuid col
                    temp_DF = intersect_DF[filt].iloc[:, 4:].mean(
                    ).reset_index().drop(columns='level_1')

                    temp_DF.columns = ['vital_state', 'drug', 'project',
                                       'mean_beta_values']

                    # ## short routine to benchmark the result
                    temp_DF_bench = temp_DF.groupby(
                        by=['project', 'vital_state']).mean()
                    median_diffs = []
                    for project in project_set:
                        alive_med = temp_DF_bench.loc[project,
                                                      :].loc['alive', :]
                        dead_med = temp_DF_bench.loc[project, :].loc['dead', :]
                        median_diffs.append(abs(alive_med - dead_med))
                    projects_diff = round(float(sum(median_diffs)), ndigits=4)

                    # we also flag the chromosome in the title, after that,
                    # increase the index chromosome for title:
                    chr_tle = DF_range.iloc[DF_table_index, 0]
                    DF_table_index += 1
                    file_name = 'boxplot_means_range_' + str(
                        i[0]) + '-' + str(i[1]) + '_' + str(projects_diff
                                                            ) + '_' + '.pdf'
                    file_name = os.path.join('boxplot_range', file_name)

                    plt.figure()
                    # print("\nintersect_DF_complete_Position:\n", temp_DF)
                    sns.set(style="ticks", palette="bright")
                    # sns.despine(offset=10, trim=True)
                    plot = sns.boxplot(x="project", y="mean_beta_values",
                                       hue="vital_state", palette=["m", "g"],
                                       data=temp_DF)
                    plt.title('range ' + str(i[0]) + ' - ' + str(i[1])
                              + '\n' + chr_tle)
                    plot.figure.savefig(file_name)
                    logger = set_logger(
                        OUTPUT_PATH, PROJECT_title, DRUGS_title)
                    file_name_DF = file_name.replace('.pdf', 'DF.tsv')
                    temp_DF.to_csv(file_name_DF, sep='\t', index=False)
                    if isinstance(PROJECT, str):
                        logger.info('create_mean_boxplot_10:\t{}'.format(
                            os.path.join(PROJECT_title, DRUGS_title, gender,
                                         met_dir, file_name)))
                        logger.info('create_mean_boxplot_10:\t{}'.format(
                            os.path.join(PROJECT_title, DRUGS_title, gender,
                                         met_dir, file_name_DF)))
                    else:
                        logger.info('create_mean_boxplot_15:\t{}'.format(
                            os.path.join(PROJECT_title, DRUGS_title, gender,
                                         met_dir, file_name)))
                        logger.info('create_mean_boxplot_15:\t{}'.format(
                            os.path.join(PROJECT_title, DRUGS_title, gender,
                                         met_dir, file_name_DF)))
                    # #### plot every position within the current range

                    # locatet in intersec_out.tsv -> intersect_DF contain all
                    # whats needed:
                    intersect_DFi = intersect_DF.set_index(
                        ['Start', 'End', 'Chromosome', 'REF'])
                    intersect_DFi.columns.names = ['vital_state', 'case_id',
                                                   'drug', 'project']
                    # DF_filt_temp = intersect_DFi['']
                    # save the start position which are within the range range
                    starts = []
                    for start in intersect_DFi.index:
                        # start grater than range start and smaller than range
                        # end?
                        if start[0] >= i[0] and start[0] <= i[1]:
                            starts.append(start[0])

                    DF_to_plot = intersect_DFi.loc[starts, :]

                    DF_to_plot = DF_to_plot.droplevel([1, 2], axis=1)
                    DF_to_plot = DF_to_plot.droplevel(
                        ['End', 'Chromosome', 'REF'])
                    DF_to_plot = DF_to_plot.unstack()

                    DF_to_plot = DF_to_plot.reset_index(
                        ['vital_state', 'project', 'Start'])
                    DF_to_plot = DF_to_plot.rename(
                        {0: 'mean_beta_value'}, axis=1)

                    DF_candidate_all_project = DF_to_plot.groupby(
                        ['project', 'vital_state',
                         'Start']).median().rename({'mean_beta_value':
                                                   'median_of_means'}, axis=1)
                    # drop Start index to directly compare the medians
                    DF_candidate_all_project = \
                        DF_candidate_all_project.reset_index([0, 2])
                    bool_DF = DF_candidate_all_project.loc[
                        'alive'].reset_index(
                    )['median_of_means'] > DF_candidate_all_project.loc[
                            'dead'].reset_index()['median_of_means']
                    # bool_DF is project specific, to test every position over
                    # all projects:
                    median_all_proj_alive = []
                    median_all_proj_dead = []
                    for project in project_set:
                        temp_DF_al = DF_candidate_all_project.loc[
                            'alive', :
                        ].set_index('project', append=True).loc[
                            (slice(None), project), :].set_index(
                                'Start').rename(
                                    {'median_of_means': 'medians_'
                                     + project + '_alive'}, axis=1)
                        temp_DF_dead = DF_candidate_all_project.loc[
                            'dead', :
                        ].set_index('project', append=True).loc[
                            (slice(None), project), :].set_index(
                                'Start').rename(
                                    {'median_of_means': 'medians_' + project
                                     + '_dead'}, axis=1)
                        median_all_proj_alive.append(temp_DF_al)
                        median_all_proj_dead.append(temp_DF_dead)
                    median_all_proj_DF_alive = pd.concat(
                        median_all_proj_alive, axis=1)
                    median_all_proj_DF_dead = pd.concat(
                        median_all_proj_dead, axis=1)
                    test_list = []
                    for alives in median_all_proj_DF_alive.columns:
                        for deads in median_all_proj_DF_dead.columns:
                            test_list.append(median_all_proj_DF_alive[alives]
                                             - median_all_proj_DF_dead[deads])
                    test_DF = pd.concat(test_list, axis=1)
                    test_DF = test_DF >= 0
                    # now concatenate all cols in one col and check number of
                    # unique values, if 1, every median of alive is higher than
                    # in deads, or v.v.
                    test_list = []
                    for test_col in test_DF.columns:
                        test_list.append(test_DF[test_col])
                    bool_DF_ap = pd.concat(test_list)
                    if bool_DF.nunique() == 1:
                        if bool_DF_ap.nunique() == 1:
                            file_name_2 = \
                                'boxplot_means_range_' + str(i[0]) + '-' + str(
                                    i[1]) + '_with_positions_CAND.pdf'
                            file_name_2 = os.path.join(
                                'boxplot_range', file_name_2)
                        else:
                            file_name_2 = \
                                'boxplot_means_range_' + str(i[0]) + '-' + str(
                                    i[1]) + '_with_positions_cand.pdf'
                            file_name_2 = os.path.join(
                                'boxplot_range', file_name_2)
                    else:
                        file_name_2 = 'boxplot_means_range_' + str(
                            i[0]) + '-' + str(i[1]) + '_with_positions.pdf'
                        file_name_2 = os.path.join(
                            'boxplot_range', file_name_2)

                    DF_to_plot['vital_project_merge'] = DF_to_plot[
                        'vital_state'] + '_' + DF_to_plot['project']
                    DF_to_plot = DF_to_plot.sort_values('vital_project_merge')
                    # file_name_2 = 'boxplot_means_range_with_positions_' +
                    # str( i[0]) + '-' + str(i[1]) + '.pdf' file_name_2 =
                    # os.path.join('boxplot_range', file_name_2)
                    plt.figure(figsize=(10, 5))
                    # len of palette:
                    # boxplot the means
                    palette_len = len(DF_to_plot.value_counts('project')) * 2
                    sns.set(
                        style="ticks",
                        palette=sns.color_palette(
                            'coolwarm_r',
                            palette_len))
                    plot = sns.boxplot(
                        x="Start",
                        y="mean_beta_value",
                        hue="vital_project_merge",
                        data=DF_to_plot)
                    plt.title('range ' + str(i[0]) + ' - ' + str(i[1])
                              + '\n' + chr_tle)
                    plot.set_xticklabels(plot.get_xticklabels(), rotation=45)
                    plt.legend(
                        title='vital state of all projects', bbox_to_anchor=(
                            1.05, 1), loc=2, borderaxespad=0.)
                    plt.tight_layout()
                    plot.figure.savefig(file_name_2)
                    plt.clf()
                    # ### plot just the median as linegraph:
                    DF_line_plot = DF_to_plot.groupby(
                        ['vital_state', 'project', 'Start']).median(
                    ).reset_index(
                        ['Start', 'vital_state', 'project']).rename(
                        {'mean_beta_value': 'median of means beta value'},
                        axis=1)
                    DF_line_plot['vital_proj'] = DF_line_plot['vital_state']\
                        + '_' + DF_line_plot['project']

                    plt.figure(figsize=(15, 5))
                    # plt.rcParams["figure.autolayout"] = True
                    # DF_line_plot['Start'] =
                    # pd.Categorical(DF_line_plot['Start'])
                    plt.figure()
                    plt.ticklabel_format(
                        style='plain', axis='x', useOffset=False)
        # style='vital_proj', markers=True
                    plot = sns.lineplot(
                        data=DF_line_plot,
                        x="Start",
                        marker='o',
                        y="median of means beta value",
                        linewidth=3,
                        hue="vital_proj")
                    plt.title('range ' + str(i[0]) + ' - ' + str(i[1])
                              + '\n' + chr_tle)
                    plt.setp(plot.get_xticklabels(), rotation=45)

                    legend = plt.legend(
                        title='vital state of all projects',
                        bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
                    for legobj in legend.legendHandles:
                        legobj.set_linewidth(3.0)
                    # plt.legend(
                        # title='vital state of all projects',
                        # bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
                    # plt.legend(
                    # title='vital state of all projects', bbox_to_anchor=(1,
                    # 1))
                    plt.tight_layout()
                    file_name_3 = file_name_2.replace('.pdf', '_linegraph.pdf')
                    plot.figure.savefig(file_name_3)

                    plt.clf()
                    # save the table to get the start positions while creating
                    # the REPORT also add the chr, s.t. no region is found on
                    # the wrong chromosome and the pvalue for later ranking on
                    # which CANDidate plots have the best p_value..:
                    # i holds the current range, acces directly through
                    # chr_tle,
                    # i[0] and i[1] (min and max)
                    # chr3    137763095       137772322       0.025473
                    # 0.047770        47      0.46099 0.41322

                    file_name_3_DF = file_name_3.replace('.pdf', '.tsv')
                    DF_line_plot = DF_line_plot.assign(chr=chr_tle)
                    # also add every information out of metilene_qval.0.05.out
                    # and the actual range out of min and max:
                    q_value = met_out_DF.loc[(chr_tle, i[0], i[1]), 'q-value']
                    m_m_d = met_out_DF.loc[
                        (chr_tle, i[0], i[1]), 'mean_methylation_difference']
                    CpGs = met_out_DF.loc[(chr_tle, i[0], i[1]), 'CpGs']
                    mean_alive = met_out_DF.loc[
                        (chr_tle, i[0], i[1]), 'mean_alive']
                    mean_dead = met_out_DF.loc[
                        (chr_tle, i[0], i[1]), 'mean_dead']
                    DF_line_plot = DF_line_plot.assign(
                        q_value=q_value,
                        mean_methylation_differences=m_m_d,
                        CpGs=CpGs,
                        mean_alive=mean_alive, mean_dead=mean_dead,
                        min_pos=str(i[0]), max_pos=str(i[1]))
                    DF_line_plot.to_csv(file_name_3_DF, sep='\t', index=False)

                    if isinstance(PROJECT, str):
                        logger.info('create_mean_boxplot_10:\t{}'.format(
                            os.path.join(PROJECT_title, DRUGS_title, gender,
                                         met_dir, file_name_2)))
                        logger.info('create_mean_boxplot_10:\t{}'.format(
                            os.path.join(PROJECT_title, DRUGS_title, gender,
                                         met_dir, file_name_3)))
                        logger.info('create_mean_boxplot_10:\t{}'.format(
                            os.path.join(PROJECT_title, DRUGS_title, gender,
                                         met_dir, file_name_3_DF)))
                    else:
                        logger.info('create_mean_boxplot_15:\t{}'.format(
                            os.path.join(PROJECT_title, DRUGS_title, gender,
                                         met_dir, file_name_2)))
                        logger.info('create_mean_boxplot_15:\t{}'.format(
                            os.path.join(PROJECT_title, DRUGS_title, gender,
                                         met_dir, file_name_3)))
                        logger.info('create_mean_boxplot_15:\t{}'.format(
                            os.path.join(PROJECT_title, DRUGS_title, gender,
                                         met_dir, file_name_3_DF)))
                    # #### end of position within ranger


# fct 11 in main
def create_table_all_projects(OUTPUT_PATH, DRUGS, PROJECT_DRUG_UUID,
                              DRUGS_title, met_dir, cutoff):
    """
    :param: OUTPUT_PATH: path for metilene pipeline outputs
    :type: OUTPUT_PATH: str
    :param: PROJECT: list of projects chosen
    :type: PROJECT: list of str
    :param: DRUGS_title: merged drug title out of multiple drugs
    :type: DRUGS_title: str

    walking in project/drugs/gender folder and aggregating those tables,
    summary_dead_alive.tsv

    vital_status;case_id;druglist;gender;project

        * vital_status
        * case_id
        * gender
        * therapeutics
        * project

    writing the new table in
    OUTPUT_PATH/PROJECT_title/DRUGS_title/gender/
    summary_dead_alive_dropped_NA_dropped_REF_sorted.tsv

    .. _function_11:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 11 to your call,
        # example:
        $ python main_metilene.py -p TCGA-CESC -d cisplatin -f 11
    """
    projects = []
    for project in PROJECT_DRUG_UUID:
        projects.append(project)
    PROJECT_title = '_'.join(sorted(map(str.upper, projects)))
    os.makedirs(os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title),
                exist_ok=True)

    print('fct 11, create_table_all_projects')
    print("\nprojects:\n", projects)

    # take first 3 cols and temp save those(Chromosome,
    # Start)

    summary_DF_list = []
    for project in projects:
        gender_list = []
        # we need to determine the gender_list in prior of loading summary
        # tables, s.t. cases are not added multiply, either male, or female, or
        # both should be loaded (not male, female and both!!!)
        for gender in ['female', 'male', 'both']:
            if not os.path.exists(os.path.join(OUTPUT_PATH, project,
                                               DRUGS_title, gender)):
                continue
            else:
                gender_list.append(gender)
        if len(gender_list) == 3:
            gender_list = ['both']
        for gender in gender_list:
            DF_dropped_REF_sorted = pd.read_csv(
                    os.path.join(
                        OUTPUT_PATH,
                        project,
                        DRUGS_title,
                        gender, 'dropped_REF_sorted.tsv'), sep='\t')
            summary_name = \
                'summary_dead_alive_dropped_NA_dropped_REF_sorted.tsv'
            summary_DF_temp = \
                pd.read_csv(
                    os.path.join(
                        OUTPUT_PATH,
                        project,
                        DRUGS_title,
                        gender, summary_name), sep='\t')
            # use the belonging REF positions as index, to join later every
            # DF on it, in case of aggregating both (female and male) the Y
            # REF positions are dropped, since females have no Y and the
            # aggregation with them can't be done, in males those positions
            # are preserved
            summary_DF_temp.set_index(
                DF_dropped_REF_sorted['Composite Element REF'],
                inplace=True)
            summary_DF_temp.set_index(
                ['Chromosome', 'Start'], append=True, inplace=True)
            summary_DF_list.append(summary_DF_temp)
    summary_DF = pd.concat(summary_DF_list, axis=1)
    # >>> for i in summary_DF_list:
    # ...     i.shape
    # ...
    # (393518, 116)
    # (393520, 22)
    # (393520, 89)
    # (387951, 109)
    # (393520, 7)
    # (393520, 12)
    # (391158, 17)
    #     concat al frames
    # we end up with different nr. of rows, concat to one DF, than select
    # gender specific, keep NA in -> the findal DF shoul have the max amount of
    # rows of all single DF's
    # (Pdb) final_summary_DF.shape
    # (393520, 358)
    # TODO how exactly manages metilene missing values?
    # #############
    # now: natsort, build the MI, select through the MI, build back the ';'
    # separated headers, write tables:
    # ###########

    # # WICHTIG, sort first Chromosome, then start: -> therefore make col
    # # 'Chromosome' categorial first, sort it then
    temp_chromosome_start_DF = summary_DF.reset_index(
        level=['Chromosome', 'Start'])
    temp_chromosome_start_DF = summary_DF.reset_index(
        level=['Chromosome', 'Start'])
    temp_chromosome_col = temp_chromosome_start_DF.loc[:, 'Chromosome']
    # prepare the categorial col, tell how the nonunique factors shall be
    # sorted:
    temp_chromosome_start_DF.loc[:, 'Chromosome'] = pd.Categorical(
        temp_chromosome_col, ordered=True, categories=natsorted(
            temp_chromosome_col.unique()))
    # do the actual sorting:
    temp_chromosome_start_DF = temp_chromosome_start_DF.sort_values(
        by=["Chromosome", "Start"], ascending=[True, True])
    # DF is correctly sorted now, set chr and start as index again:)
    summary_DF = temp_chromosome_start_DF.set_index(['Chromosome', 'Start'],
                                                    append=True)

    # build the Multiindex out of the col names:
    MI_dict = {'vital_state': [], 'UUID': [], 'DRUGS': [], 'gender': [],
               'PROJECT': []}
    col_name_list = list(summary_DF.columns)
    for col in col_name_list:
        MI_dict['vital_state'].append(col.split(';')[0])
        MI_dict['UUID'].append(col.split(';')[1])
        MI_dict['DRUGS'].append(col.split(';')[2])
        MI_dict['gender'].append(col.split(';')[3])
        MI_dict['PROJECT'].append(col.split(';')[4])
    # the multiindex can't be build directly from the dict, first create DF out
    # of dict, than use from_frame method:
    MI_index = pd.MultiIndex.from_frame(pd.DataFrame(MI_dict))
    summary_DF = summary_DF.T.set_index(MI_index).T
    # ### check if double UUID's occur:
    # temp_DF = summary_DF.T.reset_index('UUID')
    # temp_DF.loc[:, ('UUID', slice(None), slice(None),
    # slice(None))].value_counts()
    # search for col_names matching the DRUG pattern exactly:
    summary_DF = summary_DF.loc[:, summary_DF.columns.isin(DRUGS,
                                                           level='DRUGS')]
    # check for the genders
    # values, counts = np.unique(
    #     summary_DF.columns.isin(
    #         ['female'], level='gender'), return_counts=True)
    gender_list = []
    for gender in ['female', 'male']:
        values = np.unique(
            summary_DF.columns.isin([gender], level='gender'))
        if True in values:
            gender_list.append(gender)
    if len(gender_list) == 2:
        gender_list.append('both')
    logger = set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title)
    # for export, bring the header in that form:
    # Chromosome  Start   alive;aa10d6da-...;carboplatin...;female;TCGA-HNSC
    # chr1    15865   0.844660062689115
    # and also save the belonging REF in:
    # Composite Element REF
    # cg13869341
    # cg14008030
    for gender in gender_list:
        path_gen = os.path.join(
            OUTPUT_PATH, PROJECT_title, DRUGS_title, gender)
        os.makedirs(os.path.join(path_gen, met_dir), exist_ok=True)
        new_col_list = ['Chromosome', 'Start']
        # for both, write the whole table
        if gender == 'both':
            temp_REF_DF = summary_DF.reset_index(
                'Composite Element REF').loc[:, (
                    'Composite Element REF',
                    slice(None), slice(None), slice(None))]
            # drop the 4 MI levels not needed:
            temp_REF_DF = temp_REF_DF.T.reset_index([1, 2, 3, 4], drop=True).T
            temp_REF_DF.to_csv(os.path.join(
                path_gen, 'dropped_REF_sorted.tsv'), sep='\t', index=False)
            # the REF are now saved to the DF_dropped_REF_sorted.tsv
            # this index can be dropped now:
            temp_DF = summary_DF.reset_index(
                ['Chromosome', 'Start']).reset_index('Composite Element REF',
                                                     drop=True)
            for col_i in range(2, len(temp_DF.columns)):
                new_col_list.append(';'.join(temp_DF.columns[col_i]))
            temp_DF.columns = new_col_list
            summary_name = \
                'summary_dead_alive_dropped_NA_dropped_REF_sorted.tsv'
            temp_DF.to_csv(os.path.join(path_gen, summary_name),
                           sep='\t', index=False)
            logger.info(
                'create_table_for_metilene_11:\t{}'.format(
                    os.path.join(
                        PROJECT_title, DRUGS_title, gender, summary_name)))
            logger.info('create_table_for_metilene_11:\t{}'.format(
                os.path.join(PROJECT_title, DRUGS_title, gender,
                             'dropped_REF_sorted.tsv')))
            # to_csv(os.path.join(path_gen, 'summary_DF'), sep='\t')
        # or write the gender limited parts:
        else:
            # summary_DF.loc[:, (slice(None), slice(None), slice(None),
            # gender)].to_csv(os.path.join(path_gen, 'summary_DF'), sep='\t')
            temp_REF_DF = summary_DF.loc[
                :, (slice(None), slice(None), slice(None), gender)
            ].reset_index(
                'Composite Element REF').loc[:, (
                    'Composite Element REF',
                    slice(None), slice(None), slice(None))]
            # drop the 4 MI levels not needed:
            temp_REF_DF = temp_REF_DF.T.reset_index([1, 2, 3, 4], drop=True).T
            temp_REF_DF.to_csv(os.path.join(
                path_gen, 'dropped_REF_sorted.tsv'), sep='\t', index=False)
            # the REF are now saved to the DF_dropped_REF_sorted.tsv
            # this index can be dropped now:
            temp_DF = summary_DF.loc[
                :, (slice(None), slice(None), slice(None), gender)
            ].reset_index(
                ['Chromosome', 'Start']).reset_index(
                    'Composite Element REF', drop=True)
            for col_i in range(2, len(temp_DF.columns)):
                new_col_list.append(';'.join(temp_DF.columns[col_i]))
            temp_DF.columns = new_col_list
            temp_DF.to_csv(os.path.join(
                path_gen,
                'summary_dead_alive_dropped_NA_dropped_REF_sorted.tsv'
            ), sep='\t', index=False)
            summary_name = \
                'summary_dead_alive_dropped_NA_dropped_REF_sorted.tsv'
            logger.info('create_table_for_metilene_11:\t{}'.format(
                        os.path.join(
                            PROJECT_title, DRUGS_title, gender, summary_name)))
            logger.info('create_table_for_metilene_11:\t{}'.format(
                        os.path.join(
                            PROJECT_title,
                            DRUGS_title, gender, 'dropped_REF_sorted.tsv')))

    ######################################################################
    #  START aggregating the saved cases out of the single projects, which are
    # added due to a cutoff parameter
    DF_all_added_cases = pd.DataFrame()
    for project in projects:
        file_path = os.path.join(OUTPUT_PATH, project, DRUGS_title,
                                 'cutoff_cases_add_' + str(cutoff) + '.tsv')
        try:
            DF_temp = pd.read_csv(file_path, sep='\t')
        except FileNotFoundError:
            continue
        DF_all_added_cases = pd.concat([DF_temp, DF_all_added_cases])
    if not DF_all_added_cases.empty:
        DF_all_added_cases.to_csv(
            os.path.join(
                OUTPUT_PATH,
                PROJECT_title,
                DRUGS_title,
                'cutoff_cases_add_' + str(cutoff) + '.tsv'),
            sep='\t', index=False)
        logger.info('create_table_for_metilene_11:\t{}'.format(
                    os.path.join(
                        PROJECT_title,
                        DRUGS_title,
                        'cutoff_cases_add_' + str(cutoff) + '.tsv')))

    #####################################################################

# ./metilene_linux64 -m 3 -M 1000 -a alive -b dead
# /scr/dings/PEVO/Methylation/new_version/TCGA-CESC/cisplatin/
# summary_dead_alive_dropped_NA_dropped_REF_sorted.tsv
# | sort -k1,1V -k2,2n > met_output_sorted.txt
# ./metilene_output.pl -q met_output_sorted.txt -c 3 -a alive -b dead

#
# awk '{print $1,'\t',$2,'\t'$3,'\t',$4,'\t',$5,'\t',$6, '\t',$9,'\t',$10}'
# met_output_sorted.txt > met_output_for_R.txt
# call rscript with pdf outname:
# ./metilene_output.R met_output_for_R.txt R_out.pdf


def set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title):
    '''
    :param: OUTPUT_PATH: path for metilene pipeline outputs
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


def create_log_for_A(OUTPUT_PATH, PROJECT_list, PROJECT_title, DRUGS_title,
                     api_manifest):
    '''
    in case the -A option is set without the -D, the log file has to be createt
    in advance. With that, it is simultaniously checked if the files needed
    (downloaded datafiles, etc) are actually already present
    '''
    os.makedirs(os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title),
                exist_ok=True)
    logger = set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title)
    # logger.info(
    # 'log_file:\t{}'.format(
    # os.path.join(PROJECT_title, DRUGS_title, 'test_log.log')))

    for PROJECT in PROJECT_list:
        # check for the complete_summary.tsv in PROJECT/
        filename_manifest = api_manifest
        if os.path.exists(os.path.join(OUTPUT_PATH, api_manifest)):
            logger.info(f'download_GDC_manifest_1:\t{filename_manifest}')
        else:
            meta_path = os.path.join(OUTPUT_PATH, PROJECT,
                                     'complete_summary.tsv')
            print(f'could not find {meta_path}, ', end='')
            print('have you performed the download steps with the -D ', end='')
            print('option in advance?')
            os._exit(0)

        # check for 'nationwidechildrens.org_clinical_drug_cesc.txt'
        # in PROJECT/ dir:
        infix = PROJECT.replace('TCGA-', '').lower()
        drug_table_name = os.path.join(
            OUTPUT_PATH, PROJECT, 'aux_files',
            'nationwidechildrens.org_clinical_drug_' + infix + '.txt')
        if os.path.exists(drug_table_name):
            logger.info(
                'download_clinical_tables_2:\t{}'.format(
                    drug_table_name.replace(OUTPUT_PATH + os.path.sep, '')))
        else:
            print('could not find {}, have you performed the download steps\
                  with -D option in advance?'.format(drug_table_name))
            os._exit(0)

        # check for meta_info_druglist_merged.tsv,
        # meta_info_druglist_merged_drugs_combined.tsv in PROJECT/ dir:
        druglist_merged = os.path.join(OUTPUT_PATH, PROJECT,
                                       'meta_info_druglist_merged.tsv')
        druglist_combined = os.path.join(
            OUTPUT_PATH,
            PROJECT, 'meta_info_druglist_merged_drugs_combined.tsv')
        if os.path.exists(druglist_merged) and os.path.exists(
                druglist_combined):
            logger.info(
                'create_merged_metatable_3:\t{}'.format(
                    druglist_merged.replace(OUTPUT_PATH + os.path.sep, '')))
            logger.info(
                'create_merged_metatable_3:\t{}'.format(
                    druglist_combined.replace(OUTPUT_PATH + os.path.sep, '')))
        else:
            print('could not find {} or {}, have you performed the download\
                  steps with -D option in advance?'.format(
                druglist_merged, druglist_combined))
            os._exit(0)

        # check for sep_down files, first check existence of dir (OUTPUT_PATH,
        # PROJECT, 'TCGA-****_data_files'), then glob every file out of it:
        data_dir = os.path.join(OUTPUT_PATH, PROJECT, PROJECT + '_data_files')
        if os.path.exists(data_dir):
            for data_file in glob.glob(os.path.join(data_dir, '*')):
                logger.info('sep_down_data_files_4:\t{}'.format(
                    data_file.replace(OUTPUT_PATH + os.path.sep, '')))
        else:
            print('could not find {}, have you performed the download steps\
                  with -D option in advance?'.format(data_dir))
            os._exit(0)

        # check for 'complete_summary.tsv' in OUTPUT_PATH, PROJECT:
        # summary_table = os.path.join(
        #     OUTPUT_PATH, PROJECT, 'complete_summary.tsv')
        # if os.path.exists(summary_table):
        #     logger.info('create_summary_table_5:\t{}'.format(
        #         summary_table.replace(OUTPUT_PATH + os.path.sep, '')))
        # else:
        #     print('could not find {}, have you performed the download steps\
        #           with -D option in advance?'.format(
        #         summary_table.replace(OUTPUT_PATH + os.path.sep, '')))
        #     os._exit(0)


def sanity_check(OUTPUT_PATH, DRUGS):
    drugs = ','.join(DRUGS)
    OUTPUT_PATH = os.path.join(OUTPUT_PATH, drugs)
    beta_val_DF = pd.read_csv(os.path.join(
        OUTPUT_PATH,
        'summary_dead_alive_dropped_NA_dropped_REF_sorted.tsv'), sep='\t')
    ref_DF = pd.read_csv(
        os.path.join(OUTPUT_PATH, 'dropped_REF_sorted.tsv'), sep='\t')
    print("\nbeta_val_DF:\n", beta_val_DF)
    print("\nref_DF:\n", ref_DF)
    ref_DF.rename(columns={'Composite Element REF': 'REF'}, inplace=True)
    concated_DF = pd.concat([beta_val_DF, ref_DF], axis=1)
    # Chromosome and Start col not needed:
    concated_DF.drop(['Chromosome', 'Start'], axis=1, inplace=True)
    concated_DF.set_index('REF', inplace=True)
    print("\nconcated_DF.columns:\n", concated_DF.columns)
    print("\nconcated_DF:\n", concated_DF)
    # split out vital state and project name out of header, make multiindex out
    # of new lists:
    vital_list = []
    project_list = []
    for column in concated_DF.columns:
        vital_list.append(column.split(';')[0])
        project_list.append(column.split(';')[3])
    print("\nvital_list:\n", vital_list)
    index_list = [vital_list, project_list]
    concated_DF.columns = index_list
    concated_DF.columns.names = ['vital_state', 'project']
    # print("\nconcated_DF:\n", concated_DF)
    # print("\nconcated_DF.columns.names:\n", concated_DF.columns.names)
    concated_DF = concated_DF.T.reset_index()
    print("\nconcated_DF:\n", concated_DF)

    OUTPUT_PATH = os.path.join(OUTPUT_PATH, 'candidates_to_check')
    os.chdir(OUTPUT_PATH)
    file_list = os.listdir()
    # parse out the cg REF flags out of the filenames, and replace filenames
    # with that identifier:
    pattern = re.compile(r'cg\d+', re.I)
    for pos in range(0, len(file_list)):
        matches = pattern.search(file_list[pos])
        file_list[pos] = matches.group()
    # print("\nfile_list:\n", file_list)
    # now file_list just holds the ref cg... identifiers, with them, check the
    # beta_values at those postiions:

    for i in range(0, len(file_list)):
        temp_DF = concated_DF.loc[:, ['vital_state', 'project', file_list[i]]]
        plt.figure()
        sns.set(style="ticks", palette="bright")
        # sns.despine(offset=10, trim=True)
        plot = sns.boxplot(
            x="project",
            y=file_list[i],
            hue="vital_state", palette=["m", "g"], data=temp_DF)
        # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        file_name = 'boxplot_check_' + file_list[i] + '.pdf'
        file_name = os.path.join(OUTPUT_PATH, file_name)
        print('save plot in: ', file_name)
        plot.figure.savefig(file_name)
        plt.clf()


def sanity_check_2(OUTPUT_PATH, DRUGS, PROJECT_DRUG_UUID):
    '''
    taking results of single metilene runs out of project dirs and cut one
    table for intersect together
    '''
    drugs = ','.join(DRUGS)
    # concat an empty DF, check if the met_qval.0.05.bedraph files are empty:
    met_filtered_DF = pd.DataFrame()
    for PROJECT in PROJECT_DRUG_UUID:
        file_name = os.path.join(OUTPUT_PATH, PROJECT, drugs,
                                 'metilene_qval.0.05.out')
        if os.path.isfile(file_name) and os.path.getsize(file_name) != 0:
            temp_DF = pd.read_csv(file_name, header=None, sep='\t')
            print("\ntemp_DF:\n", temp_DF)
            met_filtered_DF = pd.concat([met_filtered_DF, temp_DF])
    # doing whats done with sort -V (versionsort)
    # https://stackoverflow.com/questions/29580978
    # /naturally-sorting-pandas-dataframe
    # sort met_filtered_DF chr and then position wise:
    # met_filtered_DF[0] = met_filtered_DF[0].astype('category')
    # met_filtered_DF[0].cat.reorder_categories(
    #     natsorted(set(met_filtered_DF[0])), inplace=True, ordered=True)

    # You can convert values to ordered categorical with sorted catgories by
    # natsorted and then use sort_values:
    met_filtered_DF[0] = pd.Categorical(
        met_filtered_DF[0],
        ordered=True, categories=natsorted(met_filtered_DF[0].unique()))
    met_filtered_DF.sort_values(by=[0, 1], inplace=True)
    print("\nmet_filtered_DF:\n", met_filtered_DF)
    # ##### write this concated table in the drugs dir, s.t. the intersect and
    # consecutive steps can be performed with this table:
    file_name = os.path.join(
        OUTPUT_PATH, drugs, 'metilene_qval_aggregated.tsv')
    met_filtered_DF.to_csv(file_name, sep='\t', header=False, index=False)
    print("\nfile_name:\n", file_name)


# fct 16 in main
def create_snake_config(OUTPUT_PATH, PROJECT_title, DRUGS_title, PROJECT_list,
                        DRUGS, SCRIPT_PATH, cutoff):
    '''
    :param: OUTPUT_PATH: path for metilene pipeline outputs
    :type: OUTPUT_PATH: str
    :param: PROJECT_title: concatenated str of multiple projects
    :type: PROJECT_title: str
    :param: DRUGS_title: concatenated str of multiple drugs
    :type: DRUGS_title: str
    :param: PROJECT_list: list of applied projects
    :type: PROJECT_list: list of str
    :param: DRUGS: applied drug(s)
    :type: DRUGS: list of str
    :param: SCRIPT_PATH: path to the metilene_pipeline repo
    :type: SCRIPT_PATH: str

    out of the log files in PROJECT_title/DRUGS_title/test_log.log
    parse out all outputfiles of the applied run and create
    PROJECT_title/DRUGS_title/snakemake_config.yaml


    .. _function_16:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 16 to your call,
        # example:
        $ python main_metilene.py -p TCGA-CESC -d cisplatin -f 16
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
    logger.info('config_file:\t{}'.format(conf_file_logger))

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

    print('a snakemake configuration file has been created with the following\
          projects')
    print(PROJECT_list)
    # additionally cp that file with the description in the name in the script
    # dir:
    name_copied = os.path.join(SCRIPT_PATH, 'Snakes', 'snakemake_config_'
                               + PROJECT_title + '_' + DRUGS_title + '.yaml')

    os.system('cp ' + conf_file + ' ' + name_copied)
    print('a copy of the Snakemake configuration is created in your ')
    print('SCRIPT_PATH:\n{}'.format(name_copied))


# def cutoff_table_changes(OUTPUT_PATH, PROJECT, DRUGS, DRUGS_title, met_dir,
#                          cutoff, summary_DF):
#     '''
#     merge on cases the drugs and patients tables:
#     nationwidechildrens.org_clinical_drug_cesc.txt
#     nationwidechildrens.org_clinical_patient_cesc.txt
#     information needed to get
#     patient table:
#         - last_contact_days_to
#         - death_days_to
#     -> second and third line can be omitted
#     (col 'bcr_patient_uuid' is the case_id)
#     '''
#     project = PROJECT.split('-')[1].lower()
#     # (141, 34):
#     meta_drugs_DF = pd.read_csv(
#         os.path.join(
#             OUTPUT_PATH, PROJECT,
#             'meta_info_druglist_merged_drugs_combined.tsv'), sep='\t',
#         index_col=0)
#     # (307, 2):
#     patient_DF = pd.read_csv(
#         os.path.join(
#             OUTPUT_PATH, PROJECT,
#             'nationwidechildrens.org_clinical_patient_' + project +
#             '.txt'), sep='\t', skiprows=[1, 2],
#         index_col=0).reindex(columns=['death_days_to',
#         'last_contact_days_to'])
#     patient_DF.index = patient_DF.index.str.lower()
#     patient_DF = patient_DF.reset_index().drop_duplicates(
#         'bcr_patient_uuid').set_index('bcr_patient_uuid')
#     # (141, 36):
#     merged_DF = meta_drugs_DF.join(patient_DF)
def snake_meta(PROJECT, FILE_TYPE, OUTPUT_PATH, DRUGS_title, SCRIPT_PATH,
               PROJECT_title):
    '''
    copy the respective
    SCRIPT_PATH/Snakes/meta_infos/PROJECT_title/DRUGS_title/meta_info.dat
    into the actual OUTPUT_PATH/PROJECT/ path
    '''
    # its ok to glob the DRUGS_title_cutoff value, the meta_info.dat ist not
    # affected by cutoff
    temp_dir = PROJECT_title + '_' + DRUGS_title
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
