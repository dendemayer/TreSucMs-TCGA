#!/usr/bin/env python3.7
import requests
import json
import tarfile
import os
import pandas as pd
import numpy as np
# from matplotlib import pyplot as plt
# from natsort import natsorted
import glob
# import logging
# import numpy as np  # for aranging REF at x
import re
# import requests
# import seaborn as sns
import set_logger
import subprocess


# fct 1 in main
# this will be held available in the repo:
def download_GDC_manifest(PROJECT, FILE_TYPE, OUTPUT_PATH, DRUGS_title,
                          SCRIPT_PATH, PROJECT_title, snakerun, api_manifest,
                          gtf_UUID):
    """
    - downloading the complete manifest of the 31 data release (the last
    providing unnormalized htseq counts for deseq2)
    - downloading the annotation file, this will be completed here in place
    with the long name version taken from mygene lib, s.t. the analyses part,
    where those description are needed can work offline

    gdc_manifest_20211029_data_release_31.0_active.tsv.gz
    gencode.v36.annotation.gtf.gz

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
    dest_manifest = os.path.join(OUTPUT_PATH, api_manifest)
    dest_gtf = os.path.join(OUTPUT_PATH, 'gencode.v36.annotation.gtf.gz')
    log_gtf = 'gencode.v36.annotation.gtf.gz'
    # in case the file is already present, return ,but check also here the
    # md5sum:
    md5sum_mani = '26d406d171e43f470ba338f4ec480cc0'
    md5sum_gtf = '291330bdcff1094bc4d5645de35e0871'
# https://api.gdc.cancer.gov/data/25aa497c-e615-4cb7-8751-71f744f9691f
    mani_bool = False
    gtf_bool = False
    # checking existence of manifestfile
    if os.path.exists(dest_manifest):
        # checking md5 of manifestfile
        md5sum_file = subprocess.check_output(
            ['md5sum', dest_manifest]).decode(
                'utf-8').split(' ')[0]
        if md5sum_mani == md5sum_file:
            print(f'{dest_manifest} already loaded')
            print(f'with verified md5sum of \n{md5sum_file}')
            mani_bool = True
        else:
            print(f'please remove file {dest_manifest}, file might corrupted')
            print(f'{md5sum_mani} != {md5sum_file}')
            print('exiting now')
            os._exit(0)

    # checking existence of gtf_file:
    if os.path.exists(dest_gtf):
        # checking md5 of gtf_file
        md5sum_file = subprocess.check_output(
            ['md5sum', dest_gtf]).decode(
                'utf-8').split(' ')[0]
        if md5sum_gtf == md5sum_file:
            print(f'{dest_gtf} already loaded')
            print(f'with verified md5sum of \n{md5sum_file}')
            gtf_bool = True
            # return
        else:
            print(f'please remove file {dest_gtf}, file might corrupted')
            print(f'{md5sum_gtf} != {md5sum_file}')
            print('exiting now')
            os._exit(0)

    if gtf_bool and mani_bool:
        return

    # manifest file not present or md5 wrong, download it:
    # download Data Release 31.0
    # https://docs.gdc.cancer.gov/Data/Release_Notes/gdc_manifest_20211029_data_release_31.0_active.tsv.gz
    if not mani_bool:
        url_gdc_manifest = ('https://docs.gdc.cancer.gov/Data/Release_' +
                            'Notes/' + api_manifest)
        print(f"\ngdc_manifest from:\n{url_gdc_manifest}\nis downloaded to")
        print(f"{OUTPUT_PATH}")
        data = requests.get(url_gdc_manifest)
        logger = set_logger.set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
        with open(dest_manifest, 'wb') as f:
            try:
                f.write(data.content)
                logger.info(f'download_GDC_manifest_1:\t{api_manifest}')
            except BaseException:
                print('problems while downloading the GDC manifest')
                print('is https://portal.gdc.cancer.gov/ online?')
                print('exiting the program')
                os._exit(0)

        md5sum_file = subprocess.check_output(
            ['md5sum', dest_manifest]).decode('utf-8').split(' ')[0]
        if md5sum_file == md5sum_mani:
            print(f'md5sum of\n{md5sum_file}\nis verified')

    # annotation gtf not present or md5 wrong, download it:
    if not gtf_bool:
        gtf = os.path.join('https://api.gdc.cancer.gov/data', gtf_UUID)
        print(f'downloading annotation gtf to\n{dest_gtf}')
        data = requests.get(gtf)
        logger = set_logger.set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
        with open(dest_gtf, 'wb') as f:
            try:
                f.write(data.content)
                logger.info(f'download_GDC_manifest_1:\t{log_gtf}')
            except BaseException:
                print(f'problems while downloading {dest_gtf}')
                print('is https://portal.gdc.cancer.gov/ online?')
                print('exiting the program')
                os._exit(0)

        md5sum_file = subprocess.check_output(
            ['md5sum', dest_gtf]).decode('utf-8').split(' ')[0]
        if md5sum_file == md5sum_gtf:
            print(f'md5sum of\n{md5sum_file}\nis verified')
    # adding long name to the genename with mygene library:
    # DF_annot = pd.read_csv(
    #     os.path.join(
    #         OUTPUT_PATH,
    #         'gencode.v36.annotation.gtf.gz'),
    #     sep='\t', skiprows=5, header=None)
    # DF_annot = DF_annot[(DF_annot[2] == 'gene')]

    # # filtering out the ENSG and create a new col for it:
    # print(f'filtering and adding long_name col to the annotation file:')
    # DF_annot['ENSG'] = DF_annot[8].apply(
    #     lambda x: re.search(r'ENSG\d+', x).group() if re.search(r'ENSG\d+',
    #                                                             x) else x)
    # DF_annot.set_index('ENSG', inplace=True)
    # DF_annot[8].apply(lambda x: re.search('".*"', x.split(';')[1]).group())
    # DF_annot[8].apply(
    #     lambda x: re.search('".*"', x.split(';')[1]).group().strip('"'))
    # # filtering out the gene_type (called 'type_of_gene')
    # DF_annot['type_of_gene'] = DF_annot[8].apply(
    #     lambda x: re.search('".*"', x.split(';')[1]).group().strip('"'))
    # # for a better representation in the REPORT.pdf, replace _ with -,
    # # s.t. Carriage returns are introduced within the tables
    # DF_annot['type_of_gene'] = DF_annot[
    #     'type_of_gene'].str.replace('_', '-')
    # # filtering out the 'symbol' -> gene symbol in annot file called
    # # gene_name must be named 'symbol'
    # DF_annot['symbol'] = DF_annot[8].apply(
    #     lambda x: re.search('".*"', x.split(';')[3]).group().strip('"'))
    # DF_annot['long_name'] = ''
    # mg = mygene.MyGeneInfo()
    # for ENSG in DF_annot.index:
    #     try:
    #         DF_annot.at[ENSG, 'long_name'] = mg.getgene(
    #             ENSG, fields='name')['name']
    #     except Exception as e:
    #         DF_annot.at[ENSG, 'long_name'] = DF_annot.loc[ENSG, 'symbol']
    # DF_annot.to_csv()


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
            OUTPUT_PATH, filename_manifest),
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

    logger = set_logger.set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
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
    logger = set_logger.set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
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
    manifest_file = os.path.join(SCRIPT_PATH, os.path.pardir, 'resources',
                                 'GCv36_Manifests', f'{PROJECT}.tsv')
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
    # proj_low_suf = PROJECT.replace('TCGA-', '').lower()
    vital_table_path = glob.glob(path_aux + '/*clinical_follow_up_v*')[0]
    vital_base = os.path.basename(vital_table_path)
    print(f'\nusing {vital_base} as clinical followup resource\n')
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
    # correct drugnames of all indexes:
    drug_DF['pharmaceutical_therapy_drug_name'] = drug_DF[
        'pharmaceutical_therapy_drug_name'].apply(set_logger.apply_list)
    # those nupis are from the drug table!
    not_uniq_filt = drug_DF['bcr_patient_barcode'].value_counts() > 1
    not_uniq_patient_index = not_uniq_filt[not_uniq_filt].index
    drug_DF = drug_DF.set_index('bcr_patient_barcode', drop=False)
    # out of the multiple drugs, drop duplicates, make a sorted list:

    for nupi in not_uniq_patient_index:
        sorted_drug_list = sorted(drug_DF.loc[nupi, :][
                'pharmaceutical_therapy_drug_name'].drop_duplicates(
                ).to_list())
        # correct the drugs within the list:
        # print(f'sorted_drug_list before correction:\n{sorted_drug_list}')
        # for index, value in enumerate(sorted_drug_list):
        #     sorted_drug_list[index] = apply_list(value)
        # print(f'sorted_drug_list after correction:\n{sorted_drug_list}\n')
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
    logger = set_logger.set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
    try:
        drug_DF.to_csv(drug_out_path, sep='\t', index=False)
        logger.info('create_merged_metatable_3:\t{}'.format(
            os.path.join(PROJECT, 'meta_info_druglist_merged.tsv')))
        print(f'saved meta_info_druglist_merged.tsv to {drug_out_path}')
    except Exception as e:
        print(f'Exception occured: {e}, could not save {drug_out_path}')
        print('exiting')
        os._exit(0)

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
    # in rare cases it can happen that both, death_days_to and
    # last_contact_days_to hold values, thats not Applicable, in that case take
    # the value of the either and set the other on nan, depending on vital
    # status
    complete_DF['death_days_to'] = pd.to_numeric(
        complete_DF['death_days_to'], errors='coerce')
    complete_DF['last_contact_days_to'] = pd.to_numeric(
        complete_DF['last_contact_days_to'], errors='coerce')
    # check index where both values are set to non NA
    double_value_ind = complete_DF.loc[:, [
        'last_contact_days_to', 'death_days_to']].dropna().index
    # example for HNSC:
    # vital_status  last_contact_days_to  death_days_to
    # 13         dead                1438.0         2319.0
    # 15         dead                 546.0          546.0
    # 37         dead                 548.0          548.0
    for index in double_value_ind:
        if complete_DF.loc[index, 'vital_status'] == 'dead':
            complete_DF.loc[index, 'last_contact_days_to'] = pd.NA
        else:
            complete_DF.loc[index, 'death_days_to'] = pd.NA
    try:
        complete_DF.to_csv(complete_path, sep='\t', index=False)
        logger.info('create_merged_metatable_3:\t{}'.format(
            os.path.join(PROJECT,
                         'meta_info_druglist_merged_drugs_combined.tsv')))
        print('saved meta_info_druglist_merged_drugs_combined.tsv ')
        print(f'in {complete_path}')
    except Exception as e:
        print(f'Exception occured: {e}, could not save {complete_path}')
        print('exiting')
        os._exit(0)
    # ############## keep compatibility with previous deseq routines,
    # year of birth and year of death are actually not needed and also not
    # refered to in the following steps...
    # create a merge-table (DF_3t_both_with_DRUG_combi.tsv) in the form of:
    # name                    ║ value          ║ equivalent of
    #                                         meta_info_druglist_merged

    # UUID                    ║ 00276f29-b-... ║ id
    # case_id                 ║ 6ff12a54-1-... ║ bcr_patient_uuid
    # gender                  ║ female         ║ gender
    # vital_status            ║ alive          ║ vital_status
    # drugnames               ║ cisplatin      ║ pharmaceutical_
    #                                            therapy_drug_name
    # survivaltime            ║ NaN            ║ death_days_to
    # years_to_last_follow_up ║ 3.0219178      ║ last_contact_days_to
    # #### this is omitted   -->>
    # year_of_birth           ║ 1946           ║
    # year_of_death           ║ NaN            ║
    # #### this is omitted  <<--
    # age_at_diagnosis        ║ 64.561643      ║ age_at_diagnosis
    # PROJECT                 ║ TCGA-CESC      ║ project_id

    # this is DESeq2_pipeline specific:
    if workflow_type == 'HTSeq - Counts':
        complete_DF.rename({'id': 'UUID', 'bcr_patient_uuid': 'case_id',
                            'pharmaceutical_therapy_drug_name': 'drugnames',
                            'death_days_to': 'survivaltime',
                            'last_contact_days_to': 'years_to_last_follow_up',
                            'project_id': 'PROJECT'}, axis=1, inplace=True)
        # make years out of day specificatins: -> therefore set NaN correctly.
        # that means str which cannot be converted to float are NaN
        # float(complete_DF['death_days_to'].loc[0]) -> 570.0
        # float(complete_DF['death_days_to'].loc[1]) -> *** ValueError: could
        # not convert string to float: '[Not Applicable]'

        def days_to_years(value):
            try:
                year = float(value) / 365
                return year
            except ValueError:
                return np.nan

        # make years out of last_concact_days_to col:
        complete_DF['survivaltime'] = complete_DF['survivaltime'].apply(
            days_to_years)
        complete_DF['years_to_last_follow_up'] = complete_DF[
            'years_to_last_follow_up'].apply(days_to_years)
        # save that table with name DF_3t_both_with_DRUG_combi
        log_path = os.path.join(PROJECT, 'DF_3t_both_with_DRUG_combi.tsv')
        try:
            complete_DF.to_csv(os.path.join(OUTPUT_PATH, log_path), sep='\t')
            logger.info('create_merged_metatable_3:\t{}'.format(
                os.path.join(log_path)))
        except Exception as e:
            print(f'Exception occured: {e}, could not save {log_path}')
            print('exiting')
            os._exit(0)


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

    logger = set_logger.set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
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
