import os
import pandas as pd
from yaml import safe_load
def download_help_files(OUTPUT_PATH, config_file_shared):
    """
    # creating the file names which should be requested through Snakemake
    downloading manifest file and the gtf annotation file from TCGA:
    # we are using a shared Snakefile, but refere here to the metilene specific
    # Snakefile, which again refers to the shared Snakefile via include
    # statement
    this fct returns the gtf and manifest filenames provided in the shared
    config the must not yet be present yet, the files are downloaded then via
    shared Snakefile: rule download_helpfiles:
    """
    # reading the manifest and gtf file out of the
    # SCRIPT_PATH/src/shared/config.yaml
    with open(config_file_shared, 'r') as f:
        df = pd.json_normalize((safe_load(f)))
    mani = df['manifest_file'].tolist()[0][0]
    mani_path = os.path.join(OUTPUT_PATH, 'metadata', mani)
    gtf = df['gtf.gz'].tolist()[0][0]
    gtf_path = os.path.join(OUTPUT_PATH, 'metadata', gtf)
    help_file_list = [gtf_path, mani_path]
    return help_file_list
    # if we are in this file,
    # snakemake -p --cores 40 /scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines/metadata/gdc_manifest_20211029_data_release_31.0_active.tsv.gz


def download_aux_files(OUTPUT_PATH, PROJECT, config_file_shared):
    """
    the manifest file is now downloaded vie Snakemake, provide
    ############
    create the aux filenames which should be requested to the Snakefile
    with help of the previously loaded  manifest file -> this manifest file is
    the input for this rule! and must be already loaded beforehand!
    ###########

    :param: OUTPUT_PATH: path for metilene pipeline outputs
    :type: OUTPUT_PATH: str
    :param: PROJECT: list of projects chosen
    :type: PROJECT: list of str

    with the UUID the clinical tables will be downloaded in the
    OUTPUT_PATH/PROJECT:

        * nationwidechildrens.org_clinical_patient_****.txt
        * nationwidechildrens.org_clinical_drug_****.txt
    """
    # get the UUID of all patient related meta tables out of the gdc_manifest:
    # read in the table:
    with open(config_file_shared, 'r') as f:
        df = pd.json_normalize((safe_load(f)))
    mani = df['manifest_file'].tolist()[0][0]

    mani_path = os.path.join(OUTPUT_PATH, 'metadata', mani)

    # filter on the project_id, this workds also for mutliple projects
    gdc_manifest_DF = pd.read_csv(
        mani_path, sep='\t').set_index(
            'project_id').loc[PROJECT, :].reset_index()

    # we just need the cols project_id and filename
    gdc_manifest_temp = gdc_manifest_DF[gdc_manifest_DF[
        'filename'].str.contains(
            'children.*.txt', regex=True)].loc[:, ['project_id', 'filename']]

    # create absolute paths to the filenames where they should be saved to:
    gdc_manifest_temp['filepath'] = OUTPUT_PATH + '/' + gdc_manifest_temp['project_id'] + '/aux_files/' + gdc_manifest_temp['filename']
    aux_file_list = [i for i in gdc_manifest_temp['filepath']]
    return aux_file_list

def download_data_files(OUTPUT_PATH, PROJECT, config_file_shared, file_type):
    """
    create the data filenames which should be requested to the Snakefile
    with help of the previously loaded  manifest file -> this manifest file is
    the input for this rule! and must be already loaded beforehand!

    :param: OUTPUT_PATH: path for metilene pipeline outputs
    :type: OUTPUT_PATH: str
    :param: PROJECT: list of projects chosen
    :type: PROJECT: list of str

    with the UUID the clinical tables will be downloaded in the
    OUTPUT_PATH/PROJECT:

        * nationwidechildrens.org_clinical_patient_****.txt
        * nationwidechildrens.org_clinical_drug_****.txt
    """
    # get the UUID of all patient related meta tables out of the gdc_manifest:
    # read in the table:
    with open(config_file_shared, 'r') as f:
        df = pd.json_normalize((safe_load(f)))
    mani = df['manifest_file'].tolist()[0][0]

    mani_path = os.path.join(OUTPUT_PATH, 'metadata', mani)

    # filter on the project_id, this workds also for mutliple projects
    gdc_manifest_DF = pd.read_csv(
        mani_path, sep='\t').set_index(
            'project_id').loc[PROJECT, :].reset_index()

    # we just need the cols project_id and filename
    # this time we search pipeline specific datafiles, either htseq files or
    # HumanMethylation450 files
    gdc_manifest_temp = gdc_manifest_DF[
        gdc_manifest_DF[ 'filename'].str.contains(
            f'.*{file_type}.*', regex=True)].loc[:, ['project_id', 'filename']]

    # create absolute paths to the filenames where they should be saved to:
    if file_type == 'htseq':
        exp_type = 'DESeq2'
    elif file_type == 'HumanMethylation450':
        exp_type = 'metilene'
    gdc_manifest_temp['filepath'] = OUTPUT_PATH + '/' + gdc_manifest_temp['project_id'] + f'/{exp_type}/data_files/' + gdc_manifest_temp['filename']
    aux_file_list = [i for i in gdc_manifest_temp['filepath']]
    return aux_file_list
