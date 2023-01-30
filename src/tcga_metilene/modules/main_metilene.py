# ##### main_metilene ############
'''
    within this, every parameter needed for the analysis is set

    :param: out_path: path for metilene_pipeline outputs
    :type: out_path: str
    :param: script_path: path to the metilene_pipeline repo
    :type: script_path: str
    :param: function: apply, if a specific function should be executed solely\
        (also multiple functions possible)
    :type: function: int
    :param: drugs: applied drug(s)
    :type: drugs: list of str
    :param: project: list of projects chosen
    :type: project: list of str
    :param: download_data: bool flag, whether raw data needs to be downloaded
    :type: download_data: bool
    :param: analyse_data: bool flag, whether the metilene analyses shall be\
        started
'''
from tcga_metilene.modules import methyl
from shared import download_with_api
import os
import click
import choose_therapy
import create_report
import lifeline_plots
import set_logger


# version = 'Version 1.0'
# to specify the SCRIPT_PATH with __file__ python 3.9 is mandatory, else just
# the relative path would be returned (main_metilene.py) and not the absolute
# path tho the file (here
# /homes/biertruck/gabor/phd/test_git_doc/methyl/main_metilene.py)
# the absolute path is needed to extract the path to all scripts!
SCRIPT_PATH = os.path.split(__file__)[0]
with open(os.path.join(SCRIPT_PATH, 'version.txt'), 'r') as f:
    version = f.readline().strip()


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(version)
    ctx.exit()


CWD = os.environ.get('PWD')


@click.command()
@click.option('--download_data', '-D', help='perform download and merging' +
              ' steps, without the metilene analisis', is_flag=True)
@click.option('--analyse_data', '-A',
              help='perform the metilene analysis, ' +
              'the methylation data for your projects' +
              ' must be downloaded in prior to choose this option',
              is_flag=True)
@click.option('--out_path', '-o', default=os.path.join(CWD, 'metilene_data'),
              show_default=True,
              help='path to save the result files')
# @click.option('--script_path', '-s', default=os.path.join(CWD),
#               # help='The Path, where the metilene_pipeline repo is located')
@click.option('--project', '-p', default=[], multiple=True,
              help='TCGA project to be applied. Any TCGA project can be' +
              ' chosen, like: ' +
              '-p TCGA-CESC -p TCGA-HNSC ...')
@click.option('--drugs', '-d', default=[], multiple=True, show_default=False,
              help='drug(s), like: -d drug1 -d drug2 ' +
              'or drugcombination(s), like: -d drug1,drug2')
@click.option('--function', '-f', multiple=True,
              help='running functions separately, important for the ' +
              'snakemake functionality. Look up the documentation for ' +
              'detailed description  of each function', type=int)
# @click.option('--snakerun', '-s', is_flag=True)
@click.option('--snakerun', '-s', is_flag=True, help='needed for the ' +
              'snakemake functionality, NOT to be set manually by the user')
@click.option('--cutoff', '-c',
              default=0, show_default=True, type=float,
              help='cutoff parameter: convert dead cases to alive ' +
              'cases if they outlive the cutoff parameter (in years) ' +
              'no conversion performed if default (or no) value applied')
@click.option('--version', '-v',
              help='printing out version information: {}'.format(version),
              is_flag=True, callback=print_version,
              expose_value=False, is_eager=True)
def call_with_options(out_path, project, drugs, function,
                      download_data, analyse_data, cutoff, snakerun):
    '''
    "metilene_pipeline" a tool to choose, harvest and analyse methylation data
    of the TCGA-projects with help of the metilene package.\n

    Build and activate the provided conda env:

        $ conda env create -f metilene_env.yaml

        $ conda activate metilene_pipeline

    call the script without any options to enter the interactive mode and set
    each option step by step:

        $ python main_metilene.py

    print help page:

        $ python main_metilene.py --help
    '''
    OUTPUT_PATH = out_path
    print("\nOUTPUT_PATH:\n", OUTPUT_PATH)
    # SCRIPT_PATH = script_path
    SCRIPT_PATH = os.path.split(__file__)[0]
    print("\nSCRIPT_PATH:\n", SCRIPT_PATH)
    # check if the main_metilene.py exist in the SCRIPT_PATH:
    if not os.path.isfile(os.path.join(SCRIPT_PATH, 'main_metilene.py')):
        print('there is no "main_metilene.py" in your SCRIPT_PATH ' +
              'please check your input, exiting now')
        os._exit(0)
    # print("\nSCRIPT_PATH:\n", SCRIPT_PATH)
    # the projects shall be applied in list form, afterwards the correct hash
    # pairs are handed over, the UUIDs are present here
    PROJECT_DRUG_PRESET = {
        'TCGA-CESC': '47dcf82b-de38-4921-aba0-3ff6e07d9959',
        'TCGA-LUSC': 'b5cceec9-8006-45e1-9f92-949c23d594bc',
        'TCGA-HNSC': '75d436b1-36db-470f-ae07-fec3b432336e',
        'TCGA-ESCA': 'b7f651bd-a017-4cc5-a553-dc51e801d065',
        'TCGA-BRCA': '1e6b79ff-9787-4cbe-b19d-ebabb6b43589',  # breast invasive
        # carcinoma
        'TCGA-GBM': '3b0a2c87-68cc-4c4d-8a1c-5a2d8dac2125',  # glioblastoma
        # multiforme
        'TCGA-OV': '09a4994f-d50e-423c-895e-0be2d6dd5762',  # ovarian serous
        # cystadenocarcinaoma
        'TCGA-LUAD': '59015dfd-045b-4a80-9b4a-3cc3a41707f0',  # lung
        # adenocarcinoma
        'TCGA-UCEC': '3d5ec43f-d487-42a5-b08f-2d8bacae6501',  # uterine corpus
        # endometrial carinoma
        'TCGA-KIRC': '619ada6f-7f67-4c0e-ad0a-93a44bfb14cd',  # kindney renal
        # clear cell carcinoma
        'TCGA-LGG': 'c1f6276e-6e34-4c18-b599-88ffe9ad20bc',  # brain lower
        # grade glioma
        'TCGA-THCA': 'b64997f4-6785-4b05-9cf8-ce82315420e2',  # thyroid
        # carcinoma
        'TCGA-PRAD': '8be27284-6da7-49c0-ab9e-fc408a687cef',  # prostate
        # adenocarcinoma
        'TCGA-SKCM': '8da29d73-2299-49f0-9740-d04b3793fd8d',  # skin cutaneous
        # melanoma
        'TCGA-COAD': 'dae92461-d098-40a3-997e-1a2d3e5b8dfe',  # colon
        # adenocarcinoma
        'TCGA-STAD': 'ff7cdc16-ef2c-412c-9a32-53235c04cc08',  # stomach
        # adenocarcinoma
        'TCGA-BLCA': '9014e5dd-78fa-4507-bb07-156370362990',  # bladder
        # urothelial carcinoma (Squamous Cell Neoplasms)
        'TCGA-LIHC': 'e14d0b40-f21d-4fc2-ab11-9ebf8703f86f',  # liver
        # hepatocellular carcinoma
        'TCGA-KIRP': '8a05c3dd-8fb7-454b-ae84-fb6ac336411f',  # kidney renal
        # papillary cell carcinoma
        'TCGA-SARC': 'a0595964-3b4c-4902-9025-9da3d69fd87f',  # sarcoma
        'TCGA-ESCA': 'b7f651bd-a017-4cc5-a553-dc51e801d065',  # Esophageal
        # Carcinoma (Squamous Cell Neoplasms)
        'TCGA-PAAD': '830c4482-f8e3-409b-ab3f-3dbf226b7634',  # pancreatic
        # adenocarcinoma
        'TCGA-PCPG': '8fa3c739-15e3-490d-90d1-2be157608d78',  #
        # pheochromocytoma and paraganglioma
        'TCGA-READ': '4bab4ce6-239f-468d-8bbb-1cb69532f85a',  # rectum
        # adenocarcinoma
        'TCGA-TGCT': '64130c9e-7ed0-4908-a9b4-1a1199d0d603',  # testicular
        # germcelltumors
        'TCGA-THYM': 'c3c455d1-e5f0-4315-9b1b-f261231d88a6',  # thymoma
        'TCGA-KICH': '0287f41f-1384-4ac8-a9f2-2e82ce2b11c7',  # kidney
        # chromophobe
        'TCGA-ACC': '4c259f88-b551-4e17-b78f-4c1249b288da',  # adrenochordical
        # carcinoma
        'TCGA-MESO': 'd5115b64-b1a6-4b86-ac06-c38995fa95df',  # mesothelioma
        'TCGA-UVM': 'e11025b4-3427-4023-b52b-8e3945565df8',  # uveal melanoma
        'TCGA-DLBC': 'f0884d78-864f-4aa0-a1ea-3cb281d1226a',  # lymphoid
        # neoplasm diffuse large b-cell lymphoma
        'TCGA-UCS': '1629b69f-16ce-40f4-9f40-5f3a2c1db8fb',  # uterine
        # carcinoma
        'TCGA-CHOL': '271e2937-7860-42a7-8c2b-2c478b3d7d3d'
        # cholangiocarcinoma
    }

    PROJECT_PATIENT_PRESET = {
        'TCGA-CESC': '2efbcfd6-0284-4981-820f-2d23455fe6be',
        'TCGA-LUSC': '7aaee8d4-3bea-4dfc-b392-c7aaa3e9a550',
        'TCGA-HNSC': '516e1008-c07e-4420-9b7b-dbbd219cae5e',
        'TCGA-ESCA': '916e5115-af71-427d-82de-231ac023ce41',
        'TCGA-BRCA': '8162d394-8b64-4da2-9f5b-d164c54b9608',
        'TCGA-GBM': '57683e22-a8ea-4eca-bfcf-f708cf459546',
        'TCGA-OV': '30c149ac-9ac2-4f51-88a3-68bb4afb50a9',
        'TCGA-LUAD': '42bf5eb2-bc49-45be-b18a-290f712b006c',
        'TCGA-UCEC': '23529ac8-43ca-408a-86ec-cd555f58edbd',
        'TCGA-KIRC': '8c43d640-c32d-439c-8c46-22c46e8f9ba0',
        'TCGA-LGG': 'f3a1bc62-9552-4553-b318-7d9c21d21ce7',
        'TCGA-THCA': '5c81bd45-ab70-4bcd-895c-5fe916f272d6',
        'TCGA-PRAD': '06efd272-a76f-4703-98b8-dfa751c0f019',
        'TCGA-SKCM': '58cbbc07-5ec4-47c7-9295-11ccbf7693f4',
        'TCGA-COAD': '4060482f-eedf-4959-97f1-f8b6c529c368',
        'TCGA-STAD': '0089d221-5807-47f1-a382-1e2e336df201',
        'TCGA-BLCA': 'eaa71705-960a-4abd-b5d7-f5fdc0d0c5af',
        'TCGA-LIHC': '88b30bf9-b937-4381-8684-42bbcce98fa0',
        'TCGA-KIRP': '5940e28d-75ee-4ecf-b29d-317f5ef89a05',
        'TCGA-SARC': 'ee49cde2-db41-47c9-8c41-99662f96d3c2',
        'TCGA-PAAD': '92be86dd-9f8c-4dfe-ba38-22cef2f6a9c7',
        'TCGA-PCPG': 'bbb2453e-3794-472b-a061-789513b36e61',
        'TCGA-READ': 'a4e80133-a7b1-4adc-8c69-1178c1871c78',
        'TCGA-TGCT': 'e086b524-696f-477c-87db-fa2a5ddb2a1c',
        'TCGA-THYM': '7cca5722-26cf-4ac8-a4a6-b803459f1861',
        'TCGA-KICH': '2d041399-6315-4b3f-a3f8-dd6088c452f1',
        'TCGA-ACC': '909f0719-b165-4beb-9fb3-f21aff071469',
        'TCGA-MESO': 'ddeb78fc-91e8-4e00-a89a-7c63f9bfc03d',
        'TCGA-UVM': '14bdc241-ad7e-421e-b49c-8ab5ed23a5d2',
        'TCGA-DLBC': '03fbcbed-0498-46b7-937b-51f77622776c',
        'TCGA-UCS': 'e1ca87d5-3964-4880-a3eb-30e4e2d2b26f',
        'TCGA-CHOL': '8f5ff93a-08b6-4d86-bbe2-7a02127b2f48',
    }

    # ####### get here the projects and drug through Choose_Therapy
    # depending of the terminal call, the drugs and project are set, if not,
    # request them interactively, sort both:
    if len(project) == 0:
        PROJECT = choose_therapy.Choose_project()
        PROJECT = sorted(map(str.upper, PROJECT))
    else:
        PROJECT = sorted(map(str.upper, project))
    if len(drugs) == 0:
        DRUGS = choose_therapy.Choose_drugs(SCRIPT_PATH, PROJECT)
    else:
        DRUGS = sorted(map(str.lower, drugs))

    # with the project list create the Project drug uuid hash:
    PROJECT_DRUG_UUID = {}
    for proj in PROJECT:
        if proj in PROJECT_DRUG_PRESET.keys():
            PROJECT_DRUG_UUID.update({proj: PROJECT_DRUG_PRESET[proj]})
        else:
            print('the applied project {} is not available'.format(proj))
            os._exit(0)
    PROJECT_PATIENT_UUID = {}
    for proj in PROJECT:
        if proj in PROJECT_PATIENT_PRESET.keys():
            PROJECT_PATIENT_UUID.update({proj: PROJECT_PATIENT_PRESET[proj]})

    print('PROJECT:\t', PROJECT)
    print('DRUGS:\t\t', DRUGS)
    FILE_TYPE = "dna methylation"
    # FILE_TYPE = "HTSeq - Counts"   # choose between htseq, FPKM-UQ, FPKM.txt

    # if no -A and -D specified, ask here whether to just download raw data or
    # do both download and analyse steps:
    if download_data and analyse_data:
        function = (100,)
    if download_data and not analyse_data:
        function = (1, 2, 3, 4)  # until creating table:
        # meta_info_druglist_merged_drugs_combined.tsv and download all data
    if analyse_data and not download_data:  # all the rest
        function = (5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)
    # if the first 5 fct are skipped (because those projects are
    # downloaded already), the info in the snakemake_config.yaml would miss,
    # because all those log entries are skipped also. if just -A is set,
    # we must look up all those files and write it in the config

    # whats set so far are the projects and drugs, still needed analyse options
    # and OUTPUT_PATH
    if not analyse_data and not download_data and not function:
        OUTPUT_PATH, function, analyse_data,\
            download_data = choose_therapy.Choose_path_and_option(
                OUTPUT_PATH, PROJECT, DRUGS, function, SCRIPT_PATH,
                analyse_data, download_data)

    OUTPUT_PATH = OUTPUT_PATH.rstrip(os.sep)

    ####################################
    met_opt = {'-m': '3', '-M': '1000', '-d': '0.03'}
    # met_opt = {'-m': '3', '-M': '1000'}
    met_opt_list = []
    for i in met_opt:
        met_opt_list.append(i)
        met_opt_list.append(met_opt[i])
    #############################################
    met_dir = '_'.join(met_opt_list)
    met_dir = 'metilene_' + met_dir
    # print("\nmet_dir:\n", met_dir)

    # every fct is a selfcontained step, to apply which fct is needed in
    # snakemake, enum the fct and apply the nr. to the options:
    # if min(function) >= 1 and max(function) <= 10 or 100 in function:
    DRUGS_title = '_'.join(sorted(map(str.lower, DRUGS)))
    if cutoff > 0:
        DRUGS_title = DRUGS_title + '_' + str(cutoff)
    projects = []
    for project in PROJECT_DRUG_UUID:
        projects.append(project)
    PROJECT_title = '_'.join(sorted(map(str.upper, projects)))
    PROJECT_list = projects
    api_manifest = 'gdc_manifest_20211029_data_release_31.0_active.tsv.gz'
    gtf_UUID = '25aa497c-e615-4cb7-8751-71f744f9691f'
    workflow_type = 'Liftover'
    # in case that just -A is set, go through the single project dirs, check
    # the existence of all files needed, and
    # create the logfile in the PROJECTS_title/DRUGS_title/ dir
    if analyse_data and not download_data:
        methyl.create_log_for_A(OUTPUT_PATH, PROJECT_list, PROJECT_title,
                                DRUGS_title, api_manifest)
    for PROJECT in PROJECT_DRUG_UUID:
        # log file is written to OUT/singl_proj/DRUGS_title
        os.makedirs(os.path.join(
            OUTPUT_PATH, PROJECT, DRUGS_title), exist_ok=True)
        logger = set_logger.set_logger(
            OUTPUT_PATH, PROJECT, DRUGS_title)
        logger.info('Version:\t{}'.format(version))
        if 1 in function or 100 in function:
            download_with_api.download_GDC_manifest(PROJECT, FILE_TYPE,
                                                    OUTPUT_PATH, DRUGS_title,
                                                    SCRIPT_PATH, PROJECT_title,
                                                    snakerun, api_manifest,
                                                    gtf_UUID)

        if 2 in function or 100 in function:
            download_with_api.download_clinical_tables(PROJECT, OUTPUT_PATH,
                                                       SCRIPT_PATH,
                                                       DRUGS_title,
                                                       api_manifest)

        if 3 in function or 100 in function:
            download_with_api.create_merged_metatable(OUTPUT_PATH, PROJECT,
                                                      DRUGS_title,
                                                      workflow_type,
                                                      SCRIPT_PATH)

        if 4 in function or 100 in function:
            download_with_api.sep_down_data_files(OUTPUT_PATH, PROJECT,
                                                  DRUGS_title)
        # # with the joined both_table, create the summary table
        if 5 in function or 100 in function:
            # catch here the -A option TODO -> if the summary table is already
            # present this step an all the following shall just be run if the
            # fct is explicitly called with -f 5, otherwise the analyses are
            # considedered as already performed -> continue to the next project
            methyl.create_summary_table(OUTPUT_PATH, PROJECT, DRUGS_title)
            # creates the complete summary table with every
            # drugcombination in it
            # pay attention to that the files have all the same dim!
            # DONE!!!
            # with bash:
            # for i in *; do
            # filesize=$( wc -l <"$i"); if [ "$filesize" -lt 450000 ];
            # then mv $i ../bad_data/ ; fi; done
            # temporary solution, find a better one!

            # it is supposed, that if the table_for_metilene is created, the
            # complete analyse part for single projects is performed, s.t. if
            # this file exists, the single project analyse of this project can
            # be skipped: 'summary_dead_alive.tsv' in  OUTPUT_PATH, PROJECT,
            # DRUGS_title
            # TODO, uncomment this statement when finished:
        # if os.path.exists(os.path.join(
        #         OUTPUT_PATH, PROJECT, DRUGS_title,
        #         'summary_dead_alive.tsv')):
        #     continue
        if 6 in function or 100 in function:
            methyl.create_table_for_metilene(
                OUTPUT_PATH,
                PROJECT, DRUGS, DRUGS_title, met_dir, cutoff)
            # creates the filtered table with just specific drug
            # combinations in it

        if 7 in function or 100 in function:
            methyl.provide_metilene_table(
                OUTPUT_PATH,
                DRUGS, SCRIPT_PATH, met_opt_list, met_dir, PROJECT,
                DRUGS_title)
            # call metilene and the sppl scripts on the before created
            # tables

        if 8 in function or 100 in function:
            methyl.call_bed_intersect(
                OUTPUT_PATH, DRUGS, met_opt_list, met_dir, PROJECT,
                DRUGS_title)
            # call bed intersect for the single projects

        if 9 in function or 100 in function:
            methyl.create_plots(
                OUTPUT_PATH, DRUGS, met_opt_list, met_dir, PROJECT,
                DRUGS_title)

        # boxplot for ranges
        if 10 in function or 100 in function:
            methyl.create_mean_boxplot(
                OUTPUT_PATH, DRUGS, met_opt_list, met_dir, PROJECT,
                DRUGS_title)

    # don't call multi project fcts if just one project is applied!, but create
    # the  snake config and the report file
    # the report must be created before the config, otherwise it would miss
    if len(PROJECT_DRUG_UUID) == 1:
        if 17 in function or 100 in function:
            # lifeline_plots.prepare_lifeline_plots(OUTPUT_PATH, DRUGS,
                                                  # SCRIPT_PATH,
                                                  # PROJECT_DRUG_UUID, met_dir,
                                                  # cutoff, DRUGS_title)
            # lifeline_plots.prepare_complement_lifeline_plots(OUTPUT_PATH,
                                                             # DRUGS,
                                                             # SCRIPT_PATH,
                                                             # PROJECT_DRUG_UUID,
                                                             # met_dir, cutoff,
                                                             # DRUGS_title)
            lifeline_plots.create_lifeline_plots(OUTPUT_PATH, DRUGS,
                                                 SCRIPT_PATH,
                                                 PROJECT_DRUG_UUID, met_dir,
                                                 cutoff, DRUGS_title)
            create_report.create_report_pdf(OUTPUT_PATH, DRUGS, SCRIPT_PATH,
                                            PROJECT_DRUG_UUID, met_dir,
                                            cutoff, DRUGS_title)
        if 16 in function or 100 in function:
            methyl.create_snake_config(OUTPUT_PATH, PROJECT_title,
                                       DRUGS_title, PROJECT_list, DRUGS,
                                       SCRIPT_PATH, cutoff)
            os._exit(0)

    # check here, if the cutoff is considered automatically TODO
    if 11 in function or 100 in function:
        methyl.create_table_all_projects(
            OUTPUT_PATH, DRUGS, PROJECT_DRUG_UUID, DRUGS_title, met_dir,
            cutoff)

    # call metilene for every project and for the aggregated summary:
    if 12 in function or 100 in function:
        methyl.provide_metilene_table(
            OUTPUT_PATH, DRUGS, SCRIPT_PATH, met_opt_list, met_dir,
            PROJECT_DRUG_UUID, DRUGS_title)

    # call bed intersect for every project
    if 13 in function or 100 in function:
        methyl.call_bed_intersect(
            OUTPUT_PATH, DRUGS, met_opt_list, met_dir, PROJECT_DRUG_UUID,
            DRUGS_title)

    if 14 in function or 100 in function:
        methyl.create_plots(OUTPUT_PATH, DRUGS, met_opt_list, met_dir,
                            PROJECT_DRUG_UUID, DRUGS_title)

    if 15 in function or 100 in function:
        methyl.create_mean_boxplot(
            OUTPUT_PATH, DRUGS, met_opt_list, met_dir, PROJECT_DRUG_UUID,
            DRUGS_title)

    # the report must be created before the config, otherwise it would missed
    # in the ovarall config
    if len(PROJECT_DRUG_UUID) > 1:
        if 17 in function or 100 in function:
            print(f'creating report for {PROJECT_DRUG_UUID}')
            # make the lifeline plots for the multiproject:
            lifeline_plots.prepare_lifeline_plots(OUTPUT_PATH, DRUGS,
                                                  SCRIPT_PATH,
                                                  PROJECT_DRUG_UUID, met_dir,
                                                  cutoff, DRUGS_title)
            lifeline_plots.prepare_complement_lifeline_plots(OUTPUT_PATH,
                                                             DRUGS,
                                                             SCRIPT_PATH,
                                                             PROJECT_DRUG_UUID,
                                                             met_dir, cutoff,
                                                             DRUGS_title)
            lifeline_plots.create_lifeline_plots(OUTPUT_PATH, DRUGS,
                                                 SCRIPT_PATH,
                                                 PROJECT_DRUG_UUID, met_dir,
                                                 cutoff, DRUGS_title)
            create_report.create_report_pdf(OUTPUT_PATH, DRUGS, SCRIPT_PATH,
                                            PROJECT_DRUG_UUID, met_dir, cutoff,
                                            DRUGS_title)
            # if multi projects are applied, create report in each single proj
            # dir:
            for PROJECT in PROJECT_DRUG_UUID:
                print(f'creating report for {PROJECT}')
                lifeline_plots.prepare_lifeline_plots(OUTPUT_PATH, DRUGS,
                                                    SCRIPT_PATH,
                                                    PROJECT,
                                                    met_dir, cutoff,
                                                    DRUGS_title)
                lifeline_plots.prepare_complement_lifeline_plots(OUTPUT_PATH,
                                                                DRUGS,
                                                                SCRIPT_PATH,
                                                                PROJECT,
                                                                met_dir,
                                                                cutoff,
                                                                DRUGS_title)
                lifeline_plots.create_lifeline_plots(OUTPUT_PATH, DRUGS,
                                                    SCRIPT_PATH,
                                                    PROJECT,
                                                    met_dir, cutoff,
                                                    DRUGS_title)
                create_report.create_report_pdf(OUTPUT_PATH, DRUGS,
                                                SCRIPT_PATH, PROJECT,
                                                met_dir, cutoff, DRUGS_title)

        if 16 in function or 100 in function:
            methyl.create_snake_config(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                    PROJECT_list, DRUGS, SCRIPT_PATH, cutoff)
        # ##### adding a sanitiy check for good canditates:
        #  methyl.sanity_check(OUTPUT_PATH, DRUGS)

        # ### check the single runs for each project of metilene and put them
        # together for intersect and visualisation
        # methyl.sanity_check_2(OUTPUT_PATH, DRUGS, PROJECT_DRUG_UUID)


if __name__ == '__main__':
    call_with_options()
