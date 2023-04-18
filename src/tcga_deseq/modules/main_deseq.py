import snakemake
import os
from tcga_deseq.modules import create_summary_table
from tcga_deseq.modules import create_deseq_output


"""
coming from src/shared/modules/main.py
generating the files requested through the deseq Snakefile in
src/tcga_deseq/Snakefile
"""

def entry_fct(OUTPUT_PATH, PROJECT, DRUGS, Snakemake_all_files, cutoffs,
              threshold, cores):

    SCRIPT_PATH = os.path.split(__file__)[0]

    PROJECTS = []
    if len(PROJECT) > 1:
        PROJECTS.extend(PROJECT)
        PROJECTS.append('_'.join(sorted([x.upper() for x in PROJECT])))

    DRUG_str = '_'.join(DRUGS)

    shared_workdir = os.path.join(
        os.path.split(os.path.split(SCRIPT_PATH)[0])[0], 'shared')
    Snakefile = os.path.join(os.path.split(SCRIPT_PATH)[0], 'Snakefile')

    summary_tables_list = create_summary_table.create_summary_table(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs)
    Snakemake_all_files = summary_tables_list

    deseq_output_list = create_deseq_output.create_deseq_output(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs)
    Snakemake_all_files = Snakemake_all_files + deseq_output_list


    snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
                        workdir=shared_workdir, cores=cores, forceall=False,
                        force_incomplete=True, dryrun=False, use_conda=True)

# #!/usr/bin/env python3.8
# # ##### DESEQ2 ############
# '''
#     within this, every parameter needed for the analysis is set

#     :param: out_path: path for DESeq2 pipeline outputs
#     :type: out_path: str
#     :param: script_path: path to the DESeq2_pipeline repo
#     :type: script_path: str
#     :param: function: apply, if a specific function should be executed solely\
#         (also multiple functions possible)
#     :type: function: int
#     :param: drugs: applied drug(s)
#     :type: drugs: list of str
#     :param: project: list of projects chosen
#     :type: project: list of str
#     :param: download_data: bool flag, whether raw data needs to be downloaded
#     :type: download_data: bool
#     :param: analyse_data: bool flag, whether the deseq analyses shall be\
#         started
#     :param: threshold: parameter for the lifeline plots helping for the\
#         classification of expression data
#     :type: threshold: int
# '''


# import create_matrix_new
# import lifeline_summary_test_2
# import walk_all_drug_frequency
# import choose_therapy
# import create_report
# import os
# import click
# import download_with_api

# # version = 'Version 1.0'
# SCRIPT_PATH = os.path.split(__file__)[0]
# with open(os.path.join(SCRIPT_PATH, 'version.txt'), 'r') as f:
#     version = f.readline().strip()


# def print_version(ctx, param, value):
#     if not value or ctx.resilient_parsing:
#         return
#     click.echo(version)
#     ctx.exit()
#     # check if the main_deseq.py exist in the SCRIPT_PATH:
#     # if not os.path.isfile(os.path.join(SCRIPT_PATH, 'main_deseq.py')):
#     #    # print('there is no "main_deseq.py" in your SCRIPT_PATH ' +
#     #          # 'please check your input, exiting now')
#     #    # os._exit(0)


# CWD = os.environ.get('PWD')


# @click.command()
# @click.option('--download_data', '-D',
#               help='perform download and merging steps, without the DESeq2' +
#               ' analysis', is_flag=True)
# @click.option('--analyse_data', '-A',
#               help='perform the DESeq2 analysis (the download step must' +
#               ' be completed in prior for your chosen projects)',
#               is_flag=True)
# @click.option('--out_path', '-o', default=os.path.join(CWD, 'DESeq2_data'),
#               show_default=True,
#               help='The Path, where the results are saved to')
# @click.option('--drugs', '-d', default=[],
#               multiple=True, show_default=False,
#               help='drug(s), like: -d drug1 -d drug2 ' +
#               'or drugcombination(s), like: -d drug1,drug2')
# @click.option('--project', '-p', default=[], multiple=True,
#               help='Project, that shall be applied, choose between ' +
#               'TCGA-CESC, TCGA-HNSC, TCGA-LUSC, TCGA-ESCA,' +
#               ' or combinations out of them like -p project1 -p project2')
# @click.option('--function', '-f', multiple=True,
#               help='running functions separately, important for the ' +
#               'snakemake functionality. Look up the documentation for ' +
#               'detailed description  of each function', type=int)
# @click.option('--threshold', '-t',
#               default=(0,), show_default=True, type=int, multiple=True,
#               help='threshold for a higher division of the kaplan meier plots')
# @click.option('--snakerun', '-s', is_flag=True, help='needed for the ' +
#               'snakemake functionality, NOT to be set manually by the user')
# @click.option('--cutoff', '-c',
#               default=0, show_default=True, type=float,
#               help='cutoff parameter: convert dead cases to alive ' +
#               'cases if they outlive the cutoff parameter (in years) ' +
#               'no conversion performed if default (or no) value applied')
# # @click.option('--script_path', '-s', default=os.path.split(__file__)[0],
# #              # help='your script_path is set to {}'.format(
# #                  # os.path.split(__file__)[0]))
# @click.option('--version', '-v',
#               help='printing out version information: {}'.format(version),
#               is_flag=True, callback=print_version,
#               expose_value=False, is_eager=True)
# def call_with_options(out_path, function, drugs, project,
#                       download_data, analyse_data, threshold, snakerun,
#                       cutoff):
#     '''
#     "DESeq2_pipeline" a tool to choose, harvest and analyse expression data of
#     the TCGA-projects with help of the DESeq2 R package.\n

#     Build and activate the provided conda env:

#         $ conda env create -f deseq_env.yaml

#         $ conda activate deseq_pipeline

#     Call the script without any options to enter the interactive mode and set
#     each option step by step:

#         $ python main_deseq.py

#     print help page:

#         $ python main_deseq.py --help
#     '''
#     # # if not function:
#     #     # print('\nplease provide the -D and/or the -A function\n')
#     #     # ctx = click.get_current_context()
#     #     # click.echo(ctx.get_help())
#     #     # ctx.exit()
#     # check if the OUTPU_PATH is a absolute path, with leading / or if it is a
#     # relative path without leading /, if not absolut, add the cwd in prior
#     OUTPUT_PATH = out_path.rstrip(os.path.sep)
#     if OUTPUT_PATH[1] != os.path.sep:
#         OUTPUT_PATH = os.path.join(os.getcwd(), OUTPUT_PATH)

#     # OUTPUT_PATH = os.path.join(os.getcwd(),)
#     print("\nOUTPUT_PATH:\n", OUTPUT_PATH)
#     # SCRIPT_PATH = script_path.rstrip(os.path.sep)
#     # if SCRIPT_PATH[1] != os.path.sep:
#     #    # SCRIPT_PATH = os.path.join(os.getcwd(), SCRIPT_PATH)
#     SCRIPT_PATH = os.path.split(__file__)[0]
#     print("\nSCRIPT_PATH:\n", SCRIPT_PATH)
#     PROJECT_DRUG_PRESET = {
#         'TCGA-CESC': '47dcf82b-de38-4921-aba0-3ff6e07d9959',
#         # Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma
#         'TCGA-LUSC': 'b5cceec9-8006-45e1-9f92-949c23d594bc',
#         'TCGA-HNSC': '75d436b1-36db-470f-ae07-fec3b432336e',
#         # head and neck squamous cell carcinoma
#         'TCGA-ESCA': 'b7f651bd-a017-4cc5-a553-dc51e801d065',
#         # Esophageal Carcinoma (Squamous Cell Neoplasms)
#         'TCGA-BRCA': '1e6b79ff-9787-4cbe-b19d-ebabb6b43589',
#         # breast invasive carcinoma (Squamous Cell Neoplasms)
#         'TCGA-GBM': '3b0a2c87-68cc-4c4d-8a1c-5a2d8dac2125',
#         #  glioblastoma multiforme
#         'TCGA-OV': '09a4994f-d50e-423c-895e-0be2d6dd5762',
#         #  ovarian serous cystadenocarcinaoma
#         'TCGA-LUAD': '59015dfd-045b-4a80-9b4a-3cc3a41707f0',
#         #  lung adenocarcinoma
#         'TCGA-UCEC': '3d5ec43f-d487-42a5-b08f-2d8bacae6501',
#         #  uterine corpus endometrial carinoma
#         'TCGA-KIRC': '619ada6f-7f67-4c0e-ad0a-93a44bfb14cd',
#         #  kindney renal clear cell carcinoma
#         'TCGA-LGG': 'c1f6276e-6e34-4c18-b599-88ffe9ad20bc',
#         #  brain lower grade glioma
#         'TCGA-THCA': 'b64997f4-6785-4b05-9cf8-ce82315420e2',
#         #  thyroid carcinoma
#         'TCGA-PRAD': '8be27284-6da7-49c0-ab9e-fc408a687cef',
#         #  prostate adenocarcinoma
#         'TCGA-SKCM': '8da29d73-2299-49f0-9740-d04b3793fd8d',
#         #  skin cutaneous melanoma
#         'TCGA-COAD': 'dae92461-d098-40a3-997e-1a2d3e5b8dfe',
#         #  colon adenocarcinoma
#         'TCGA-STAD': 'ff7cdc16-ef2c-412c-9a32-53235c04cc08',
#         #  stomach adenocarcinoma
#         'TCGA-BLCA': '9014e5dd-78fa-4507-bb07-156370362990',
#         #  bladder urothelial carcinoma (Squamous Cell Neoplasms)
#         'TCGA-LIHC': 'e14d0b40-f21d-4fc2-ab11-9ebf8703f86f',
#         #  liver hepatocellular carcinoma
#         'TCGA-KIRP': '8a05c3dd-8fb7-454b-ae84-fb6ac336411f',
#         #  kidney renal papillary cell carcinoma
#         'TCGA-SARC': 'a0595964-3b4c-4902-9025-9da3d69fd87f',
#         #  sarcoma
#         'TCGA-PAAD': '830c4482-f8e3-409b-ab3f-3dbf226b7634',
#         #  pancreatic adenocarcinoma
#         'TCGA-PCPG': '8fa3c739-15e3-490d-90d1-2be157608d78',
#         #  pheochromocytoma and paraganglioma
#         'TCGA-READ': '4bab4ce6-239f-468d-8bbb-1cb69532f85a',
#         #  rectum adenocarcinoma
#         'TCGA-TGCT': '64130c9e-7ed0-4908-a9b4-1a1199d0d603',
#         #  testicular germcelltumors
#         'TCGA-THYM': 'c3c455d1-e5f0-4315-9b1b-f261231d88a6',
#         #  thymoma
#         'TCGA-KICH': '0287f41f-1384-4ac8-a9f2-2e82ce2b11c7',
#         #  kidney chromophobe
#         'TCGA-ACC': '4c259f88-b551-4e17-b78f-4c1249b288da',
#         #  adrenochordical carcinoma
#         'TCGA-MESO': 'd5115b64-b1a6-4b86-ac06-c38995fa95df',
#         #  mesothelioma
#         'TCGA-UVM': 'e11025b4-3427-4023-b52b-8e3945565df8',
#         #  uveal melanoma
#         'TCGA-DLBC': 'f0884d78-864f-4aa0-a1ea-3cb281d1226a',
#         #  lymphoid neoplasm diffuse large b-cell lymphoma
#         'TCGA-UCS': '1629b69f-16ce-40f4-9f40-5f3a2c1db8fb',
#         #  uterine carcinoma
#         'TCGA-CHOL': '271e2937-7860-42a7-8c2b-2c478b3d7d3d',
#         #  cholangiocarcinoma
#         }

#     PROJECT_PATIENT_PRESET = {
#         'TCGA-CESC': '2efbcfd6-0284-4981-820f-2d23455fe6be',
#         'TCGA-LUSC': '7aaee8d4-3bea-4dfc-b392-c7aaa3e9a550',
#         'TCGA-HNSC': '516e1008-c07e-4420-9b7b-dbbd219cae5e',
#         'TCGA-ESCA': '916e5115-af71-427d-82de-231ac023ce41',
#         'TCGA-BRCA': '8162d394-8b64-4da2-9f5b-d164c54b9608',
#         'TCGA-GBM': '57683e22-a8ea-4eca-bfcf-f708cf459546',
#         'TCGA-OV': '30c149ac-9ac2-4f51-88a3-68bb4afb50a9',
#         'TCGA-LUAD': '42bf5eb2-bc49-45be-b18a-290f712b006c',
#         'TCGA-UCEC': '23529ac8-43ca-408a-86ec-cd555f58edbd',
#         'TCGA-KIRC': '8c43d640-c32d-439c-8c46-22c46e8f9ba0',
#         'TCGA-LGG': 'f3a1bc62-9552-4553-b318-7d9c21d21ce7',
#         'TCGA-THCA': '5c81bd45-ab70-4bcd-895c-5fe916f272d6',
#         'TCGA-PRAD': '06efd272-a76f-4703-98b8-dfa751c0f019',
#         'TCGA-SKCM': '58cbbc07-5ec4-47c7-9295-11ccbf7693f4',
#         'TCGA-COAD': '4060482f-eedf-4959-97f1-f8b6c529c368',
#         'TCGA-STAD': '0089d221-5807-47f1-a382-1e2e336df201',
#         'TCGA-BLCA': 'eaa71705-960a-4abd-b5d7-f5fdc0d0c5af',
#         'TCGA-LIHC': '88b30bf9-b937-4381-8684-42bbcce98fa0',
#         'TCGA-KIRP': '5940e28d-75ee-4ecf-b29d-317f5ef89a05',
#         'TCGA-SARC': 'ee49cde2-db41-47c9-8c41-99662f96d3c2',
#         'TCGA-PAAD': '92be86dd-9f8c-4dfe-ba38-22cef2f6a9c7',
#         'TCGA-PCPG': 'bbb2453e-3794-472b-a061-789513b36e61',
#         'TCGA-READ': 'a4e80133-a7b1-4adc-8c69-1178c1871c78',
#         'TCGA-TGCT': 'e086b524-696f-477c-87db-fa2a5ddb2a1c',
#         'TCGA-THYM': '7cca5722-26cf-4ac8-a4a6-b803459f1861',
#         'TCGA-KICH': '2d041399-6315-4b3f-a3f8-dd6088c452f1',
#         'TCGA-ACC': '909f0719-b165-4beb-9fb3-f21aff071469',
#         'TCGA-MESO': 'ddeb78fc-91e8-4e00-a89a-7c63f9bfc03d',
#         'TCGA-UVM': '14bdc241-ad7e-421e-b49c-8ab5ed23a5d2',
#         'TCGA-DLBC': '03fbcbed-0498-46b7-937b-51f77622776c',
#         'TCGA-UCS': 'e1ca87d5-3964-4880-a3eb-30e4e2d2b26f',
#         'TCGA-CHOL': '8f5ff93a-08b6-4d86-bbe2-7a02127b2f48',
#     }
#     # ####### get here the projects and drug through Choose_Therapy
#     # depending of the terminal call, the drugs and project are set, if not,
#     # request them interactively, sort both:
#     if len(project) == 0:
#         PROJECT = choose_therapy.Choose_project()
#         PROJECT = sorted(map(str.upper, PROJECT))
#     else:
#         PROJECT = sorted(map(str.upper, project))
#     if len(drugs) == 0:
#         DRUGS = choose_therapy.Choose_drugs(SCRIPT_PATH, PROJECT)
#     else:
#         DRUGS = sorted(map(str.lower, drugs))
#         # DRUGS = ','.join(sorted(DRUGS))

#     # with the project list create the Project DRUG uuid hash:
#     PROJECT_DRUG_UUID = {}
#     for proj in PROJECT:
#         if proj in PROJECT_DRUG_PRESET.keys():
#             PROJECT_DRUG_UUID.update({proj: PROJECT_DRUG_PRESET[proj]})
#         else:
#             print('the applied project {} is not available'.format(proj))
#             os._exit(0)
#     # with the project list create the Project PATIENT uuid hash:
#     PROJECT_PATIENT_UUID = {}
#     for proj in PROJECT:
#         if proj in PROJECT_PATIENT_PRESET.keys():
#             PROJECT_PATIENT_UUID.update({proj: PROJECT_PATIENT_PRESET[proj]})

#     print('PROJECT:\t', PROJECT)
#     print('DRUGS:\t\t', DRUGS)
#     FILE_TYPE = "HTSeq - Counts"   # choose between htseq, FPKM-UQ, FPKM.txt

#     # if no -A and -D specified, ask here whether to just download raw data or
#     # do both download and analyse steps:
#     if download_data and analyse_data:
#         function = (100,)
#     if download_data and not analyse_data:
#         function = (1, 2, 3, 4, 5)
#     if analyse_data and not download_data:
#         function = (7, 8, 9, 10, 11, 12, 13, 14, 15)  # 11 -> lifelines
#         # 12 -> walk all drug # 13 -> create report  14 ->
#         # drug_frequency_all_single projects
#     # TODO: ask for analysis and download steps

#     # whats set so far are the projects and drugs, still needed analyse options
#     # and OUTPUT_PATH
#     if not analyse_data and not download_data and not function:
#         OUTPUT_PATH, function, analyse_data,\
#             download_data, threshold = choose_therapy.Choose_path_and_option(
#                 OUTPUT_PATH, PROJECT, DRUGS, function, SCRIPT_PATH,
#                 analyse_data, download_data, threshold)
#     threshold = tuple(set(threshold))

#     OUTPUT_PATH = OUTPUT_PATH.rstrip(os.path.sep)
#     if OUTPUT_PATH[1] != os.path.sep:
#         OUTPUT_PATH = os.path.join(os.getcwd(), OUTPUT_PATH)
#     # SCRIPT_PATH = script_path.rstrip(os.path.sep)
#     # if SCRIPT_PATH[1] != os.path.sep:
#         # SCRIPT_PATH = os.path.join(os.getcwd(), SCRIPT_PATH)

#     projects = []
#     for project in PROJECT_DRUG_UUID:
#         projects.append(project)
#     PROJECT_list = projects
#     DRUGS_title = '_'.join(sorted(map(str.lower, DRUGS)))
#     if cutoff > 0:
#         DRUGS_title = DRUGS_title + '_' + str(cutoff)

#     PROJECT_title = '_'.join(sorted(map(str.upper, PROJECT_list)))
#     workflow_type = 'HTSeq - Counts'
#     api_manifest = 'gdc_manifest_20211029_data_release_31.0_active.tsv.gz'
#     gtf_UUID = '25aa497c-e615-4cb7-8751-71f744f9691f'

#     # every fct is a selfcontained step, to apply which fct is needed in
#     # snakemake, enum the fct and apply the nr. to the options:
#     # if min(function) >= 1 and max(function) <= 11 or 100 in function:
#     # if min(function) >= 1 and max(function) <= 11 or 100 in function:

#     # in case that just -A is set, go through the single project dirs, check
#     # the existence of all files needed, and
#     # create the logfile in the PROJECTS_title/DRUGS_title/ dir
#     # TODO -> check that create_log_for_A functionality
#     # if analyse_data and not download_data:
#     #     create_matrix_new.create_log_for_A(
#     #         OUTPUT_PATH, PROJECT_list, PROJECT_title, DRUGS_title)
#     for PROJECT in PROJECT_DRUG_UUID:
#         os.makedirs(os.path.join(
#             OUTPUT_PATH, PROJECT, DRUGS_title), exist_ok=True)
#         logger = create_matrix_new.set_logger(
#             OUTPUT_PATH, PROJECT, DRUGS_title)
#         logger.info('Version:\t{}'.format(version))
#         if 1 in function or 100 in function:
#             download_with_api.download_GDC_manifest(
#                 PROJECT, FILE_TYPE, OUTPUT_PATH, DRUGS_title, SCRIPT_PATH,
#                 PROJECT_title, snakerun, api_manifest, gtf_UUID)

#         if 2 in function or 100 in function:
#             download_with_api.download_clinical_tables(PROJECT, OUTPUT_PATH,
#                                                        SCRIPT_PATH,
#                                                        DRUGS_title,
#                                                        api_manifest)

#         if 3 in function or 100 in function:
#             download_with_api.create_merged_metatable(OUTPUT_PATH, PROJECT,
#                                                       DRUGS_title,
#                                                       workflow_type,
#                                                       SCRIPT_PATH)

#         if 4 in function or 100 in function:
#             download_with_api.sep_down_data_files(OUTPUT_PATH, PROJECT,
#                                                   DRUGS_title)
#         # # with the joined both_table, create the summary table
#         if 5 in function or 100 in function:
#             # catch here the -A option TODO -> if the summary table is already
#             # present this step an all the following shall just be run if the
#             # fct is explicitly called with -f 5, otherwise the analyses are
#             # considedered as already performed -> continue to the next project
#             create_matrix_new.create_summary_table(OUTPUT_PATH, PROJECT,
#                                                    DRUGS_title, drugs)
#             # creates the complete summary table with every
#             # drugcombination in it
#             # pay attention to that the files have all the same dim!
#             # DONE!!!
#             # with bash:
#             # for i in *; do
#             # filesize=$( wc -l <"$i"); if [ "$filesize" -lt 450000 ];
#             # then mv $i ../bad_data/ ; fi; done
#             # temporary solution, find a better one!

#             # it is supposed, that if the table_for_metilene is created, the
#             # complete analyse part for single projects is performed, s.t. if
#             # this file exists, the single project analyse of this project can
#             # be skipped: 'summary_dead_alive.tsv' in  OUTPUT_PATH, PROJECT,
#             # DRUGS_title
#             # TODO, uncomment this statement when finished:
#         # if os.path.exists(os.path.join(
#         #         OUTPUT_PATH, PROJECT, DRUGS_title,
#         #         'summary_dead_alive.tsv')):
#         #     continue

#         # a new table is created with add col drugnames (drugcombinations)
#         # DF_3t_both_with_DRUG_combi.tsv
#         # aggregation of same case_id's and their drugcombinations, make
#         # those lists as sets in case after correcting the drugnames
#         # duplicates arise order the entries, drug combis in field

#         # drugnames are comma seperated filter case_id of interest and
#         # exchange case_id in summary with vital status, call DESeq2

#         # the last file created in the -A part (single project)is the
#         # DESeq2_out_DRUG_combi*/results_statistics.tsv file in
#         # OUTPUT/project/DRUGS_title, if this file exists in the available
#         # DESeq2_out_DRUG_combi* dirs, the deseq runs are already completed for
#         # the single projects, continue to the next single project:

#         # continue_flag = False
#         # for dir_iter in glob.glob(os.path.join(
#         #         OUTPUT_PATH, PROJECT, DRUGS_title,
#         #         'DESeq2_out_DRUG_combi*')):
#         #     if os.path.exists(os.path.join(dir_iter,
#         #                                    'results_statistics.tsv')):
#         #         continue_flag = True
#         #     else:
#         #         continue_flag = False
#         #         # if one 'results_statistic.tsv' is missing, the deseq
#         #         analysis
#         #         # has to be started all over again, we can break here
#         #         break
#         # # TODO uncomment this, when finished cutoff stuff
#         # if continue_flag:
#         #     continue

#         if 7 in function or 100 in function:
#             create_matrix_new.provide_DESeq2_table(PROJECT, OUTPUT_PATH,
#                                                    DRUGS, SCRIPT_PATH, logger,
#                                                    cutoff, DRUGS_title)
#         if 8 in function or 100 in function:
#             create_matrix_new.create_statistics_from_DESeq2_tables(
#                 OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT, logger, DRUGS_title)

#     # if min(function) >= 11 and max(function) <= 15:

#     # dont call that fct 9 and 10 if just one project is applied!
#     os.makedirs(os.path.join(
#         OUTPUT_PATH, PROJECT_title, DRUGS_title), exist_ok=True)
#     logger = create_matrix_new.set_logger(
#         OUTPUT_PATH, PROJECT_title, DRUGS_title)
#     logger.info('Version:\t{}'.format(version))
#     if 9 in function or 100 in function:
#         if len(PROJECT_DRUG_UUID) == 1:
#             pass
#             # os._exit(0)
#         else:
#             # aggregate the summary tables of every PROJECT for a single
#             # multifactorial DESeq2 run:
#             create_matrix_new.provide_DESeq2_table(
#                 PROJECT_DRUG_UUID, OUTPUT_PATH, DRUGS, SCRIPT_PATH, logger,
#                 cutoff, DRUGS_title)

#     if 10 in function or 100 in function:
#         if len(PROJECT_DRUG_UUID) == 1:
#             pass
#         else:
#             create_matrix_new.create_statistics_from_DESeq2_tables(
#                 OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID, logger,
#                 DRUGS_title)

#     if 11 in function or 100 in function:
#         for threshold_iter in threshold:
#             lifeline_summary_test_2.lifelines_ENSG(
#                 OUTPUT_PATH, PROJECT_DRUG_UUID, threshold_iter, DRUGS_title)
#         # in case new thresholds are added after the complete run with -A,
#         # those
#         # outputs have to be added right away to the snake config:
#         create_matrix_new.create_snake_config(
#             OUTPUT_PATH, PROJECT_title, DRUGS_title, PROJECT_list, DRUGS,
#             SCRIPT_PATH, cutoff)

#     if 12 in function or 100 in function:
#         walk_all_drug_frequency.drug_frequency(
#             OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID, DRUGS_title)
#     if 13 in function or 100 in function:
#         create_report.create_report_pdf(
#             OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID, threshold,
#             cutoff, DRUGS_title)
#         # if multi projects are applied, create report in each single proj dir:
#         if len(PROJECT_DRUG_UUID) > 1:
#             for PROJECT in PROJECT_DRUG_UUID:
#                 create_report.create_report_pdf(
#                     OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT, threshold,
#                     cutoff, DRUGS_title)
#     if 14 in function or 100 in function:
#         walk_all_drug_frequency.drug_frequency_all_single_projects(
#             OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID, DRUGS_title)
#     if 15 in function or 100 in function:
#         create_matrix_new.create_snake_config(
#             OUTPUT_PATH, PROJECT_title, DRUGS_title, PROJECT_list, DRUGS,
#             SCRIPT_PATH, cutoff)


# if __name__ == '__main__':
#     call_with_options()
