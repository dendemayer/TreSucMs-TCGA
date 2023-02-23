import snakemake
import os
from tcga_metilene.modules import create_summary_table

def entry_fct(OUTPUT_PATH, PROJECT, DRUGS, Snakemake_all_files, cutoffs, threshold, cores):
    SCRIPT_PATH = os.path.split(__file__)[0]
    metilene_snake_workdir = os.path.join(
        os.path.split(os.path.split(SCRIPT_PATH)[0])[0], 'tcga_metilene')

    shared_workdir = os.path.join(
        os.path.split(os.path.split(SCRIPT_PATH)[0])[0], 'shared')
    Snakefile = os.path.join(os.path.split(SCRIPT_PATH)[0], 'Snakefile')
    config_file = os.path.join(os.path.split(SCRIPT_PATH)[0], 'config.yaml')

    summary_tables = create_summary_table.return_summary_tables(
        OUTPUT_PATH, PROJECT, DRUGS, cutoffs)

    Snakemake_all_files = Snakemake_all_files + summary_tables
    snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files, workdir=shared_workdir, cores=cores, forceall=False, force_incomplete=True, dryrun=False, debug=False, use_conda=True)
# # ##### main_metilene ############

# from tcga_metilene.modules import methyl
# from shared import download_with_api
# import os
# import click
# from shared import choose_therapy
# # import choose_therapy
# # import create_report
# # import lifeline_plots
# # import set_logger

# '''
#     within this, every parameter needed for the analysis is set

#     :param: out_path: path for metilene_pipeline outputs
#     :type: out_path: str
#     :param: script_path: path to the metilene_pipeline repo
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
#     :param: analyse_data: bool flag, whether the metilene analyses shall be\
#         started
# '''

# # version = 'Version 1.0'
# # to specify the SCRIPT_PATH with __file__ python 3.9 is mandatory, else just
# # the relative path would be returned (main_metilene.py) and not the absolute
# # path tho the file (here
# # /homes/biertruck/gabor/phd/test_git_doc/methyl/main_metilene.py)
# # the absolute path is needed to extract the path to all scripts!
# SCRIPT_PATH = os.path.split(__file__)[0]
# with open(os.path.join(SCRIPT_PATH, 'version.txt'), 'r') as f:
#     version = f.readline().strip()


# def print_version(ctx, param, value):
#     if not value or ctx.resilient_parsing:
#         return
#     click.echo(version)
#     ctx.exit()


# CWD = os.environ.get('PWD')


# @click.command()
# @click.option('--out_path', '-o', default=os.path.join(CWD, 'metilene_data'),
#               show_default=True,
#               help='path to save the result files')
# @click.option('--project', '-p', default=[], multiple=True,
#               help='TCGA project to be applied. Any TCGA project can be' +
#               ' chosen, like: ' +
#               '-p TCGA-CESC -p TCGA-HNSC ...')
# @click.option('--drugs', '-d', default=[], multiple=True, show_default=False,
#               help='drug(s), like: -d drug1 -d drug2 ' +
#               'or drugcombination(s), like: -d drug1,drug2')
# @click.option('--version', '-v',
#               help='printing out version information: {}'.format(version),
#               is_flag=True, callback=print_version,
#               expose_value=False, is_eager=True)
# def call_with_options(out_path, project, drugs):
#     '''
#     "metilene_pipeline" a tool to choose, harvest and analyse methylation data
#     of the TCGA-projects with help of the metilene package.\n

#     Build and activate the provided conda env:

#         $ conda env create -f metilene_env.yaml

#         $ conda activate metilene_pipeline

#     call the script without any options to enter the interactive mode and set
#     each option step by step:

#         $ python main_metilene.py

#     print help page:

#         $ python main_metilene.py --help
#     '''
#     OUTPUT_PATH = out_path
#     print("\nOUTPUT_PATH:\n", OUTPUT_PATH)
#     # SCRIPT_PATH = script_path
#     SCRIPT_PATH = os.path.split(__file__)[0]
#     print("\nSCRIPT_PATH:\n", SCRIPT_PATH)
#     # check if the main_metilene.py exist in the SCRIPT_PATH:
#     if not os.path.isfile(os.path.join(SCRIPT_PATH, 'main_metilene.py')):
#         print('there is no "main_metilene.py" in your SCRIPT_PATH ' +
#               'please check your input, exiting now')
#         os._exit(0)
#     # print("\nSCRIPT_PATH:\n", SCRIPT_PATH)
#     # the projects shall be applied in list form, afterwards the correct hash
#     # pairs are handed over, the UUIDs are present here
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

#     print('PROJECT:\t', PROJECT)
#     print('DRUGS:\t\t', DRUGS)

#     # if no -A and -D specified, ask here whether to just download raw data or
#     # do both download and analyse steps:
#     if download_data and analyse_data:
#         function = (100,)
#     if download_data and not analyse_data:
#         function = (1, 2, 3, 4)  # until creating table:
#         # meta_info_druglist_merged_drugs_combined.tsv and download all data
#     if analyse_data and not download_data:  # all the rest
#         function = (5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)
#     # if the first 5 fct are skipped (because those projects are
#     # downloaded already), the info in the snakemake_config.yaml would miss,
#     # because all those log entries are skipped also. if just -A is set,
#     # we must look up all those files and write it in the config

#     # whats set so far are the projects and drugs, still needed analyse options
#     # and OUTPUT_PATH
#     if not analyse_data and not download_data and not function:
#         OUTPUT_PATH, function, analyse_data,\
#             download_data = choose_therapy.Choose_path_and_option(
#                 OUTPUT_PATH, PROJECT, DRUGS, function, SCRIPT_PATH,
#                 analyse_data, download_data)

#     OUTPUT_PATH = OUTPUT_PATH.rstrip(os.sep)

#     ####################################
#     # met_opt = {'-m': '3', '-M': '1000', '-d': '0.03'}
#     # met_opt = {'-m': '3', '-M': '1000'}
#     # met_opt_list = []
#     # for i in met_opt:
#     #     met_opt_list.append(i)
#     #     met_opt_list.append(met_opt[i])
#     #############################################
#     # met_dir = '_'.join(met_opt_list)
#     # met_dir = 'metilene_' + met_dir
#     # print("\nmet_dir:\n", met_dir)

#     # every fct is a selfcontained step, to apply which fct is needed in
#     # snakemake, enum the fct and apply the nr. to the options:
#     # if min(function) >= 1 and max(function) <= 10 or 100 in function:
#     DRUGS_title = '_'.join(sorted(map(str.lower, DRUGS)))
#     if cutoff > 0:
#         DRUGS_title = DRUGS_title + '_' + str(cutoff)
#     projects = []
#     for project in PROJECT_DRUG_UUID:
#         projects.append(project)
#     PROJECT_title = '_'.join(sorted(map(str.upper, projects)))
#     PROJECT_list = projects
#     api_manifest = 'gdc_manifest_20211029_data_release_31.0_active.tsv.gz'
#     gtf_UUID = '25aa497c-e615-4cb7-8751-71f744f9691f'
#     workflow_type = 'Liftover'
#     # in case that just -A is set, go through the single project dirs, check
#     # the existence of all files needed, and
#     # create the logfile in the PROJECTS_title/DRUGS_title/ dir
#     if analyse_data and not download_data:
#         methyl.create_log_for_A(OUTPUT_PATH, PROJECT_list, PROJECT_title,
#                                 DRUGS_title, api_manifest)
#     for PROJECT in PROJECT_DRUG_UUID:
#         # log file is written to OUT/singl_proj/DRUGS_title
#         os.makedirs(os.path.join(
#             OUTPUT_PATH, PROJECT, DRUGS_title), exist_ok=True)
#         logger = set_logger.set_logger(
#             OUTPUT_PATH, PROJECT, DRUGS_title)
#         logger.info('Version:\t{}'.format(version))
#         if 1 in function or 100 in function:
#             download_with_api.download_GDC_manifest(PROJECT, FILE_TYPE,
#                                                     OUTPUT_PATH, DRUGS_title,
#                                                     SCRIPT_PATH, PROJECT_title,
#                                                     snakerun, api_manifest,
#                                                     gtf_UUID)

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
#             methyl.create_summary_table(OUTPUT_PATH, PROJECT, DRUGS_title)
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
#         if 6 in function or 100 in function:
#             methyl.create_table_for_metilene(
#                 OUTPUT_PATH,
#                 PROJECT, DRUGS, DRUGS_title, met_dir, cutoff)
#             # creates the filtered table with just specific drug
#             # combinations in it

#         if 7 in function or 100 in function:
#             methyl.provide_metilene_table(
#                 OUTPUT_PATH,
#                 DRUGS, SCRIPT_PATH, met_opt_list, met_dir, PROJECT,
#                 DRUGS_title)
#             # call metilene and the sppl scripts on the before created
#             # tables

#         if 8 in function or 100 in function:
#             methyl.call_bed_intersect(
#                 OUTPUT_PATH, DRUGS, met_opt_list, met_dir, PROJECT,
#                 DRUGS_title)
#             # call bed intersect for the single projects

#         if 9 in function or 100 in function:
#             methyl.create_plots(
#                 OUTPUT_PATH, DRUGS, met_opt_list, met_dir, PROJECT,
#                 DRUGS_title)

#         # boxplot for ranges
#         if 10 in function or 100 in function:
#             methyl.create_mean_boxplot(
#                 OUTPUT_PATH, DRUGS, met_opt_list, met_dir, PROJECT,
#                 DRUGS_title)

#     # don't call multi project fcts if just one project is applied!, but create
#     # the  snake config and the report file
#     # the report must be created before the config, otherwise it would miss
#     if len(PROJECT_DRUG_UUID) == 1:
#         if 17 in function or 100 in function:
#             # lifeline_plots.prepare_lifeline_plots(OUTPUT_PATH, DRUGS,
#                                                   # SCRIPT_PATH,
#                                                   # PROJECT_DRUG_UUID, met_dir,
#                                                   # cutoff, DRUGS_title)
#             # lifeline_plots.prepare_complement_lifeline_plots(OUTPUT_PATH,
#                                                              # DRUGS,
#                                                              # SCRIPT_PATH,
#                                                              # PROJECT_DRUG_UUID,
#                                                              # met_dir, cutoff,
#                                                              # DRUGS_title)
#             lifeline_plots.create_lifeline_plots(OUTPUT_PATH, DRUGS,
#                                                  SCRIPT_PATH,
#                                                  PROJECT_DRUG_UUID, met_dir,
#                                                  cutoff, DRUGS_title)
#             create_report.create_report_pdf(OUTPUT_PATH, DRUGS, SCRIPT_PATH,
#                                             PROJECT_DRUG_UUID, met_dir,
#                                             cutoff, DRUGS_title)
#         if 16 in function or 100 in function:
#             methyl.create_snake_config(OUTPUT_PATH, PROJECT_title,
#                                        DRUGS_title, PROJECT_list, DRUGS,
#                                        SCRIPT_PATH, cutoff)
#             os._exit(0)

#     # check here, if the cutoff is considered automatically TODO
#     if 11 in function or 100 in function:
#         methyl.create_table_all_projects(
#             OUTPUT_PATH, DRUGS, PROJECT_DRUG_UUID, DRUGS_title, met_dir,
#             cutoff)

#     # call metilene for every project and for the aggregated summary:
#     if 12 in function or 100 in function:
#         methyl.provide_metilene_table(
#             OUTPUT_PATH, DRUGS, SCRIPT_PATH, met_opt_list, met_dir,
#             PROJECT_DRUG_UUID, DRUGS_title)

#     # call bed intersect for every project
#     if 13 in function or 100 in function:
#         methyl.call_bed_intersect(
#             OUTPUT_PATH, DRUGS, met_opt_list, met_dir, PROJECT_DRUG_UUID,
#             DRUGS_title)

#     if 14 in function or 100 in function:
#         methyl.create_plots(OUTPUT_PATH, DRUGS, met_opt_list, met_dir,
#                             PROJECT_DRUG_UUID, DRUGS_title)

#     if 15 in function or 100 in function:
#         methyl.create_mean_boxplot(
#             OUTPUT_PATH, DRUGS, met_opt_list, met_dir, PROJECT_DRUG_UUID,
#             DRUGS_title)

#     # the report must be created before the config, otherwise it would missed
#     # in the ovarall config
#     if len(PROJECT_DRUG_UUID) > 1:
#         if 17 in function or 100 in function:
#             print(f'creating report for {PROJECT_DRUG_UUID}')
#             # make the lifeline plots for the multiproject:
#             lifeline_plots.prepare_lifeline_plots(OUTPUT_PATH, DRUGS,
#                                                   SCRIPT_PATH,
#                                                   PROJECT_DRUG_UUID, met_dir,
#                                                   cutoff, DRUGS_title)
#             lifeline_plots.prepare_complement_lifeline_plots(OUTPUT_PATH,
#                                                              DRUGS,
#                                                              SCRIPT_PATH,
#                                                              PROJECT_DRUG_UUID,
#                                                              met_dir, cutoff,
#                                                              DRUGS_title)
#             lifeline_plots.create_lifeline_plots(OUTPUT_PATH, DRUGS,
#                                                  SCRIPT_PATH,
#                                                  PROJECT_DRUG_UUID, met_dir,
#                                                  cutoff, DRUGS_title)
#             create_report.create_report_pdf(OUTPUT_PATH, DRUGS, SCRIPT_PATH,
#                                             PROJECT_DRUG_UUID, met_dir, cutoff,
#                                             DRUGS_title)
#             # if multi projects are applied, create report in each single proj
#             # dir:
#             for PROJECT in PROJECT_DRUG_UUID:
#                 print(f'creating report for {PROJECT}')
#                 lifeline_plots.prepare_lifeline_plots(OUTPUT_PATH, DRUGS,
#                                                     SCRIPT_PATH,
#                                                     PROJECT,
#                                                     met_dir, cutoff,
#                                                     DRUGS_title)
#                 lifeline_plots.prepare_complement_lifeline_plots(OUTPUT_PATH,
#                                                                 DRUGS,
#                                                                 SCRIPT_PATH,
#                                                                 PROJECT,
#                                                                 met_dir,
#                                                                 cutoff,
#                                                                 DRUGS_title)
#                 lifeline_plots.create_lifeline_plots(OUTPUT_PATH, DRUGS,
#                                                     SCRIPT_PATH,
#                                                     PROJECT,
#                                                     met_dir, cutoff,
#                                                     DRUGS_title)
#                 create_report.create_report_pdf(OUTPUT_PATH, DRUGS,
#                                                 SCRIPT_PATH, PROJECT,
#                                                 met_dir, cutoff, DRUGS_title)

#         if 16 in function or 100 in function:
#             methyl.create_snake_config(OUTPUT_PATH, PROJECT_title, DRUGS_title,
#                                     PROJECT_list, DRUGS, SCRIPT_PATH, cutoff)
#         # ##### adding a sanitiy check for good canditates:
#         #  methyl.sanity_check(OUTPUT_PATH, DRUGS)

#         # ### check the single runs for each project of metilene and put them
#         # together for intersect and visualisation
#         # methyl.sanity_check_2(OUTPUT_PATH, DRUGS, PROJECT_DRUG_UUID)


# if __name__ == '__main__':
#     call_with_options()
