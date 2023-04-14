
# import os
# import click
# from shared.modules import choose_therapy
# from shared.modules import download_with_api
# import snakemake

# SCRIPT_PATH = os.path.split(__file__)[0]
# with open(os.path.join(SCRIPT_PATH, 'version.txt'), 'r') as f:
#     version = f.readline().strip()


# def print_version(ctx, param, value):
#     if not value or ctx.resilient_parsing:
#         return
#     click.echo(version)
#     ctx.exit()


# HOME = os.getenv('HOME')


# @click.command()
# @click.option('--out_path', '-o', default=os.path.join(HOME, 'TCGA-pipelines'),
#               show_default=True,
#               help='path to save the result files')
# @click.option('--project', '-p', default=[], multiple=True,
#               help='TCGA project to be applied. Any TCGA project can be' +
#               ' chosen, like: ' +
#               '-p TCGA-CESC -p TCGA-HNSC ...')
# @click.option('--drugs', '-d', default=[], multiple=True, show_default=False,
#               help='drug(s), like: -d drug1 -d drug2 ' +
#               'or drugcombination(s), like: -d drug1,drug2')
# @click.option('--cores', '-c', default=1, multiple=False, show_default=True,
#               type=int, help='number of cores provided to snakemake',
#               required=False)
# @click.option('--version', '-v',
#               help='printing out version information: {}'.format(version),
#               is_flag=True, callback=print_version,
#               expose_value=False, is_eager=True)
# def call_with_options(out_path, project, drugs, cores):
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

#     shared_workdir = os.path.join(
#         os.path.split(os.path.split(SCRIPT_PATH)[0])[0], 'shared')
#     Snakemake_all_files = []

#     Snakefile_path = os.path.join(os.path.split(SCRIPT_PATH)[0], 'Snakefile')
#     config_file_shared = os.path.join(shared_workdir, 'config.yaml')

#     # help files for both pipelines, like:
#     # OUTPUT_PATH/metadata/gdc_manifest_20211029_data_release_31...,
#     # gencode.v36.annotation
#     help_file_list = download_with_api.download_GDC_manifest(OUTPUT_PATH,
#                                                              config_file_shared)
#     Snakemake_all_files = Snakemake_all_files + help_file_list

#     # auxfiles for both pipelines:
#     # OUTPUT_PATH/PROJECT/aux_files/nationwidechildrens.....
#     aux_file_list = download_with_api.download_aux_files(OUTPUT_PATH, PROJECT,
#                                                          config_file_shared)

#     Snakemake_all_files = Snakemake_all_files + aux_file_list

#     # Datafiles: OUTPUT_PATH/PROJECT/Diffexpression/PROJECT_data_files/...
#     data_file_list = download_with_api.download_data_files(OUTPUT_PATH, PROJECT,
#                                                           config_file_shared,
#                                                           'HumanMethylation450')

#     Snakemake_all_files = Snakemake_all_files + data_file_list
#     print('running snakemake with\n')
#     print(f'Snakefile_path:\t{Snakefile_path}')
#     print(f'shared_workdir:\t{shared_workdir}')
#     print(f'Snakemake_all_files:\t{Snakemake_all_files}')
#     print(f'cores:\t, {cores}')

#     snakemake.snakemake(snakefile=Snakefile_path, targets=Snakemake_all_files,
#                         workdir=shared_workdir, cores=cores, forceall=False,
#                         force_incomplete=True, dryrun=True)