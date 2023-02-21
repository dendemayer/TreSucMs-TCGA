import os
import click
from shared.modules import choose_therapy
from shared.modules import download_with_api
from tcga_metilene.modules import main_metilene
from tcga_deseq.modules import main_deseq
import snakemake

SCRIPT_PATH = os.path.split(__file__)[0]
with open(os.path.join(SCRIPT_PATH, 'version.txt'), 'r') as f:
    version = f.readline().strip()


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(version)
    ctx.exit()


HOME = os.getenv('HOME')


@click.command()
@click.option('--out_path', '-o', default=os.path.join(HOME, 'TCGA-pipelines'),
              show_default=True,
              help='path to save the result files')
@click.option('--project', '-p', default=[], multiple=True,
              help='TCGA project to be applied. Any TCGA project can be' +
              ' chosen, like: ' +
              '-p TCGA-CESC -p TCGA-HNSC ...')
@click.option('--drugs', '-d', default=[], multiple=True, show_default=False,
              help='drug(s), like: -d drug1 -d drug2 ' +
              'or drugcombination(s), like: -d drug1,drug2')
@click.option('--cores', '-c', default=1, multiple=False, show_default=True,
              type=int, help='number of cores provided to snakemake',
              required=False)
@click.option('--cutoff', '-C', default=0, multiple=False, show_default=True,
              type=int, help='Cut-off parameter',
              required=False)
@click.option('--threshold', '-t', default=0, multiple=False, show_default=True,
              type=int, help='threshold parameter',
              required=False)
@click.option('--execute', '-e', default=['Deseq2', 'metilene'], multiple=True,
              show_default=True, help='choose which pipeline shall be\
              executed')
@click.option('--version', '-v',
              help='printing out version information: {}'.format(version),
              is_flag=True, callback=print_version,
              expose_value=False, is_eager=True)
def call_with_options(out_path, project, drugs, cores, execute, cutoff,
                      threshold):
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
    print("\nOUTPUT_PATH:\t\t", OUTPUT_PATH)
    # SCRIPT_PATH = script_path
    SCRIPT_PATH = os.path.split(__file__)[0]
    print("SCRIPT_PATH:\t\t", SCRIPT_PATH)
    # check if the main_metilene.py exist in the SCRIPT_PATH:
    print("PIPELINES executed:\t", execute)
    if len(project) == 0:
        PROJECT = choose_therapy.Choose_project()
        PROJECT = sorted(map(str.upper, PROJECT))
    else:
        PROJECT = sorted(map(str.upper, project))
    if len(drugs) == 0:
        DRUGS = choose_therapy.Choose_drugs(SCRIPT_PATH, PROJECT)
        DRUGS = sorted(map(str.upper, DRUGS))
    else:
        DRUGS = sorted(map(str.lower, drugs))

    print('PROJECT:\t\t', PROJECT)
    print('DRUGS:\t\t\t', DRUGS)
    print(f'cores:\t\t\t{cores}')
    print(f'cutoff:\t\t\t{cutoff}')
    print(f'threshold:\t\t{threshold}')

    shared_workdir = os.path.join(
        os.path.split(os.path.split(SCRIPT_PATH)[0])[0], 'shared')
    Snakemake_all_files = []

    Snakefile = os.path.join(os.path.split(SCRIPT_PATH)[0], 'Snakefile')
    config_file_shared = os.path.join(shared_workdir, 'config.yaml')

    # help files for both pipelines, like:
    # OUTPUT_PATH/metadata/gdc_manifest_20211029_data_release_31...,
    # gencode.v36.annotation
    help_file_list = download_with_api.download_GDC_manifest(
        OUTPUT_PATH, config_file_shared)
    Snakemake_all_files = Snakemake_all_files + help_file_list
    # once we have to call snakemake in prior, s.t. the manifest file is
    # present on which all the following selections are done on, make sure that
    # here the dryrun flag is not set to False
    snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
                        workdir=shared_workdir, cores=cores, forceall=False,
                        force_incomplete=True, dryrun=False)

    # auxfiles for both pipelines:
    # OUTPUT_PATH/PROJECT/aux_files/nationwidechildrens.....
    aux_file_list = download_with_api.download_aux_files(OUTPUT_PATH, PROJECT,
                                                         config_file_shared)

    Snakemake_all_files = Snakemake_all_files + aux_file_list

    # translate here the applied pipeline which shall be executet:
    # Datafiles: OUTPUT_PATH/PROJECT/Diffexpression/PROJECT_data_files/...
    def map_execute(pipeline):
        if pipeline == 'Deseq2':
            return 'htseq'
        elif pipeline == 'metilene':
            return 'HumanMethylation450'

    data_file_list = []
    file_types = map(map_execute, execute)
    for file_type in file_types:
        data_file_list = (data_file_list +
                          download_with_api.download_data_files(
                              OUTPUT_PATH, PROJECT, config_file_shared,
                              file_type))

    Snakemake_all_files = Snakemake_all_files + data_file_list

    merged_tables_list = []
    # applied to the snakemake rule merge_meta_tables:, it also creates the
    # meta_info_druglist_merged_drugs_combined.tsv tables
    for proj in PROJECT:
        for pipeline in execute:
            merged_tables_list.append(os.path.join(
                OUTPUT_PATH, proj, pipeline,'merged_meta_files',
                'merged_meta_tables.tsv'))

    # merged_tables_list = [os.path.join(
    #     OUTPUT_PATH, x,
    #     'merged_meta_files', 'merged_meta_tables.tsv') for x in PROJECT]
    # add the merged and processed metatables: specific on pipeline:
    Snakemake_all_files = Snakemake_all_files + merged_tables_list

    print('running snakemake with\n')
    print(f'Snakefile:\t{Snakefile}')
    print(f'shared_workdir:\t{shared_workdir}')
    # print(f'Snakemake_all_files:\t{Snakemake_all_files}')

    # also add the multi proj meta_info_druglist_merged_drugs_combined.tsv
    # which is just the concatenation of the single proj pendants:
    projects = '_'.join(PROJECT)
    merged_drugs_combined_list = []
    for pipeline in execute:
        merged_drugs_combined_list.append(os.path.join(
            OUTPUT_PATH, projects, pipeline,'merged_meta_files',
            'meta_info_druglist_merged_drugs_combined.tsv'))

    # Snakemake_all_files = Snakemake_all_files + merged_drugs_combined_list

    snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
                        workdir=shared_workdir, cores=cores, forceall=False,
                        force_incomplete=True, dryrun=False)
    #########################################################################
    #########################################################################
    # from here the shared modules and Snakemake scripts are getting pipeline
    # specific, hand over all outputfiles created so far and enter the pipeline
    # specific main files:
    # main_metilene.entry_fct(OUTPUT_PATH, PROJECT, DRUGS, Snakemake_all_files,
    #                         cutoff, threshold, cores)


