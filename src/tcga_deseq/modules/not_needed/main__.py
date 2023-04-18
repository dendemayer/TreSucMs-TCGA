import os
import click
from shared.modules import choose_therapy
from shared.modules import download_with_api


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
@click.option('--version', '-v',
              help='printing out version information: {}'.format(version),
              is_flag=True, callback=print_version,
              expose_value=False, is_eager=True)
def call_with_options(out_path, project, drugs, cores):
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
    if len(project) == 0:
        PROJECT = choose_therapy.Choose_project()
        PROJECT = sorted(map(str.upper, PROJECT))
    else:
        PROJECT = sorted(map(str.upper, project))
    if len(drugs) == 0:
        DRUGS = choose_therapy.Choose_drugs(SCRIPT_PATH, PROJECT)
    else:
        DRUGS = sorted(map(str.lower, drugs))

    print('PROJECT:\t', PROJECT)
    print('DRUGS:\t\t', DRUGS)

    # Snakefile_path = os.path.join(SCRIPT_PATH, os.path.pardir, 'Snakefile')
    Snakefile_path = os.path.join(os.path.split(SCRIPT_PATH)[0], 'Snakefile')
    workdir = os.path.split(Snakefile_path)[0]
    # downloading the manifest file and the gtf file
    mani = 'gdc_manifest_20211029_data_release_31.0_active.tsv.gz'
    mani_path = os.path.join(OUTPUT_PATH, 'metadata', mani)
    # download_with_api.download_GDC_manifest(Snakefile_path, mani_path,
    # workdir,
    #                                         cores)
    gtf_path = os.path.join(OUTPUT_PATH, 'metadata',
                            'gencode.v36.annotation.gtf.gz')
    # important to hand over the workdir, the CWD can be anywhere, the workdir
    # is derived from the scriptpath
    help_file_list = [gtf_path, mani_path]
    download_with_api.download_GDC_manifest(
        Snakefile_path, help_file_list, workdir, cores)

    # if no -A and -D specified, ask here whether to just download raw data or
    # do both download and analyse steps:
    # if download_data and analyse_data:
    #     function = (100,)
    # if download_data and not analyse_data:
    #     function = (1, 2, 3, 4)  # until creating table:
    #     # meta_info_druglist_merged_drugs_combined.tsv and download all data
    # if analyse_data and not download_data:  # all the rest
    #     function = (5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)
    # # if the first 5 fct are skipped (because those projects are
    # # downloaded already), the info in the snakemake_config.yaml would miss,
    # # because all those log entries are skipped also. if just -A is set,
    # # we must look up all those files and write it in the config

    # whats set so far are the projects and drugs, still needed analyse options
    # and OUTPUT_PATH
    # if not analyse_data and not download_data and not function:
    #     OUTPUT_PATH, function, analyse_data,\
    #         download_data = choose_therapy.Choose_path_and_option(
    #             OUTPUT_PATH, PROJECT, DRUGS, function, SCRIPT_PATH,
    #             analyse_data, download_data)

    # OUTPUT_PATH = OUTPUT_PATH.rstrip(os.sep)
