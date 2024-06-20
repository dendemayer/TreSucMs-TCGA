import os
import click
from shared.modules import choose_therapy
from shared.modules import download_with_api
from tcga_metilene.modules import main_metilene
from tcga_deseq.modules import main_deseq
import snakemake
from itertools import compress
# import re

SCRIPT_PATH = os.path.split(__file__)[0]
with open(os.path.join(SCRIPT_PATH, 'version.py'), 'r') as f:
    version = f.readline().strip()

pipeline_list = ['DESeq2', 'metilene']


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(version)
    ctx.exit()


HOME = os.getenv('HOME')


@click.command()
@click.option('--out_path', '-o', default=os.path.join(HOME, 'TREMSUCS'),
              show_default=True,
              help='path to save the result files')
@click.option('--project', '-p', default=[], multiple=True,
              help='TCGA project(s) to be applied. Any TCGA project can be' +
              ' chosen, like: ' +
              '-p TCGA-CESC -p TCGA-HNSC ...')
@click.option('--drugs', '-d', default=[], multiple=True, show_default=False,
              help='drug(s), like: -d drug1 -d drug2 ' +
              'or drugcombination(s), like: -d drug1,drug2')
@click.option('--cores', '-c', default=1, multiple=False, show_default=True,
              type=int, help='number of cores provided to snakemake',
              required=False)
@click.option('--cutoff', '-C', default=[0], multiple=True, show_default=True,
              type=float, help='Cut-off parameter',
              required=False)
@click.option('--threshold', '-t', default=[0], multiple=True, show_default=True,
              type=float, help='threshold parameter',
              required=False)
@click.option('--execute', '-e', default=pipeline_list, multiple=True,
              show_default=True, help='choose which pipeline shall be executed')
@click.option('--dryrun', '-N', default=False, multiple=False,
              show_default=True, is_flag=True, help='snakemake dryrun',
              required=False)
@click.option('--download', '-D', default=False, multiple=False,
              show_default=True, is_flag=True, help='''if set, just download raw
              and meta data for given projects and analysis types, revise them,
              link them, but do not run any analysis''',
              required=False)
@click.option('--unlock', '-u', default=False, multiple=False,
              show_default=True, is_flag=True, help='''in case the analysis
              crashs, snakemake locks the output directory, run with -u to
              unlock, then repeat the analysis''',
              required=False)
# this one is hidden, debug purposes
@click.option('--report', '-r', default=False, multiple=False,
              show_default=True, is_flag=True, help='just create a report',
              required=False, hidden=True)
@click.option('--version', '-v',
              help='printing out version information: {}'.format(version),
              is_flag=True, callback=print_version,
              expose_value=False, is_eager=True)
def call_with_options(out_path, project, drugs, cores, execute, cutoff,
                      threshold, dryrun, report, download, unlock):
    '''
    "TREMSUCS" a tool to choose, harvest and analyse expression and methylation data
    of the TCGA-projects for revealing Biomarkers which indicate threapy
    specific treatment success predictions.

    Calling the pipeline without any argument starts the interactive mode to
    help setting all needed parameters for the analysis.
    '''
    OUTPUT_PATH = out_path
    # print("\nOUTPUT_PATH:\t\t", OUTPUT_PATH)
    # SCRIPT_PATH = script_path
    SCRIPT_PATH = os.path.split(__file__)[0]
    # print("SCRIPT_PATH:\t\t", SCRIPT_PATH)
    # make sure that the pipelines to execute also exist, every entry must be
    # present in : ['DESeq2', 'metilene']
    temp_check = [True if i not in pipeline_list else False for i in execute]
    if True in temp_check:
        print('\nyou misspelled a pipeline name, make sure the ', end='')
        print(f'-e option set is within the set of {pipeline_list}, ', end='')
        print('wrong pipeline name: ', end='')
        print(f'{list(compress(execute, temp_check))}, ', end='')
        print('exiting now')
        os._exit(0)
    execute = sorted(list(execute))
    project = [i.strip() for i in project]
    interacitve_proj = False
    interacitve_drugs = False
    if len(project) == 0:
        PROJECT = choose_therapy.Choose_project()
        interacitve_proj = True
    else:
        PROJECT = sorted(map(str.upper, project))
    if len(drugs) == 0:
        DRUGS = choose_therapy.Choose_drugs(PROJECT)
        interacitve_drugs = True
    else:
        DRUGS = sorted(map(str.lower, drugs))
    # in case PROJECT and DRUGS are set via interactive mode, also ask every
    # other parameter and if the default value is ok:
    # set all parameters in interactive mode, also ask for Cutoff, threshold,
    # cores, output path, pipelines to execute
    if interacitve_proj and interacitve_drugs:
        # just ask those parameters wich deviate from the default value:
        if OUTPUT_PATH == os.path.join(HOME, 'TREMSUCS'):
            OUTPUT_PATH = choose_therapy.update_parameters(OUTPUT_PATH, 'OUTPUT_PATH')
        if cores == 1:
            cores = choose_therapy.update_parameters(cores, 'cores')
        if execute == ['DESeq2', 'metilene']:
            execute = choose_therapy.update_parameters(execute, 'pipelines')
        if cutoff == (0.0,):
            cutoff = choose_therapy.update_parameters(cutoff, 'cutoff')
        if threshold == (0.0,):
            threshold = choose_therapy.update_parameters(threshold, 'threshold')

    cutoffs = list(cutoff)
    for index, cutoff in enumerate(cutoffs):
        if cutoff % 1 == 0:
            cutoffs[index] = round(cutoff)
    if 0 not in cutoffs:
        cutoffs.append(0)
    cutoffs = sorted(list(set(cutoffs)))

    thresholds = list(threshold)
    for index, threshold in enumerate(thresholds):
        if threshold % 1 == 0:
            thresholds[index] = round(threshold)
    if 0 not in thresholds:
        thresholds.append(0)
    threshold = sorted(list(set(thresholds)))

    thresh_list = [f'threshold_{str(i)}' for i in threshold]
    thresh_str = '_'.join(thresh_list)
    print("SCRIPT_PATH:\t\t", SCRIPT_PATH)
    print('OUTPUT_PATH:\t\t', OUTPUT_PATH)
    print('PROJECT:\t\t', PROJECT)
    print('DRUGS:\t\t\t', DRUGS)
    print('pipelines executed:\t', execute)
    print(f'cores:\t\t\t{cores}')
    print(f'cutoff:\t\t\t{cutoffs}')
    print(f'threshold:\t\t{threshold}')
    print('the command to issue this analysis would be:\n')
    print(f'TREMSUCS {"-p " + " -p ".join(PROJECT)} {"-d " + " -d ".join(DRUGS)} {"-C " + "-C ".join([str(i) + " " for i in cutoffs])}{"-t " + "-t ".join([str(i) + " " for i in threshold])}-o {OUTPUT_PATH} -c {cores}\n')

    # if the parameters were set by interactive mode, ask here one last time if
    # the analysis shall be started:
    if interacitve_proj and interacitve_drugs:
        while True:
            print('press ENTER to start or q to quit:')
            start_or_not = input()
            if start_or_not == 'q':
                os._exit(0)
            if start_or_not == '':
                break
            else:
                continue

    shared_scriptdir = os.path.join(
        os.path.split(os.path.split(SCRIPT_PATH)[0])[0], 'shared')

    Snakemake_all_files = []

    Snakefile = os.path.join(os.path.split(SCRIPT_PATH)[0], 'Snakefile')
    config_file_shared = os.path.join(shared_scriptdir, 'config.yaml')

    # help files for both pipelines, like:
    # OUTPUT_PATH/metadata/gdc_manifest_20211029_data_release_31...,
    # gencode.v36.annotation
    help_file_list = download_with_api.download_help_files(
        OUTPUT_PATH, config_file_shared)
    Snakemake_all_files = Snakemake_all_files + help_file_list
    projects = '_'.join(PROJECT)

    def return_type(pipeline):
        if pipeline == "DESeq2":
            return "norm_count"
        elif pipeline == "metilene":
            return "beta_vals"
    types = [return_type(i) for i in execute]
    drug_str = '_'.join(DRUGS)

    # count_type = []  # the count type is later specified within either deseq or
    count_type = {'metilene': ['beta_vals'], 'DESeq2': ['norm_count']}
    # metilene main module, must be defined here already since it is also set in
    # the shared Snakefile
    config = {'thresh': thresh_str, 'thresh_list': thresh_list, 'pipelines':
              execute, 'projects_str': projects, 'cutoffs': cutoffs, 'types':
              types, 'OUTPUT_PATH': OUTPUT_PATH, 'drug_str': drug_str,
              'count_type': count_type}
    # once we have to call snakemake in prior, s.t. the manifest file is
    # present on which all the following selections are done on, make sure that
    # here the dryrun flag is not set to True
    # # TODO uncomment this !!!
    if not report:
        workflow = snakemake.snakemake(snakefile=Snakefile,
                                       targets=Snakemake_all_files,
                                       rerun_triggers='mtime',
                                       workdir=OUTPUT_PATH, cores=cores,
                                       forceall=False, force_incomplete=True,
                                       dryrun=dryrun, use_conda=True,
                                       configfiles=[config_file_shared],
                                       config=config, unlock=unlock)
        if unlock:
            print(f'The output dir:\n{OUTPUT_PATH}\nis unlocked now')
            os._exit(0)
        if not workflow:
            print('snakemake execution failed, exiting now')
            os._exit(0)
    # TODO uncomment this !!!

    # auxfiles for both pipelines:
    # OUTPUT_PATH/PROJECT/aux_files/nationwidechildrens.....
    aux_file_list = download_with_api.download_aux_files(OUTPUT_PATH, PROJECT,
                                                         config_file_shared)

    Snakemake_all_files = Snakemake_all_files + aux_file_list

    def map_execute(pipeline):
        """
        translate here the applied pipeline  which shall be executet to its datatype:
        Datafiles: OUTPUT_PATH/PROJECT/analyse_type/data_files/...

        :param str pipeline: either "DESeq2" or "metilene"
        """
        if pipeline == 'DESeq2':
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

    print('running snakemake with\n')
    print(f'Snakefile:\t{Snakefile}')
    print(f'shared_scriptdir:\t{shared_scriptdir}')

    # also add the multi proj meta_info_druglist_merged_drugs_combined.tsv
    # which is just the concatenation of the single proj pendants:
    # by that those singl proj meta tables are created aswell
    merged_drugs_combined_list = []

    cutoffs_str = [f'cutoff_{str(i)}' for i in cutoffs]
    for pipeline in execute:
        for cutoff in cutoffs_str:
            merged_drugs_combined_list.append(os.path.join(
                OUTPUT_PATH, projects, pipeline, 'merged_meta_files', cutoff,
                'meta_info_druglist_merged_drugs_combined.tsv'))

    Snakemake_all_files = Snakemake_all_files + merged_drugs_combined_list

    # ########################################################################
    # # TODO uncomment this !!!
    # ########################################################################
    # # rerun_trggers='mtime' IMPORTANT -> the needed input data is generated via
    # # input function in shared snakefile for rule merge_meta_tables:
    # # without this option this rule would be ran everytime, since everything
    # # following is based on those aux file, every rule would be triggered
    # # then again
    if not report:
        workflow = snakemake.snakemake(snakefile=Snakefile,
                                       targets=Snakemake_all_files,
                                       workdir=OUTPUT_PATH, cores=cores,
                                       forceall=False, force_incomplete=True,
                                       dryrun=dryrun, use_conda=True,
                                       rerun_triggers='mtime',
                                       configfiles=[config_file_shared],
                                       config=config)
        if not workflow:
            print('snakemake execution failed, exiting now')
            os._exit(0)
    ########################################################################
    # TODO uncomment this !!!
    ########################################################################

    # to download every meta table :
    # run this command and exit here with os._exit(0)
    # temp=("-p TCGA-CESC" "-p TCGA-HNSC" "-p TCGA-LUSC" "-p TCGA-ESCA" "-p TCGA-BRCA" "-p TCGA-GBM" "-p TCGA-OV" "-p TCGA-LUAD" "-p TCGA-UCEC" "-p TCGA-KIRC" "-p TCGA-LGG" "-p TCGA-THCA" "-p TCGA-PRAD" "-p TCGA-SKCM" "-p TCGA-COAD" "-p TCGA-STAD" "-p TCGA-BLCA" "-p TCGA-LIHC" "-p TCGA-KIRP" "-p TCGA-SARC" "-p TCGA-PAAD" "-p TCGA-PCPG" "-p TCGA-READ" "-p TCGA-TGCT" "-p TCGA-THYM" "-p TCGA-KICH" "-p TCGA-ACC" "-p TCGA-MESO" "-p TCGA-UVM" "-p TCGA-DLBC" "-p TCGA-UCS" "-p TCGA-CHOL")
    # for i in ${temp[@]}; do echo tcga_pipelines -p $i -d cisplatin -o /scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_7 -c 40; done
    # if the -D option is set, we can quit right here:
    if download:
        os._exit(0)

    # from here the shared modules and Snakemake scripts are getting pipeline
    # specific, hand over all outputfiles requested so far and enter the
    # pipeline specific main files:
    Snakemake_report_met = []
    Snakemake_report_des = []
    if 'metilene' in execute:
        Snakemake_report_met = main_metilene.entry_fct(OUTPUT_PATH, PROJECT,
                                                       DRUGS,
                                                       Snakemake_all_files,
                                                       cutoffs, threshold,
                                                       cores, 'metilene',
                                                       config_file_shared,
                                                       config, dryrun,
                                                       cutoffs_str, report)
    if 'DESeq2' in execute:
        Snakemake_report_des = main_deseq.entry_fct(OUTPUT_PATH, PROJECT, DRUGS,
                                                    Snakemake_all_files,
                                                    threshold, cores, 'DESeq2',
                                                    config_file_shared, config,
                                                    dryrun, cutoffs_str, report)

    Snakemake_report_files = Snakemake_report_des + Snakemake_report_met

    # # one or both pipelines are finished here, final aggregation over both
    # pipelines:

    # the final majority vote file, aggregated over all pipelines:
    # depending on the pipeline chosen:

    # TODO also include threshold into majortiy vote!!
    major_file = os.path.join(
        OUTPUT_PATH, projects, '_'.join(execute),
        '_'.join(DRUGS),'_'.join([f'threshold_{str(i)}' for i in threshold]), 'final_majority_vote.tsv.gz')
    major_file_pdf = os.path.join(
        OUTPUT_PATH, projects, '_'.join(execute),
        '_'.join(DRUGS), '_'.join([f'threshold_{str(i)}' for i in threshold]), 'final_majority_vote_pipeline_project_final.pdf')

    p_val_prod_sum = []
    if 'DESeq2' in execute:
        p_val_prod_sum.extend([os.path.join(OUTPUT_PATH, projects, 'DESeq2', 'DESeq2_output', '_'.join(DRUGS), 'female_male', '-'.join(cutoffs_str), '_'.join([f'threshold_{str(i)}' for i in threshold]), 'DESeq2-norm_count_p_prod_sum.pdf')])
    if 'metilene' in execute:
        p_val_prod_sum.extend([os.path.join(OUTPUT_PATH, projects, 'metilene', 'metilene_output', '_'.join(DRUGS), 'female_male', '-'.join(cutoffs_str), '_'.join([f'threshold_{str(i)}' for i in threshold]), 'metilene-beta_vals_p_prod_sum.pdf')])

    shared_pipeline_files = [major_file, major_file_pdf] + p_val_prod_sum
    Snakemake_report_files = Snakemake_report_files + shared_pipeline_files + p_val_prod_sum
    Snakemake_report_files.sort()

    if not report:
        workflow = snakemake.snakemake(snakefile=Snakefile,
                                       targets=shared_pipeline_files,
                                       workdir=OUTPUT_PATH, cores=cores,
                                       forceall=False, force_incomplete=True,
                                       dryrun=dryrun, use_conda=True,
                                       configfiles=[config_file_shared],
                                       rerun_triggers='mtime', config=config)
        if not workflow:
            print('snakemake execution failed, exiting now')
            os._exit(0)

    report_file = os.path.join(OUTPUT_PATH, projects, '_'.join(execute),
                               '_'.join(DRUGS), '_'.join([f'threshold_{str(i)}' for i in threshold]), 'report.html')

    ###########################################################################
    #                            final REPORT creation                        #
    ###########################################################################
    # sort the report files according to their subcategory:

    # also with dryrun set, the report file would be created, catch that
    # beforehand:
    # TODO adjust the
    # /homes/biertruck/gabor/phd/test_git_doc/TREMSUCS/src/shared/report_src/workflow.rst
    # s.t. the pipeline call is included there:
    # before creating the report, adjust the rst file, path is directly given
    # within the shared Snakefile
    rst_file = os.path.join(os.path.dirname(__file__), os.path.pardir, 'report_src', 'workflow.rst' )
    with open(rst_file, 'w') as f:
        f.write('| TREMSUCS TCGA final report\n')
        f.write('| Issued CLI call:\n\n')
        temp_str = f'| TREMSUCS {"-p " + " -p ".join(PROJECT)}\n| {"-d " + " -d ".join(DRUGS)}\n| {"-C " + "-C ".join([str(i) + " " for i in cutoffs])}\n| {"-t " + "-t ".join([str(i) + " " for i in threshold])}\n| -o {OUTPUT_PATH}\n| -c {cores}'
        f.write(temp_str)

    if not dryrun:
        workflow = snakemake.snakemake(snakefile=Snakefile,
                                       targets=Snakemake_report_files,
                                       workdir=OUTPUT_PATH, cores=cores,
                                       forceall=False, force_incomplete=True,
                                       use_conda=True,
                                       configfiles=[config_file_shared],
                                       rerun_triggers='mtime', config=config,
                                       report=report_file,
                                       report_stylesheet=os.path.join(
                                           os.path.split(__file__)[0],
                                           os.pardir, "report_src",
                                           "custom-stylesheet.css"))

        if not workflow:
            print('snakemake execution failed, exiting now')
            os._exit(0)
