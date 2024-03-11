import snakemake
import os
from tcga_deseq.modules import create_summary_table
from tcga_deseq.modules import create_deseq_output
from tcga_deseq.modules import create_deseq_lifeline_plots
from tcga_deseq.modules import gzip_counts_all_cases
from tcga_deseq.modules import aggregate_lifelines_all
# from shared.modules import aggregate_lifelines_all


"""
coming from src/shared/modules/main.py
generating the files requested through the deseq Snakefile in
src/tcga_deseq/Snakefile
"""


def entry_fct(OUTPUT_PATH, PROJECT, DRUGS, Snakemake_all_files,
              threshold, cores, pipeline, config_file_shared, config, dryrun,
              cutoffs_str, report):

    """
    the different count types are created for the heatmaps, not needed for the lifeline creation:
    """
    count_types = ['norm_count']
    # config['count_type'] = count_types
    SCRIPT_PATH = os.path.split(__file__)[0]

    # thresh_list = [f'threshold_{str(i)}' for i in threshold]

    PROJECTS = []
    if len(PROJECT) > 1:
        PROJECTS.extend(PROJECT)
        PROJECTS.append('_'.join(sorted([x.upper() for x in PROJECT])))
    else:
        PROJECTS = PROJECT

    DRUG_str = '_'.join(sorted(DRUGS))

    Snakefile = os.path.join(os.path.split(SCRIPT_PATH)[0], os.path.pardir, 'shared', 'Snakefile')

    summary_tables_list = create_summary_table.create_summary_table(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs_str)  # ['cutoff_0', 'cutoff_5', 'cutoff_8']
    Snakemake_all_files = summary_tables_list

    deseq_output_list = create_deseq_output.create_deseq_output(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs_str)
    Snakemake_all_files = Snakemake_all_files + deseq_output_list

    # ##########################################################################
    # the DESeq2 results are created now, based on them the lifelineplots
    # created, ENSG is within the filenames
    # ##########################################################################
    # # # TODO uncomment this
    if not report:
        workflow = snakemake.snakemake(snakefile=Snakefile,
                                    targets=Snakemake_all_files,
                                    workdir=OUTPUT_PATH, cores=cores,
                                    forceall=False, force_incomplete=True,
                                    dryrun=dryrun, use_conda=True,
                                    printshellcmds=True, quiet=False,
                                    unlock=False, rerun_triggers='mtime',
                                    configfiles=[config_file_shared],
                                    config=config)
        if not workflow:
            print('snakemake execution failed, exiting now')
            os._exit(0)
    # # TODO uncomment this
    ###########################################################################
    ###########################################################################

    deseq_lifeline_list = create_deseq_lifeline_plots.create_lifeline_plots(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs_str, threshold, count_types)
    Snakemake_all_files = Snakemake_all_files + deseq_lifeline_list

    gz_count_files = gzip_counts_all_cases.create_gz_counts(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs_str, count_types)
    Snakemake_all_files = Snakemake_all_files + gz_count_files

    deseq_lifeline_validation_list = create_deseq_lifeline_plots.create_lifeline_plots_validation(deseq_lifeline_list)
    Snakemake_all_files = Snakemake_all_files + deseq_lifeline_validation_list

    # aggregate all outputs, we need to compare the first diff expressions
    # plottet:
    # DESeq2_log2f_DECREASE_norm_ENSG00000000005_lifeline.tsv
    # with the both validation outputs:
    # DESeq2_log2f_DECREASE_norm_ENSG00000000005_lifeline_DOWN_val.tsv
    # DESeq2_log2f_DECREASE_norm_ENSG00000000005_lifeline_UP_val.tsv

    ###########################################################################
    #     here the lifelineplots are created, both, base and validation plots #
    ###########################################################################
    # ## TODO uncomment this !!!!
    if not report:
        workflow = snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
                            workdir=OUTPUT_PATH, cores=cores, forceall=False,
                            force_incomplete=True, dryrun=dryrun, use_conda=True,
                            rerun_triggers='mtime', printshellcmds=True,
                            quiet=False, configfiles=[config_file_shared],
                                    config=config)
        if not workflow:
            print('snakemake run failed, exiting now')
            os._exit(0)
    # ## TODO uncomment this !!!!
    ###############################################################################
    ## the previous snakemake runs must be completed before requesting the next
    # aggregated files. all lifeline base and validation plots must be
    # completed before globbing them out of the threshold dirs and aggregating
    # to one table on which the ranking of potential biomarkers is performed
    #                                                                             #
    ###############################################################################
    # DESeq2_lifelines_aggregated.tsv.gz:
    aggregate_lifelines_list = aggregate_lifelines_all.aggregate_lifeline_plots(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs_str, threshold)

    Snakemake_all_files = aggregate_lifelines_list

    # DESeq2_lifelines_evaluated.pdf:
    evaluate_lifelines_list = [ i.replace(f'{pipeline}_lifelines_aggregated.tsv.gz', f'{pipeline}_lifelines_evaluated.tsv.gz') for i in aggregate_lifelines_list]
    evaluate_lifelines_list = evaluate_lifelines_list + [ i.replace(f'{pipeline}_lifelines_aggregated.tsv.gz', f'{pipeline}_lifelines_evaluated.pdf') for i in aggregate_lifelines_list]
    eval_new = []
    # just request the evaluated norm counts, s.t. they are added to the report
    for count_type in count_types:
        for eval_file in evaluate_lifelines_list:
            eval_new.append(eval_file.replace('evaluated', f'evaluated-{count_type}'))

    Snakemake_all_files = Snakemake_all_files + eval_new
    Snakemake_report_files = eval_new

    # DESeq2_plot_diffs.pdf:
    plot_diffs_all = []
    plot_types = ['base_plot', 'UP_validation', 'DOWN_validation']
    for plot_type in plot_types:
        for i in aggregate_lifelines_list:
            plot_diffs_all.append(i.replace('lifelines_aggregated.tsv.gz', f'plot_diffs_{plot_type}.pdf'))
    plot_diffs_new = []
    # also here just request the norm counts
    for plot in plot_diffs_all:
        for count_type in count_types:
            plot_diffs_new.append(plot.replace('.pdf', f'-{count_type}.pdf'))

    Snakemake_all_files = Snakemake_all_files + plot_diffs_new
    plot_diffs_eval_all = [i.replace('plot_diffs', 'plot_eval_diffs') for i in plot_diffs_new]
    Snakemake_all_files = Snakemake_all_files + plot_diffs_eval_all

    patients_overview = []
    heatmap_merged = []
    for project in PROJECTS:
        for cutoff in cutoffs_str:
            patients_overview.append(os.path.join(OUTPUT_PATH, project, pipeline, 'merged_meta_files', cutoff, DRUG_str, 'meta_info_druglist_merged_drugs_combined_final.pdf'))
            for gender in ['male', 'female', 'female_male']:
                heatmap_merged.append(os.path.join(OUTPUT_PATH, project, pipeline, pipeline + '_output', DRUG_str, gender, cutoff, 'DESeq2_heatmap_merged.pdf'))
    Snakemake_report_files = Snakemake_report_files + heatmap_merged + patients_overview

    Snakemake_all_files = Snakemake_all_files + patients_overview
    Snakemake_all_files = Snakemake_all_files + heatmap_merged
    merged_diffs_list = [i.replace('lifelines_aggregated.tsv.gz', 'plot_aggr+eval_diffs_merged.pdf') for i in aggregate_lifelines_list]
    Snakemake_all_files = Snakemake_all_files + merged_diffs_list
    Snakemake_report_files = Snakemake_report_files + merged_diffs_list
    # TODO
    if not report:
        workflow = snakemake.snakemake(snakefile=Snakefile, targets= Snakemake_all_files,
                            workdir=OUTPUT_PATH, cores=cores, forceall=False,
                            force_incomplete=True, dryrun=dryrun, use_conda=True,
                            rerun_triggers='mtime', printshellcmds=True,
                            quiet=False, verbose=False,
                                    configfiles=[config_file_shared],
                                    config=config)

        if not workflow:
            print('snakemake run failed, exiting now')
            os._exit(0)
    return Snakemake_report_files
