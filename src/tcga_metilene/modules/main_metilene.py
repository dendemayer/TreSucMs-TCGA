import snakemake
import os
# from tcga_metilene.modules import create_summary_table
# from tcga_metilene.modules import create_metilene_out
from tcga_metilene.modules import bed_intersect_metilene
from tcga_metilene.modules import create_lifeline_plot
from tcga_metilene.modules import create_lifeline_plots_validation
from shared.modules import aggregate_lifelines_all

"""
coming from src/shared/modules/main.py
generating the files requested through the metilene Snakefile in
src/tcga_metilene/Snakefile
"""


def entry_fct(OUTPUT_PATH, PROJECT, DRUGS, Snakemake_all_files, cutoffs,
              threshold, cores, pipeline, config_file_shared, config, dryrun,
              cutoffs_str, report):

    SCRIPT_PATH = os.path.split(__file__)[0]

    cutoffs = [f'cutoff_{str(i)}' for i in cutoffs]

    # config['count_type'] = ['beta_vals']
    # an extra parameter can be handed over to the snakefile via the config
    # method in he snakemake call
    # # to be able to aggregate over the thresholds, txt file must be written
    # # into file. , this is then read by Snakemake
    # thresh_file = os.path.join(os.path.split(SCRIPT_PATH)[0], 'dynamic',
    #                            'thresholds.txt')
    # thresh_list = [f'threshold_{str(i)}' for i in threshold]

    PROJECTS = []
    if len(PROJECT) > 1:
        PROJECTS.extend(PROJECT)
        PROJECTS.append('_'.join(sorted([x.upper() for x in PROJECT])))
    else:
        PROJECTS = PROJECT

    DRUG_str = '_'.join(sorted(DRUGS))

    # shared_workdir = os.path.join(
    #     os.path.split(os.path.split(SCRIPT_PATH)[0])[0], 'shared')
    # Snakefile = os.path.join(os.path.split(SCRIPT_PATH)[0], 'Snakefile')
    Snakefile = os.path.join(os.path.split(SCRIPT_PATH)[0], os.path.pardir,
                             'shared', 'Snakefile')

    metilene_intersect_tables = \
        bed_intersect_metilene.return_bed_interesect_metilene_files(
            OUTPUT_PATH, PROJECTS, DRUGS, cutoffs)
    Snakemake_all_files = Snakemake_all_files + metilene_intersect_tables

    # IMPORTANT the intersect tables define the DMRs which are requested in the
    # following, they must be completed here!
    # TODO uncomment this !!!
    if not report:
        workflow = snakemake.snakemake(snakefile=Snakefile,
                                       targets=Snakemake_all_files,
                                       workdir=OUTPUT_PATH, cores=cores,
                                       forceall=False, rerun_triggers='mtime',
                                       force_incomplete=True, dryrun=dryrun,
                                       use_conda=True,
                                       configfiles=[config_file_shared],
                                       config=config)
        if not workflow:
            print('snakemake execution failed, exiting now')
            os._exit(0)
    # TODO uncomment this !!!

    # ## up to this point, every download, preprocessing and metilene analyse
    # steps including intersection of DMR with the original beta_value input
    # tables is performed, every postprocessing step depends on
    # metilene_intersect.tsv output tables, they serve as starting input for
    # the third snakemake running instance
    # ## make sure that this snakemake process is completed before starting
    # ## consecutive processes

    # next to do, plot the regions found with its betavalues
    # input is the metilene_intersect.tsv, f.e.:
    # /scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines/TCGA-HNSC/metilene/
    # metilene_output/carboplatin,paclitaxel_cisplatin/female/cutoff_0/
    # metilene_intersect.tsv

    # plot for each range found in the intersect.tsv a regions plot

    metilene_plots = bed_intersect_metilene.return_plot_DMR_regions_plot(
        metilene_intersect_tables)

    # Start Snakemake_all_files new over, since the rest must be created at
    # this point
    Snakemake_all_files = metilene_plots

    # the lifeline regression is the first step where the threshold is invoked,
    # for every so far created metilene out, create a lifeline for each
    # threshold (the applied metilene_plots list already contains the cutoff
    # permutation...)
    lifeline_plots = create_lifeline_plot.create_lifeline_plots(
        metilene_plots, threshold)

    Snakemake_all_files = Snakemake_all_files + lifeline_plots

    # validation plots for the found DMRs:
    validation_plots = create_lifeline_plots_validation.validation_plots(
        lifeline_plots)

    Snakemake_all_files = Snakemake_all_files + validation_plots

    # merge the DMR related plots threshold specific lifeline plots and not
    # threshold specific DMR box and lineplots:

    merged_plots = bed_intersect_metilene.return_DMR_merge_plots(OUTPUT_PATH,
                                                                 PROJECTS,
                                                                 DRUGS,
                                                                 cutoffs,
                                                                 threshold)

    Snakemake_all_files = Snakemake_all_files + merged_plots
    # TODO
    if not report:
        workflow = snakemake.snakemake(snakefile=Snakefile,
                                       targets=Snakemake_all_files,
                                       workdir=OUTPUT_PATH, cores=cores,
                                       forceall=False, force_incomplete=True,
                                       dryrun=dryrun, use_conda=True,
                                       printshellcmds=True,
                                       rerun_triggers='mtime',
                                       configfiles=[config_file_shared],
                                       config=config, unlock=False)
        if not workflow:
            print('snakemake execution failed, exiting now')
            os._exit(0)
    # TODO

    aggregate_lifelines_list = \
        aggregate_lifelines_all.aggregate_lifeline_plots(OUTPUT_PATH, PROJECTS,
                                                         DRUG_str, cutoffs,
                                                         threshold, pipeline)
    Snakemake_all_files = Snakemake_all_files + aggregate_lifelines_list

    evaluate_lifelines_list = [i.replace(f'{pipeline}_lifelines_aggregated.tsv.gz', f'{pipeline}_lifelines_evaluated.tsv.gz') for i in aggregate_lifelines_list]
    evaluate_lifelines_list = evaluate_lifelines_list + [i.replace(f'{pipeline}_lifelines_aggregated.tsv.gz', f'{pipeline}_lifelines_evaluated.pdf') for i in aggregate_lifelines_list]
    eval_new = [i.replace('evaluated', 'evaluated-beta_vals') for i in evaluate_lifelines_list]

    Snakemake_all_files = Snakemake_all_files + eval_new
    Snakemake_report_files = eval_new

    plot_diffs_all = []
    plot_types = ['base_plot', 'UP_validation', 'DOWN_validation']
    for plot_type in plot_types:
        for i in aggregate_lifelines_list:
            plot_diffs_all.append(i.replace('lifelines_aggregated.tsv.gz', f'plot_diffs_{plot_type}-beta_vals.pdf'))

    Snakemake_all_files = Snakemake_all_files + plot_diffs_all
    # make the diffs for the evaluated features:
    plot_diffs_eval_all = [i.replace('plot_diffs', 'plot_eval_diffs') for i in plot_diffs_all]
    Snakemake_all_files = Snakemake_all_files + plot_diffs_eval_all

    # /scr/palinca/gabor/TCGA-pipeline_5/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/cutoff_0/meta_info_druglist_merged_drugs_combined_final.pdf
    DMR_intersect_merged = []  # they contain the metilene_intersect_[boxplot|heatmaps|lineplot|violinplot]_beta_value_{DMR}.pdf plots
    patients_overview = []
    for project in PROJECTS:
        for cutoff in cutoffs_str:
            patients_overview.append(os.path.join(OUTPUT_PATH, project, pipeline, 'merged_meta_files', cutoff, DRUG_str, 'meta_info_druglist_merged_drugs_combined_final.pdf'))
            for gender in ['male', 'female', 'female_male']:
                DMR_intersect_merged.append(os.path.join(OUTPUT_PATH, project, pipeline, pipeline + '_output', DRUG_str, gender, cutoff, 'metilene_intersect_beta_value_merged.pdf'))

    Snakemake_report_files = Snakemake_report_files + DMR_intersect_merged + patients_overview
    Snakemake_all_files = Snakemake_all_files + patients_overview + DMR_intersect_merged
    merged_diffs_list = [i.replace('lifelines_aggregated.tsv.gz', 'plot_aggr+eval_diffs_merged.pdf') for i in aggregate_lifelines_list]
    Snakemake_all_files = Snakemake_all_files + merged_diffs_list
    Snakemake_report_files = Snakemake_report_files + merged_diffs_list
    # TODO
    if not report:
        workflow = snakemake.snakemake(snakefile=Snakefile,
                                       targets=Snakemake_all_files,
                                       workdir=OUTPUT_PATH, cores=cores,
                                       forceall=False, force_incomplete=True,
                                       dryrun=dryrun, use_conda=True,
                                       printshellcmds=True,
                                       rerun_triggers='mtime',
                                       configfiles=[config_file_shared],
                                       config=config)
        if not workflow:
            print('snakemake execution failed, exiting now')
            os._exit(0)
    return Snakemake_report_files
    # # TODO
# # ##### main_metilene ############
