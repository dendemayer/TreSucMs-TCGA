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
              threshold, cores, pipeline, config_file_shared):

    SCRIPT_PATH = os.path.split(__file__)[0]

    # an extra parameter can be handed over to the snakefile via the config
    # method in he snakemake call
    # # to be able to aggregate over the thresholds, txt file must be written
    # # into file. , this is then read by Snakemake
    # thresh_file = os.path.join(os.path.split(SCRIPT_PATH)[0], 'dynamic',
    #                            'thresholds.txt')
    thresh_list = [f'threshold_{str(i)}' for i in threshold]

    # with open(thresh_file, 'w') as w:
    #     for i in thresh_list:
    #         w.write(f'{i}\n')
    thresh_str = '_'.join(thresh_list)

    PROJECTS = []
    if len(PROJECT) > 1:
        PROJECTS.extend(PROJECT)
        PROJECTS.append('_'.join(sorted([x.upper() for x in PROJECT])))
    else:
        PROJECTS = PROJECT

    DRUG_str = '_'.join(sorted(DRUGS))

    # shared_workdir = os.path.join(
    #     os.path.split(os.path.split(SCRIPT_PATH)[0])[0], 'shared')
    Snakefile = os.path.join(os.path.split(SCRIPT_PATH)[0], 'Snakefile')

    metilene_intersect_tables = \
        bed_intersect_metilene.return_bed_interesect_metilene_files(
            OUTPUT_PATH, PROJECTS, DRUGS, cutoffs)
    Snakemake_all_files = Snakemake_all_files + metilene_intersect_tables

    # IMPORTANT the intersect tables define the DMRs which are requested in the
    # following, they must be completed here!
    # TODO uncomment this !!!
    workflow = snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
                         workdir=OUTPUT_PATH, cores=cores, forceall=False,
                         force_incomplete=True, dryrun=False, use_conda=True, config={'thresh': thresh_str, 'thresh_list': thresh_list}, configfiles=[config_file_shared])
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
    workflow = snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
                         workdir=OUTPUT_PATH, cores=cores, forceall=False,
                         force_incomplete=True, dryrun=False, use_conda=True,
                         printshellcmds=True,  rerun_triggers='mtime', config={'thresh': thresh_str, 'thresh_list': thresh_list}, configfiles=[config_file_shared])
    #TODO
    if not workflow:
        print('snakemake execution failed, exiting now')
        os._exit(0)

    aggregate_lifelines_list = \
        aggregate_lifelines_all.aggregate_lifeline_plots(OUTPUT_PATH, PROJECTS,
                                                         DRUG_str, cutoffs,
                                                         threshold, pipeline)
    Snakemake_all_files = Snakemake_all_files + aggregate_lifelines_list

    evaluate_lifelines_list = [ i.replace(f'{pipeline}_lifelines_aggregated.tsv.gz', f'{pipeline}_lifelines_evaluated.tsv.gz') for i in aggregate_lifelines_list]
    evaluate_lifelines_list = evaluate_lifelines_list + [ i.replace(f'{pipeline}_lifelines_aggregated.tsv.gz', f'{pipeline}_lifelines_evaluated.pdf') for i in aggregate_lifelines_list]
    eval_new = [i.replace('evaluated', 'evaluated-beta_vals') for i in  evaluate_lifelines_list]

    Snakemake_all_files = Snakemake_all_files + eval_new

    # plot_diffs_all = [ i.replace(f'{pipeline}_lifelines_aggregated.tsv.gz', f'{pipeline}_plot_diffs_base.pdf') for i in aggregate_lifelines_list]
    # plot_diffs_all = plot_diffs_all + [ i.replace(f'{pipeline}_lifelines_aggregated.tsv.gz', f'{pipeline}_plot_diffs.tsv.gz') for i in aggregate_lifelines_list]


    plot_diffs_all = []
    plot_types = ['base_plot', 'UP_validation', 'DOWN_validation']
    for plot_type in plot_types:
        for i in aggregate_lifelines_list:
            plot_diffs_all.append(i.replace('lifelines_aggregated.tsv.gz', f'plot_diffs_{plot_type}-beta_vals.pdf'))

    Snakemake_all_files = Snakemake_all_files + plot_diffs_all
    # make the diffs for the evaluated features:
    plot_diffs_eval_all = [i.replace('plot_diffs', 'plot_eval_diffs') for i in plot_diffs_all]
    Snakemake_all_files = Snakemake_all_files + plot_diffs_eval_all


    #report_file = os.path.join(OUTPUT_PATH, PROJECTS[-1], pipeline, f'{pipeline}_output', DRUG_str, 'report.html')
    # TODO
    workflow = snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
                        workdir=OUTPUT_PATH, cores=cores, forceall=False,
                        force_incomplete=True, dryrun=False, use_conda=True,
                                   printshellcmds=True,  rerun_triggers='mtime', config={'thresh': thresh_str, 'thresh_list': thresh_list}, configfiles=[config_file_shared])
                                   #printshellcmds=True,  rerun_triggers='mtime', config={'thresh': thresh_str, 'thresh_list': thresh_list}, report=report_file)
    # TODO
    if not workflow:
        print('snakemake execution failed, exiting now')
        os._exit(0)
# # ##### main_metilene ############
