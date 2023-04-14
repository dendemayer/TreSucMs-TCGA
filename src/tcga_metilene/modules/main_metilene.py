import snakemake
import os
from tcga_metilene.modules import create_summary_table
from tcga_metilene.modules import create_metilene_out
from tcga_metilene.modules import bed_intersect_metilene
from tcga_metilene.modules import create_lifeline_plot
from tcga_metilene.modules import create_lifeline_plots_validation

"""
coming from src/shared/modules/main.py
generating the files requested through the metilene Snakefile in
src/tcga_metilene/Snakefile
"""

def entry_fct(OUTPUT_PATH, PROJECT, DRUGS, Snakemake_all_files, cutoffs,
              threshold, cores):

    SCRIPT_PATH = os.path.split(__file__)[0]

    PROJECTS = []
    if len(PROJECT) > 1:
        PROJECTS.extend(PROJECT)
        PROJECTS.append('_'.join(sorted([x.upper() for x in PROJECT])))

    shared_workdir = os.path.join(
        os.path.split(os.path.split(SCRIPT_PATH)[0])[0], 'shared')
    Snakefile = os.path.join(os.path.split(SCRIPT_PATH)[0], 'Snakefile')

    summary_tables = create_summary_table.return_summary_tables(
        OUTPUT_PATH, PROJECTS, DRUGS, cutoffs)

    Snakemake_all_files = Snakemake_all_files + summary_tables

    metilene_out_tables = \
        create_metilene_out.return_metilene_tables(
            OUTPUT_PATH, PROJECTS, DRUGS, cutoffs)

    Snakemake_all_files = Snakemake_all_files + metilene_out_tables

    metilene_intersect_tables = \
        bed_intersect_metilene.return_bed_interesect_metilene_files(
            OUTPUT_PATH, PROJECTS, DRUGS, cutoffs)

    Snakemake_all_files = Snakemake_all_files + metilene_intersect_tables

    # ## up to this point, every download, preprocessing and metilene analyse
    # steps including intersection of DMR with the original beta_value input
    # tables is performed, every postprocessing step depends on
    # metilene_intersect.tsv output tables, they serve as starting input for
    # the third snakemake running instance
    ### make sure that this snakemake process is completed before starting
    ### consecutive processes
    # TODO uncomment this !!!
    # snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
    #                     workdir=shared_workdir, cores=cores, forceall=False,
    #                     force_incomplete=True, dryrun=False, use_conda=True)
    # TODO uncomment this !!!

    merged_plots = bed_intersect_metilene.return_DMR_merge_plots(OUTPUT_PATH, PROJECTS, DRUGS, cutoffs, threshold)

    # rerun_triggers='mtime'

    # next to do, plot the regions found with its betavalues
    # input is the metilene_intersect.tsv, f.e.:
    # /scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines/TCGA-HNSC/metilene/metilene_output/carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_intersect.tsv
    # plot for each range found in the intersect.tsv a regions plot

    # metilene_plots = bed_intersect_metilene.return_plot_DMR_regions_plot(
        # metilene_intersect_tables)

    # Start Snakemake_all_files new over, since the rest must be created at
    # this point and the DAG resolvement is faster than
    # Snakemake_all_files = metilene_plots

    # the lifeline regression is the first step where the threshold is invoked,
    # for every so far created metilene out, create a lifeline for each
    # threshold (the applied metilene_plots list already contains the cutoff
    # permutation...)
    # lifeline_plots = create_lifeline_plot.create_lifeline_plots(
        # metilene_plots, threshold)

    # Snakemake_all_files = Snakemake_all_files + lifeline_plots

    # validation plots for the found DMRs:
    # validation_plots = create_lifeline_plots_validation.validation_plots(lifeline_plots)

    # Snakemake_all_files = Snakemake_all_files + validation_plots

    # merge the DMR related plots threshold specific lifeline plots and not
    # threshold specific DMR box and lineplots:
    # merged_plots = [os.path.join(j, 'metilene_merged_lifeline_plot.pdf') for j in list(set([os.path.split(i)[0] for i in lifeline_plots]))]
    # merged_plots = merged_plots + [os.path.join(j, 'metilene_merged_boxplot_beta_value.pdf') for j in list(set([os.path.split(os.path.split(i)[0])[0] for i in lifeline_plots]))]
    # merged_plots = merged_plots + [os.path.join(j, 'metilene_merged_lineplot_median_beta_value.pdf') for j in list(set([os.path.split(os.path.split(i)[0])[0] for i in lifeline_plots]))]
    # merged_plots

    Snakemake_all_files = Snakemake_all_files + merged_plots

    # TODO
    snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
                        workdir=shared_workdir, cores=cores, forceall=False,
                        force_incomplete=True, dryrun=False, use_conda=True,
                        printshellcmds=True,  rerun_triggers='mtime')
    # TODO

# # ##### main_metilene ############
