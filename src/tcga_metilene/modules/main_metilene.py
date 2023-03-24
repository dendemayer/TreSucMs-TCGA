import snakemake
import os
from tcga_metilene.modules import create_summary_table
from tcga_metilene.modules import create_metilene_out
from tcga_metilene.modules import bed_intersect_metilene
from tcga_metilene.modules import lifeline_tables

"""
coming from src/shared/modules/main.py
generating the files requested through the metilene Snakefile in
src/tcga_metilene/Snakefile
"""

def entry_fct(OUTPUT_PATH, PROJECT, DRUGS, Snakemake_all_files, cutoffs,
              threshold, cores):
    SCRIPT_PATH = os.path.split(__file__)[0]
    # metilene_snake_workdir = os.path.join(
    #     os.path.split(os.path.split(SCRIPT_PATH)[0])[0], 'tcga_metilene')

    shared_workdir = os.path.join(
        os.path.split(os.path.split(SCRIPT_PATH)[0])[0], 'shared')
    Snakefile = os.path.join(os.path.split(SCRIPT_PATH)[0], 'Snakefile')
    # config_file = os.path.join(os.path.split(SCRIPT_PATH)[0], 'config.yaml')

    summary_tables = create_summary_table.return_summary_tables(
        OUTPUT_PATH, PROJECT, DRUGS, cutoffs)

    Snakemake_all_files = Snakemake_all_files + summary_tables

    metilene_out_tables = \
        create_metilene_out.return_metilene_tables(
            OUTPUT_PATH, PROJECT, DRUGS, cutoffs)

    Snakemake_all_files = Snakemake_all_files + metilene_out_tables

    metilene_intersect_tables = \
        bed_intersect_metilene.return_bed_interesect_metilene_files(
            OUTPUT_PATH, PROJECT, DRUGS, cutoffs)

    Snakemake_all_files = Snakemake_all_files + metilene_intersect_tables

    # ## up to this point, every download, preprocessing and metilene analyse
    # steps including intersection of DMR with the original beta_value input
    # tables is performed, every postprocessing step depends on
    # metilene_intersect.tsv output tables, they serve as starting input for
    # the third snakemake running instance
    ### make sure that this snakemake process is completed before starting
    ### consecutive processes
    # TODO uncomment this !!!
    snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
                        workdir=shared_workdir, cores=cores, forceall=False,
                        force_incomplete=True, dryrun=False, use_conda=True)
    # TODO uncomment this !!!

    # forcerun=['/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines/TCGA-HNSC/metilene/metilene_output/carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_intersect.tsv'])
    # rerun_triggers='mtime'

    # next to do, plot the regions found with its betavalues
    # input is the metilene_intersect.tsv, f.e.:
    # /scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines/TCGA-HNSC/metilene/metilene_output/carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_intersect.tsv
    # plot for each range found in the intersect.tsv a regions plot

    metilene_plots = bed_intersect_metilene.return_plot_DMR_regions_plot(
        metilene_intersect_tables)
    # Start Snakemake_all_files new over, since the rest must be created at
    # this point and the DAG resolvement is faster than
    Snakemake_all_files = metilene_plots
    snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
                        workdir=shared_workdir, cores=cores, forceall=False,
                        force_incomplete=True, dryrun=False, debug=False,
                        use_conda=True)
    # lifeline_tables = lifeline_tables

# # ##### main_metilene ############
