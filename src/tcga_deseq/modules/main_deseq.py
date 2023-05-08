import snakemake
import os
from tcga_deseq.modules import create_summary_table
from tcga_deseq.modules import create_deseq_output
from tcga_deseq.modules import create_deseq_lifeline_plots
from tcga_deseq.modules import gzip_counts_all_cases


"""
coming from src/shared/modules/main.py
generating the files requested through the deseq Snakefile in
src/tcga_deseq/Snakefile
"""

def entry_fct(OUTPUT_PATH, PROJECT, DRUGS, Snakemake_all_files, cutoffs,
              threshold, cores):

    SCRIPT_PATH = os.path.split(__file__)[0]

    PROJECTS = []
    if len(PROJECT) > 1:
        PROJECTS.extend(PROJECT)
        PROJECTS.append('_'.join(sorted([x.upper() for x in PROJECT])))
    else:
        PROJECTS = PROJECT

    DRUG_str = '_'.join(DRUGS)

    shared_workdir = os.path.join(
        os.path.split(os.path.split(SCRIPT_PATH)[0])[0], 'shared')
    Snakefile = os.path.join(os.path.split(SCRIPT_PATH)[0], 'Snakefile')

    summary_tables_list = create_summary_table.create_summary_table(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs)
    Snakemake_all_files = summary_tables_list

    deseq_output_list = create_deseq_output.create_deseq_output(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs)
    Snakemake_all_files = Snakemake_all_files + deseq_output_list

    ###########################################################################
    # the DESeq2 results are created now, based on them the lifelineplots
    # created, ENSG is within the filenames
    ###########################################################################
    # TODO uncomment this
    snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
                        workdir=shared_workdir, cores=cores, forceall=False,
                        force_incomplete=True, dryrun=True, use_conda=True)
    # TODO uncomment this
    ###########################################################################
    # the DESeq2 results are created now, based on them the lifelineplots
    # created, ENSG is within the filenames
    ###########################################################################

    deseq_lifeline_list = create_deseq_lifeline_plots.create_lifeline_plots(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs, threshold)
    Snakemake_all_files = Snakemake_all_files + deseq_lifeline_list

    gz_count_files = gzip_counts_all_cases.create_gz_counts(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs)
    Snakemake_all_files = Snakemake_all_files + gz_count_files

    deseq_lifeline_validation_list = create_deseq_lifeline_plots.create_lifeline_plots_validation(deseq_lifeline_list)
    Snakemake_all_files = Snakemake_all_files + deseq_lifeline_validation_list

    snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
                        workdir=shared_workdir, cores=cores, forceall=False,
                        force_incomplete=True, dryrun=True, use_conda=True, rerun_triggers='mtime', printshellcmds=True)
