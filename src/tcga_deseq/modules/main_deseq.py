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


def entry_fct(OUTPUT_PATH, PROJECT, DRUGS, Snakemake_all_files, cutoffs,
              threshold, cores, pipeline):

    SCRIPT_PATH = os.path.split(__file__)[0]

    thresh_list = [f'threshold_{str(i)}' for i in threshold]
    thresh_str = '_'.join(thresh_list)

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
    # snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
    #                     workdir=shared_workdir, cores=cores, forceall=False,
    #                     force_incomplete=True, dryrun=False, use_conda=True, printshellcmds=True, quiet=False, unlock=False)
    # TODO uncomment this
    ###########################################################################
    ###########################################################################

    deseq_lifeline_list = create_deseq_lifeline_plots.create_lifeline_plots(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs, threshold)
    Snakemake_all_files = Snakemake_all_files + deseq_lifeline_list

    gz_count_files = gzip_counts_all_cases.create_gz_counts(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs)
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
    ### TODO uncomment this !!!!
    # workflow = snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
    #                     workdir=shared_workdir, cores=cores, forceall=False,
    #                     force_incomplete=True, dryrun=False, use_conda=True,
    #                     rerun_triggers='mtime', printshellcmds=True,
    #                     quiet=False, config={'thresh': thresh_str, 'thresh_list': thresh_list})
    # ### TODO uncomment this !!!!
    # if not workflow:
    #     print('snakemake run failed, exiting now')
    #     os._exit(0)
    ###############################################################################
    ## the previous snakemake runs must be completed before requesting the next
    # aggregated files. all lifeline base and validation plots must be
    # completed before globbing them out of the threshold dirs and aggregating
    # to one table on which the ranking of potential biomarkers is performed
    #                                                                             #
    ###############################################################################
    # DESeq2_lifelines_aggregated.tsv.gz:
    aggregate_lifelines_list = aggregate_lifelines_all.aggregate_lifeline_plots(OUTPUT_PATH, PROJECTS, DRUG_str, cutoffs, threshold)

    Snakemake_all_files = aggregate_lifelines_list

    # DESeq2_lifelines_evaluated.pdf:
    evaluate_lifelines_list = [ i.replace(f'{pipeline}_lifelines_aggregated.tsv.gz', f'{pipeline}_lifelines_evaluated.tsv.gz') for i in aggregate_lifelines_list]
    evaluate_lifelines_list = evaluate_lifelines_list + [ i.replace(f'{pipeline}_lifelines_aggregated.tsv.gz', f'{pipeline}_lifelines_evaluated.pdf') for i in aggregate_lifelines_list]

    Snakemake_all_files = Snakemake_all_files + evaluate_lifelines_list

    # DESeq2_plot_diffs.pdf:
    plot_diffs_all = []
    plot_types = ['base_plot', 'UP_validation', 'DOWN_validation']
    for plot_type in plot_types:
        for i in aggregate_lifelines_list:
            plot_diffs_all.append(i.replace('lifelines_aggregated.tsv.gz', f'plot_diffs_{plot_type}.pdf'))

    Snakemake_all_files = Snakemake_all_files + plot_diffs_all

    workflow = snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files,
                        workdir=shared_workdir, cores=cores, forceall=False,
                        force_incomplete=True, dryrun=False, use_conda=True,
                        rerun_triggers='mtime', printshellcmds=True,
                        quiet=False, verbose=False, config={'thresh': thresh_str, 'thresh_list': thresh_list})

    if not workflow:
        print('snakemake run failed, exiting now')
        os._exit(0)
