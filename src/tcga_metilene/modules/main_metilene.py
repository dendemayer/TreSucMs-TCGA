import snakemake
import os
from tcga_metilene.modules import create_summary_table
from tcga_metilene.modules import create_metilene_out

def entry_fct(OUTPUT_PATH, PROJECT, DRUGS, Snakemake_all_files, cutoffs, threshold, cores):
    SCRIPT_PATH = os.path.split(__file__)[0]
    metilene_snake_workdir = os.path.join(
        os.path.split(os.path.split(SCRIPT_PATH)[0])[0], 'tcga_metilene')

    shared_workdir = os.path.join(
        os.path.split(os.path.split(SCRIPT_PATH)[0])[0], 'shared')
    Snakefile = os.path.join(os.path.split(SCRIPT_PATH)[0], 'Snakefile')
    config_file = os.path.join(os.path.split(SCRIPT_PATH)[0], 'config.yaml')

    summary_tables = create_summary_table.return_summary_tables(
        OUTPUT_PATH, PROJECT, DRUGS, cutoffs)

    Snakemake_all_files = Snakemake_all_files + summary_tables

    metilene_out_tables = create_metilene_out.return_metilene_tables(OUTPUT_PATH, PROJECT, DRUGS, cutoffs)

    Snakemake_all_files = Snakemake_all_files + metilene_out_tables
    snakemake.snakemake(snakefile=Snakefile, targets=Snakemake_all_files, workdir=shared_workdir, cores=cores, forceall=False, force_incomplete=True, dryrun=False, debug=False, use_conda=True)
# # ##### main_metilene ############
