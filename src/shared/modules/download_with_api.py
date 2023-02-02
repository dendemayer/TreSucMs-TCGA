import snakemake
def download_GDC_manifest(Snakefile_path, file_to_download, workdir):
    """
    downloading manifest file and the gtf annotation file from TCGA:
    """
    print(f'running snakemake with {Snakefile_path}, {workdir} , {file_to_download}')
    breakpoint()
    snakemake.snakemake(snakefile=Snakefile_path, cores=20,
                        targets=[file_to_download], workdir=workdir,
                        forceall=False)
    # snakemake -p --cores 40 /scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines/metadata/gdc_manifest_20211029_data_release_31.0_active.tsv.gz
