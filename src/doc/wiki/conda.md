# install conda:


wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh

mamba install -y -c conda-forge r-base pandas matplotlib lifelines pandoc r-sjmisc requests; mamba install -y -c bioconda bioconductor-deseq2 bioconductor-pcaexplorer snakemake mygene; mamba install -y -c pkgs/main click seaborn




mamba install -c bioconda bioconductor-deseq2=1.30.1
mamba install -c bioconda bioconductor-pcaexplorer=2.16.0

conda install mamba
mamba install r-basemamba install r-base=4.0.3
mamba install -c bioconda bioconductor-deseq2=1.30.0

mamba install bioconductor-pcaexplorer=2.16.0
mamba install -c conda-forge pandas=1.2.1
mamba install -c conda-forge matplotlib=3.3.4
mamba install -c conda-forge lifelines=0.25.9
mamba install -c pkgs/main click=7.1.2
    conda env export > deseq_env.yaml
mamba install -c conda-forge r-sjmisc=2.8.6
    conda env export > deseq_env_2.yaml
mamba install -c pkgs/main seaborn=0.11.1
    conda env export > deseq_env_3.yaml
mamba install -c conda-forge visidata=2.2.1
mamba install -c conda-forge sphinx=3.5.0
mamba install -c conda-forge pandoc=2.11.4
    conda env export > deseq_env_4.yaml
%%mamba install -c pkgs/main reportlab=3.5.60
mamba install -c bioconda snakemake=6.1.1
mamba install -c bioconda mygene=3.2.2
