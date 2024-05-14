# TreSucMs README:  
## got to the documentation:  
[Documentation](https://dendemayer.github.io/TreSucMs-TCGA/#)  
## an example Report html file can be downloaded here:  
[example Report download link](https://media.githubusercontent.com/media/dendemayer/TreSucMs-TCGA/main/suppl/report.html?download=true).  
This report file has a size of about 300 MB.

# installing from github.com:
```bash
$ git clone https://github.com/dendemayer/TreSucMs-TCGA.git
$ cd TreSucMs-TCGA
$ pip install .
```

# Start the pipeline interactively or via CLI:

- To start the analysis with help of the interactive mode, call the pipeline
  without any argument:  

```bash
$ TreSucMs
```

- Calling the help or the manual page:  

```bash
$ TreSucMs --help
$ man TreSucMs
```

- example configuration with CLI:  

```bash
TreSucMs -p TCGA-HNSC -p TCGA-CESC -p TCGA-LUSC -d cisplatin -d carboplatin,paclitaxel -d carboplatin -o TreSucM -c 40 -e metilene -t 5 -t 10 -C 8 -C 5 -e DESeq2
```
- content of help page:  

```
Usage: TreSucMs [OPTIONS]  
  
  "TreSucMs" a tool to choose, harvest and analyse expression and methylation  
  data of the TCGA-projects for revealing Biomarkers which indicate treatment  
  success predictions.  
  
  Calling the pipeline without any argument starts the interactive mode to  
  help setting all needed parameters for the analysis.  
  
Options:  
  -o, --out_path TEXT    path to save the result files  [default:  
                         /homes/biertruck/gabor/TreSucMs]  
  -p, --project TEXT     TCGA project(s) to be applied. Any TCGA project can  
                         be chosen, like: -p TCGA-CESC -p TCGA-HNSC ...  
  -d, --drugs TEXT       drug(s), like: -d drug1 -d drug2 or  
                         drugcombination(s), like: -d drug1,drug2  
  -c, --cores INTEGER    number of cores provided to snakemake  [default: 1]  
  -C, --cutoff FLOAT     Cut-off parameter  [default: 0]  
  -t, --threshold FLOAT  threshold parameter  [default: 0]  
  -e, --execute TEXT     choose which pipeline shall be executed  [default:  
                         DESeq2, metilene]  
  -N, --dryrun           snakemake dryrun  
  -D, --download         if set, just download raw and meta data for given  
                         projects and analysis types, revise them, link them,  
                         but do not run any analysis  
  -v, --version          printing out version information: Version 1.0  
  --help                 Show this message and exit.  
```
