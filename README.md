Usage: TreMSuc [OPTIONS]

```bash
TreMSuc -p TCGA-HNSC -p TCGA-CESC -p TCGA-LUSC -d cisplatin -d carboplatin,paclitaxel -d carboplatin -o /scr/palinca/gabor/TCGA-pipeline -c 40 -e metilene -t 5 -t 10 -C
 8 -C 5 -e DESeq2
```

  "TreMSuc" a tool to choose, harvest and analyse  data of
  the TCGA-projects with help of the metilene and DESeq2 package to find
  biomarkers which can have impact on treatment success.

  install the package:
    $ git clone git@github.com:dendemayer/tcga_piplines.git
    $ cd tcga_piplines
    $ pip install .

  call the script without any options to enter the interactive mode and set
  each option step by step:

      $ TreMSuc

  print help page:

      $ TreMSuc --help

Options:  

  -D, --download_data

      perform download and merging steps, without the metilene analysis  

  -A, --analyse_data

      perform the metilene analysis  (the download step must be completed in
      prior for your chosen projects)
    
  -o, --out_path TEXT
  
      The Path, where the results are saved to  
      [default: ./metilene_data]  

  -p, --project TEXT
  
      Project, that shall be applied, choose between  
      TCGA- CESC, TCGA-HNSC, TCGA-LUSC, TCGA-ESCA, or combinations  
      out of them like -p project1 -p project2  
  -d, --drugs TEXT
  
      drug(s), like: -d drug1 -d drug2 or drugcombination(s),  
      like: -d drug1,drug2  


  -f, --function INTEGER
  
      running functions seperately lookup the for detailed  
      description of each function, lookup the documentation  
      (link at end of page)

  -v, --version
  
      printing out version information: Version 1.0
  --help
  
      Show this message and exit.  
  
For detailed explanations please refer to the [documentation](doc/_build/html/index.html)

