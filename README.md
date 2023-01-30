Usage: python main_metilene.py [OPTIONS]

  "metilene_pipeline" a tool to choose, harvest and analyse methylation data of
  the TCGA-projects with help of the metilene package.

  build and activate the provided conda env:

      $ conda env create -f metilene_env.yaml

      $ conda activate metilene_pipeline

  call the script without any options to enter the interactive mode and set
  each option step by step:

      $ python main_metilene.py

  print help page:

      $ python main_metilene.py --help

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

