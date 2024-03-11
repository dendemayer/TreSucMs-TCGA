# testversion with automatic api download:

## TODO overview on start of report:
* give an overview of the data invoked: in table form, like:

Project   |drugs      |female | male
----------|-----------|-------|-----
TCGA-CESC |cisplatin  |4      |1 
TCGA-CESC |carboplatin|5      |1
TCGA-HNSC |cisplatin  |3      |5 
TCGA-HNSC |carboplatin|7      |9

commit 578d41f05ab9d97f2f2968f7c3f3caf1dc20bc2f (HEAD -> master, origin/master)
Author: gabor <gabor@bioinf.uni-leipzig.de>
Date:   Wed Jun 2 21:49:53 2021 +0200

    new deseq_env.yaml without prefix line, with mygene library
    
....

## DONE include the walk all drug frequencies into report, put the at the end of the report 

> Author: gabor <gabor@bioinf.uni-leipzig.de>
> Date:   Tue Apr 27 14:03:35 2021 +0200
> 
>     included automated snakemake workflow creation
> 
> commit 63475484c173fefee7d5cbdf1c7ac64f93493879
> Author: gabor <gabor@bioinf.uni-leipzig.de>
> Date:   Sat Apr 17 10:11:00 2021 +0200
> 
>     before changing snake again
> 
> commit 00ce45df8d59b702856cb44a4664e01f87233f57
> Author: gabor <gabor@bioinf.uni-leipzig.de>
> Date:   Mon Apr 12 12:28:11 2021 +0200
> 
>     added REPORT.pdf



## do it like so:
- since we also include combinations like drug1,drug2, do not join the drug
  arrays anymore, change that in every function where it is performed
- if no projects are applied via terminal call, ask them interactively, as well
  as ask the projects if they are not specified 
- the script dir is default set to cwd, no specifying is needed anymore

- the complement stuff is in future just the complete oposite of the selected
  DF
>commit 4320a20f80ec76d80d04d5a5e75daf918376423d
>commit 4320a20f80ec76d80d04d5a5e75daf918376423d
>Author: gabor <gabor@bioinf.uni-leipzig.de>
>Date:   Wed Feb 24 04:05:53 2021 -0500
>
>    version with background included in DESeq2 approach, lifeline with treshold finnishe, everythings working so far
>
__TODO__:
roll back the integration of the background directly in the DESeq2 approach,
instead run those background in an extra run in dir:
OUTPUT_PATH, DRUGS, 'DESeq2_out', project_title + 'background'

__DONE__:
pca plots back to pdf, svg is not working there

__DONE__:
change the DESeq2 plot output from pdf to svg to dynamically create md
reportfiles and infoke some plots and then convert to single html file

> commit 47b0f2d30b278c447ad8d20d7f143c336cdcb180
> Author: gabor <gabor@bioinf.uni-leipzig.de>
> Date:   Mon Feb 22 05:27:10 2021 -0500

    finnished new directory management for lifeline plots

__DONE__:  
how to handle save locations of the lifelineplots, until now, only multiproject
lifelines are created add single project lifeline plots

__DONE__  
rename the ALL_PROJECTS_summary_dead_alive_NOT_reduced.tsv in the drug dirs in the correct
aggregation of the summarized projects, for example: TCGA-CESC_TCGA-HNSC_summary_dead_alive_N...

> commit 6911d747007c11037848f5f98ff92fb904bcdc3c
> Author: gabor <gabor@bioinf.uni-leipzig.de>
> Date:   Tue Feb 16 08:48:03 2021 -0500
> 
> changing heatmap order and dirs of the multiproject summary tables

__DONE__   
>Error in DESeqDataSet(se, design = design, ignoreRank) : 
  >design contains one or more variables with all samples having the same value,
  >remove these variables from the design
  
__DONE__:  
R deseq with applying the background  

- the creation of the drug filtered summary tables is performed in fct 7  
- first provide table ALL_PROJECTS_summary_dead_alive_reduced_INFO.tsv  
with additional cisplatin row (and within the factors TRUE and FALSE)  
->   edit fct 8 s.t. the cisplatin true or false is also handed over
> fct 8 in main
>def provide_DESeq2_table_all_PROJECTS(PROJECT_DRUG_UUID, OUTPUT_PATH, DRUGS,
>                                      SCRIPT_DIR):

> commit 5394be09497818be527061d5f75a9a6062dcdfa1 (HEAD -> master, origin/master)
> Author: gabor <gabor@bioinf.uni-leipzig.de>
> Date:   Thu Feb 4 11:34:52 2021 -0500
> 
>     implementation of lifelines finnished with treshold and all 3 counttables

__TODO__:  
the raw counts ALL_PROJECTS_summary_dead_alive_reduced.tsv have the col numbers
in the table when created, delete that when the table is created and take care
of the following steps this line has to be deleted for proper working of the
lifelineplots!!


> commit 914ce275bd1e9372c2e460110d7279af27f6cf28 (HEAD -> master, origin/master)
> Author: gabor <gabor@bioinf.uni-leipzig.de>
> Date:   Mon Feb 1 15:55:00 2021 +0100
>
>     fixed the sphinx doc stuff(run make html in /doc)

the thresholds are going to be implemented anyway, create new dirs for every
percentage given, also add tables to every kaplan meier plot created



> commit 7a4fbbd2c66d4ccb30d77e44e7a2300eeecff109 (HEAD -> master, origin/master)
> Author: gabor <gabor@bioinf.uni-leipzig.de>
> Date:   Wed Jan 27 14:35:18 2021 +0100
>
>       added dir Lifeline_plots, starting to implement treshold for kaplan meier separation

treshhold is not very revealing.... 
in this particular approach (CESC with HNSC) it might be useful to exclude the
males from the HNSC trail, although gender can be added in DESeq2 (is this automatically
considered as batch effect?) 
ddsMF <- DESeqDataSetFromMatrix(countData = count_data, colData = vital_cancer, design= ~ cancer + gender + condition )
it might be also usefull in the metilene approach...
