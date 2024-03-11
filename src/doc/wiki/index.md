# DESeq2 pipeline  

## variables of the result table:

The first column, __baseMean__, is a just the average of the normalized count
values, divided by the size factors, taken over all samples in the
DESeqDataSet. The remaining four columns refer to a specific contrast, namely
the comparison of the trt (alive) level over the untrt (dead) level for the
factor variable dex. We will find out below how to obtain other contrasts.

The column __log2FoldChange__ is the __effect size estimate__. It tells us how much the
gene’s expression seems to have changed due to treatment with dexamethasone in
comparison to untreated samples. This value is reported on a logarithmic scale
to base 2: for example, a log2 fold change of 1.5 means that the gene’s
expression is increased by a multiplicative factor of 21.5≈2.8221.5≈2.82.

Of course, this estimate has an uncertainty associated with it, which is
available in the column __lfcSE__, the standard error estimate for the log2 fold
change estimate. We can also express the uncertainty of a particular effect
size estimate as the result of a statistical test. The purpose of a test for
differential expression is to test whether the data provides sufficient
evidence to conclude that this value is really different from zero. DESeq2
performs for each gene a hypothesis test to see whether evidence is sufficient
to decide against the null hypothesis that there is zero effect of the
treatment on the gene and that the observed difference between treatment and
control was merely caused by experimental variability (i.e., the type of
variability that you can expect between different samples in the same treatment
group). As usual in statistics, the result of this test is reported as a p
value, and it is found in the column pvalue. Remember that a p value indicates
the probability that a fold change as strong as the observed one, or even
stronger, would be seen under the situation described by the null hypothesis.




TODO -> change the DESeq2 script to handle the htseq unnormalized type!
DONE -> the htseq functionality serves just the reading in of the tables, this is already accomplished 

TODO -> ->  modules/create_report.py l 740 the pdflatex part is not working with conda installation...


## table flow see [meta_table_schema](../../../metilene_pipeline/doc/wiki/meta_table_schema.md):

- overview of headers of al meta tables of CESC projects from TCGA
- link primary tumor samples with exact filenames:
- __?? still up to date??__ the direct info is available in
  table:nationwidechildrens.org_biospecimen_sample_cesc col:sample_type see
[meta_table_schema](../../../metilene_pipeline/doc/wiki/meta_table_schema.md)
  

! LaTeX Error: File `article.cls' not found.
! LaTeX Error: File `booktabs.sty' not found.
! LaTeX Error: File `l3backend-pdftex.def' not found.
pdfTeX warning: pdflatex (file pdftex.map): cannot open font map file
] (./template_latex_all_invoked_cases.aux)
kpathsea: Running mktexpk --mfmode / --bdpi 600 --mag 1+0/600 --dpi 600 cmr10
mktexpk: Cannot find mktex.opt; check your installation.
kpathsea: Appending font creation commands to missfont.log.
 )
!pdfTeX error: pdflatex (file cmr10): Font cmr10 at 600 not found
 ==> Fatal error occurred, no output PDF file produced!
if  /usr/bin/pdflatex invoked, it runs through..


download functions 1 to 4 are working properly like in metilene pipeline 
create the summary table based on whats in the new
meta_info_druglist_merged_drugs_combined.tsv file

- [doc_progress_deseq2](doc_progress_deseq2)
- [conda](conda) install deseq pipeline
- [DESeq2](DESeq2) multifactorial run stuff from the vignette (https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#multi-factor-designs)
- git [Readme](../test_git_doc/9_create_matrix_with_doc/Readme) shall be the same like default help page
* [environment](environment) stuff

# cutoff problem when cutoffs are added with -A and then snakerun

# preserve database for snakerun:

* the meta_info.dat out of fct 1 must be preserved:  

associated_entities.0.case_id	cases.0.demographic.gender	cases.0.demographic.vital_status	cases.0.demographic.year_of_birth	cases.0.demographic.year_of_death	cases.0.diagnoses.0.age_at_diagnosis	cases.0.diagnoses.0.days_to_last_follow_up	id  
3d190f07-2f79-4ee6-9f63-0b563dc05c74	female	Alive	1969		15365	1221	0bf72111-cbfe-4d06-9db1-d04353eff0b7  
929505f7-0270-475e-a9d6-00bf38dc78ba	female	Alive	1979		12775	287	c90c4e59-6c5a-4599-aaf4-ace45a97a27d  
55add854-e63f-4335-bea9-23c3eb1e766a	female	Alive	1970		15918	90	a234dc97-644b-4ebc-8375-861d9b82516a  

also the file uuid's which are actually downloaded must be perserved, they are
saved in the  

    print("\nfile_uuid_list:\n", file_uuid_list)  
of fct. 3

file uuids are all the same in:  
  * MANIFEST.txt id col
  * TCGA-CESC_data_files
  * meta_info.dat id col

nevertheless, the meta_info.dat file must be preserved because of possible
changes in the meta data (survivaltime and such..)

-> this goes into the Snakes dir in
SCRIPT_PATH/meta_infos/TCGA-CESC_TCGA-HNSC_TCGA-LUSC_carboplatin_carboplatin,paclitaxel_cisplatin/
TCGA-CESC_meta_info.dat
TCGA-HNSC_meta_info.dat
TCGA-LUSC_meta_info.dat
