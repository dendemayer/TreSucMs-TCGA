#!/usr/bin/env Rscript
suppressPackageStartupMessages({
require(sjmisc)  # for the str_contains() fct
require(DESeq2)
require(dplyr)
require(apeglm)
require(pcaExplorer)
require(ggplot2)
require(pheatmap)
require(RColorBrewer)
})
args <- commandArgs(trailingOnly = TRUE)

cat("arguments:\n")
j=1
for (i in args){
    cat(j, ": ", i, "\n")
    j <- j + 1
}


### test files:



#count_inname <- file.path('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-HNSC/DESeq2/DESeq2_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_2/summary_for_DESeq2.tsv')
#info_inname  <- file.path('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-HNSC/DESeq2/DESeq2_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_2/summary_for_DESeq2_INFO.tsv')
#OUTPUT_PATH  <- file.path('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_2')
#PROJECTS <- 'TCGA-HNSC'

#count_inname <- file.path('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-LUSC/DESeq2/DESeq2_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_2/summary_for_DESeq2.tsv')
#info_inname  <- file.path('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-LUSC/DESeq2/DESeq2_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_2/summary_for_DESeq2_INFO.tsv')
#OUTPUT_PATH  <- file.path('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_2')
#PROJECTS <- 'TCGA-LUSC'

#count_inname <- file.path('/scr/dings/PEVO/NEW_downloads_3/DEseq_env_test/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/carboplatin_carboplatin,paclitaxel_cisplatin/DRUG_combi_summary_dead_alive_short.tsv')
#info_inname <- file.path('/scr/dings/PEVO/NEW_downloads_3/DEseq_env_test/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/carboplatin_carboplatin,paclitaxel_cisplatin/DRUG_combi_summary_dead_alive_INFO.tsv')
#OUTPUT_PATH  <- file.path('/scr/dings/PEVO/NEW_downloads_3/DEseq_env_test/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/carboplatin_carboplatin,paclitaxel_cisplatin/DESeq2_out_DRUG_combi')
#PROJECTS <- 'TCGA-CESC_TCGA-HNSC_TCGA_LUSC'


#count_inname <- file.path('/scr/dings/PEVO/NEW_downloads_3/DEseq_31_4/TCGA-CESC/carboplatin_carboplatin,paclitaxel_cisplatin/DRUG_combi_female_summary_dead_alive.tsv')
#info_inname <- file.path('/scr/dings/PEVO/NEW_downloads_3/DEseq_31_4/TCGA-CESC/carboplatin_carboplatin,paclitaxel_cisplatin/DRUG_combi_female_summary_dead_alive_INFO.tsv')
#OUTPUT_PATH  <- file.path('/scr/dings/PEVO/NEW_downloads_3/DEseq_31_4/TCGA-CESC/carboplatin_carboplatin,paclitaxel_cisplatin/DESeq2_out_DRUG_combi_female')
#PROJECTS <- 'TCGA-CESC'

#count_inname <- file.path('/scr/dings/PEVO/NEW_downloads_3/DEseq_31_6/TCGA-CESC/carboplatin_carboplatin,paclitaxel_cisplatin/DRUG_combi_complement_summary_dead_alive.tsv')
#info_inname <- file.path('/scr/dings/PEVO/NEW_downloads_3/DEseq_31_6/TCGA-CESC/carboplatin_carboplatin,paclitaxel_cisplatin/DRUG_combi_complement_summary_dead_alive_INFO.tsv')
#OUTPUT_PATH  <- file.path('/scr/dings/PEVO/NEW_downloads_3/DEseq_31_6/TCGA-CESC/carboplatin_carboplatin,paclitaxel_cisplatin/DESeq2_out_DRUG_combi_complement')
#PROJECTS <- 'TCGA-CESC'


#count_inname <- file.path('/scr/dings/PEVO/NEW_downloads_3/DEseq_31_6/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/carboplatin_carboplatin,paclitaxel_cisplatin/DRUG_combi_complement_summary_dead_alive.tsv')
#info_inname <- file.path('/scr/dings/PEVO/NEW_downloads_3/DEseq_31_6/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/carboplatin_carboplatin,paclitaxel_cisplatin/DRUG_combi_complement_summary_dead_alive_INFO.tsv')
#OUTPUT_PATH  <- file.path('/scr/dings/PEVO/NEW_downloads_3/DEseq_31_6/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/carboplatin_carboplatin,paclitaxel_cisplatin/DESeq2_out_DRUG_combi_complement')
#PROJECTS <- 'TCGA-CESC_TCGA-HNSC_TCGA-LUSC'

#####

#count_data <- as.matrix(read.table(count_inname,header=T,row.names=1))
## the default sep value for read.delim is '\t' and handles better the spaces
## which occur in some drug name cols
#info_data <- as.matrix(read.delim(info_inname,row.names=1))
##colnames(count_data) = info_data[5,]

#temp_name = paste0(count_inname, '_test.tsv')
#write.table(info_data, file=temp_name, quote=FALSE, sep='\t', col.names = NA)
##### END TEST

# ### ommit this block if you want to test the script with some test data>>>
count_inname <- args[1]  # count_table
info_inname <- args[2]  # info_table
OUTPUT_PATH  <- args[3]  # deseq <- output
PROJECTS <- args[4]
# ### ommit this block if you want to test the script with some test data<<<

count_data <- as.matrix(read.delim(count_inname,header=T,row.names=1))
# the default sep value for read.delim is '\t' and handles better the spaces
# which occur in some drug name cols
info_data <- as.matrix(read.delim(info_inname,row.names=1))
# put the right case_id to the colnames
colnames(count_data) <- info_data[5,]

vital_status <- as.vector(info_data[1,])
vital_status <- data.frame(vital_status)
colnames(vital_status)<-c("vital_state")
vital_status$vital_state <- as.factor(vital_status$vital_state)
#str(vital_status)
#head(vital_status)
#dim(unique(vital_status))[1]
tmp_cancer_column  <- as.vector(info_data[3,])
#info_data
#tmp_cancer_column
tmp_cancer_column  <- data.frame(tmp_cancer_column)
colnames(tmp_cancer_column) <- c("cancer")
tmp_cancer_column$cancer <- as.factor(tmp_cancer_column$cancer)
#head(tmp_cancer_column)
#dim(unique(tmp_cancer_column))[1]

tmp_gender_column  <- as.vector(info_data[2,])
tmp_gender_column  <- data.frame(tmp_gender_column)
colnames(tmp_gender_column) <- c("gender")
tmp_gender_column$gender <- as.factor(tmp_gender_column$gender)
#head(tmp_gender_column)
#dim(unique(tmp_gender_column))[1]

# including the drugname
drugnames_column <- as.vector(info_data[4,])
drugnames_column  <- data.frame(drugnames_column)
colnames(drugnames_column) <-c("drugnames")
drugnames_column$drugnames <- as.factor(drugnames_column$drugnames)
#head(drugnames_column)
#dim(unique(drugnames_column))[1]
# We create a copy of the DESeqDataSet,
# so that we can rerun the analysis using a multi-factor design.

#TODO Error in DESeqDataSet(se, design = design, ignoreRank) : design contains
#one or more variables with all samples having the same value, remove these
#variables from the design vital_cancer = cbind(vital_status,
#tmp_gender_column, tmp_cancer_column, drugnames_column) Check every col if there are
#different values, otherwise, the gender (CESC all female) and projects (in
#single runs) can be omitted  dim(unique(vital_status))[1] -> not allowed to be
#1
#vital_cancer = cbind(tmp_gender_column, tmp_cancer_column, vital_status)
vital_cancer <- cbind(tmp_gender_column, tmp_cancer_column, drugnames_column,
                     vital_status)
#cat('factors available:\n')
#head(vital_cancer)
get_col_flags  <- function()  # get colnames of ambigous cols
{
    #print('available cols of count DataFrame: ')
    #print(colnames(vital_cancer))
    col_flags <- vector()
    for (col_name in colnames(vital_cancer))
    {
        if (dim(unique(vital_cancer[col_name]))[1] > 1)
            {   col_flags <- c(col_flags, col_name)
                #print(paste0('adding following to DESeq2 design:'))
                #print(paste0(col_name))
            }
    }
    return(col_flags)
}

# get the right cols for the DESeq design, of all cols available choose the
# ones wich have ambiguos values(unique value count higher 1)
col_flags = get_col_flags()
col_flags
# possible combinations: 
# gender, cancer, treated, vital_state
# gender, cancer, vital_state
# gender, treated, vital_state
# cancer, treated, vital_state
# gender, vital_state
# cancer, vital_state
# treated, vital_state
# vital_state

# EXCLUGING treated: possible combinations: 
# gender, cancer,  vital_state
# gender, vital_state
# cancer, vital_state
# vital_state
# also check length of vector, otherwise test for example c(gender vital_state) would be
# truth in actually c(gender cancer vital_state)

# including drugnames:
# possible combinations: 
#   gender, cancer, drugnames, vital_state
#   gender, cancer, vital_state
#   gender, drugnames, vital_state
#   cancer, drugnames, vital_state
#   gender, vital_state
#   cancer, vital_state
#   drugnames, vital_state
#   vital_state
get_ddsMF  <- function(col_flags){ # (create the deseq dataset depending on col_flags)
if (str_contains(col_flags,c("gender", "cancer", "drugnames", "vital_state"), logic= "and"))
    {
        cat('\nDESeq2 Multifactor design with: \n \t gender, cancer, drugnames, vital_state\n\n')
    ddsMF <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = vital_cancer,
                                design= ~ gender + cancer + drugnames + vital_state)
    } else if (str_contains(col_flags,c("gender", "cancer", "vital_state"), logic= "and"))
    {
        cat('\nDESeq2 Multifactor design with: \n \t gender, cancer, vital_state\n\n')
    ddsMF <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = vital_cancer,
                                design= ~ gender + cancer + vital_state)
    } else if (length(col_flags)==2 & str_contains(col_flags,c("gender", "drugnames", "vital_state"), logic= "and")) {
        cat('\nDESeq2 Multifactor design with: \n \t gender, drugnames, vital_state\n\n')
        ddsMF <- DESeqDataSetFromMatrix(countData = count_data,
                                    colData = vital_cancer,
                                    design= ~ gender + drugnames + vital_state ) # colnames of
    } else if (length(col_flags)==2 & str_contains(col_flags,c("cancer", "drugnames", "vital_state"), logic= "and")) {
        cat('\nDESeq2 Multifactor design with: \n \t cancer, drugnames, vital_state\n\n')
        ddsMF <- DESeqDataSetFromMatrix(countData = count_data,
                                    colData = vital_cancer,
                                    design= ~ cancer + drugnames + vital_state ) # colnames of
    } else if (length(col_flags)==2 & str_contains(col_flags,c("gender", "vital_state"), logic= "and")) {
        cat('\nDESeq2 Multifactor design with: \n \t gender, vital_state\n\n')
        ddsMF <- DESeqDataSetFromMatrix(countData = count_data,
                                    colData = vital_cancer,
                                    design= ~ gender + vital_state )
    } else if (length(col_flags)==2 & str_contains(col_flags,c("cancer", "vital_state"), logic= "and")) {
        cat('\nDESeq2 Multifactor design with: \n \t cancer, vital_state\n\n')
        ddsMF <- DESeqDataSetFromMatrix(countData = count_data,
                                    colData = vital_cancer,
                                    design= ~ cancer + vital_state )
    } else if (length(col_flags)==2 & str_contains(col_flags,c("drugnames", "vital_state"), logic= "and")) {
        cat('\nDESeq2 Multifactor design with: \n \tdrugnames, vital_state\n\n')
        ddsMF <- DESeqDataSetFromMatrix(countData = count_data,
                                        colData = vital_cancer,
                                        design= ~ drugnames + vital_state )
    } else {
        cat('\nDESeq2 Singlefactor design with: \n \t vital_state\n\n')
        ddsMF <- DESeqDataSetFromMatrix(countData = count_data,
                                        colData = vital_cancer,
                                        design= ~ vital_state ) # colnames of
    }
}


ddsMF <- get_ddsMF(col_flags)
#table vital_cancer
#head(ddsMF) 

colData(ddsMF)

# We can account for the different types of sequencing,
# and get a clearer picture of the differences attributable to the treatment.
# As condition is the variable of interest, we put it at the end of the formula.
# Thus the results function will by default pull the condition
# results unless contrast or name arguments are specified.

#design(ddsMF) <- formula(~ cancer + condition)

# Pre-filtering.
keep = rowSums(counts(ddsMF)) >= 10
#ddsMF_bak = ddsMF
ddsMF <- ddsMF[keep,]
#### -> this step leads to the pheatmap error:
#> pheatmap(count_data[select_decrease_log2Fold,], cluster_rows=TRUE, show_rownames=TRUE,show_colnames = FALSE, cluster_cols=FALSE, annotation_col = df, scale = "row", main=PROJECTS, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList) Fehler in hclust(d, method = method) :
  #NA/NaN/Inf in externem Funktionsaufruf (arg 10)
#Error in hclust(d, method = method) :
#  NA/NaN/Inf in foreign function call (arg 10)
#Calls: pheatmap -> cluster_mat -> hclust

### begin change
# with introducing compelement analysis this error occurs:
    #fitting model and testing 2108 rows did not converge in beta, labelled in
    #mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest
# suggested solution from: https://support.bioconductor.org/p/65091/:
    #dds <- estimateSizeFactors(dds)
    #nc <- counts(dds, normalized=TRUE)
    #filter <- rowSums(nc >= 10) >= 2
    #dds <- dds[filter,]
#-> in contrast to the problematic pre filtering from vignette here
#estimateSizeFactors and adidtional >= 2 for the filter is suggested

#ddsMF <- estimateSizeFactors(ddsMF)
#nc  <- counts(ddsMF, normalized=TRUE)
#filter  <- rowSums(nc >= 10) >= 4
#ddsMF <- ddsMF[filter,]

##Some more options, try increasing the 2 sample requirement to 3 or 4.
##To deal with any remaining rows, you can either omit them from the results step: ddsClean <- dds[which(mcols(dds)$betaConv),], because these are typically genes with very small counts and little power, or you can increase the maximum iterations with the following code.  __Instead__ of running DESeq(), run:
##mcols(DESeq_ddsMF)
### if thats not sufficient use
#ddsMF <- estimateSizeFactors(ddsMF)
#ddsMF <- estimateDispersions(ddsMF)
##ddsMF <- nbinomWaldTest(ddsMF, maxit=500)
## 
#cat('new maxit with 1000:\n')
#DESeq_ddsMF <- nbinomWaldTest(ddsMF, maxit = 1000)
##> DESeq_ddsMF <- nbinomWaldTest(ddsMF, maxit=500) 1794 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest
#DESeq_ddsMF <- DESeq_ddsMF[which(mcols(DESeq_ddsMF)$betaConv), ]

### end change

# Then we can run DESeq: this code is substituted by the above one...

DESeq_ddsMF <- DESeq(ddsMF)
#assay(DESeq_ddsMF)
#rowData(DESeq_ddsMF)
#colData(DESeq_ddsMF)

### shall we include contrast thing?:
#dds <- DESeq(dds)
#res <- results(dds)
#res <- results(dds, name="condition_treated_vs_untreated")
#res <- results(dds, contrast=c("condition","treated","untreated"))
### NO because we are interested in the vital_state_dead_vs_alive condition
### anyway and this is considered in the design:
#design(ddsMF) <- formula(~ cancer + condition)

#resultsNames(DESeq_ddsMF)
#[1] "Intercept"
#[2] "gender_male_vs_female"
#[3] "cancer_TCGA.HNSC_vs_TCGA.CESC"
#[4] "cancer_TCGA.LUSC_vs_TCGA.CESC"
#[5] "drugnames_carboplatin.paclitaxel_vs_carboplatin"
#[6] "drugnames_cisplatin_vs_carboplatin"
#[7] "vital_state_dead_vs_alive"

results_ddsMF <- results(DESeq_ddsMF)

# some rows lack the padj value, omit NAs:
results_ddsMF <- na.omit(results_ddsMF)
# also rows with padj values greater than 0.05 can be omitted
# that can be done in the -f 8 and -f10 step:
# results_DECREASE_10_lgfch_ENSG_and_Gene_info.tsv on that basis the tables and
# lifeline plots are created
#TODO: this also triggers the phetamap error...
#results_ddsMF <- results_ddsMF[results_ddsMF$padj < 0.05, ]

# using the LFC shrinkage as folowing results:
#results_ddsMF <- lfcShrink(DESeq_ddsMF, coef="vital_state_dead_vs_alive", type="apeglm")
#results_ddsMF <- lfcShrink(DESeq_ddsMF, coef="vital_state_dead_vs_alive", type="normal")

#head(results_ddsMF)
results_ddsMF_sort_pv <- results_ddsMF[order(results_ddsMF$padj),]
 ####### writing result table ordered pvalue
results_outname = file.path(OUTPUT_PATH, paste('DESeq2_results.tsv', sep=''))
cat('writing DESeq2 MF results in:\t', results_outname, '\n')
write.table(as.data.frame(results_ddsMF_sort_pv), file = results_outname, quote = FALSE, sep='\t', col.names=NA)

 ####### saving the summary of the results:
results_summary = file.path(OUTPUT_PATH, paste('DESeq2_results_summary.tsv', sep=''))
cat('writing DESeq2 MF results SUMMARY in:\t', results_summary, '\n')
sink(file = results_summary)
summary(results_ddsMF)
sink(file = NULL)


### START Log fold change shrinkage for visualization and ranking 

#Shrinkage of effect size (LFC estimates) is useful for visualization and
#ranking of genes. To shrink the LFC, we pass the dds object to the function
#lfcShrink. Below we specify to use the apeglm method for effect size shrinkage
#(Zhu, Ibrahim, and Love 2018), which improves on the previous estimator.

#We provide the dds object and the name or number of the coefficient we want to
#shrink, where the number refers to the order of the coefficient as it appears
#in resultsNames(dds).
#resLFC <- lfcShrink(DESeq_ddsMF, coef="vital_state_dead_vs_alive", type="apeglm")
#resLFC
### END Log fold change shrinkage for visualization and ranking 

### plot counts:

#[7] ""
#[7] "vital_state_dead_vs_alive"
#plotCounts(ddsMF, gene=which.min(results_ddsMF$padj), intgroup="vital_statee")

d <- plotCounts(ddsMF, gene=which.min(results_ddsMF$padj), intgroup="vital_state", returnData=TRUE)
#library("ggplot2")
pdf(file.path(OUTPUT_PATH, paste('DESeq2_gene_counts.pdf', sep='')))
ggplot(d, aes(x=vital_state, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) + labs(title = paste("gene counts for ", PROJECTS))
dev.off()


# # ####### writing the results 'dot'plot.pdf
#### 
### this is now the following:
#It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
resLFC <- lfcShrink(DESeq_ddsMF, coef="vital_state_dead_vs_alive", type="apeglm")
results_outname_pdf = file.path(OUTPUT_PATH, paste('DESeq2_results.pdf', sep=''))
cat("save plot in:\t", results_outname_pdf, '\n')
pdf(results_outname_pdf)
#plotMA(results_ddsMF, ylim=c(-3,3), main=PROJECTS)
main_tile = paste(PROJECTS, ', LFC shrinkage with apeglm method')
plotMA(resLFC, ylim=c(-3,3), main=main_tile)
dev.off()

# ###############writing normalized count table (original count table for DESeq
# input is the input summary)
table_counts_normalized_MF <- counts(DESeq_ddsMF, normalized=TRUE)
table_counts_raw_MF <- counts(DESeq_ddsMF)
#head(table_counts_normalized_MF)
normalized_outname = file.path(OUTPUT_PATH, paste('DESeq2_normalized_counts.tsv', sep=''))
cat("writing normalized MF counts in:\t", normalized_outname, "\n")
write.table(table_counts_normalized_MF, file = normalized_outname, quote = FALSE, sep = "\t", col.names=FALSE)


# ###############writing normtransform counts table:

table_counts_norm_transform <- assay(normTransform(DESeq_ddsMF))
normalized_outname = file.path(OUTPUT_PATH, paste('DESeq2_norm_transform_counts.tsv', sep=''))
cat("writing norm_transform MF counts in:\t", normalized_outname, "\n")
write.table(table_counts_norm_transform, file = normalized_outname, quote = FALSE, sep = "\t", col.names=FALSE)


# ###########HEATMAP with normalised, ntd, vsd values and scaled rows:
# 
# #############################################
#library("pheatmap")
ntd <- normTransform(DESeq_ddsMF)
# (to make use in heatmaps use assay(ntd))
vsd <- vst(DESeq_ddsMF, blind=FALSE) # variance stabilizing transformation
#rld <- rlog(DESeq_ddsMF, blind=FALSE) # regularized log transformation
df <- as.data.frame(colData(DESeq_ddsMF)[,col_flags]) # use as annotation col
#head(df)

# needed s.t. the color scale is not based on negative value, instead starts
# at zero
#library(RColorBrewer)

# ###### sort plot log2fold wise for normalized score:

# ################################# start of INCREASE norm logfch heatmaps: -> usage of the **table_counts_normalized_MF** counts

select_increase_log2Fold  <- with(results_ddsMF,  order(log2FoldChange))
results_ddsMF_temp <- results_ddsMF[select_increase_log2Fold,]
results_ddsMF_temp <- results_ddsMF_temp[results_ddsMF_temp$padj < 0.05,]
if (dim(results_ddsMF_temp)[1] < 60) {
    select_size <- dim(results_ddsMF_temp)[1]
} else {
    select_size  <- 60
}

row_names_inc <- row.names(results_ddsMF_temp[1:select_size,])
    
# the breaksList should be updated for every count matrix
# limit the range of breaks to the selection of ENSG which shall be plotted:
#old: [1]       0.0  941414.3 1882828.6 2824243.0 3765657.3 4707071.6 5648485.9 
breaksList <- seq(min(table_counts_normalized_MF[row_names_inc, ]), max(table_counts_normalized_MF[row_names_inc, ]), by = max(table_counts_normalized_MF[row_names_inc, ]) / 10000)
# new: [1]      0.0  90464.2 180928.4 271392.6 361856.8 452321.0 542785.2
# with those rownames, access the count table:

#cat('results_ddsMF[select_increase_log2Fold,] normalized')
#head(results_ddsMF[select_increase_log2Fold,])
heatmap_outname_pdf=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fINCREASE_norm.pdf")
cat('writing DESeq2 normTransform in:\t', heatmap_outname_pdf, '\n')

pdf(heatmap_outname_pdf)
pheatmap(table_counts_normalized_MF[row_names_inc,], show_rownames=TRUE,show_colnames = FALSE, cluster_cols=FALSE, cluster_rows=FALSE, annotation_col = df, main=paste(PROJECTS, "\nnormalized counts\nincreasing log2foldChange", sep=" ") , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)
#pheatmap(table_counts_normalized_MF[select_increase_log2Fold,], show_rownames=TRUE,show_colnames = FALSE, cluster_cols=FALSE, annotation_col = df, main=PROJECTS)
dev.off()
# #### writing the belonging count table and result table of selected size:
heatmap_outname_tsv=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fINCREASE_norm_counts.tsv")
write.table(as.data.frame(table_counts_normalized_MF[row_names_inc,]), file = heatmap_outname_tsv, quote = FALSE, sep = "\t", col.names=TRUE)

heatmap_outname_tsv=file.path(OUTPUT_PATH, "DESeq2_log2fINCREASE_result.tsv")
write.table(as.data.frame(results_ddsMF[row_names_inc,]), file = heatmap_outname_tsv, quote = FALSE, sep = "\t", col.names=TRUE)
# ########################################### end of INCREASE norm logfch heatmaps

# ###########################################start of DECREASE norm
#select_decrease_log2Fold  <- order(results_ddsMF$log2FoldChange, decreasing=TRUE)[1:select_size]
#select_decrease_log2Fold  <- with(results_ddsMF, order(-log2FoldChange,padj, pvalue))[1:select_size]
select_decrease_log2Fold  <- with(results_ddsMF,  order(-log2FoldChange))
results_ddsMF_temp <- results_ddsMF[select_decrease_log2Fold,]
results_ddsMF_temp <- results_ddsMF_temp[results_ddsMF_temp$padj < 0.05,]
row_names_dec <- row.names(results_ddsMF_temp[1:select_size,])

breaksList <- seq(min(table_counts_normalized_MF[row_names_dec, ]), max(table_counts_normalized_MF[row_names_dec, ]), by = max(table_counts_normalized_MF[row_names_dec, ]) / 10000)
#cat('results_ddsMF[select_decrease_log2Fold,] normalized')
#head(results_ddsMF[select_decrease_log2Fold,])
heatmap_outname_pdf=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fDECREASE_norm.pdf")
cat('writing DESeq2 normTransform in:\t', heatmap_outname_pdf, '\n')
pdf(heatmap_outname_pdf)
# #The assay function is used to extract the matrix of normalized values.
pheatmap(table_counts_normalized_MF[row_names_dec,], show_rownames=TRUE,show_colnames = FALSE, cluster_cols=FALSE, cluster_rows=FALSE, annotation_col = df, main=paste(PROJECTS, "\nnormalized counts\ndecreasing log2foldChange", sep=" "), color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)
#pheatmap(table_counts_normalized_MF[row_names_dec,], cluster_rows=TRUE, show_rownames=TRUE,show_colnames = FALSE, cluster_cols=FALSE, annotation_col = df, main=PROJECTS, scale="row")
dev.off()
# #### writing the belonging count table and result table of selected size:
heatmap_outname_tsv=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fDECREASE_norm_counts.tsv")
write.table(as.data.frame(table_counts_normalized_MF[row_names_dec,]), file = heatmap_outname_tsv, quote = FALSE, sep = "\t", col.names=TRUE)
heatmap_outname_tsv=file.path(OUTPUT_PATH, "DESeq2_log2fDECREASE_result.tsv")
write.table(as.data.frame(results_ddsMF[row_names_dec,]), file = heatmap_outname_tsv, quote = FALSE, sep = "\t", col.names=TRUE)
# ###########################################end of DECREASE norm


## #############################################start of INCREASE nt:

breaksList <- seq(min(assay(ntd)[row_names_inc,]), max(assay(ntd)[row_names_inc,]), by = max(assay(ntd)[row_names_inc,]) / 10000)

heatmap_outname_tsv=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fINCREASE_nt_counts.tsv")
write.table(as.data.frame(assay(ntd)[row_names_inc,]), file = heatmap_outname_tsv, quote = FALSE, sep = "\t", col.names=TRUE)

heatmap_outname_pdf=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fINCREASE_nt.pdf")
cat('writing DESeq2 normTransform Increase in:\t', heatmap_outname_pdf, '\n')
pdf(heatmap_outname_pdf)
# #The assay function is used to extract the matrix of normalized values.
pheatmap(assay(ntd)[row_names_inc,], show_rownames=TRUE,show_colnames = FALSE, cluster_cols=FALSE, cluster_rows=FALSE, annotation_col = df, main=paste(PROJECTS, "\nnormtransformed counts\nincreasing log2foldChange", sep=" "), color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)
dev.off()

heatmap_outname_tsv=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fDECREASE_norm_counts.tsv")
write.table(as.data.frame(table_counts_normalized_MF[row_names_dec,]), file = heatmap_outname_tsv, quote = FALSE, sep = "\t", col.names=TRUE)
## #############################################end of INCREASE nt:

## ###########################################start of DECREASE nt

breaksList <- seq(min(assay(ntd)[row_names_dec,]), max(assay(ntd)[row_names_dec,]), by = max(assay(ntd)[row_names_dec,]) / 10000)
#select_decrease_log2Fold <- with(results_ddsMF, order(-log2FoldChange, padj, pvalue))[1:select_size]

#cat('results_ddsMF[select_decrease_log2Fold,] nt')
#head(results_ddsMF[select_decrease_log2Fold,])
heatmap_outname_pdf=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fDECREASE_nt.pdf")
cat('writing DESeq2 normTransform in:\t', heatmap_outname_pdf, '\n')
pdf(heatmap_outname_pdf)
# #The assay function is used to extract the matrix of normalized values.
pheatmap(assay(ntd)[row_names_dec,], show_rownames=TRUE,show_colnames = FALSE, cluster_cols=FALSE, cluster_rows=FALSE, annotation_col = df, main=paste(PROJECTS, "\nnormtransformed counts\ndecreasing log2foldChange", sep=" "), color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)
dev.off()

heatmap_outname_tsv=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fDECREASE_nt_counts.tsv")
write.table(as.data.frame(assay(ntd)[row_names_dec,]), file = heatmap_outname_tsv, quote = FALSE, sep = "\t", col.names=TRUE)
## ###########################################end of DECREASE nt

##############################################start of INCREASE raw counts

breaksList <- seq(min(count_data[row_names_inc,]), max(count_data[row_names_inc,]), by = max(count_data[row_names_inc,]) / 10000)

heatmap_raw_counts_pdf=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fINCREASE_raw_counts.pdf")
pdf(heatmap_raw_counts_pdf)
pheatmap(count_data[row_names_inc,],  show_rownames=TRUE,show_colnames = FALSE, cluster_cols=FALSE, cluster_rows=FALSE, annotation_col = df, main = paste(PROJECTS, "\nraw counts\nincreasing log2foldChange", sep=" "), color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)
dev.off()
heatmap_outname_tsv=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fINCREASE_raw_counts.tsv")
write.table(as.data.frame(count_data[row_names_inc,]), file = heatmap_outname_tsv, quote = FALSE, sep = "\t", col.names=TRUE)
##################################################end of INCREASE raw counts

##################################################start of DECREASE raw counts
breaksList <- seq(min(count_data[row_names_dec,]), max(count_data[row_names_dec,]), by = max(count_data[row_names_dec,]) / 10000)

heatmap_raw_counts_pdf=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fDECREASE_raw_counts.pdf")
pdf(heatmap_raw_counts_pdf)
pheatmap(count_data[row_names_dec,], show_rownames=TRUE,show_colnames = FALSE, cluster_cols=FALSE, cluster_rows=FALSE, annotation_col = df, main= paste(PROJECTS, "\nraw counts\ndecreasing log2foldChange" , sep=" "), color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)
dev.off()

heatmap_outname_tsv=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fDECREASE_raw_counts.tsv")
write.table(count_data[row_names_dec,], file = heatmap_outname_tsv, quote = FALSE, sep = "\t", col.names=TRUE)
##################################################end of DECREASE raw counts


##############################################start of INCREASE vsd counts

breaksList <- seq(min(assay(vsd)[row_names_inc,]), max(assay(vsd)[row_names_inc,]), by = max(assay(vsd)[row_names_inc,]) / 10000)

heatmap_vsd_counts_pdf=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fINCREASE_vsd_counts.pdf")
pdf(heatmap_vsd_counts_pdf)
pheatmap(assay(vsd)[row_names_inc,],  show_rownames=TRUE,show_colnames = FALSE, cluster_cols=FALSE, cluster_rows=FALSE, annotation_col = df, main = paste(PROJECTS, "\nvsd counts\nincreasing log2foldChange", sep=" "), color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)
dev.off()
heatmap_outname_tsv=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fINCREASE_vsd_counts.tsv")
write.table(assay(vsd)[row_names_inc,], file = heatmap_outname_tsv, quote = FALSE, sep = "\t", col.names=TRUE)
##################################################end of INCREASE vsd counts

##################################################start of DECREASE vsd counts
breaksList <- seq(min(assay(vsd)[row_names_dec,]), max(assay(vsd)[row_names_dec,]), by = max(assay(vsd)[row_names_dec,]) / 10000)

heatmap_vsd_counts_pdf=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fDECREASE_vsd_counts.pdf")
pdf(heatmap_vsd_counts_pdf)
pheatmap(assay(vsd)[row_names_dec,], show_rownames=TRUE,show_colnames = FALSE, cluster_cols=FALSE, cluster_rows=FALSE, annotation_col = df, main= paste(PROJECTS, "\nvsd counts\ndecreasing log2foldChange" , sep=" "), color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)
dev.off()

heatmap_outname_tsv=file.path(OUTPUT_PATH, "DESeq2_heatmap_log2fDECREASE_vsd_counts.tsv")
write.table(assay(vsd)[row_names_dec,], file = heatmap_outname_tsv, quote = FALSE, sep = "\t", col.names=TRUE)
##################################################end of DECREASE vsd counts



# plot the raw count data sorted by logfold:

 #adding PCA (principal component analysis)
#rlog() may take a long time with 50 or more samples,
#vst() is a much faster transformation

#library("pcaExplorer")
#rlt <- DESeq2::rlogTransformation(ddsMF)
#rlt <- vst(DESeq_ddsMF, blind=FALSE) # variance stabilizing transformation
#rlt <- rlogTransformation(ddsMF)
rlt <- vst(DESeq_ddsMF)

#pcaobj <- prcomp(t(assay(rlt)))
pcaobj <- prcomp((assay(rlt)))
head(pcaobj, 1)

groups <- colData(DESeq_ddsMF)$vital_state
groups <- factor(groups,levels=unique(groups))
cols <- scales::hue_pal()(2)[groups]

pca_filename_pdf = file.path(OUTPUT_PATH, paste('DESeq2_pca.pdf', sep=''))
pca_filename_tsv = file.path(OUTPUT_PATH, paste('DESeq2_pca.tsv', sep=''))
pdf(pca_filename_pdf)
genespca(rlt,ntop=100,arrowColors=cols,groupNames=groups)
dev.off()
genespca_DF <- genespca(rlt,ntop=100,arrowColors=cols,groupNames=groups, returnData = TRUE)
write.table(as.data.frame(genespca_DF), file = pca_filename_tsv, quote = FALSE, sep='\t', col.names=NA)

pca_filename_groups_pdf = file.path(OUTPUT_PATH, paste('DESeq2_pca_groups.pdf', sep=''))
pca_filename_groups_tsv = file.path(OUTPUT_PATH, paste('DESeq2_pca_groups.tsv', sep=''))
pdf(pca_filename_groups_pdf)
pcaplot(rlt, intgroup = "vital_state", ntop = 500, title = PROJECTS, pcX = 1, pcY = 2, text_labels = FALSE, point_size = 3, ellipse = TRUE, ellipse.prob = 0.95)
dev.off()
pca_DF <- pcaplot(rlt, intgroup = "vital_state", ntop = 500, returnData = TRUE, title = PROJECTS, pcX = 1, pcY = 2, text_labels = TRUE, point_size = 3, ellipse = TRUE, ellipse.prob = 0.95)
write.table(as.data.frame(pca_DF), file = pca_filename_groups_tsv, quote = FALSE, sep='\t', col.names=NA)

pca_filename_test_pdf = file.path(OUTPUT_PATH, paste('DESeq2_pca_plot.pdf', sep=''))
pca_filename_test_tsv = file.path(OUTPUT_PATH, paste('DESeq2_pca_plot.tsv', sep=''))
pdf(pca_filename_test_pdf)
plotPCA(rlt, intgroup = "vital_state")
dev.off()
plotPCA_DF  <- plotPCA(rlt, intgroup = "vital_state", returnData = TRUE)
write.table(as.data.frame(plotPCA_DF), file = pca_filename_test_tsv, quote = FALSE, sep='\t', col.names=NA)
