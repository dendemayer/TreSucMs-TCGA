#!/usr/bin/env Rscript
suppressPackageStartupMessages({
require(sjmisc)  # for the str_contains() fct
require(DESeq2)
})
args <- commandArgs(trailingOnly = TRUE)

cat("arguments:\n")
j=1
for (i in args){
    cat(j, ": ", i, "\n")
    j <- j + 1
}

#count_inname <-file.path('/scr/palinca/gabor/TCGA-pipeline/TCGA-BRCA/DESeq2/DESeq2_input_table/adriamycin,cytoxan,paclitaxel,tamoxifen_cyclophosphamide,doxorubicin,paclitaxel,tamoxifen_tamoxifen/female/cutoff_10/summary_for_counts.tsv')
#info_inname <- file.path('/scr/palinca/gabor/TCGA-pipeline/TCGA-BRCA/DESeq2/DESeq2_input_table/adriamycin,cytoxan,paclitaxel,tamoxifen_cyclophosphamide,doxorubicin,paclitaxel,tamoxifen_tamoxifen/female/cutoff_10/summary_for_counts_INFO.tsv') 
#OUTPUT_PATH <- file.path('/scr/palinca/gabor/TCGA-pipeline/TCGA-BRCA/DESeq2/DESeq2_output/adriamycin,cytoxan,paclitaxel,tamoxifen_cyclophosphamide,doxorubicin,paclitaxel,tamoxifen_tamoxifen/female/cutoff_10')
#PROJECTS <-   'TCGA-BRCA'

### test files:
#count_inname <-file.path('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC/DESeq2/DESeq2_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/summary_for_counts.tsv')
#info_inname <- file.path('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC/DESeq2/DESeq2_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/summary_for_counts_INFO.tsv')
#OUTPUT_PATH <- file.path('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0')
#PROJECTS <-              'TCGA-CESC'

#count_inname <-file.path('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/summary_for_counts.tsv')
#info_inname <- file.path('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/summary_for_counts_INFO.tsv')
#OUTPUT_PATH <- file.path('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0')
#PROJECTS <- 'TCGA-CESC_TCGA-HNSC_TCGA-LUSC'
##### END TEST

# ### ommit this block if you want to test the script with some test data>>>
count_inname <- args[1]  # count_table
info_inname <- args[2]  # info_table
OUTPUT_PATH  <- args[3]  # deseq <- output
PROJECTS <- args[4]
# ### ommit this block if you want to test the script with some test data<<<

count_data <- as.matrix(read.delim(count_inname,header=T,row.names=1))
info_data <- as.matrix(read.delim(info_inname,row.names=1))
colnames(count_data) <- info_data[5,]

vital_status <- as.vector(info_data[1,])
vital_status <- data.frame(vital_status)
colnames(vital_status)<-c("vital_state")
vital_status$vital_state <- as.factor(vital_status$vital_state)
tmp_cancer_column  <- as.vector(info_data[3,])
tmp_cancer_column  <- data.frame(tmp_cancer_column)
colnames(tmp_cancer_column) <- c("cancer")
tmp_cancer_column$cancer <- as.factor(tmp_cancer_column$cancer)

tmp_gender_column  <- as.vector(info_data[2,])
tmp_gender_column  <- data.frame(tmp_gender_column)
colnames(tmp_gender_column) <- c("gender")
tmp_gender_column$gender <- as.factor(tmp_gender_column$gender)

# including the drugname
drugnames_column <- as.vector(info_data[4,])
drugnames_column  <- data.frame(drugnames_column)
colnames(drugnames_column) <-c("drugnames")
drugnames_column$drugnames <- as.factor(drugnames_column$drugnames)
# We create a copy of the DESeqDataSet,
# so that we can rerun the analysis using a multi-factor design.

vital_cancer <- cbind(tmp_gender_column, tmp_cancer_column, drugnames_column,
                     vital_status)
get_col_flags  <- function() { # get colnames of ambigous cols
    #print('available cols of count DataFrame: ')
    #print(colnames(vital_cancer))
    col_flags <- vector()
    for (col_name in colnames(vital_cancer)) {
        if (dim(unique(vital_cancer[col_name]))[1] > 1) {
            col_flags <- c(col_flags, col_name)
            }
    }
    return(col_flags)
}

# get the right cols for the DESeq design, of all cols available choose the
# ones wich have ambiguos values(unique value count higher 1)
col_flags <- get_col_flags()
col_flags

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


#ddsMF <- get_ddsMF(col_flags)

# since we just need the normalisation of the count data, no fancy design is
# needed, just take the vital state as design
ddsMF  <- DESeqDataSetFromMatrix(countData = count_data, colData = vital_cancer, design= ~ vital_state )

colData(ddsMF)
# Pre-filtering.
#keep = rowSums(counts(ddsMF)) >= 10
#ddsMF <- ddsMF[keep,]

# the output should be:
#/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/DESeq2_nt_counts_all_cases.tsv
#/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/DESeq2_raw_counts_all_cases.tsv
#/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/DESeq2_norm_counts_all_cases.tsv
#/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-CESC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/DESeq2_vsd_counts_all_cases.tsv
# no DE analysis necessary, just create the counts:
raw_counts <- counts(ddsMF)
write.table(raw_counts, file=file.path(OUTPUT_PATH, "DESeq2_raw_counts_all_cases.tsv"), sep="\t", quote=F, col.names=NA)

ddsMF <- estimateSizeFactors(ddsMF)
sizeFactors(ddsMF)

normalized_counts <- counts(ddsMF, normalized=TRUE)
write.table(normalized_counts, file=file.path(OUTPUT_PATH, "DESeq2_norm_counts_all_cases.tsv"), sep="\t", quote=F, col.names=NA)


ntd <- as.data.frame(assay(normTransform(ddsMF)))
write.table(ntd, file=file.path(OUTPUT_PATH, "DESeq2_nt_counts_all_cases.tsv"), sep="\t", quote=F, col.names=NA)

vsd <- as.data.frame(assay(vst(ddsMF, blind=FALSE))) # variance stabilizing transformation
write.table(vsd, file=file.path(OUTPUT_PATH, "DESeq2_vsd_counts_all_cases.tsv"), sep="\t", quote=F, col.names=NA)
