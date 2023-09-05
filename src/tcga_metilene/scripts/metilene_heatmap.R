suppressPackageStartupMessages({
require(pheatmap)
require(RColorBrewer)
})
args <- commandArgs(trailingOnly = TRUE)

temp_args <- c('beta_val_inname', 'header_temp', 'pdf_heatmap_out', 'PROJECTS', 'gender', 'drug_combi', 'cutoff', 'DMR', 'q_value', 'mmd', 'CpGs')
cat("arguments:\n")
j <- 1
for (i in args){
    cat(temp_args[j],  "<-", i, "\n")
    j <- j + 1
}

# test set:

#beta_val_inname <- '/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/metilene_intersect_heatmaps_beta_value_chr7_83648624_83648822_input_temp.tsv'
#header_temp <- '/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/metilene_intersect_heatmaps_beta_value_chr7_83648624_83648822_header_temp.tsv'
#pdf_heatmap_out <- '/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_5/metilene_intersect_heatmaps_beta_value_chr7_83648624_83648822.pdf'
#PROJECTS <- 'TCGA-CESC'
#gender <- 'female'
#drug_combi <- 'carboplatin_carboplatin,paclitaxel_cisplatin'
#cutoff <- '5'
#DMR <- 'chr7_83648624_83648822'
#q_value <- '0.00064942'
#mmd <- '0.199076'
#CpGs <- '5' 

beta_val_inname <- args[1]  # beta_table
header_temp <- args[2]  # info_table
pdf_heatmap_out <- args[3] # file path to pdf
PROJECTS  <- args[4] # formatted header string for the heatmap plot, including PROJECTs, drugs, gender, cutoff, DMR
gender  <- args[5]
drug_combi <- args[6]
cutoff  <- args[7]
DMR  <- args[8]
q_value  <- args[9]
mmd  <-  args[10]
CpGs  <-  args[11]


#beta_val_inname <- file.path('/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_intersect_temp.tsv')
#header_temp  <- file.path('/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_intersect_temp_info.tsv')
beta_val_temp <- read.table(beta_val_inname,header=T,row.names=1)
### the default sep value for read.delim is '\t' and handles better the spaces
### which occur in some drug name cols
info_data <- t(read.table(header_temp,row.names=1))

df_info_data <-  as.data.frame(info_data)
pdf(pdf_heatmap_out)
pheatmap(as.data.frame(beta_val_temp), show_rownames=TRUE,show_colnames = FALSE, cluster_cols=FALSE, cluster_rows=FALSE, main=paste(PROJECTS,"beta values", "\n",gender, "\n", drug_combi, "\nDMR:", DMR, "\nCutoff:", cutoff, "; mmd:", mmd, "; q-value:", q_value, "\n#CpGs:", CpGs, sep <- " "), annotation_col = df_info_data )
dev.off()
