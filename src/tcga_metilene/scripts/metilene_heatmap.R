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


#count_inname <- file.path('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-LUSC/DESeq2/DESeq2_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_2/summary_for_DESeq2.tsv')
#info_inname  <- file.path('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-LUSC/DESeq2/DESeq2_input_table/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_2/summary_for_DESeq2_INFO.tsv')
#OUTPUT_PATH  <- file.path('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_3/TCGA-LUSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_2')
#PROJECTS <- 'TCGA-LUSC'
metilene_intersect <- '/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/metilene_intersect.tsv'
read.table(metilene_intersect, nrows=5)

thresh_pos = '/scr/palinca/gabor/TCGA-pipeline_2/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female_male/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr10_122879357_122879744.tsv'

pdf(heatmap_outname_pdf)
pheatmap(table_counts_normalized_MF[row_names_inc,], show_rownames=TRUE,show_colnames = FALSE, cluster_cols=FALSE, cluster_rows=FALSE, annotation_col = df, main=paste(PROJECTS, "\nnormalized counts\nincreasing log2foldChange", sep=" ") , color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)
#pheatmap(table_counts_normalized_MF[select_increase_log2Fold,], show_rownames=TRUE,show_colnames = FALSE, cluster_cols=FALSE, annotation_col = df, main=PROJECTS)
dev.off()
