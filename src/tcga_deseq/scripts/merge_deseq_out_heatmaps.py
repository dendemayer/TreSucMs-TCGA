import PyPDF2
import os
import sys

sys.stderr = sys.stdout = open(snakemake.log[0], 'w')

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

print('# snakemake params:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.params.items()]

heatmaps = snakemake.input.heatmaps
heatmaps_merged = snakemake.output.heatmaps_merged
###############################################################################
#                             test_set                                        #
###############################################################################
# # snakemake inputs:
# heatmaps = ["/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/DESeq2_heatmap_log2fINCREASE_nt_counts.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DEn,paclitaxel_cisplatin/male/cutoff_5/DESeq2_heatmap_log2fDECREASE_nt_counts.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/DESeq2_heatmap_log2fINCREASE_TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/DESeq2_heatmap_log2fDECREASE_raw_counts.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboptoff_5/DESeq2_heatmap_log2fINCREASE_norm_counts.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/DESeq2_heatmap_log2fDECREASE_norm_counts.pdf", "/scr/palincaESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/DESeq2_heatmap_log2fINCREASE_vsd_counts.pdf", "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/DECREASE_vsd_counts.pdf"]
# # snakemake output:
# heatmaps_merged = "/scr/palinca/gabor/TCGA-pipeline_7/TCGA-HNSC/DESeq2/DESeq2_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_5/DESeq2_heatmap_merged.pdf"
# # snakemake wildcards:
# output_path = "/scr/palinca/gabor/TCGA-pipeline_7"
# project = "TCGA-HNSC"
# drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
# gender = "male"
# cutoff = "cutoff_5"
###############################################################################
#                             test_set                                        #
###############################################################################

heatmaps = [i for i in heatmaps if os.path.getsize(i) > 1]
if len(heatmaps) == 0:
    open(heatmaps_merged, 'a').close()
    os._exit(0)
merger = PyPDF2.PdfMerger()
[merger.append(pdf) for pdf in heatmaps]
merger.write(heatmaps_merged)
merger.close()
