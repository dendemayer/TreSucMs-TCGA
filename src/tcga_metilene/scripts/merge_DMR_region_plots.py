# import subprocess as sbp
import os
import glob
import sys
import re
from PyPDF2 import PdfMerger

"""
important, don't invoke pdf's which are to small:
"""
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print('# snakemake inputs:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

print('# snakemake output:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

print('# snakemake wildcards:')
[ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]
All_plots_merged = snakemake.output.All_plots_merged

"""
if the input list has just len 1, that means that the intersect table is empty,
no DMRs could be found, write an empty All_plots_merged.pdf file and exit:
# (Pdb) snakemake.input
# ['/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_4/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/male/cutoff_0/metilene_intersect.tsv']
"""

if len(snakemake.input) == 1:
    print(f'no merge with qpdf possible, writing the empty {All_plots_merged}')
    open(All_plots_merged, 'a').close()
    os._exit(0)



# we have:
# "_boxplot_beta_value" plots
# "_lineplot_median_beta_value" plots
# beyond the threshold dir:
# "_lifeline_plot_"

search_path =  os.path.split(os.path.split(All_plots_merged)[0])[0]
boxplots =  [pdf_file for pdf_file in glob.glob(os.path.join(search_path, '*_boxplot_beta_value*pdf')) if os.lstat(pdf_file)[6] > 0]
# include DMRs only if also a boxplot is available, meaning, leaving out empty
# plots
DMRs = [re.search('chr.*', os.path.split(i)[1].replace('.pdf', '')).group() for i in  boxplots if re.search('chr.*', os.path.split(i)[1].replace('.pdf', ''))]
lineplots = [pdf_file for pdf_file in  glob.glob(os.path.join(search_path, '*_lineplot_median_beta_value*pdf')) if os.lstat(pdf_file)[6] > 0]
search_path_lifeplots = os.path.split(All_plots_merged)[0]
lifeline_list = []
for dmr in DMRs:
    # lifeline_list.append([os.path.split(i)[1] for i in glob.glob(os.path.join(search_path_lifeplots, f'*{dmr}*pdf'))])
    lifeline_list.append([i for i in glob.glob(os.path.join(search_path_lifeplots, f'*{dmr}*pdf'))])
# pdfs_to_merge_str = " ".join(pdfs_to_merge)
pdfs_to_merge = []
for i,j,k in zip(boxplots, lineplots, lifeline_list):
    pdfs_to_merge = pdfs_to_merge + [i] + [j] + k

# last check, whether all plots are not empty (the lifeline plots are not
# checked yet
pdfs_to_merge = [i for i in pdfs_to_merge if os.lstat(i)[6] > 0]

# with pdf merger just the pdfs_to_merge list is needed:

merger = PdfMerger()

for pdf in pdfs_to_merge:
    merger.append(pdf)

merger.write(All_plots_merged)
merger.close()


# # pdfs_to_merge_str = " ".join(pdfs_to_merge)
# final_list = ["qpdf", "--empty",  "--pages"] + pdfs_to_merge +  ["--", All_plots_merged]
# final_list = ["pdfcombine"] + pdfs_to_merge + ["-o", All_plots_merged]
# final_str = " ".join(final_list)
# # sbp.call(args=final_list)
# print(f'calling qpdf with {final_str}')
# # trys = 3
# process = sbp.Popen(final_list, stdout=sys.stdout, stderr=sys.stderr)
# exit_code=process.wait()
# if exit_code != 0:
#     raise Exception(f"pdfcombine merge failed with exit code {exit_code}")
#     # while trys > 0:
#     #     print(f'{exit_code}, trying {trys} more')
#     #     process = sbp.Popen(final_list, stdout=sys.stdout, stderr=sys.stderr)
#     #     exit_code=process.wait()
#     #     if exit_code == 0:
#     #         print('qpdf succeeded')
#     #         os._exit(0)
#     #     else:
#     #         trys -= 1
#     #     # raise Exception(f"qpdf merge failed with exit code {exit_code}")

# # sbp.Popen(args=["qpdf --empty --pages", {pdfs_to_merge_str} -- {All_plots_merged}", shell=True, stderr=sys.stderr, stdout=sys.stdout)
