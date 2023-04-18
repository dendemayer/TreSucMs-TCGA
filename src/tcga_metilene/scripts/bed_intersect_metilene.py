import pandas as pd
from pybedtools import BedTool
import os
from natsort import natsort_keygen
import sys

"""
DF_met_out: (in system without header)
chr     start           end             p_value         mean_methylation_difference #CpGs   mean g1 mean g2
chr7    95396298        95397846        8.6644e-05      0.061329                    25      0.36234 0.30101
chr19   57769432        57770082        0.0076932       0.114709                    9       0.29953 0.18483
chr20   58850459        58852155        5.2322e-11      0.050191                    45      0.62156 0.56942

# for every region we want to report the locations we can annotate:
DF_data:

       Composite Element REF  Beta_value Chromosome      Start        End  ...                                          Gene_Type                                      Transcript_ID                               Position_to_TSS                CGI_Coordinate Feature_Type
0                 cg00000029    0.708892      chr16   53434200   53434201  ...       protein_coding;protein_coding;protein_coding  ENST00000262133.9;ENST00000544405.5;ENST000005...                                -221;-1420;222   CGI:chr16:53434489-53435297      N_Shore
1                 cg00000108         NaN       chr3   37417715   37417716  ...  lincRNA;lincRNA;lincRNA;lincRNA;lincRNA;lincRN...  ENST00000328376.8;ENST00000332506.6;ENST000004...  18552;18552;6505;31445;18143;447;18552;18552    CGI:chr3:37451927-37453047            .
-> add here a col in the form of chrxx_start_end

(Pdb) DF_met_input:
       Chromosome     Start  ...  alive;4da121da-9193-4ac5-8206-03cb1f7c2c14;cisplatin;female;TCGA-CESC  alive;e89d0cb1-6ec2-43ee-a6d9-b922b680af34;cisplatin;female;TCGA-CESC
0            chr1     15865  ...                                           0.916271                                                               0.862795
1            chr1     18827  ...                                           0.642255                                                               0.658495
those are the beta values, we want to link them for later kaplan meier plots


3 infos needed, the summary table:
    /scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines/TCGA-CESC/metilene/merged_meta_files/meta_info_druglist_merged_drugs_combined.tsv
the metilene out result:
/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines/TCGA-LIHC/metilene/metilene_output/carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_complement_qval.0.05.out
and the metilene input table with which the metilene output was created with:
'/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines/TCGA-CESC/metilene/metilene_input_table/carboplatin,paclitaxel_cisplatin/male/cutoff_0/summary_for_metilene_complement.tsv'
"""
sys.stdout = sys.stderr = open(snakemake.log[0], "w")

OUTPUT_PATH = snakemake.wildcards.output_path
# meta_info = snakemake.input.meta_info
qval_out = snakemake.input.qval_out
metilene_input = snakemake.input.metilene_input  # this is the beta value summary from all patients applied
intersected_out = snakemake.output.out_file
PROJECT = snakemake.wildcards.project


# data_mid = 'metilene/data_files/'
try:
    # DF_summary = pd.read_table(meta_info)
    DF_met_out = pd.read_table(qval_out, header=None)
    DF_met_input = pd.read_table(metilene_input, na_values='.').dropna(how='any')
    # if the DF_met_out DF is empty, nothing can be intersected, independent of
    # the content of the metilene input (metilene_input) table, at this point, an
    # empty DF out can be written as result of this step
    if DF_met_out.empty or DF_met_out.shape == (1, 1):
        pd.DataFrame().to_csv(intersected_out)
        os._exit(0)
except Exception as e:
    print(f'cought exception {e}')
    # print('writing empty files:')
    # breakpoint()
    # empty_DF = pd.DataFrame()
    os._exit(0)



DF_met_input['End'] = DF_met_input['Start'] + 1
DF_met_input = pd.concat([DF_met_input.loc[:, ['Chromosome', 'Start', 'End']], DF_met_input.drop(['Chromosome', 'Start', 'End'], axis=1)], axis=1).sort_values(['Chromosome', 'Start'], key=natsort_keygen())


#######
# first plot the distribution of betavalues of the regions found
# DF_met_input (betavalues) intersected with DF_met_out:
met_inp_bed = BedTool.from_dataframe(DF_met_input)
met_out_bed = BedTool.from_dataframe(DF_met_out)
met_inp_out_int = met_inp_bed.intersect(met_out_bed)
DF_met_inp_out = met_inp_out_int.to_dataframe(names=DF_met_input.columns)
# now we want to know to which region each position belongs, s.t. the plots can
# be then filtered on that col:
def return_region_entry(row):
    position_bed = BedTool.from_dataframe(pd.Series({'chr': row.name[0], 'start': row.name[1], 'end': row.name[2]}).to_frame().T)
    # return the region out of met_out, in which the applied pos was found in:
    return '_'.join([str(i) for i in met_out_bed.intersect(position_bed, wa=True).to_dataframe().iloc[0, [0,1,2]].to_list()])


# make a multiindex of the vital_status;case_id;PROJECT;DRUGS header, s.t. in
# can be read on that basis by the following pandas methods
# pd.read_table('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines/TCGA-CESC/metilene/metilene_output/carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_intersect.tsv', header=[0,1,2,3,4], index_col=[0,1,2])
DF_met_inp_out.set_index(['Chromosome', 'Start', 'End'], inplace=True)
col_t = [tuple(x) for x in [i.split(';') for i in DF_met_inp_out.columns]]
# MI_start.extend(col_t)
MI = pd.MultiIndex.from_tuples(col_t, names=('vital_status', 'case_id', 'drugs', 'gender', 'projects'))
DF_met_inp_out.columns = MI

DF_met_inp_out['region'] = DF_met_inp_out.apply(return_region_entry, axis=1)
DF_met_inp_out = DF_met_inp_out.set_index('region', append=True)

# also regions are reportet, which do not have any beta value at all, they
# would cause errors while plotting the DMRs
# region chr19_58228367_58228578
# '/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_2/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin_paclitaxel/female_male/cutoff_0/metilene_complement_intersect.tsv'
# range_ = 'chr19_58228367_58228578'
# DF_met_inp_out.loc[(slice(None), slice(None), slice(None), range_), :]
# problem solved with dropping not present betavalues whilst reading in the
# input for metilene:
# DF_met_input = pd.read_table(metilene_input, na_values='.').dropna(how='all')
# if intersected_out=='/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_2/TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin_paclitaxel/female_male/cutoff_0/metilene_complement_intersect.tsv':
    # breakpoint()
DF_met_inp_out.to_csv(intersected_out, sep='\t')

# to read it correctly with the MI use header and index_col:
# pd.read_table('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines/TCGA-CESC/metilene/metilene_output/carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_intersect.tsv', header=[0,1,2,3,4], index_col=[0,1,2])



# annot_bed = BedTool.from_dataframe(DF_data_bed)
# annot_met_int = annot_bed.intersect(met_out_bed)
# met_annot_int = met_out_bed.intersect(annot_bed)
# DF_annot_met = BedTool.to_dataframe(annot_met_int)
# DF_met_annot = BedTool.to_dataframe(met_annot_int)
# breakpoint()

# sort the DFs, and make first 3 cols chr, start, stop


# for every DMR in DF_met add the available infos out of the annotation and to
# the case
# DF_empty = pd.read_table('/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines/TCGA-CESC/metilene/metilene_output/carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_out_sorted_complement.tsv')
