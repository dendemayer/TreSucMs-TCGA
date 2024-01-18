import os
import pandas as pd
import sys
import re

"""
what would be expect: which group has a longer lifeexpectation in the base lifeline plot:

                cases in therapy:                  |            all cases:
                (base_plot)                        |           (validation_plots)
                                                   |  - UP regulated in therapy have a higher expectation
                                                   |   (injects - UP regulated in not in therapy have lower life expectation)
                                                   |/
            UP regulated - high expectation          |
      (injects DOWN regulated - low expectation)   |\
                                                   |  - DOWN regulated in therapy have higher expectation
                                                   |  (injects - DOWN regulated not in therapy have lower expectation)
                                                   |
                                                   |  - UP regulated in therapy have a higher expectation
                                                   |/ (injects - UP regulated not in therapy have a lower expectation)
            UP regulated - low expectation         |
(injects DOWN regulated - higher expectation)      |\
                                                   |  - DOWN regulated in therapy have higher expectation
                                                   |  ( injects - DOWN regulated in therapy have lower expectation)

- within the base lifeline plot the UP regulated are marked depending on
whether they have low or high expectation,
testing on:
 UP_regulated_high_expectation   --> UP_regulated, in_therapy have higher expectation    --> DOWN_regulated
         ---small p_value---                       ---small p_value---                       ---high p_value---
 DOWN_regulated_high_expectation -->  DOWN_regulated, in therapy have higher expactatoin --> UP_regulated

TODO:
add info which group has the higher life expancy to each of the 3 plot types


| plot_type       | fst_life_mean             | scnd_life_mean                | CMP                      |
| -               | -                         | -                             | -                        |
| base_plot       | UP_life_expectancy_mean   | DOWN_life_expectancy_mean     | CMP_UP-DOWN              |
| UP_validation   | UP_in_therapy_life_mean   | UP_not_in_therapy_life_mean   | CMP_in or not in therapy |
| DOWN_validation | DOWN_in_therapy_life_mean | DOWN_not_in_therapy_life_mean | CMP_in or not in therapy |

base_plot: UP is higher -> expectation would be that in UP_validation the UP_in_therapy is also higher
base_plot: DOWN is higher -> expectation would be that in DOWN_validation the DOWN_in_therapy is also higher

create CMP col which gives the higher value each plot type -> either UP or DOWN ->

for base_plot that would be comparison between UP and DOWN lifeexpancy,
for UP_validation, that would be comparison between IN therapy or NOT IN therapy
for DOWN_validation, that would be comparison between IN therapy or NOT IN therapy


a base plot is labeled an UP plot if the first life mean is higher than the scnd life mean (and vice verca) (fst is UP)
the UP_validation plot is labeled UP if fst_life_mean is higher than the scnd life mean (and vice verca)  (fst is in therapy)
the DOWN_validation plot is labeled DOWN if the fst_life_mean is higher than the scnd life mean (and vice verca) (fst is in therapy)

also include the threshold info.
"""

###############################################################################
#                              snakemake inputs                               #
###############################################################################
if "snakemake" in dir():
    sys.stderr = sys.stdout = open(snakemake.log[0], 'w')

    print('# snakemake inputs:')
    [ print(f'{i[0]} = "{i[1]}"') for i in snakemake.input.items()]

    print('# snakemake output:')
    [ print(f'{i[0]} = "{i[1]}"') for i in snakemake.output.items()]

    print('# snakemake wildcards:')
    [ print(f'{i[0]} = "{i[1]}"') for i in snakemake.wildcards.items()]

    print('# snakemake params:')
    [ print(f'{i[0]} = "{i[1]}"') for i in snakemake.params.items()]

    lifeline_all = snakemake.input.lifeline_all
    metilene_out = snakemake.input.metilene_out
    metilene_lifeline_aggregated = snakemake.output[0]
else:
    # snakemake inputs:
    lifeline_all = ["/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chrX_72305564_72308502_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr8_140238482_140238568_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr20_58850459_58852587_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr6_3848399_3850498.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr7_95396094_95397226_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chrX_72305564_72308502.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr20_58850459_58852587_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chrX_16711802_16712814_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr7_83648624_83648822_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr6_26183279_26184017_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr6_26183279_26184017_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr14_85529381_85531409.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chrX_16711802_16712814.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr14_85529381_85531409_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr6_3848399_3850498_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chrX_72305564_72308502_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr7_95396094_95397226_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chrX_16711802_16712814_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr20_58850459_58852587.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr7_95396094_95397226.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr6_3848399_3850498_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr8_140238482_140238568.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr6_26183279_26184017.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr8_140238482_140238568_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr7_83648624_83648822.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr14_85529381_85531409_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr7_83648624_83648822_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chrX_72305564_72308502_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr8_140238482_140238568_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr20_58850459_58852587_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr6_3848399_3850498.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr7_95396094_95397226_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chrX_72305564_72308502.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr20_58850459_58852587_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chrX_16711802_16712814_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr7_83648624_83648822_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr6_26183279_26184017_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr6_26183279_26184017_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr14_85529381_85531409.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chrX_16711802_16712814.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr14_85529381_85531409_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr6_3848399_3850498_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chrX_72305564_72308502_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr7_95396094_95397226_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chrX_16711802_16712814_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr20_58850459_58852587.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr7_95396094_95397226.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr6_3848399_3850498_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr8_140238482_140238568.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr6_26183279_26184017.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr8_140238482_140238568_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr7_83648624_83648822.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr14_85529381_85531409_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_5/metilene_intersect_lifeline_plot_chr7_83648624_83648822_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chrX_72305564_72308502_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr8_140238482_140238568_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr20_58850459_58852587_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr6_3848399_3850498.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr7_95396094_95397226_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chrX_72305564_72308502.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr20_58850459_58852587_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chrX_16711802_16712814_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr7_83648624_83648822_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr6_26183279_26184017_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr6_26183279_26184017_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr14_85529381_85531409.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chrX_16711802_16712814.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr14_85529381_85531409_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr6_3848399_3850498_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chrX_72305564_72308502_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr7_95396094_95397226_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chrX_16711802_16712814_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr20_58850459_58852587.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr7_95396094_95397226.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr6_3848399_3850498_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr8_140238482_140238568.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr6_26183279_26184017.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr8_140238482_140238568_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr7_83648624_83648822.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr14_85529381_85531409_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_10/metilene_intersect_lifeline_plot_chr7_83648624_83648822_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chrX_72305564_72308502_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr8_140238482_140238568_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr20_58850459_58852587_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr6_3848399_3850498.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr7_95396094_95397226_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chrX_72305564_72308502.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr20_58850459_58852587_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chrX_16711802_16712814_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr7_83648624_83648822_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr6_26183279_26184017_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr6_26183279_26184017_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr14_85529381_85531409.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chrX_16711802_16712814.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr14_85529381_85531409_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr6_3848399_3850498_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chrX_72305564_72308502_UP_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr7_95396094_95397226_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chrX_16711802_16712814_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr20_58850459_58852587.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr7_95396094_95397226.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr6_3848399_3850498_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr8_140238482_140238568.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr6_26183279_26184017.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr8_140238482_140238568_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr7_83648624_83648822.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr14_85529381_85531409_DOWN_val.tsv", "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_20/metilene_intersect_lifeline_plot_chr7_83648624_83648822_UP_val.tsv"]
    metilene_out = "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_qval.0.05.out"
    script_file = "/homes/biertruck/gabor/phd/test_git_doc/tcga_piplines/src/shared/../tcga_metilene/scripts/aggregate_lifelines_all.py"
    # snakemake output:
    metilene_lifeline_aggregated = "/scr/palinca/gabor/TCGA-pipeline_8/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0_threshold_5_threshold_10_threshold_20/metilene_lifelines_aggregated.tsv.gz"
    # snakemake wildcards:
    output_path = "/scr/palinca/gabor/TCGA-pipeline_8"
    project = "TCGA-CESC_TCGA-HNSC_TCGA-LUSC"
    drug_combi = "carboplatin_carboplatin,paclitaxel_cisplatin"
    gender = "female"
    cutoff = "cutoff_0"

###############################################################################
#                                   content                                   #
###############################################################################


if len(lifeline_all) == 0:
    print(f'saving empty table in {metilene_lifeline_aggregated}')
    pd.DataFrame().to_csv(metilene_lifeline_aggregated, sep='\t', compression='gzip')
    os._exit(0)

# check on empty DFs, do not integrate them:
# DF_aggr_list = [pd.read_table(i).loc[0, ['count_type', 'p_value', 'ENSG', 'plot_type','threshold', 'fst_life_mean', 'scnd_life_mean']] for i in lifeline_all if not pd.read_table(i).empty]
# ##########
# first = '/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr20_58851489_58852155_DOWN_val.tsv'
# second = '/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chr7_27112963_27116443_UP_val.tsv'
# third = '/scr/palinca/gabor/TCGA-pipeline/TCGA-CESC/metilene/metilene_output/carboplatin_carboplatin,paclitaxel_cisplatin/female/cutoff_0/threshold_0/metilene_intersect_lifeline_plot_chrX_16711440_16712814.tsv'
# temp = [first, second, third]

# (Pdb) pd.read_table(first).columns
# Index(['case_id', 'group', 'vital_status', 'drugs', 'gender', 'projects',
#        'beta_values', 'T', 'E', 'in_therapy', 'p_value', 'start', 'chr',
#        'threshold', 'fst_life_mean', 'scnd_life_mean', 'plot_type', 'ENSG',
#        'gene_type', 'gene_status', 'gene_name'],
#       dtype='object')
# (Pdb) pd.read_table(second).columns
# Index(['case_id', 'group', 'vital_status', 'drugs', 'gender', 'projects',
#        'beta_values', 'T', 'E', 'in_therapy', 'p_value', 'start', 'chr',
#        'threshold', 'fst_life_mean', 'scnd_life_mean', 'plot_type', 'ENSG',
#        'gene_type', 'gene_status', 'gene_name'],
#       dtype='object')
# (Pdb) pd.read_table(third).columns
# Index(['Unnamed: 0', 'case_id', 'drugs', 'gender', 'projects', 'UP_or_DOWN',
#        'beta_value', 'vital_status', 'survivaltime', 'years_to_last_follow_up',
#        'T', 'E', 'median', 'DMR', 'p_value', 'start', 'chr', 'threshold',
#        'fst_life_mean', 'scnd_life_mean', 'plot_type', 'ENSG', 'gene_type',
#        'gene_status', 'gene_name'],
#       dtype='object')
# ###########
# DF_aggr_list = [pd.concat([pd.read_table(i).loc[0, ['count_type', 'p_value', 'ENSG', 'plot_type','threshold', 'fst_life_mean', 'scnd_life_mean']],  pd.Series([i], index=['file_path'])]) for i in lifeline_all if not pd.read_table(i).empty]
DF_aggr_list = [pd.concat([pd.read_table(i).loc[0, ['p_value', 'ENSG', 'plot_type','threshold', 'fst_life_mean', 'scnd_life_mean', 'start', 'chr']],  pd.Series([i], index=['file_path'])]) for i in lifeline_all if not pd.read_table(i).empty]

# in LUSC (just one dead case, no kaplan meier estimation possible) it can
# occur that the lifeline tables are all empty, if that occurs, write an empty
# table:
if len(DF_aggr_list) == 0:
    pd.DataFrame().to_csv(metilene_lifeline_aggregated, sep='\t')
    os._exit(0)
DF_aggr = pd.concat(DF_aggr_list, axis=1).T.drop_duplicates()  # from where do such duplicates arise, check that
DF_aggr['threshold'] = DF_aggr['threshold'].astype('int16')

DF_aggr['chr_start'] = DF_aggr['chr'] + '_' + DF_aggr['start'].astype('str')

# ENSGs_list = DF_aggr['ENSG'].value_counts().index.tolist()
chr_start_list = DF_aggr['chr_start'].value_counts().index.tolist()
# DF_aggr = DF_aggr.set_index(['ENSG', 'plot_type', 'threshold']).sort_index()
DF_aggr = DF_aggr.set_index(['chr_start', 'plot_type', 'threshold']).sort_index()
thresholds = DF_aggr.index.to_frame()['threshold'].value_counts().index.tolist()
# check whether we have all 3 sets available needed for the evaluation, the
# base plot, and the up and down plot, otherwise delete those ENSG - count
# combinations:
ENSGs_to_delete = []
# DF_aggr = DF_aggr.reset_index().drop_duplicates(subset=['ENSG', 'plot_type', 'threshold', 'p_value', 'fst_life_mean', 'scnd_life_mean']).set_index(['ENSG', 'plot_type', 'threshold'])
DF_aggr = DF_aggr.reset_index().drop_duplicates(subset=['chr_start', 'plot_type', 'threshold', 'p_value', 'fst_life_mean', 'scnd_life_mean']).set_index(['chr_start', 'plot_type', 'threshold'])
# ENSGs_to_delete = []
chr_start_to_delete = []
for thresh in thresholds:
    for start in chr_start_list:
        try:
            if len(DF_aggr.loc[(start, slice(None), thresh),:]) != 3:
                chr_start_to_delete.append((start, thresh))
        except KeyError:
            continue
    # [ENSGs_to_delete.append((ensg, count, thresh)) for ensg in ENSGs_list if len(DF_aggr.loc[(ensg,count, slice(None), thresh),:]) != 3]

# delete the ENSG, count_type and Thresh combinations which do not hold all 3
# needed plots:
if len(chr_start_to_delete) != 0:
    MI = pd.MultiIndex.from_tuples(chr_start_to_delete)
    DF_aggr = DF_aggr.reset_index('plot_type').drop(MI).reset_index().set_index(['chr_start', 'plot_type', 'threshold'])


"""
plot categorisation is done, for all plottypes, now give a rank fct::

UP:   baseplot p_val down  -> UP_val_plot p_val_down   -> DOWN_val_down
            p_val          +    p_val                 +   not of interest
DOWN: baseplot p_val_down  -> UP_val_plot p_val_high   -> DOWN_val_down p_val_down
            p_val          +    not of interest       +   p_val
-> factors that can be included are:
p_values of interest, but also the diff of the life_means ? what would be
interesting is, whether or not an treatment can cause a shorter lifeexpancy

a base plot is labeled an UP plot if the first life mean is higher than the scnd life mean (and vice verca) (fst is UP)
the UP_validation plot is labeled UP if fst_life_mean is higher than the scnd life mean (and vice verca)  (fst is in therapy)
the DOWN_validation plot is labeled DOWN if the fst_life_mean is higher than the scnd life mean (and vice verca) (fst is in therapy)
"""
DF_aggr['CMP'] = 'CMP'

for plot_type in ['base_plot', 'UP_validation']:
    try:
        UP_bool = (DF_aggr.loc[(slice(None), plot_type, slice(None)), 'fst_life_mean'] > DF_aggr.loc[(slice(None), plot_type, slice(None)), 'scnd_life_mean'])
    except KeyError:
        continue
    UP_bool_index = UP_bool[UP_bool].index
    DOWN_bool_index = UP_bool[~UP_bool].index
    DF_aggr.loc[UP_bool_index,'CMP'] = 'UP'
    DF_aggr.loc[DOWN_bool_index,'CMP']= 'DOWN'


# the DOWN_validation is handled differently from the both before, if the fst_life_mean is higher, set them to DOWN
for plot_type in ['DOWN_validation']:
    try:
        DOWN_bool = (DF_aggr.loc[(slice(None), plot_type, slice(None)), 'fst_life_mean'] > DF_aggr.loc[(slice(None), plot_type, slice(None)), 'scnd_life_mean'])
    except KeyError:
        continue
    DOWN_bool_index = DOWN_bool[DOWN_bool].index
    UP_bool_index = DOWN_bool[~DOWN_bool].index
    DF_aggr.loc[DOWN_bool_index,'CMP'] = 'DOWN'
    DF_aggr.loc[UP_bool_index,'CMP']= 'UP'

# plot categorisation is done, for all plottypes, now give a rank fct::

# UP:   baseplot p_val down  -> UP_val_plot p_val_down   -> DOWN_val_down
            # p_val          +    p_val                 +   not of interest
# DOWN: baseplot p_val_down  -> UP_val_plot p_val_high   -> DOWN_val_down p_val_down
            # p_val          +    not of interest       +   p_val
# -> factors that can be included are:
# p_values of interest, but also the diff of the life_means ? what would be
# interesting is, whether or not an treatment can cause a shorter lifeexpancy


DF_aggr = DF_aggr.set_index('CMP', append=True)
# DF_aggr['p_sum'] = 0
# DF_aggr['p_prod'] = 0
DF_aggr['scored'] = False

# depending on UP or DOWN in base_plot, we score either UP or DOWN in
# validation, so out of 3 available plots, we just score 2 of them. the ones
# scored get the scored col set to True
for chr_start in DF_aggr.reset_index()['chr_start'].value_counts().index.to_list():
    for UP_or_DOWN in ['UP', 'DOWN']:
        for thresh in thresholds:
            try:
                if len(DF_aggr.loc[(chr_start, ['base_plot', f'{UP_or_DOWN}_validation'], thresh, UP_or_DOWN),:]) == 2:
                    # sum_ =  DF_aggr.loc[(chr_start, ['base_plot', f'{UP_or_DOWN}_validation'], thresh, UP_or_DOWN), 'p_value'].sum()
                    # DF_aggr.loc[(chr_start, ['base_plot', f'{UP_or_DOWN}_validation'], thresh, UP_or_DOWN), 'p_sum'] = sum_
                    # prod =  DF_aggr.loc[(chr_start, ['base_plot', f'{UP_or_DOWN}_validation'], thresh, UP_or_DOWN), 'p_value'].prod()
                    # DF_aggr.loc[(chr_start, ['base_plot', f'{UP_or_DOWN}_validation'], thresh, UP_or_DOWN), 'p_prod'] = prod
                    DF_aggr.loc[(chr_start, ['base_plot', f'{UP_or_DOWN}_validation'], thresh, UP_or_DOWN), 'scored'] = True
                else:
                    continue
            except Exception as e:
                # print(e)
                continue

# drop every p_sum with value of 0
# DF_aggr = DF_aggr[DF_aggr['p_sum'] != 0] # maybe those 3er sets are still
# needed, do that filtering in the evaluation step

DF_met = pd.read_table(metilene_out, header=None)
DF_met.columns = ['chr', 'start', 'stop', 'q-value', 'mean_methylation_difference', '#CpGs', 'mean_alive', 'mean_dead']
DF_met['DMR'] = DF_met['chr'] + '_' + DF_met['start'].astype('str') + '_' + DF_met['stop'].astype('str')

# get the DMRs out of the absolute filepath:
DF_aggr['DMR'] = DF_aggr['file_path'].apply(lambda x: re.search('chr\w*_\d*_\d*', x).group(0).strip('_'))
# add met out information based on shared DMR col:
DF_aggr = DF_aggr.reset_index().merge(DF_met.iloc[:,-6:], how='left', on='DMR')
# for not so obvious reasons, there are nans in mean_dead, this can be fixed
#  mmd = mean_alive - mean_dead | + mean_dead
#  mmd + mean_dead = mean_alive | - mmd
#  mean_dead = mean_alive - mmd
# now, recalculate the mean_dead based on mean_methylation_difference an
DF_aggr['mean_dead'] = DF_aggr['mean_alive'] - DF_aggr['mean_methylation_difference']
# also include the life_mean_diff
DF_aggr['life_mean_diff'] = DF_aggr['fst_life_mean'] - DF_aggr['scnd_life_mean']
# reaarange the cols:
DF_aggr = DF_aggr.loc[:, ['chr_start', 'ENSG', 'plot_type', 'threshold', 'CMP', 'p_value', 'fst_life_mean', 'scnd_life_mean','life_mean_diff', 'file_path', 'scored', 'DMR', 'q-value', '#CpGs', 'mean_alive', 'mean_dead', 'mean_methylation_difference']]
DF_aggr = DF_aggr.rename({'p_value': 'p_value_life'}, axis=1)
print(f'saving {metilene_lifeline_aggregated}')
DF_aggr.to_csv(metilene_lifeline_aggregated, sep='\t', index=None)
