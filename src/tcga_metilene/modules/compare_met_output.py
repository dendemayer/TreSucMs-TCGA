#!/usr/bin/env python3.7
"""
comparing:
/scr/dings/PEVO/NEW_downloads_3/metilene_api_31_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/carboplatin_carboplatin,paclitaxel_cisplatin/both/metilene_-m_3_-M_1000_-d_0.03/metilene_qval.0.05.out
/scr/dings/PEVO/NEW_downloads_3/metilene_env_test_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/carboplatin_carboplatin,paclitaxel_cisplatin/both/metilene_-m_3_-M_1000_-d_0.03/metilene_qval.0.05.out
"""
# from venn import venn
import sys
import matplotlib.pyplot as plt
import pandas as pd
import os

breakpoint()
new_met = "/scr/dings/PEVO/NEW_downloads_3/metilene_api_31_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/carboplatin_carboplatin,paclitaxel_cisplatin/both/metilene_-m_3_-M_1000_-d_0.03/metilene_qval.0.05.out"
old_met = "/scr/dings/PEVO/NEW_downloads_3/metilene_env_test_3/TCGA-CESC_TCGA-HNSC_TCGA-LUSC/carboplatin_carboplatin,paclitaxel_cisplatin/both/metilene_-m_3_-M_1000_-d_0.03/metilene_qval.0.05.out"
DF_old = pd.read_csv(old_met, sep='\t')
DF_new = pd.read_csv(old_met, sep='\t')
# venn(set_dict, fmt='{percentage:.1f}%')


