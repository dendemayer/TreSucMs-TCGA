import os


def return_summary_tables(OUTPUT_PATH, PROJECT, DRUGS, cutoffs):
    """
    on the basis of the flags, DRUGS and PROJECT the summary_tables for
    metilene are returned with this fct:
    {OUTPUT_PATH}/TCGA-CESC/metilene/metilene_input_table/carboplatn,paclitaxel_cisplatin/female/cutoff_0/threshold_0
    {OUTPUT_PATH}/TCGA-CESC/metilene/metilene_input_table/carboplatn,paclitaxel_cisplatin/male/cutoff_0/threshold_0
    {OUTPUT_PATH}/TCGA-CESC/metilene/metilene_input_table/carboplatn,paclitaxel_cisplatin/femlale_male/cutoff_0/threshold_0
    {OUTPUT_PATH}/TCGA-HNSC/metilene/metilene_input_table/carboplatn,paclitaxel_cisplatin/female/cutoff_0/threshold_0
    {OUTPUT_PATH}/TCGA-HNSC/metilene/metilene_input_table/carboplatn,paclitaxel_cisplatin/male/cutoff_0/threshold_0
    {OUTPUT_PATH}/TCGA-HNSC/metilene/metilene_input_table/carboplatn,paclitaxel_cisplatin/femlale_male/cutoff_0/threshold_0
    {OUTPUT_PATH}/TCGA-CESC_TCGA-HNSC/metilene/metilene_input_table/carboplatn,paclitaxel_cisplatin/female/cutoff_0/threshold_0
    {OUTPUT_PATH}/TCGA-CESC_TCGA-HNSC/metilene/metilene_input_table/carboplatn,paclitaxel_cisplatin/male/cutoff_0/threshold_0
    {OUTPUT_PATH}/TCGA-CESC_TCGA-HNSC/metilene/metilene_input_table/carboplatn,paclitaxel_cisplatin/femlale_male/cutoff_0/threshold_0
    """
    drugs = '_'.join(DRUGS)
    gender_list = ['female', 'male', 'female_male']
    if len(PROJECT) > 1:
        projects = '_'.join(PROJECT)
        PROJECT.append(projects)
    # cutoff = 'cutoff_' + str(cutoff)
    # threshold = 'threshold_' + str(threshold)

    temp_list = []
    for project in PROJECT:  # PROJECT contains alread the single proj and the
        # combination out of them (TCGA-CESC_TCGA-HNSC)
        for gender in gender_list:
            for cutoff in cutoffs:
                cutoff = 'cutoff_' + str(cutoff)
                temp_list.append(
                    os.path.join(
                        OUTPUT_PATH, project, 'metilene', 'metilene_input_table',
                        drugs, gender, cutoff,
                        'summary_for_metilene.tsv.gz'))

    return temp_list
