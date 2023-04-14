import os
import copy


def return_metilene_tables(OUTPUT_PATH, PROJECTS, DRUGS, cutoffs):
    """
    the 3 parameters dont have to be set explisitly, take an inclusive
    paramter set and filter later on those cols

    metilene -m 3 -M 1000 -d 0.03 -a alive -b dead summary_for_metilene.tsv | sort -V -k1,1 -k2,2n > metilene_out_sorted-m3_-M1000-d_0.03_.tsv
    perl ~/phd/test_git_doc/tcga_piplines/src/tcga_metilene/resources/metilene_output.pl -q metilene_out_sorted-m3_-M1000-d_0.03.tsv -a 'alive' -b 'dead' -c 3 -d 0.03 -o metilene_qval.0.05-m3_-M1000-d_0.03.tsv
    """
    drugs = '_'.join(DRUGS)
    gender_list = ['female', 'male', 'female_male']
    # PROJECT_temp = copy.deepcopy(PROJECT)
    # if len(PROJECT) > 1:
    #     projects = '_'.join(PROJECT_temp)
    #     PROJECT_temp.append(projects)
    temp_list_met = []
    temp_list_met_filtered = []
    for project in PROJECTS:  # PROJECTS contains alread the single proj and the
        # combination out of them (TCGA-CESC_TCGA-HNSC)
        for gender in gender_list:
            for cutoff in cutoffs:
                cutoff = 'cutoff_' + str(cutoff)
                temp_list_met.append(
                    os.path.join(
                        OUTPUT_PATH, project, 'metilene',
                        'metilene_output', drugs, gender, cutoff,
                        'metilene_out_sorted.tsv'))
                        # f'metilene_out_sorted-m_{str(m)}-M_{str(M)}-d_{str(d)}.tsv'))
                temp_list_met_filtered.append(  # metilene-m_3-M_1000-d_0.03_qval.0.05.out
                    os.path.join(
                        OUTPUT_PATH, project, 'metilene',
                        'metilene_output', drugs, gender, cutoff,
                        'metilene_qval.0.05.out'))
                        # f'metilene-m_{str(m)}-M_{str(M)}-d_{str(d)}_qval.0.05.out'))
    return temp_list_met + temp_list_met_filtered
