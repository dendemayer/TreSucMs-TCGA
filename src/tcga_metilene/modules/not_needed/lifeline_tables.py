import pandas as pd
import os
import glob
import re
# from methyl import set_logger
# from lifelines.datasets import load_waltons
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from lifelines import KaplanMeierFitter
from decimal import Decimal
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
import sys

"""
follow lifeline succession from patients with detected DMR
- needed, metadata:f.e.:
    - TCGA-pipelines_2/TCGA-CESC_TCGA-HNSC/metilene/merged_meta_files/meta_info_druglist_merged_drugs_combined.tsv
    - cols:
        # survivaltime       | years_to_last_follow_up 'bcr_patient_uuid',
- the metilene intersect out table
    - TCGA-pipelines_2/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_complement_intersect.tsv
    # needed infos are vital_status (alive dead),

# use the logrank_test at every position within a DMR to select the best
# resolution between up and down :
If the p-value of the log-rank test is less than a chosen significance level
(such as 0.05), it suggests that there is a statistically significant
difference in survival between the groups, and the Kaplan-Meier estimate can be
considered meaningful.
"""
breakpoint()

meta_table = sys.argv[1]
met_table = sys.argv[2]
DMR = 'chr1_78045659_78046029'
out_lifeline_pdf = '/scr/dings/PEVO/NEW_downloads_3/TCGA-pipelines_2/TCGA-CESC_TCGA-HNSC/metilene/metilene_output/carboplatin,paclitaxel_cisplatin/female/cutoff_0/metilene_complement_intersect_lifeline_chr1_78045659_78046029.pdf'

DF_meta = pd.read_table(meta_table, usecols=['bcr_patient_uuid', 'survivaltime', 'years_to_last_follow_up', 'vital_status'])
DF_metilene = pd.read_table(met_table, header=[0,1,2,3,4], index_col=[0,1,2,3])
# limit the DF_metilene to the recent range:
DF_metilene = DF_metilene.loc[(slice(None), slice(None), slice(None), DMR), :]

# get the used cases out of the meta_table:
used_cases = [ i[1] for i in DF_metilene.columns]
DF_meta = DF_meta.set_index('bcr_patient_uuid').loc[used_cases, :]
DF_meta['T'] = DF_meta['survivaltime'].combine_first(DF_meta['years_to_last_follow_up'])

T = DF_meta['T']
E = DF_meta['vital_status'] == 'dead' #  E is a either boolean or binary array
# representing whether the “death” was observed or not (alternatively an
# individual can be censored)

# the groups are representing the cases which are either higher or lower
# methylated than the median of all cases within the DMR, therefore calculate
# the kmf at every postition and take the one with the highest
starts = [i[1] for i in DF_metilene.index]
# pvalue
pvalue_list = []
for start in starts:
    median = DF_metilene.loc[(slice(None), start, slice(None), slice(None)), :].median(axis=1).values[0]
    DF_temp = DF_metilene.loc[(slice(None), start, slice(None), slice(None)), :].apply(lambda x: 'UP' if float(x) > median else 'DOWN' )
    ix = (DF_temp == 'UP').values
    results = logrank_test(T[ix], E[~ix], T[ix], E[~ix])
    p_value = results.p_value
    pvalue_list.append(p_value)

# take the lowest p_value,

p_value_str = f'p_value = {Decimal(str(min(pvalue_list))):.2e}'
plot_index = pvalue_list.index(min(pvalue_list))
median = DF_metilene.iloc[plot_index, :].median()
DF_plot = DF_metilene.iloc[plot_index, :].apply(lambda x: 'UP' if float(x) > median else 'DOWN' )
DF_plot.name = 'UP_DOWN'

# kmf.fit(T, event_observed=E)  # or, more succinctly, kmf.fit(T, E)

# groups = DF_plot['UP_DOWN']
ix = (DF_plot == 'UP').values
fig, ax = plt.subplots(figsize=(8, 6))

kmf_UP = KaplanMeierFitter()
kmf_UP.fit(T[ix], E[ix], label='UP')
ax = kmf_UP.plot_survival_function(ax=ax)

kmf_DOWN = KaplanMeierFitter()
kmf_DOWN.fit(T[~ix], E[~ix], label='DOWN')
kmf_DOWN.plot_survival_function(ax=ax)

add_at_risk_counts(kmf_UP, kmf_DOWN, ax=ax)
ax.set_title(f'{p_value_str} \nDMR: {": ".join(DMR.split("_")[:2]) + "-" + DMR.split("_")[2]}\nStart: {str(starts[plot_index])}')

plt.tight_layout()
print(f'saving: {out_lifeline_pdf}')
plt.savefig(out_lifeline_pdf)
plt.close()


# to compare two survivalfunctions and get CoxPHFitter is actually
# better than logrank_test, see: https://discourse.datamethods.org/t/
# when-is-log-rank-preferred-over-univariable-cox-regression/2344
# but in our case we dont want to invoke the age (a covariate var is needed)
# since the KaplanMeier depends just on survivaltime[T] and event [E]
# from lifelines import *  # lifelines has builtin parametric models. for


def prepare_complement_lifeline_plots(OUTPUT_PATH,
                                      DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID,
                                      met_dir, cutoff, DRUGS_title):
    '''
    create the complement lifeline_table like:
    lifeline_cg00414041_166530702-166530887_chr6_166530795_CAND.tsv
    for patients not treated with DRUGS,
    lifeline_cg00414041_166530702-166530887_chr6_166530795_CAND_complement.tsv
    iss then visualized by fct: def create_lifeline_plots(
    through globbing:
    glob.glob(os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title, gender,
    met_dir, 'lifeline_out', '*.tsv')):
    -> do not include previously created compelement tsv!
    '''
    print(f'preparing complement lifeline_plots for {PROJECT_DRUG_UUID}')
    if isinstance(PROJECT_DRUG_UUID, dict):  # a multi project as dict
        projects = []
        for project in PROJECT_DRUG_UUID:
            projects.append(project)
        PROJECT_title = '_'.join(sorted(map(str.upper, projects)))
    else:  # single projects are str and can be set directly
        PROJECT_title = PROJECT_DRUG_UUID
        projects = [PROJECT_DRUG_UUID]

    logger = set_logger.set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title)
    gender_list_temp = ['both', 'female', 'male']
    gender_list = []
    for gender in gender_list_temp:
        if os.path.exists(os.path.join(OUTPUT_PATH, PROJECT_title,
                                       DRUGS_title, gender)):
            gender_list.append(gender)
    # ##### BEGIN copying and concatenating all meta infos from single projects
    # # if a multi project is applied
    # to the combined projects:, f.e.
    # TCGA-CESC/meta_info_druglist_merged_drugs_combined.tsv
    # TCGA-HNSC/meta_info_druglist_merged_drugs_combined.tsv
    # to TCGA-CESC_TCGA-HNSC/meta_info_druglist_merged_drugs_combined.tsv
    if isinstance(PROJECT_DRUG_UUID, str):
        path_to_meta = os.path.join(
            OUTPUT_PATH, PROJECT_title,
            'meta_info_druglist_merged_drugs_combined.tsv')
        meta_DF = pd.read_csv(
            path_to_meta, sep='\t',
            na_values=['[Not Applicable]', '[Not Available]',
                       '[not available]', '[unknown]'])
    else:
        first = True
        for project in projects:
            if first:
                path_to_meta = os.path.join(
                    OUTPUT_PATH, project,
                    'meta_info_druglist_merged_drugs_combined.tsv')
                meta_DF = pd.read_csv(
                    path_to_meta, sep='\t',
                    na_values=['[Not Applicable]', '[Not Available]',
                               '[not available]', '[unknown]'])
                first = False
            else:
                path_to_meta = os.path.join(
                    OUTPUT_PATH, project,
                    'meta_info_druglist_merged_drugs_combined.tsv')
                meta_DF = pd.concat(
                    [meta_DF, pd.read_csv(
                        path_to_meta, sep='\t',
                        na_values=['[Not Applicable]', '[Not Available]',
                                   '[not available]', '[unknown]'])])

    path_to_meta = os.path.join(OUTPUT_PATH, PROJECT_title,
                                'meta_info_druglist_merged_drugs_combined.tsv')
    log_name = os.path.join(PROJECT_title,
                            'meta_info_druglist_merged_drugs_combined.tsv')
    try:
        meta_DF.to_csv(path_to_meta, sep='\t', index=False)
    except Exception as e:
        print('The following error occured while trying to write the\
                file: {}: {}'.format(path_to_meta, e))
    else:
        logger.info('REPORT_17:\t{}'.format(log_name))
    # ##### END copying and concatenating all meta infos from single projects
    # if a multiproject is applied

    # ######## BEGIN plot content of this summarized table like how many cases
    # are dead alive per project...
    # cols of interest:
    # project_id , birth_days_to, gender, vital_status, last_contact_days_to,
    # death_days_to
    meta_plot_DF = meta_DF.loc[
        :, [
            "project_id", "gender", "vital_status", "birth_days_to",
            "last_contact_days_to", "death_days_to"]]
    # before transforming days into years, make sure to convert str to numeric
    # where possible and set nan where not possible:
    # errors : {'ignore', 'raise', 'coerce'}, default 'raise'
    #     - If 'raise', then invalid parsing will raise an exception.
    #     - If 'coerce', then invalid parsing will be set as NaN.
    #     - If 'ignore', then invalid parsing will return the input.
    # meta_plot_DF['birth_days_to'] = pd.to_numeric(
    #     meta_plot_DF['birth_days_to'], errors='coerce')
    # # than do the actual change:
    # meta_plot_DF['birth_days_to'] = meta_plot_DF[
    #     'birth_days_to'].apply(lambda x: x/-365)
    for i in ['last_contact_days_to', 'death_days_to']:
        meta_plot_DF[i] = pd.to_numeric(meta_plot_DF[i], errors='coerce')
        meta_plot_DF[i] = meta_plot_DF[i].apply(lambda x: x/365)
    # now meta_plot['last_contact_days_to'] and meta_plot['death_days_to']
    # holds the YEAR
    # TODO, save those plots ...
    # plt.scatter( meta_plot_DF[ 'last_contact_days_to'],
    # meta_plot_DF['death_days_to'], s=meta_plot_DF['birth_days_to'])

    # sns.boxplot( x="project_id", y="birth_days_to", hue="vital_status", palette=["m", "g"], data=meta_plot_DF)
    # sns.set_theme(style="ticks", palette="pastel")
    # sns.despine(offset=10, trim=True)
    # plt.show()

    # sns.countplot(x="project_id", hue="gender", data=meta_plot_DF)
    # sns.catplot( x="project_id", hue="gender", col="vital_status",
    #             data=meta_plot_DF, kind="count", height=4, aspect=.7)
    # sns.catplot( x="project_id", hue="vital_status", col="gender",
    #             data=meta_plot_DF, kind="count")

    # ######## END plot content of this summarized table like how many cases
    # are dead alive per project...

    # ######## BEGIN collecting case_id to survivaltimes:
    DF_surv_time = meta_DF.loc[:, ['death_days_to', 'last_contact_days_to',
                                   'pharmaceutical_therapy_drug_name',
                                   'vital_status', 'gender', 'case_ids']]
    DF_surv_time.rename({'death_days_to': 'years_to_death',
                         'last_contact_days_to': 'years_to_last_followup',
                         'case_ids': 'case_id'},
                        inplace=True, axis=1)
    DF_surv_time['years_to_death'] = pd.to_numeric(
        DF_surv_time['years_to_death'], errors='coerce')
    DF_surv_time['years_to_last_followup'] = pd.to_numeric(
        DF_surv_time['years_to_last_followup'], errors='coerce')
    DF_surv_time['years_to_last_followup'] = DF_surv_time[
        'years_to_last_followup'] / 365
    DF_surv_time['years_to_death'] = DF_surv_time[
        'years_to_death'] / 365
    DF_surv_time.set_index('case_id', inplace=True)
    # create the T col with merging the both cols together:
    DF_surv_time['T'] = DF_surv_time[
            'years_to_last_followup'].combine_first(
                DF_surv_time['years_to_death'])

    # ######## END collecting case_id to survivaltimes:

    # ######### BEGIN get the case_id's which belong to the complement DRUG
    # selection:
    # TCGA-CESC_TCGA-HNSC/
    # carboplatin_carboplatin,paclitaxel_cisplatin/
    # both/summary_dead_alive_dropped_NA_dropped_REF_sorted_complement.tsv
    case_DF = pd.DataFrame()
    for gender in gender_list:
        case_path = os.path.join(
            OUTPUT_PATH, PROJECT_title, DRUGS_title, gender,
            'summary_dead_alive_dropped_NA_dropped_REF_sorted_complement.tsv')
        case_list = pd.read_csv(
            case_path, sep='\t', nrows=1).columns.to_list()[2:]
        case_list = [i.split(';')[1] for i in case_list]
        case_DF_temp = pd.DataFrame({'bcr_patient_uuid': case_list})
        case_DF_temp['gender'] = gender
        case_DF = pd.concat([case_DF, case_DF_temp])

    # add the beta-value filename associated with the case_id:
    case_DF = case_DF.merge(
        meta_DF.loc[:, [
            'bcr_patient_uuid', 'filename', 'project_id']],
        on='bcr_patient_uuid')

    def add_full_path(project_id_file_name):
        project_id = project_id_file_name['project_id']
        filename = project_id_file_name['filename']
        path = os.path.join(
            OUTPUT_PATH, project_id, project_id + '_data_files', filename)
        return glob.glob(path)[0]
    case_DF['full_path'] = case_DF.loc[:, ['project_id', 'filename']].apply(
        add_full_path, axis=1)

    # ######### END get the case_id which belong to the complement DRUG

    # #### BEGIN glob the cg flags needed:
    # ### glob every already created lifeline candidate table and parse
    # cg00414041:
    # for every gender in PROJECT_title glob the lifelinetables:
    # the median is given in those tables plus the cg flag, with that the up
    # and down grouping must be build per case_id:
    # -> first parse every cg, this is already presented in the filename:
    # carefull not to read in already created complement tables leading to:
    # lifeline_cg00691999_54226195-54232316_
    # chr4_54227844_CAND_complement_complement.tsv
    life_table_list = []
    life_table_list_gender = []
    for gender in gender_list:
        path_to_life = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                    gender, met_dir, 'lifeline_out')
        temp_list = glob.glob(f'{path_to_life}/*.tsv')
        # exclude files containing complement, those are already created
        # beforehand:
        elements_to_be_rm = []
        for i in temp_list:
            try:
                re.search(r'.*complement.*', i).group()
                elements_to_be_rm.append(i)
            except Exception as e:
                # print(f'continuing because of {e}')
                continue
        temp_list = list(set(temp_list) - set(elements_to_be_rm))
        life_table_list.extend(temp_list)
        life_table_list_gender.extend(np.repeat(gender, len(temp_list)))

    life_table_list_complement = [
        i.replace('.tsv', '_complement.tsv') for i in life_table_list]
    # build a list of equal length which represents the cg flag for every
    # table:
    cg_list = [re.search(r'cg\d*', i).group() for i in life_table_list]
    median_list = [pd.read_csv(i, sep='\t', nrows=1)['median'].values[0] for i
                   in life_table_list]
    DMR_list = [pd.read_csv(i, sep='\t', nrows=1)['DMR'].values[0] for i in
                life_table_list]
    Start_list = [pd.read_csv(i, sep='\t', nrows=1)['Start'].values[0] for i in
                  life_table_list]
    chr_list = [pd.read_csv(i, sep='\t', nrows=1)[
        'chr'].values[0] for i in life_table_list]
    # for life_table in life_table_list:
    #     pd.read_csv(life_table, sep='\t'['case_id'])

    life_table_DF = pd.DataFrame({'life_table': life_table_list,
                                  'life_table_list_complement':
                                  life_table_list_complement,
                                  'cg_flag': cg_list,
                                  'median': median_list,
                                  'DMR': DMR_list,
                                  'Start': Start_list,
                                  'chr': chr_list,
                                  'gender': life_table_list_gender})

    # for each table in life_table_DF, we need the up or down info for every
    # case in case_DF and for every cg_ flag within the life_table_DF
    # gender_list = life_table_DF['gender'].values.tolist()
    cases_to_check = case_DF['bcr_patient_uuid'].value_counts().index.to_list()
    # for life_table in life_table_list_complement:
    #     case_up_down.update({life_table: []})
    #     get the right file to the case_id:
    case_up_down = {}
    # case_beta_val = {}
    for case_id in cases_to_check:
        case_up_down.update({case_id: []})
        # case_beta_val.update({case_id: []})
        # refer to the cg_flag and median out of life_table:
        methyl_path = case_DF.set_index('bcr_patient_uuid').drop_duplicates(
            'filename').loc[case_id]['full_path']
        # duplicates can occur, because the same case can be present in the
        # f.e. male dir but also in the both dir
        beta_val_DF = pd.read_csv(
            methyl_path, sep='\t', na_values='NaN', usecols=[
                'Composite Element REF', 'Beta_value'],
            index_col='Composite Element REF')
        for cg, median in zip(cg_list, median_list):
            beta_val = beta_val_DF.loc[cg, 'Beta_value']
            # keep that file open once, and select out all cg of interest:
            # beta_val = pd.read_csv(
            #     methyl_path, sep='\t', na_values='NaN', usecols=[
            #         'Composite Element REF', 'Beta_value'],
            #     index_col='Composite Element REF').loc[cg, 'Beta_value']
            if beta_val > median:
                group = 'up'
            else:
                group = 'down'
            case_up_down[case_id].append([group, beta_val, cg])
    temp_DF = pd.DataFrame(case_up_down).T
    # no multi_index needed, the cg_name is sufficient to access via
    # life_table_DF infos like median, Start, chr, DMR
    # multi_index = pd.MultiIndex.from_arrays([cg_list, median_list],
    # names=('cg_name', 'median'))
    temp_DF.columns = cg_list
    case_DF = case_DF.set_index('bcr_patient_uuid')
    case_DF = case_DF.join(temp_DF)
    life_table_DF = life_table_DF.set_index('life_table_list_complement')
    # add the drugs to the case_DF
    temp_DF = meta_DF.set_index(
        'bcr_patient_uuid')['pharmaceutical_therapy_drug_name']
    case_DF = case_DF.join(temp_DF)
    case_DF.rename({'pharmaceutical_therapy_drug_name': 'drug_name'},
                   inplace=True, axis=1)
    life_complement_list = life_table_DF.index.tolist()
    gender_list = life_table_DF['gender'].to_list()
    for life_table in life_complement_list:
        start_DF = pd.DataFrame()
        for case_id in set(case_DF.index.to_list()):
            drug = DF_surv_time.loc[case_id,
                                    'pharmaceutical_therapy_drug_name']
            gender = meta_DF.set_index('bcr_patient_uuid').loc[case_id,
                                                               'gender']
            project = meta_DF.set_index('bcr_patient_uuid').loc[case_id,
                                                                'project_id']
            vital_state = DF_surv_time.loc[case_id, 'vital_status']
            T = DF_surv_time.loc[case_id, 'T']
            cg = life_table_DF.loc[life_table, 'cg_flag']
            if isinstance(case_DF.loc[case_id, cg], pd.Series):
                group = case_DF.loc[case_id, cg].iloc[0][0]
            elif isinstance(case_DF.loc[case_id, cg], list):
                group = case_DF.loc[case_id, cg][0]
                # (Pdb) case_DF.loc[case_id, cg]
                # ['down', 0.418406321868661, 'cg19140375']
            else:
                group = case_DF.loc[case_id, cg].iloc[0, 0][0]
            chrom = life_table_DF.loc[life_table, 'chr']
            start = life_table_DF.loc[life_table, 'Start']
            dmr = life_table_DF.loc[life_table, 'DMR']
            median = life_table_DF.loc[life_table, 'median']
            temp_DF = pd.DataFrame({'case_id': [case_id], 'drug': [drug],
                                    'gender': [gender], 'Project': [project],
                                    'vital_state': [vital_state],
                                    'T': [T], 'group': [group], 'chr': [chrom],
                                    'cg_name': [cg], 'Start': [start],
                                    'DMR': [dmr], 'median': [median]})
            start_DF = pd.concat([start_DF, temp_DF])
        log_name = life_table.replace(os.path.join(
            OUTPUT_PATH, os.path.sep), '')
        try:
            start_DF.to_csv(life_table, index=False, sep='\t')
        except Exception as e:
            print(f'The following error occured while trying to write the\
                    file: {life_table}: {e}')
        else:
            logger.info(f'REPORT_17:\t{log_name}')

    # case_id specific cols:
    # -> case_id  drug gender Project vital_state T group
    #
    # life_line_table specific cols: (the can be added after case cols are
    # build)
    # ->  chr cg_name Start DMR median

    # for every life_table, write the case_id (cases_to_check)

    # # puzzle the DF newly together:

    # we need:(lifeline_cg00414041_166530702-166530887
    # _chr6_166530795_CAND_complement.tsv)
    # _case_id_drug____gender__Project_vital_state_____T_______group___chr_____cg_name_Start___DMR_____median
    # _91310434-2a99-4e1d-81e6-3d58d097089f____cisplatin_______female__TCGA-CESC_______alive___3.0465753424657533______up______chr6____cg00414041______166530795_______166530702-166530887_____0.555592401269022
    # _3143822c-66c4-4b58-8d29-e92d065ff049____cisplatin_______female__TCGA-CESC_______dead____3.249315068493151_______up______chr6____cg00414041______166530795_______166530702-166530887_____0.555592401269022


def prepare_lifeline_plots(OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID,
                           met_dir, cutoff, DRUGS_title):
    # we search in OUTPUT_PATH, PROJECT, DRUGS, [genders],
    # intersect_header.
    # in which genes are we interested in? check all the cand CAND flags of
    # belonging genders...
    # we have several positions within a range, each position corresponds to a
    # gene name
    # out of a table like:
    # boxplot_means_range_170296969-170297418_with_positions_CAND_
    # linegraph.tsv
    # how to evaluate which position/gene to lifeline plot?
    # -> for every position within a range add up the median of means beta
    # value from one vital_state
    # compare those vital states, the one with the highest differences is going
    # to be plottet as lifelines
    # with position found, go into the beta value table and plot those values
    # there for the right postion -> intersect_header.tsv, intersect_reads.tsv

    # populate the cand s.t. hash:
    print(f'preparing lifeline_plots for {PROJECT_DRUG_UUID}')

    if isinstance(PROJECT_DRUG_UUID, dict):  # a multi project as dict
        projects = []
        for project in PROJECT_DRUG_UUID:
            projects.append(project)
        PROJECT_title = '_'.join(sorted(map(str.upper, projects)))
    else:  # single projects are str and can be set directly
        PROJECT_title = PROJECT_DRUG_UUID
        projects = [PROJECT_DRUG_UUID]

    gender_list_temp = ['both', 'female', 'male']
    gender_list = []
    for gender in gender_list_temp:
        if os.path.exists(os.path.join(OUTPUT_PATH, PROJECT_title,
                                       DRUGS_title, gender)):
            gender_list.append(gender)
    # #######################################################################
    # begin collecting case_id to survivaltimes:
    # survivaltime to the case_id:, can be found in the
    # single project dirs in meta_info_druglist_merged_drugs_combined.tsv
    # cols:
    # death_days_to | cases.0.diagnoses.0.days_to_last_follow_up | id
    DF_list = []
    for i in projects:
        path_to_table = os.path.join(
            OUTPUT_PATH, i, 'meta_info_druglist_merged_drugs_combined.tsv')
        temp_DF = pd.read_csv(
            path_to_table, sep='\t', usecols=[
                'death_days_to', 'last_contact_days_to',
                'case_ids'], na_values='[Not Applicable]')
        temp_DF.rename({'death_days_to': 'years_to_death',
                        'last_contact_days_to':
                        'years_to_last_followup', 'case_ids': 'case_id'},
                       axis=1, inplace=True)
        DF_list.append(temp_DF)
    DF_surv_time = pd.concat(DF_list).reset_index(drop=True)
    # in case string values which deviate from [Not Applicable], cast the cols
    # to numeric and if that raises an error, put NaN there:
    # in example str values like '1' can be casted to numeric 1, but '[Not
    # Applicable]' can not be castet to numeric values and are casted to NA
    # (float or np... or whatever)
    DF_surv_time['years_to_death'] = pd.to_numeric(
        DF_surv_time['years_to_death'], errors='coerce')
    DF_surv_time['years_to_last_followup'] = pd.to_numeric(
        DF_surv_time['years_to_last_followup'], errors='coerce')
    DF_surv_time['years_to_last_followup'] = DF_surv_time[
        'years_to_last_followup'] / 365
    DF_surv_time['years_to_death'] = DF_surv_time[
        'years_to_death'] / 365
    DF_surv_time.set_index('case_id', inplace=True)
    # end collecting case_id to survivaltimes:
    # #######################################################################

    # #######################################################################
    # begin candidate ranges , check every postion within the DMR range, how
    # different are the mean of medians at every position there?
    cand_list = ['CAND', 'cand']
    # save the coresponding range according to the chosen postion out of a DMR
    # final dict form:
    # {(Start, chr, dmr) : {diff: None, cand_type: None}}
    final_dict = {}
    for gender in gender_list:
        for cand in cand_list:
            glob_path = os.path.join(
                OUTPUT_PATH, PROJECT_title, DRUGS_title, gender, met_dir,
                'boxplot_range', '*' + cand + '*_linegraph.tsv')
            files_found = glob.glob(glob_path)
            for table in files_found:
                # also parse out the DMR out of the filename table:
                pattern = re.compile(r'\d*-\d*')
                dmr = pattern.search(os.path.basename(table)).group(0)
                # ##
                DF_means_range = pd.read_csv(table, sep='\t')
                pos_list = DF_means_range['Start'].unique()
                # reindex to choose dependent on vital_state and Start
                DF_means_range.set_index(['vital_state', 'Start'],
                                         inplace=True)
                # DONE how to aviod the warning:
                # PerformanceWarning: indexing past lexsort depth may impact performance. [lifeline_plots.py:409]
                # PerformanceWarning: indexing past lexsort depth may impact performance. [lifeline_plots.py:412]
                # see modules/lexsort_problem.py, sort_index() is sufficient
                # DF_means_range.index.is_lexsorted()
                # <stdin>:1: FutureWarning: MultiIndex.is_lexsorted is deprecated as a public function, users should use MultiIndex.is_monotonic_increasing instead.
                DF_means_range = DF_means_range.sort_index()
                # DF_means_range.index.is_monotonic_increasing
                # True
                diff_list = []
                for position in pos_list:
                    alive_sum = DF_means_range.loc[
                        ('alive', position), :][
                                'median of means beta value'].sum()
                    dead_sum = DF_means_range.loc[
                        ('dead', position), :][
                            'median of means beta value'].sum()
                    diff = abs(alive_sum - dead_sum)
                    diff_list.append(diff)
                    # update the final_dict:
                chromosome = DF_means_range['chr'].unique()[0]
                # to relate the right position to the chosen diff, get the
                # index position of the max diff_list:
                pos_ind = diff_list.index(max(diff_list))
                final_dict.update(
                    {(chromosome, pos_list[pos_ind], dmr, gender):
                        {'diff': max(diff_list), 'cand_type': cand}})
    DF_cand_final = pd.DataFrame().from_dict(final_dict, orient='index')
    DF_cand_final.index.name = ('chr', 'Start', 'DMR', 'gender')

    # end candidate ranges , check every postion within the DMR range, how
    # #######################################################################

    # #######################################################################
    # ##  begin connect the cand positions found to the beta_values to the
    # respective position:
    # by going through intersect_header, we limit the set automatically
    # to the right drug selection

    DF_list = []
    for gender in gender_list:
        path_to_header = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                      gender, 'intersect_header.tsv')
        DF_header = pd.read_csv(path_to_header, sep='\t', header=None)
        # header structure:
        # Chromosome      Start   End     REF
        # alive;f5d90810-8c4b-4ced-9360-383e42cb54b3;cisplatin;female;TCGA-CESC
        DF_header.dropna(inplace=True, axis=1)
        DF_header = DF_header.loc[:, range(4, DF_header.shape[1])]
        DF_header_list = DF_header.values.tolist()[0]
        multi_dict = {'vital_state': [], 'case_id': [], 'drug': [], 'gender':
                      [], 'Project': []}
        for i in DF_header_list:
            temp = i.split(';')
            multi_dict['vital_state'].append(temp[0])
            multi_dict['case_id'].append(temp[1])
            multi_dict['drug'].append(temp[2])
            multi_dict['gender'].append(temp[3])
            multi_dict['Project'].append(temp[4])

        # pd.MultiIndex.from_arrays(, names=()
        multi_ind = pd.MultiIndex.from_arrays([multi_dict['vital_state'],
                                               multi_dict['case_id'],
                                               multi_dict['drug'],
                                               multi_dict['gender'],
                                               multi_dict['Project']],
                                              names=('vital_state', 'case_id',
                                                     'drug', 'gender',
                                                     'Project'))
        # ################# BEGIN just read in rows which are actually needed,
        # therefore read in the first 2 cols, chr and start, then intersect
        # with chr start of found CAND positions, then skip all rows not needed
        # while reading in

        path_to_reads = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                     gender, 'intersect_reads.tsv')

        DF_chr_pos = pd.read_csv(path_to_reads, sep='\t', header=None,
                                 usecols=[0, 1], engine='c')
        DF_chr_pos['index'] = DF_chr_pos.index
        # append the cols as index to existing index, the numeric index is then
        # needed to skip rows
        DF_chr_pos.set_index([0, 1], inplace=True)
        # just need first 2 index pos from candidates, chr and start
        index_from_cand = DF_cand_final.index.droplevel([2, 3])
        # access with cand pos in complete reads df to gather needed row index:
        # turn it to list right away
        index_needed = DF_chr_pos.loc[index_from_cand, :]['index'].to_list()
        # ################# END just read in rows which are actually needed,

        DF_reads = pd.read_csv(path_to_reads, sep='\t', header=None,
                               index_col=[0, 1, 2, 3], engine='c',
                               skiprows=lambda x: x not in index_needed)

        # DF_reads = pd.read_csv(path_to_reads, sep='\t', header=None,
        #                        index_col=[0, 1, 2, 3], engine='c',
        #                        skiprows= index_needed)
        # make first 4 cols of DF_reads as index, then apply 4 level multiindex
        # cols.
        # DF_reads.set_index([0, 1, 2, 3], inplace=True)
        DF_reads.index.names = ['chr', 'Start', 'Stop', 'cg_name']
        DF_reads.columns = multi_ind
        # now intersect with the intersect_reads.tsv
        # complete DF_reads is finished ,
        # now acces the values of interest, out of
        # DF_cand_final:
        # check what kind of cands we actually have found:
        cand_list = DF_cand_final['cand_type'].unique()
        for cand in cand_list:
            DF_cand_final_temp = DF_cand_final[
                DF_cand_final['cand_type'] == cand]
            for index in DF_cand_final_temp.index:
                # important: preserve the DMR in index[2] !
                index_new = (index[0], index[1], slice(None), slice(None))
                # with index_new access directly chr and pos out of the reads
                median = DF_reads.loc[index_new, :].median(axis=1)
                DF_row_complete = DF_reads.loc[index_new, :].T
                DF_reads_up = DF_row_complete[
                    DF_row_complete > median].dropna()
                DF_reads_up['up_or_down'] = 'up'
                DF_reads_down = DF_row_complete[
                    DF_row_complete < median].dropna()
                DF_reads_down['up_or_down'] = 'down'
                DF_final = pd.concat([DF_reads_up, DF_reads_down])
                dmr = index[2]
                # !! don't add gender out of DF_final to the DF_list take
                # gender from the gender_list iterator!!
                # also add the medium on which the up and down classification
                # depends on:
                DF_final[('median', '', '', '')] = median.to_list()[0]
                DF_list.append([cand, DF_final, dmr, gender])

    # #######################################################################
    # ##  end connect the cand positions found to the beta_values to the
    # respective position:
    # DF_list -> holds type of cand and the belonging DF, merge the needed
    # survivaltimes
    # DF_surv_time -> holds case_id to survivaltimes
    # DF_list is a array with array of type: i[cand, DF_final, dmr, gender]
    # DONE! avoid warning:
    # FutureWarning: merging between different levels is deprecated and will be
    # removed in a future version. (4 levels on the left, 1 on the right)
    # [lifeline_plots.py:539]
    for i in DF_list:
        beta_temp = i[1]
        # temp save multi indexes:
        beta_index = beta_temp.index
        beta_columns = beta_temp.columns
        # make a numeric cols to join with DF_surv_time
        beta_temp.columns = range(beta_columns.size)
        # drop everyting out of row index except case_id, on which the join is
        # performed
        beta_temp = beta_temp.reset_index(['vital_state', 'drug', 'gender',
                                           'Project'], drop=True)
        DF_joined = beta_temp.join(DF_surv_time, on='case_id')
        # recreate the multiindexes, both in beta temp and the DF_surv_time
        beta_temp.index = beta_index
        beta_temp.columns = beta_columns
        DF_joined.rename({0: beta_columns[0], 1: beta_columns[1], 2:
                          beta_columns[2]}, axis=1, inplace=True)
        DF_joined.index = beta_index
        DF_joined['cand_type'] = i[0]
        DF_life = pd.DataFrame()
        DF_life['T'] = DF_joined[
            'years_to_last_followup'].combine_first(
                DF_joined['years_to_death'])
        DF_life['group'] = DF_joined[('up_or_down', '', '', '')]
        DF_life.reset_index('vital_state', inplace=True)
        # key: ('chr6', 170297418, 170297418, 'cg26575057')
        DF_life['chr'] = DF_joined.columns[0][0]
        chromosome = DF_joined.columns[0][0]
        DF_life['cg_name'] = DF_joined.columns[0][3]
        cg_name = DF_joined.columns[0][3]
        cand = DF_joined['cand_type'].unique()[0]
        DF_life['Start'] = DF_joined.columns[0][1]
        start = DF_joined.columns[0][1]
        DF_life['DMR'] = i[2]
        DMR = i[2]
        DF_life['median'] = beta_temp['median'].unique()[0]
        DF_life.reset_index(['case_id', 'drug', 'gender', 'Project'],
                            inplace=True)
        gender = i[3]
        # if len(DF_life['gender'].unique()) == 2:
        #     gender = 'both'
        # else:
        #     gender = DF_life['gender'].unique()[0]
        file_path = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                 gender, met_dir, 'lifeline_out')
        os.makedirs(file_path, exist_ok=True)
        # use the cg_name to flag the filename:
        file_name = 'lifeline_' + cg_name + '_' + DMR + '_' + chromosome +\
            '_' + str(start) + '_' + cand + '.tsv'
        try:
            DF_life.to_csv(os.path.join(file_path, file_name), sep='\t',
                           index=False)
        except Exception as e:
            print('The following error occured while trying to write the\
                  file: {}: {}'.format(os.path.join(file_path, file_name), e))
        else:
            logger = set_logger.set_logger(
                OUTPUT_PATH, PROJECT_title, DRUGS_title)
            logger.info('REPORT_17:\t{}'.format(
                os.path.join(PROJECT_title, DRUGS_title, gender, met_dir,
                             'lifeline_out', file_name)))


def create_lifeline_plots(OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID,
                          met_dir, cutoff, DRUGS_title):

    print(f'creating lifeline plots for {PROJECT_DRUG_UUID}')

    if isinstance(PROJECT_DRUG_UUID, dict):  # a multi project as dict
        projects = []
        for project in PROJECT_DRUG_UUID:
            projects.append(project)
        PROJECT_title = '_'.join(sorted(map(str.upper, projects)))
    else:
        PROJECT_title = PROJECT_DRUG_UUID
        projects = [PROJECT_DRUG_UUID]

    gender_list_temp = ['both', 'female', 'male']
    gender_list = []
    for gender in gender_list_temp:
        if os.path.exists(os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                       gender, met_dir, 'lifeline_out')):
            gender_list.append(gender)
    DF_list = []
    file_name_list = []
    for gender in gender_list:
        for i in glob.glob(os.path.join(OUTPUT_PATH, PROJECT_title,
                                        DRUGS_title, gender, met_dir,
                                        'lifeline_out', '*.tsv')):
            # IMPORTANT, glob also the com
            # check here whether there is actually anything in the DF to
            # analyse:
            DF_list.append(pd.read_csv(i, sep='\t'))
            file_name_list.append(i)

    def create_survival_plot(DF_temp, file_name):

        # (Pdb) ax = kmf_UP.fit(T[ix], E[ix],
        # label='up').plot_survival_function(ax=ax)
        # *** TypeError: NaNs were detected in the dataset. Try using pd.isnull
        # to find the problematic values.
        # in case T contains NaN-> make them to 0
        index_to_set_null = DF_temp[DF_temp['T'].isnull()].index
        DF_temp.loc[index_to_set_null, 'T'] = 0
        # DF_temp[DF_temp['T'].isnull() ]
        # DF_temp['T'][DF_temp['T'].isnull()] = 0
        # for DF_temp in DF_list:
        # specifiy the gender dependent on gender col:
        if len(DF_temp['gender'].unique()) == 2:
            gender = 'both'
        else:
            gender = DF_temp['gender'].unique()[0]
        T = DF_temp['T']
        # "death" event observed -> make boolean col and True equals death,
        # event is observed:
        E = DF_temp['vital_state'] == 'dead'

        groups = DF_temp['group']

        ix = (groups == 'up')
        # ax = plt.subplot()

        # kmf = KaplanMeierFitter(label="waltons_data")
        # kmf.fit(waltons['T'], waltons['E'])
        # kmf.plot()
        # sns.set(palette=sns.color_palette('coolwarm_r'))
        sns.set_theme(style='white', palette=sns.color_palette('deep'))
        # sns.set_theme(style="ticks", palette=["b", "r"])
        # sns.set_theme(palette=["cornflowerblue", "sienna"])
        kmf_DOWN = KaplanMeierFitter()
        ax = kmf_DOWN.fit(T[~ix], E[~ix],
                          label='down').plot_survival_function()

        kmf_UP = KaplanMeierFitter()
        try:
            ax = kmf_UP.fit(T[ix], E[ix], label='up').plot_survival_function(
                ax=ax)
        except Exception as e:
            print(f'exception occured: {e}')
            breakpoint()
        results = logrank_test(T[ix], E[~ix], T[ix], E[~ix])
        p_value = f'p_value = {Decimal(str(results.p_value)):.2e}'

        cg_name = DF_temp['cg_name'].unique()[0]
        # start = DF_temp['Start'].unique()[0]
        dmr = DF_temp['DMR'].unique()[0]
        # chromosome = DF_temp['chr'].unique()[0]
        title_temp = ('survival function ' + cg_name + ' DMR: ' + dmr + '\n' +
                      p_value)
        ax.set_title(title_temp)
        path = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title, gender,
                            met_dir, 'lifeline_out')
        # file_temp = 'lifeline_' + cg_name + '_' + dmr + '_' + chromosome +\
        #     '_' + str(start) + '.pdf'
        file_temp = file_name.replace('tsv', 'pdf')
        file_log = os.path.join(PROJECT_title,
                                DRUGS_title,
                                gender, met_dir, 'lifeline_out', file_temp)
        file_temp = os.path.join(path, file_temp)

        add_at_risk_counts(kmf_UP, kmf_DOWN, ax=ax, rows_to_show=['At risk'])
        plt.tight_layout()

        try:
            plt.savefig(file_temp)
            breakpoint()
        except Exception as e:
            print('The following error occured while trying to write the\
                  file: {}: {}'.format(file_temp, e))
        else:
            logger = set_logger.set_logger(
                OUTPUT_PATH, PROJECT_title, DRUGS_title)
            logger.info('REPORT_17:\t{}'.format(file_log))
        plt.clf()
        plt.cla()

    for DF, file_name in zip(DF_list, file_name_list):
        # TODO
        # ________________________________case_id_______drug__gender____Project_vital_state__________T_group___chr_____cg_name______Start__________________DMR____median
        # 0____f5d90810-8c4b-4ced-9360-383e42cb54b3__cisplatin__female__TCGA-CESC_______alive___1.268493____up__chr5__cg02343823__150904857__150904853-150905038__0.307809
        # 1____da30a845-c4d3-4c78-b8b0-210239224f8f__cisplatin__female__TCGA-CESC_______alive__12.980822____up__chr5__cg02343823__150904857__150904853-150905038__0.307809
        # TODO: get every other case, not treatet with those drugs in this form
        # and crate the complement survivalplot
        # needed tables: complete_summary.tsv,
        create_survival_plot(DF, file_name)

    # T = instanz.lifeline_df['T']
    # E = instanz.lifeline_df['E']
    # kmf = KaplanMeierFitter()
    # kmf.fit(T, event_observed=E)  # or, more succinctly, kmf.fit(T, E)
    # # After calling the fit() method, we have access to new properties
    # like
    # # survival_function_ and methods like plot(). The latter is a wrapper
    # # around Panda’s internal plotting library.
    # kmf.survival_function_
    # kmf.cumulative_density_
    # ax = plt.subplots()
    # ax = kmf.plot_survival_function()  # or just kmf.plot()
    # title_temp = 'survival_function_every_case'
    # ax.set_title(title_temp)
    # path = os.path.join(
    #     OUTPUT_PATH, instanz.PROJECT_title, 'Lifeline_plots_general')
    # file_temp = 'lifelines_survival_every_case.pdf'
    # file_temp = os.path.join(path, file_temp)

    # plt.savefig(file_temp)
    # log_file = file_temp.replace(OUTPUT_PATH + os.path.sep, '')
    # logger = create_matrix_new.set_logger(
    #     OUTPUT_PATH, instanz.PROJECT_title, instanz.DRUGS_title)
    # logger.info('lifeline_plots_11:\t{}'.format(log_file))
    # plt.clf()
    # plt.cla()
    # ax = plt.subplots()
    # ax = kmf.plot_cumulative_density()
    # title_temp = 'cumulative_density_every_case'
    # ax.set_title(title_temp)

    # file_temp = 'lifelines_cumulative_density_every_case.pdf'
    # file_temp = os.path.join(path, file_temp)

    # plt.savefig(file_temp)
    # log_file = file_temp.replace(OUTPUT_PATH + os.path.sep, '')
    # logger.info('lifeline_plots_11:\t{}'.format(log_file))
    # plt.clf()
    # plt.cla()
    # # ############################################################
