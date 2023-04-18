import pandas as pd
import os
import glob
# from pathlib import path
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal
# from seaborn import regression
import create_matrix_new
# from lifelines import *  # lifelines has builtin parametric models. for
from lifelines import KaplanMeierFitter
from lifelines import WeibullFitter
from lifelines import ExponentialFitter
from lifelines import LogNormalFitter
from lifelines import LogLogisticFitter
from lifelines import PiecewiseExponentialFitter
from lifelines import GeneralizedGammaFitter
from lifelines import SplineFitter
from lifelines.plotting import add_at_risk_counts
# from lifelines.datasets import load_rossi
# from lifelines.utils import median_survival_times
# from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test
# to compare two survivalfunctions and get CoxPHFitter is actually
# better than logrank_test, see: https://discourse.datamethods.org/t/
# when-is-log-rank-preferred-over-univariable-cox-regression/2344
# but in our case we dont want to invoke the age (a covariate var is needed)
# since the KaplanMeier depends just on survivaltime[T] and event [E]
# from lifelines import CoxPHFitter
# from lifelines.datasets import load_regression_dataset
import mygene
# example, weibull, log-normal, log-logistic, and more.
# from lifelines.utils import median_survival_times
# a useful summary stat is the median survivaltime, which represents when 50%
# of the population has died:

# from lifelines import aalenadditivefitter  # for survival regression
# from lifelines import weibullaftfitter  # for survival regression
# from lifelines.datasets import load_regression_dataset

# create a class which returns data like:
# - result_table   (df)
# - count_matrix (df with mi out of countmatrix_info)
# - ensg_list


class Lifeline_plot:
    """
    returning the ENSG on which the gene specific KaplanMeier plots are created
    """
    def __init__(self, OUTPUT_PATH, PROJECT_title, DRUGS_title, path_to_result,
                 prefix, threshold, multi_project):
        self.OUTPUT_PATH = OUTPUT_PATH
        self.PROJECT_title = PROJECT_title
        # same as PROJECT_title, depending on single or multi project
        self.DRUGS_title = DRUGS_title
        self.path_to_result = path_to_result
        self.prefix = prefix
        # we need the flag after DESeq2_out_* to relate to the right count
        # table
        self.threshold = threshold
        self.DF_res = Lifeline_plot.result_DF(self)
        self.multi_project = multi_project

    def result_DF(self):
        path = os.path.join(self.path_to_result, 'DESeq2_results.tsv')
        DF_res = pd.read_csv(path, sep='\t', index_col=0)
        return DF_res

    def create_ENSG_list(self):
        DF_res_sorted_IN = self.DF_res.sort_values(
            by='log2FoldChange').dropna()
        DF_res_sorted_IN = DF_res_sorted_IN[DF_res_sorted_IN['pvalue'] < 0.05]
        # print("\nDF_res_sorted_IN, filtered pvalue:\n", DF_res_sorted_IN)
        DF_res_sorted_IN = DF_res_sorted_IN[DF_res_sorted_IN['padj'] < 0.05]
        # print("\nDF_res_sorted_IN, filtered padj:\n", DF_res_sorted_IN)
        index_INCREASE = DF_res_sorted_IN.nsmallest(
            n=10, columns='log2FoldChange').index
        index_DECREASE = DF_res_sorted_IN.nlargest(
            n=10, columns='log2FoldChange').index
        # print('ENSG list for PROJECT_title: ', self.PROJECT_title)
        ENSG_list = {'INCREASE': list(index_INCREASE),
                     'DECREASE': list(index_DECREASE)}
        # save those 10 highest and 10 smallest ENSG, together with all info
        # out of the deseq result table plus mygene stuff
        return ENSG_list


# fct 11 in main
def lifelines_ENSG(OUTPUT_PATH, PROJECT_DRUG_UUID, threshold, DRUGS_title):
    r"""
    :param: OUTPUT_PATH: path for DESeq2 pipeline outputs
    :type: OUTPUT_PATH: str
    :param: DRUGS: applied drug(s)
    :type: DRUGS: list of strings
    :param: PROJECT_DRUG_UUID: hash holding a project to the UUID of the\
    belonging drugtable
    :type: PROJECT_DRUG_UUID: dict
    :param: threshold: parameter for the lifeline plots helping for the\
    classification of expression data
    :type: threshold: int

    | script creates the plots:
    | lifelines\_cumulative\_density.svg,
    | lifelines\_multiple\_groups.svg,
    | lifelines\_parametric\_models\_2.svg,
    | lifelines\_survival\_fct.svg,
    | lifelines\_table.tsv in OUTPUT\_PATH

    - with the different expression ENSGs (normalized counts given from\
        deseq, median from them and separated in UP and DOWN) \
        groups plottet with lifeline

    **needed**:

    - ALL\_PROJECTS\_summary\_dead\_alive\_reduced\_INFO.tsv as

    **count\_DF\_MI** (just the info for building the multiindex for\
    **count\_DF**):

    .. _multiindex_table:
    .. table:: example multiindex table

        =============   ====================================
        variable        value
        =============   ====================================
        vital\_status   alive
        gender          female
        case\_id        6ff12a54-10da-4941-bfea-7b66e19b4be9
        PROJECT         TCGA-CESC
        =============   ====================================

    - DESeq2\_MF\_normalized\_counts\_reduced.tsv as

    .. _ENSG_table:
    .. table:: ENSG table, **count_DF** (normalized counts with help of DESeq2)

        =============== ===============
        variable        value
        =============== ===============
        ENSG00000000003 4109.0073147311
        ENSG00000000005 0
        =============== ===============

    - OUTPUT\_PATH/DRUGS/DESeq2\_MF\_results\_reduced.tsv as **DF\_res**

    creates the **ENSG\_list** on the basis of the DESeq\_results in
    dependence of the resulttables, we get the 10 highest and 10 lowest
    logfoldchange

    --> here we sort log2fold change wise and merge then the most different
    INCREASING and DECREASING cases with the
    TCGA-\*/DF\_3t\_both\_with\_DRUG\_combi.tsv, in this table we have the
    fields:

    UUID case\_id gender, vital\_status, drugnames, survivaltime,
    years\_to\_last\_follow\_up, year\_of\_birth, year\_of\_death,
    age\_at\_diagnosis, PROJECT,

    [the survivaltime is created ot of the
    nationwidechildrens.org_clinical_patient_cesc.txt table, col days_to_death,
    this is performed in fct. 6, correct_drugs()]

    with it we create the table for lifeline
    "death" event observed -> make boolean col and True equals death, event
    is observed:

    .. _lifeline_table:
    .. table:: example table for lifeline input

        ======= ================== ======== ===================================
        index    T	                E       case_id
        ======= ================== ======== ===================================
        0       5.197260273972604   True    9ffa79fa-d2d8-48e1-8fd6-4b020ecf357
        1       0.8602739726027397  False   0de19185-3517-4e30-925b-7eb1f5079ec
        ======= ================== ======== ===================================

    - the up and down separation depends on the median of the normalized \
    count matrix

    - setting the threshold based on the median of the logfoldchange, delete \
        out 10 % around it if -f 11 is applied

    .. we have to decide if theres a table over several projects or if just
    .. one project abailable...

    .. _function_11:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 11 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 11

        # you can try out different thresholds, new directorys are
        # created therefore:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 11 -t 10

    """
    # DRUGS_title = '_'.join(sorted(map(str.lower, DRUGS)))
    project_list = []
    instanz_list = []
    # since new directory structure, no need to distinguishe between multi and
    # single project, deseq results are in out, proj, DRUGS_title, deseq2_out*,
    # but there can be multiple deseq2_out dirs, create an instanz also for
    # every result table...
    for project in PROJECT_DRUG_UUID:
        for path_to_result in glob.glob(os.path.join(OUTPUT_PATH,
                                                     project,
                                                     DRUGS_title,
                                                     'DESeq2_out_*')):

            pattern = os.path.join(OUTPUT_PATH, project, DRUGS_title,
                                   'DESeq2_out_')
            count_table_prefix = path_to_result.replace(pattern, '')
            instanz = Lifeline_plot(OUTPUT_PATH, project, DRUGS_title,
                                    path_to_result, count_table_prefix,
                                    threshold, False)
            ENSG_list = instanz.create_ENSG_list()
            instanz.ENSG_list = ENSG_list
            instanz_list.append(instanz)
        project_list.append(project)
    # check whether multiple projects are applied:

    if len(project_list) > 1:
        PROJECT_title = "_".join(sorted(project_list))
        for path_to_result in glob.glob(os.path.join(OUTPUT_PATH,
                                                     PROJECT_title,
                                                     DRUGS_title,
                                                     'DESeq2_out_*')):
            pattern = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                   'DESeq2_out_')
            count_table_prefix = path_to_result.replace(pattern, '')
            instanz = Lifeline_plot(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                    path_to_result, count_table_prefix,
                                    threshold, True)
            ENSG_list = instanz.create_ENSG_list()
            instanz.ENSG_list = ENSG_list
            instanz_list.append(instanz)

    # #### until now we have Lifeline_plot instances list of CESC and HNSC and
    # CESC_HNSC with:
    # DRUGS_title, OUTPUT_PATH, PROJECT_title, threshold, multi_project, DF_res
    # and ENSG_LIST now

    #########################################################################
    # extract the survivaltimes out of the DRUG_combi tables
    # they are located in the OUTPUT_PATH/single_project names:
    # add the tables projectDF_3t_both_with_DRUG_combi.tsv for single projects
    # the tables for the merged projects have to be merged
    for instanz in instanz_list:
        if not instanz.multi_project:
            path = os.path.join(instanz.OUTPUT_PATH, instanz.PROJECT_title,
                                'DF_3t_both_with_DRUG_combi.tsv')
            instanz.DF_3t_both_DRUG_combi = pd.read_csv(path, sep="\t")
    # merge the cases from both projects and apply them to the multi_projects
    # instances
    # gather the 3t_both of all single project, then set the merged table to
    # the multi_project instanz:
    DF_3t_both_merged = []
    for instanz in instanz_list:
        if not instanz.multi_project:
            DF_3t_both_merged.append(instanz.DF_3t_both_DRUG_combi)
        # print("\ninstanz.PROJECT_title:\n", instanz.PROJECT_title)
    for instanz in instanz_list:
        if instanz.multi_project:
            instanz.DF_3t_both_DRUG_combi = pd.concat(DF_3t_both_merged)
    # for instanz in instanz_list:
        # print('PROJECT_title: {} \n DF: \n {} '.format(
            # instanz.PROJECT_title, instanz.DF_3t_both_DRUG_combi))

    # #################################
    # everything here can be done to every instance:
    # merge the years_to_last_follow_up and survivaltime to one col named T

    def create_lifeline_df(instanz):
        DF = instanz.DF_3t_both_DRUG_combi
        DF['T'] = DF['years_to_last_follow_up'].combine_first(
            DF['survivaltime'])
        # make a bool vector out of vital_status:
        DF['E'] = DF['vital_status'] == 'alive'
        # the age at diagnosis is preserved for bin overview of applied cases,
        # in REPORT.pdf
        lifeline_df = DF.loc[:, ['T', 'E', 'case_id', 'age_at_diagnosis',
                                 'drugnames']]
        lifeline_df.dropna(inplace=True)  # delete empty rows
        lifeline_df.drop_duplicates(subset='case_id', inplace=True)
        # if values are to small, an error is raised, add a small value:
        filt = lifeline_df['T'] < 0.000001
        lifeline_df.loc[filt, 'T'] = lifeline_df.loc[filt, 'T'] + 0.000001
        # and negatives occur also, make the abs of every value:
        # ##################save tables:
        lifeline_df['T'] = abs(lifeline_df['T'])
        # if instanz.multi_project:
        #     path = os.path.join(
        #         OUTPUT_PATH, 'Lifeline_plots_general_' + PROJECT_title)
        #     os.makedirs(path, exist_ok=True)
        #     lifeline_df.to_csv(os.path.join(
        #         path, 'lifelines_table.tsv'), sep='\t')
        # else:  # gather the project name out of the PROJECT col
        path = os.path.join(OUTPUT_PATH, instanz.PROJECT_title,
                            'Lifeline_plots_general')
        os.makedirs(path, exist_ok=True)
        lifeline_df.to_csv(os.path.join(
            path, 'lifelines_table.tsv'), sep='\t', index=False)
        log_path = path.replace(OUTPUT_PATH + os.path.sep, '')
        logger = create_matrix_new.set_logger(OUTPUT_PATH,
                                              instanz.PROJECT_title,
                                              DRUGS_title)
        logger.info('lifeline_plots_11:\t{}'.format(
            os.path.join(log_path, 'lifelines_table.tsv')))
        # ##############################################################
        instanz.lifeline_df = lifeline_df
        # ##############################################################

    def create_general_plots(instanz):
        T = instanz.lifeline_df['T']
        E = instanz.lifeline_df['E']
        kmf = KaplanMeierFitter()
        kmf.fit(T, event_observed=E)  # or, more succinctly, kmf.fit(T, E)
        # After calling the fit() method, we have access to new properties like
        # survival_function_ and methods like plot(). The latter is a wrapper
        # around Panda’s internal plotting library.
        kmf.survival_function_
        kmf.cumulative_density_
        ax = plt.subplots()
        ax = kmf.plot_survival_function()  # or just kmf.plot()
        title_temp = 'survival_function_every_case'
        ax.set_title(title_temp)
        path = os.path.join(
            OUTPUT_PATH, instanz.PROJECT_title, 'Lifeline_plots_general')
        file_temp = 'lifelines_survival_every_case.pdf'
        file_temp = os.path.join(path, file_temp)

        plt.savefig(file_temp)
        log_file = file_temp.replace(OUTPUT_PATH + os.path.sep, '')
        logger = create_matrix_new.set_logger(
            OUTPUT_PATH, instanz.PROJECT_title, instanz.DRUGS_title)
        logger.info('lifeline_plots_11:\t{}'.format(log_file))
        plt.clf()
        plt.cla()
        ax = plt.subplots()
        ax = kmf.plot_cumulative_density()
        title_temp = 'cumulative_density_every_case'
        ax.set_title(title_temp)

        file_temp = 'lifelines_cumulative_density_every_case.pdf'
        file_temp = os.path.join(path, file_temp)

        plt.savefig(file_temp)
        log_file = file_temp.replace(OUTPUT_PATH + os.path.sep, '')
        logger.info('lifeline_plots_11:\t{}'.format(log_file))
        plt.clf()
        plt.cla()
        # ############################################################

        # ############################################################
        # start of different fitters:
        # Instead of the Kaplan-Meier estimator, you may be interested
        # in a parametric
        # model. lifelines has builtin parametric models. For example,
        # Weibull, Log-Normal, Log-Logistic, and more.
        fig, axes = plt.subplots(3, 3, figsize=(13.5, 7.5))
        kmf = KaplanMeierFitter().fit(T, E, label='KaplanMeierFitter')
        wbf = WeibullFitter().fit(T, E, label='WeibullFitter')
        exf = ExponentialFitter().fit(T, E, label='ExponentialFitter')
        lnf = LogNormalFitter().fit(T, E, label='LogNormalFitter')
        llf = LogLogisticFitter().fit(T, E, label='LogLogisticFitter')
        pwf = PiecewiseExponentialFitter([40, 60]).fit(
            T, E, label='PiecewiseExponentialFitter')
        ggf = GeneralizedGammaFitter().fit(
            T, E, label='GeneralizedGammaFitter')
        sf = SplineFitter(
            np.percentile(T.loc[E.astype(bool)], [0, 50, 100])).fit(
                T, E, label='SplineFitter')
        wbf.plot_survival_function(ax=axes[0][0])
        exf.plot_survival_function(ax=axes[0][1])
        lnf.plot_survival_function(ax=axes[0][2])
        kmf.plot_survival_function(ax=axes[1][0])
        llf.plot_survival_function(ax=axes[1][1])
        pwf.plot_survival_function(ax=axes[1][2])
        ggf.plot_survival_function(ax=axes[2][0])
        sf.plot_survival_function(ax=axes[2][1])

        ax.set_title('lifelines survival parametric models_all_cases')

        file_temp = 'lifelines_survival_parametric_models_all_cases.pdf'
        file_temp = os.path.join(path, file_temp)

        plt.savefig(file_temp)
        log_file = file_temp.replace(OUTPUT_PATH + os.path.sep, '')
        logger.info('lifeline_plots_11:\t{}'.format(log_file))
        plt.clf()
        plt.cla()

    for instanz in instanz_list:
        create_lifeline_df(instanz)

    for instanz in instanz_list:
        try:
            create_general_plots(instanz)
        except Exception as e:
            print('skipping because of {}'.format(e))
            continue

    # DONE add here the complement count_tables:
    def create_count_DF(instanz):
        # info file for MI creation, multi_project:
        count_DF_MI = pd.MultiIndex.from_frame(
            pd.read_csv(os.path.join(
                instanz.OUTPUT_PATH, instanz.PROJECT_title,
                instanz.DRUGS_title,
                instanz.prefix + '_summary_dead_alive_INFO.tsv'), sep='\t').T)
        count_DF_MI.set_names(count_DF_MI[0], inplace=True)
        count_DF_MI = count_DF_MI[1:len(count_DF_MI)]
        # ##### add raw count table to instance:
        count_DF = pd.read_csv(
            os.path.join(
                instanz.OUTPUT_PATH, instanz.PROJECT_title,
                instanz.DRUGS_title, instanz.prefix +
                '_summary_dead_alive.tsv'),
            sep='\t', index_col=0)
        count_DF = count_DF.T.set_index(count_DF_MI)
        count_DF = count_DF.T
        count_DF.index.set_names('ENSG', inplace=True)
        temp_path = os.path.join(
            instanz.OUTPUT_PATH, instanz.PROJECT_title, instanz.DRUGS_title,
            instanz.prefix + '_Lifeline_raw_counts')
        instanz.count_path_DF_hash.update({temp_path: count_DF})

        # #######add normalised counttables:
        count_DF = pd.read_csv(
            os.path.join(
                instanz.OUTPUT_PATH,
                instanz.PROJECT_title, instanz.DRUGS_title,
                'DESeq2_out_' + instanz.prefix,
                'DESeq2_normalized_counts.tsv'),
            sep='\t', index_col=0).T.set_index(count_DF_MI).T
        count_DF.index.set_names('ENSG', inplace=True)
        temp_path = os.path.join(
            instanz.OUTPUT_PATH, instanz.PROJECT_title, instanz.DRUGS_title,
            instanz.prefix + '_Lifeline_norm_counts')
        instanz.count_path_DF_hash.update({temp_path: count_DF})

        # #######add nt counttables:
        count_DF = pd.read_csv(
            os.path.join(
                instanz.OUTPUT_PATH,
                instanz.PROJECT_title, instanz.DRUGS_title,
                'DESeq2_out_' + instanz.prefix,
                'DESeq2_norm_transform_counts.tsv'),
            sep='\t', index_col=0).T.set_index(count_DF_MI).T
        count_DF.index.set_names('ENSG', inplace=True)
        temp_path = os.path.join(
            instanz.OUTPUT_PATH, instanz.PROJECT_title, instanz.DRUGS_title,
            instanz.prefix + '_Lifeline_nt_counts')
        instanz.count_path_DF_hash.update({temp_path: count_DF})
        return instanz

    for instanz in instanz_list:
        instanz.count_path_DF_hash = {}
        # also include the complement count tables to later plot thos lifelines
        # for complement in ['', '_complement']:
        instanz = create_count_DF(instanz)

    # # ######################################################
    # # gene specific plot
    # # theses plots show differences in expression rate(UP
    # # and DOWN of a specific ENSG, based on median of normalised counts)
    # # if threshold is applied, those plots are saved in a tresh specific dir
    # # ######################################################

    for instanz in instanz_list:
        for ENSG_key in instanz.ENSG_list:
            # ENSG_key either 'INCREASE' or 'DECREASE'
            for ENSG in instanz.ENSG_list[ENSG_key]:
                for path in instanz.count_path_DF_hash:
                    # print("\nENSG_key:\n", ENSG_key)
                    # print("\nENSG:\n", ENSG)
                    # for count_DF, life_path in instanz.count_DF,
                    # instanz.life_path:
                    # count_DF is filtered on the DRUGS_TITLE selection
                    count_DF = instanz.count_path_DF_hash[path]
                    count_PATH = path
                    lifeline_df = instanz.lifeline_df
                    # if ENSG == 'ENSG00000178522':
                    # print("\ncount_DF:\n", count_DF)
                    # print("\ncount_DF.loc[", ENSG, ", :]\n",
                    # count_DF.loc[ENSG, :])
                    # # for evaluating the median, drop zeros, then get median:
                    # median = count_DF.loc[ENSG, (count_DF.loc[ENSG, :] !=
                    #                              0)].median()
                    # zeros are better kept:
                    median = count_DF.loc[ENSG, :].median()
                    # print("\nmedian:\n", median)
                    # print("\ncount_DF.loc[ENSG, :].min:\n",
                    #     count_DF.loc[ENSG, :].min())
                    # print("\ncount_DF.loc[ENSG, :].max:\n",
                    #     count_DF.loc[ENSG, :].max())

                    # # insert an threshold:
                    # if a threshold !=0 percent is applied, integrate it here,
                    # # create then new dirs
                    # # keep the actual values from the counts UP_DF_temp
                    if threshold != 0:
                        median_tresh_up = median * (1 + (threshold/2/100))
                        # filt = (count_DF.loc[ENSG, :] > (median * 1.1))
                        filt = (count_DF.loc[ENSG, :] > (median_tresh_up))
                        UP_DF = count_DF.loc[ENSG, filt]
                        UP_DF.loc[:] = 'UP'
                        median_tresh_down = median * (1 - (threshold/2/100))
                        filt = (count_DF.loc[ENSG, :] < (median_tresh_down))
                        DOWN_DF = count_DF.loc[ENSG, filt]
                        DOWN_DF.loc[:] = 'DOWN'
                        median_col_name = str(median_tresh_down
                                              ) + ' - ' + str(median_tresh_up)
                    else:
                        filt = (count_DF.loc[ENSG, :] > median)
                        UP_DF = count_DF.loc[ENSG, filt]
                        UP_DF.loc[:] = 'UP'
                        DOWN_DF = count_DF.loc[ENSG, ~filt]
                        DOWN_DF.loc[:] = 'DOWN'
                        median_col_name = str(median)

                    UP_DOWN_DF = pd.concat([UP_DF, DOWN_DF])
                    UP_DOWN_DF = UP_DOWN_DF.reset_index(2)
                    # setting the ENSG index to
                    # # ENSG index to col

                    # # #######################################################
                    # # df_lifeline_ENSG is the ENSG specific table for the
                    # lifeline plots
                    # # #######################################################

                    # this is also the filtering on drugs!
                    df_lifeline_ENSG = pd.merge(
                        lifeline_df, UP_DOWN_DF, on='case_id')
                    # # access every information given in the count_DF for the
                    # specific # ENSG -> those are the indexes, they can be
                    # selected directly, # drop the col multiindex levels, s.t.
                    # count_DF and # df_lifeline_ENSG can be merged
                    temp_df = count_DF.loc[ENSG, :].reset_index()
                    df_lifeline_ENSG_info = pd.merge(
                        df_lifeline_ENSG, temp_df, on='case_id')
                    # ensg duplicates arose, make them unique:
                    df_lifeline_ENSG_info.drop_duplicates(subset=['case_id'],
                                                          inplace=True)
                    df_lifeline_ENSG_info['median'] = median_col_name

                    T = df_lifeline_ENSG_info['T']
                    E = df_lifeline_ENSG_info['E']

                    # # T is an array of durations, E is a either boolean or
                    # binary # array representing whether the “death” was
                    # observed or not # (alternatively an individual can be
                    # censored). We will fit a # Kaplan Meier model to this,
                    # implemented as:
                    kmf = KaplanMeierFitter()
                    try:
                        kmf.fit(T, event_observed=E)
                        # or, more succinctly, kmf.fit(T, E)
                    except ValueError:
                        print('ValueError: Empty array/Series passed in.\
                              continuing')
                        continue

                    # # After calling the fit() method, we have access to new
                    # properties # like survival_function_ and methods like
                    # plot(). The latter is a
                    # # wrapper around
                    # # Panda’s internal plotting library.

                    # # ###### TODO is that needed?
                    # # # By specifying the timeline keyword argument in fit(),
                    # # # we can change how the
                    # # # above models are indexed:
                    # # kmf.fit(T, E, timeline=range(0, 100, 2))
                    # # kmf.survival_function_ # index is now the same
                    # asrange(0, 100, 2) # kmf.confidence_interval_  # index is
                    # the same as range(0, 100, 2)

                    # # # A useful summary stat is the median survival time,
                    # # # which represents when 50%
                    # # # of the population has died:
                    # # median_ = kmf.median_survival_time_
                    # # median_confidence_interval_ = median_survival_times(
                    # # kmf.confidence_interval_)
                    # # #######################################################

                    # # Multiple groups
                    # # check here if UP and DOWN values are present:
                    groups = df_lifeline_ENSG[ENSG]
                    if not groups.isin(['UP']).any() or not \
                            groups.isin(['DOWN']).any():
                        print('no groups for comparison available', end='')
                        print(' continuing ({})'.format(ENSG))
                        # print("\ngroups:\n", groups)
                        continue
                    # print("\ngroups:\n", groups)
                    ix = (groups == 'UP')
                    # TODO specify the timeline argument in kmf.fit()
                    # ax = plt.subplot()

                    ###########################################################
                    # including the complements:
                    # plt.cla()
                    # plt.clf()
                    # data = load_dd()
                    # kmf_temp = KaplanMeierFitter()
                    # T_ = data["duration"]
                    # E_ = data["observed"]
                    # kmf_temp.fit(T_, E_)
                    # kmf_temp.survival_function_.plot()
                    # median_ci =
                    # median_survival_times(kmf.confidence_interval_)
                    # ###measure if two survival_functions are different ->
                    # begin logrank_test
                    results = logrank_test(T[ix], E[~ix], T[ix], E[~ix])
                    # results.print_summary()
                    p_value = f'p_value = {Decimal(str(results.p_value)):.2e}'

                    # end logrank_test

                    # ### CoxPHFitter example:
                    # cph = CoxPHFitter()
                    # # cph.fit(
                    #     # df_lifeline_ENSG_info.loc[
                    #         # :, ['T', 'E', 'age_at_diagnosis']], 'T', 'E')
                    # cph.print_summary()
                    # cph.plot_partial_effects_on_outcome(
                    # covariates='age_at_diagnosis', values=[
                    # 30, 40, 50, 60, 70], cmap='coolwarm')
                    ###########################################################

                    kmf_DOWN = KaplanMeierFitter()
                    ax = kmf_DOWN.fit(
                        T[~ix], E[~ix], label='DOWN').plot_survival_function()
                    # ax = kmf.plot()
                    kmf_UP = KaplanMeierFitter()
                    ax = kmf_UP.fit(
                        T[ix], E[ix], label='UP').plot_survival_function(ax=ax)
                    # ax = kmf.plot(ax=ax)
                    # ax = kmf_UP.plot_survival_function(ax=ax)
                    # add_at_risk_counts(kmf_DOWN, kmf_UP, ax=ax)
                    add_at_risk_counts(kmf_UP, kmf_DOWN, ax=ax,
                                       rows_to_show=['At risk'])

                    # ENSG holds the ENSG identifier, here the genenames can be
                    # added:
                    # take the symbol out of the gtf annotation, the long name
                    # out of the mygene library
                    annot_DF = pd.read_csv(
                        os.path.join(
                            instanz.path_to_result,
                            'DESeq2_results_with_annotation.tsv'), sep='\t',
                        index_col='ENSG')
                    # the symbol is within the used annotation, it must be
                    # present, no try block necessary.
                    gene_symbol = annot_DF.loc[ENSG, 'symbol']
                    mg = mygene.MyGeneInfo()
                    # if the long gene name can be found, not None is returned:
                    try:
                        gene_name = mg.getgene(ENSG, fields='name')['name']
                    except Exception as e:
                        print(f'Exception: {e}')
                        print('No long gene name available for ', end='')
                        print(f'{ENSG} in mygene library ', end='')
                        print(f'using the gene_symbol instead: {gene_symbol}')
                        gene_name = gene_symbol
                    ax.set_title('multiple_groups (UP, DOWN)' + ENSG + '_' +
                                 ENSG_key + '\n' + 'genename: ' + gene_name +
                                 '  symbol: ' + gene_symbol + '\n' +
                                 'median: ' + str(median) + '  ' + p_value)
                    file_temp_start = 'lifelines_multiple_groups_' + ENSG + \
                        '_' + ENSG_key + '_' + gene_symbol + '.pdf'
                    tresh_dir = "threshold_" + str(threshold)
                    os.makedirs(os.path.join(count_PATH, tresh_dir),
                                exist_ok=True)
                    file_temp = os.path.join(count_PATH, tresh_dir,
                                             file_temp_start)
                    plt.tight_layout()
                    plt.savefig(file_temp)
                    plt.clf()
                    log_file = file_temp.replace(OUTPUT_PATH + os.path.sep, '')
                    logger = create_matrix_new.set_logger(
                        OUTPUT_PATH,
                        instanz.PROJECT_title, instanz.DRUGS_title)
                    logger.info('lifeline_plots_11:\t{}'.format(log_file))
                    # # also save the belonging table to that plot with infos
                    # from # count_table, cut away .svg set .tsv instead:
                    df_lifeline_ENSG_info.to_csv(
                        os.path.splitext(file_temp)[0] + '.tsv', sep='\t')
                    log_file = os.path.splitext(file_temp)[0] + '.tsv'
                    log_file = log_file.replace(OUTPUT_PATH + os.path.sep, '')
                    logger.info('lifeline_plots_11:\t{}'.format(log_file))
                    # print("\nfile_temp:\n", file_temp)
                    # print("\nos.path.splitext(file_temp)[0] + '.tsv':\n",
                    #     os.path.splitext(file_temp)[0] + '.tsv')
                    # # # provife here also the tables, merge them with the
                    # specific # values
                    plt.clf()
                    plt.cla()

        # single = ['summary_dead_alive_', 'deseq2_normalized_counts_',
                # 'deseq2_norm_transform_counts_']
        # multi = [PROJECT_title + '_summary_dead_alive_',
                # 'deseq2_normalized_counts_',
                # 'deseq2_norm_transform_counts_']
        # suffix = ['raw_counts',  'norm_counts', 'nt_counts']
        # apply an instanz to the fuct, create 3 DF within instanz, depending
        # on if its single or multi project...

        # the raw count tables are finnished
        # ### add all 3 count tables to all instances:
        # combinations based on the resulttables of CESC, HNSC, CESC-HNSC,
        # f.e.:
        # - CESC: OUT, PROJECT_title, DRUG, DESEQ_out, -->
        # DESeq2_results_reduced.tsv the ensg list is based on that table and
        # ENSG specific lifelineplots
        # are created with the 3 counttables:
        #
        # ####################### example for singleproject:
        #  ***RAW***:
        # - CESC:  OUT, PROJECT_title, DRUG,  -->
        # summary_dead_alive_reduced_INFO.tsv -->
        # summary_dead_alive_reduced.tsv
        # ***NORM_COUNTS****:
        # - CESC: OUT, PROJECT_title, DRUG, DESEQ_out, -->
        # DESeq2_normalized_counts_reduced.tsv ***NT-COUNTS***
        # - CESC: OUT, PROJECT_title, DRUG, DESEQ_out, -->
        # DESeq2_norm_transform_counts_reduced.tsv
        # ####################### example for multiproject:
        # ***RAW***:
        # -CESC_HNSC:  OUT, DRUG, --> _summary_dead_alive_reduced_INFO.tsv
        #                               --> project,
        #                               _summary_dead_alive_reduced.tsv
        # ***NORM-COUNTS***
        # -CESC_HNSC: OUT, DRUG, DESEQ_out, project, -->
        # DESeq2_normalized_counts_reduced.tsv
        # *** NT-COUNTS***
        # - CESC_HNSC: OUT, DRUG, DESEQ_out, project, -->
        # DESeq2_norm_transform_counts_reduced.tsv
