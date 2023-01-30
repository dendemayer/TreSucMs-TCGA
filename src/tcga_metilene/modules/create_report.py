import subprocess       # to call pandoc
import os
import re
import glob
import matplotlib.pyplot as plt  # the included frequency heatmap is also
# created here
from methyl import set_logger
# import re
import numpy as np
import pandas as pd
import seaborn as sns
from PyPDF2 import PdfFileReader, PdfFileWriter


# fct 17 in main
def create_report_pdf(OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID,
                      met_dir, cutoff, DRUGS_title):
    """
    :param: OUTPUT_PATH: path for metilene pipeline outputs
    :type: OUTPUT_PATH: str
    :param: DRUGS: applied drug(s)
    :type: DRUGS: list of str
    :param: SCRIPT_PATH: path to the metilene_pipeline repo
    :type: SCRIPT_PATH: str
    :param: PROJECT_DRUG_UUID: hash holding a project to the UUID of the\
    belonging drugtable
    :type: PROJECT_DRUG_UUID: dict

    creating a report file with all visual outputs created in an analysis,
    saved at OUTPUT_PATH/PROJECT_title/DRUGS_title/met_dir/REPORT.pdf

    .. _function_17:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 17 to your call,
        # example:
        $ python main_metilene.py -p TCGA-CESC -d cisplatin -f 17

    """
    print(f'creating report for {PROJECT_DRUG_UUID}')
    if isinstance(PROJECT_DRUG_UUID, dict):  # a multi project as dict
        projects = []
        for project in PROJECT_DRUG_UUID:
            projects.append(project)
        PROJECT_title = '_'.join(sorted(map(str.upper, projects)))
    else:
        PROJECT_title = PROJECT_DRUG_UUID
        projects = PROJECT_DRUG_UUID

    gender_list_temp = ['both', 'female', 'male']
    gender_list = []
    for gender in gender_list_temp:
        if os.path.exists(os.path.join(OUTPUT_PATH, PROJECT_title,
                                       DRUGS_title, gender)):
            gender_list.append(gender)
    gender_str = '[' + '|'.join(gender_list) + ']'
    pdf_file = os.path.join(
        OUTPUT_PATH, PROJECT_title, DRUGS_title, 'REPORT.pdf')
    md_file = os.path.join(
        OUTPUT_PATH, PROJECT_title, DRUGS_title, 'REPORT.md')
    with open(md_file, 'w') as f:
        f.write('# Report file for the metilene_pipeline  \n\n')
        f.write('## Query:  \n\n')
        f.write('- PROJECT:\t')
        if isinstance(PROJECT_DRUG_UUID, dict):  # a multi project as
            # dict
            f.write('__' + '; '.join(projects) + '__  \n')
        else:
            f.write('__' + PROJECT_title + '__  \n')
        # f.write('__' + '; '.join(projects) + '__  \n')
        f.write('- DRUGS:\t')
        f.write('__' + '; '.join(DRUGS) + '__  \n')
        f.write('- OUTPUT_PATH:\t')
        f.write('__' + os.path.basename(OUTPUT_PATH) + '__  \n')
        f.write('- metilene results:\t')
        f.write('__' + '/ '.join(os.path.join(
            os.path.basename(OUTPUT_PATH),
            PROJECT_title, DRUGS_title, gender_str, met_dir).split(
                '/')) + '__  \n')
        f.write('- SCRIPT_PATH:\t')
        f.write('__' + os.path.basename(SCRIPT_PATH) + '__  \n')

        name_copied = os.path.join(
            'Snakes', 'snakemake_config_'
            + PROJECT_title + '_' + DRUGS_title + '.yaml')
        f.write('- your configfile within the SCRIPT_PATH:\n\n   ')
        f.write(name_copied)
        f.write('   \n\n')
        if cutoff > 0.0:
            f.write('- cutoff Parameter:\t')
            f.write('__')
            if cutoff % 1 == 0:
                f.write(str(int(cutoff)))
                cut_str = str(int(cutoff))
            else:
                f.write(str(cutoff))
                cut_str = str(cutoff)
            f.write('__  \n')
        f.write('   \n\n')

        f.write(r'\pagebreak')
        f.write('  \n\n')

        # ###################### START Overview of applied cases: #############
        # just print 1 table, in case both genders appear (gender_list =
        # ['female', 'male', 'both']), just print both, if just one gender,
        # print that gender:
        if len(gender_list) == 1:
            gender = gender_list[0]
        else:
            gender = 'both'
        DF_invoked_cases = create_overview_table(
            OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID,
            DRUGS_title, PROJECT_title, gender)
        # ####### with the created DF,
        # ######
        # return(total_cases, complete_DF)
        # total_count = case_count_list[0]
        total_count = len(DF_invoked_cases)
        # DF_info_table = DF_invoked_cases
        DF_invoked_cases['age_at_diagnosis'] = DF_invoked_cases[
            'age_at_diagnosis'].dropna().astype('int32')
        DF_invoked_cases_fc = DF_invoked_cases.loc[
            :, ['vital_status',
                'gender', 'project', 'drugnames']].value_counts()
        alive_count = DF_invoked_cases['vital_status'].value_counts()[
            'alive']
        dead_count = DF_invoked_cases['vital_status'].value_counts()['dead']
        # f.write('## {} gender:  \n\n'.format(gender))

        tex_file_name = 'latex_all_invoked_cases.tex'

        pdf_name = create_pdf_from_DF(
            OUTPUT_PATH, SCRIPT_PATH, DRUGS_title,
            PROJECT_title, DF_invoked_cases_fc, tex_file_name, '# of cases')

        f.write('# Overview of your applied cases to the analysis  \n\n')
        f.write(
            '![factors of your applied cases](' + pdf_name + ')  \n\n')
        f.write('__total cases: ' + str(total_count) + '__  \n')
        f.write('__alive cases: ' + str(alive_count) + '__  \n')
        f.write('__dead cases: ' + str(dead_count) + '__  \n\n')
        if len(gender_list) > 1:
            female_count = DF_invoked_cases['gender'].value_counts()[
                'female']
            male_count = DF_invoked_cases['gender'].value_counts()['male']
            f.write('__female cases: ' + str(female_count) + '__  \n')
            f.write('__male cases: ' + str(male_count) + '__  \n\n')

            # ### start of cutoff added cases:
        if cutoff > 0:
            path = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                'cutoff_cases_add_' + str(cutoff) + '.tsv')
            try:
                cutoff_add_DF = pd.read_csv(path, sep='\t')
                # ### the same frequency table like the DF_invoked_cases, but
                # just
                # the cutoff dependent added cases:
                DF_freq_cutoff_add = cutoff_add_DF.T.loc[
                    ['vital_status', 'gender', 'PROJECT', 'drugnames'],
                    :].T.value_counts()
                tex_file_name = 'latex_cases_cutoff_add.tex'
                # create the pdf and include it in the report:
                pdf_name = create_pdf_from_DF(
                    OUTPUT_PATH,
                    SCRIPT_PATH,
                    DRUGS_title,
                    PROJECT_title,
                    DF_freq_cutoff_add, tex_file_name, '# of cases')
                f.write(r'\pagebreak')
                f.write('\n')
                f.write('## Cases added due to the cutoff of ' +
                        cut_str + ' years   \n\n')
                f.write('While having a vital status of dead, those cases ')
                f.write('outlived the cutoff limit and are ')
                f.write('therefore taken as alive cases into account to your ')
                f.write('analyses.  \n\n')
                f.write(
                    '![factors of added cases](' +
                    pdf_name + ')  \n\n')

                f.write('__total cases: ' + str(len(cutoff_add_DF)) + '__  \n')
                gender_total = {}
                for gender in cutoff_add_DF['gender'].value_counts().index:
                    gender_total.update(
                        {gender: str(
                            cutoff_add_DF['gender'].value_counts()[gender])})
                for gender in gender_total:
                    f.write(
                        '__' + gender + ' cases: ' + gender_total[
                            gender] + '__  \n')
                f.write('\n')
                # ### include a table showing the survivaltimes:
                DF_surv_time = cutoff_add_DF.loc[
                    :, ['gender',
                        'PROJECT',
                        'drugnames',
                        'survivaltime']].sort_values('survivaltime')
                DF_surv_time = DF_surv_time.round(1)
                DF_surv_time.rename(
                    {'survivaltime': 'survivaltime in years'},
                    axis=1, inplace=True)
                tex_file_name = 'latex_cutoff_add_surv_time.tex'
                pdf_name = create_pdf_from_DF(
                    OUTPUT_PATH,
                    SCRIPT_PATH,
                    DRUGS_title,
                    PROJECT_title,
                    DF_surv_time, tex_file_name)
                f.write('## Survivaltimes of the added cases:   \n\n')
                f.write(
                    '![survival_times of added cases](' +
                    pdf_name + ')  \n\n')
            except FileNotFoundError:
                print('not able to find {}'.format(path))

        # #################END Overview of applied cases: ###################

        # #######################barplot age at diagnose #####################
        # palette=sns.color_palette('coolwarm_r',6))
        # sns.set(palette=sns.color_palette('deep'))
        # sns.set_theme(style='whitegrid',palette=sns.color_palette(
        # 'coolwarm_r')[-2:])
        # plot_kws={'alpha':1
        sns.set_theme(
            style='whitegrid', palette=sns.color_palette('deep'))
        # sns.set_theme(palette=["cornflowerblue", "sienna"])
        plot = sns.displot(data=DF_invoked_cases, x='age_at_diagnosis',
                           col='gender', kde=True)
        plot_path = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                 gender, 'displot.pdf')
        plot.savefig(plot_path)
        plt.clf()
        plt.close()

        f.write(r'\pagebreak')
        f.write('  \n\n')
        f.write('## Distribution of age at diagnosis, all cases:  \n\n')
        f.write('![dis_plot](' + plot_path + ') \
                {} in \n\n **{}**'.format(
            os.path.basename(
                plot_path), '/ '.join(
                            os.path.split(
                                plot_path)[0].replace(
                                    OUTPUT_PATH + os.path.sep, '').split(
                                        '/'))))
        f.write('\n\n')
        f.write(r'\pagebreak')

        # ####################################################################
        # ## include the q-values sorted DMR's:

        for gender in gender_list:
            include_met_table(OUTPUT_PATH, PROJECT_title, DRUGS_title, gender,
                              met_dir, f)

        # ### populate the {candidate(int):gender} dict, for both the CAND and
        # cand candidates
        # gets: OUTPUT_PATH, PROJECT_title, DRUGS_title, gender, met_dir,
        # gender_list
        # return glob_path, CAND_gender_dict
        cand_hash, CAND_gender_dict = populate_can_gen_dict(OUTPUT_PATH,
                                                            PROJECT_title,
                                                            DRUGS_title,
                                                            gender_list,
                                                            met_dir)
        for gender in gender_list:
            f.write('  \n\n')
            f.write(r'\pagebreak')
            f.write('  \n\n')
            search_flag = 'CAND'
            # the function is just allowed to call if the gender also really
            # appears in the glob_path, this could be checked via the keys of
            # the CAND_gender_dict
            if gender in [*cand_hash[search_flag]]:
                print_plots(OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID,
                            DRUGS_title, PROJECT_title, search_flag, met_dir,
                            gender, f, cand_hash[search_flag])
            # ############################# end CANDidate plots:

            # #######################################
            # ### begin of cand candidates -> good results:

            f.write('  \n\n')
            f.write(r'\pagebreak')
            f.write('  \n\n')

            # in case this is a single project report, we have no good
            # candidates, since every found candidate satisfies the condition
            # of best candidate:
            if isinstance(projects, list):
                search_flag = 'cand'
                if gender in [*cand_hash[search_flag]]:
                    print_plots(OUTPUT_PATH, DRUGS, SCRIPT_PATH,
                                PROJECT_DRUG_UUID, DRUGS_title, PROJECT_title,
                                search_flag, met_dir, gender, f,
                                cand_hash[search_flag])

        # #### start of original metilene visual output:
        f.write('# original visual metilene outputs:  \n\n')
        for gender in gender_list:
            f.write('\n\n')
            f.write('## metilene results, {} gender:  \n\n'.format(gender))
            f.write('### plots located in:\n  {}  \n\n'.format(
                '/ '.join(os.path.join(os.path.basename(
                    OUTPUT_PATH), PROJECT_title, DRUGS_title, gender,
                    met_dir).split('/'))))
            for path in glob.glob(
                os.path.join(
                    OUTPUT_PATH, PROJECT_title,
                    DRUGS_title, gender, met_dir, 'metilene_qval.0.05_*.pdf')):
                f.write('![met_result](' + path + ') **{}**  \n\n'.format(
                    os.path.basename(path)))
                f.write('  \n\n')
                f.write(r'\pagebreak')
                f.write('  \n\n')
        # #### end of original metilene visual output:

            # ############################# end cand- candidate plots:
            # end of gender specific stuff belonging to metilene results
            # #######################################
        #  add drugfrequencies according to the projects, but not restricted to
        #  the chosen drugs, to give an overview of which drugs may also be
        #  added to the analysis if, based on the complete table in ScriptPath:
        #  DRUG_combi_frequency_all.tsv:
        # TODO check whether those are the ones from the api 31 version!
        logger = set_logger.set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title)
        complete_frequency_DF = pd.read_csv(os.path.join(
            SCRIPT_PATH, os.path.pardir, 'resources',
            'DRUG_combi_frequency_all.tsv.gz'), sep='\t')
        complete_frequency_DF = complete_frequency_DF.set_index('Unnamed: 0')
        complete_frequency_DF.index.name = 'Drugs'
        complete_frequency_DF = complete_frequency_DF.loc[:, projects].dropna(
            how='all')

        # distinguish between DF (multiproject) and Series:
        if len(complete_frequency_DF.shape) == 2:
            # a DF can be to large for plot, filter out every row with 1 as
            # highest value
            complete_frequency_DF = complete_frequency_DF[
                complete_frequency_DF.max(axis=1) > 1]
            fig, ax = plt.subplots(figsize=(14, 14))
            sns.heatmap(
                complete_frequency_DF,
                xticklabels=complete_frequency_DF.columns,
                ax=ax,
                yticklabels=1,
                annot=True,
                fmt='g').set_title(
                'DRUG_combi_frequency for\n{}\n'.format(PROJECT_title),
                fontsize=18)
            # cbar_kws={'label': 'number of Drugs'}
            # .set_xticklabels(ax.get_xticklabels(), rotation=90)
            plt.xlabel('Projects', fontsize=18)
            plt.ylabel('Drugs', fontsize=18)
            plt.tick_params(which='major', labelsize=12, labelbottom=False,
                            bottom=False, top=False, labeltop=True)
            # plt.xticks(rotation=90)
            plt.tight_layout()
            file_name = os.path.join(
                OUTPUT_PATH,
                PROJECT_title, DRUGS_title,
                'DRUG_combi_frequency_heatmap.pdf')
            plt.savefig(file_name)
            plt.clf()
            plt.close()
            log_file = file_name.replace(OUTPUT_PATH + os.path.sep, '')
            logger.info('REPORT_17:\t{}'.format(log_file))
        else:
            # in case of a Series, a barplot is created:
            complete_frequency_DF = complete_frequency_DF[
                complete_frequency_DF > 1]
            ax = complete_frequency_DF.plot.barh(
                figsize=(10, 8),
                legend=True, title='frequency counts drugnames')
            plt.tight_layout()
            fig = ax.get_figure()
            file_name = os.path.join(
                OUTPUT_PATH,
                PROJECT_title, DRUGS_title,
                'DRUG_combi_frequency.pdf')
            fig.savefig(os.path.join(file_name))
            log_file = file_name.replace(OUTPUT_PATH + os.path.sep, '')
            logger.info('REPORT_17:\t{}'.format(log_file))
            plt.close()

        # #### now include the just created heatmap:
        f.write('## Drugfrequencies:   \n\n')
        f.write('### Every drug/drugcombination available according to your\
                choice of projects, NOT restricted to your choice of\
                therapeutics :   \n\n')
        f.write('* to keep the plot well-arranged,\
                rows with a maximum value of 1 are filtered out   \n\n')
        f.write('![frequency](' + file_name + ') \
                {} for {}'.format(os.path.basename(file_name),
                                  PROJECT_title))
        f.write('\n\n')

        f.close()

        # log the md file:
        log_file = md_file.replace(OUTPUT_PATH + os.path.sep, '')
        logger.info('REPORT_17:\t{}'.format(log_file))
        # the missing iftex.sty is located in the SCRIPT_PATH, cd to it s.t. it
        # can be found by latex
        os.chdir(os.path.join(SCRIPT_PATH, os.path.pardir, 'resources',
                              'iftex'))
        pandoc_sequence = [
            'pandoc', md_file, '-f', 'markdown-implicit_figures',
            '-t', 'latex', '-o', pdf_file]
        print('calling pandoc with:\n{}'.format(pandoc_sequence))
        # in case the pdf file exists already, delete it in prior:
        if os.path.isfile(pdf_file):
            subprocess.check_call(['rm', pdf_file])
        try:
            subprocess.check_call(pandoc_sequence)
            # log the pdf
            log_file = pdf_file.replace(OUTPUT_PATH + os.path.sep, '')
            logger.info('REPORT_17:\t{}'.format(log_file))
        except subprocess.CalledProcessError:
            breakpoint()
            print('\nfailed to create the report pdf, ')
            print('do you have installed latex?\n')
            return


def include_met_table(OUTPUT_PATH, PROJECT_title, DRUGS_title, gender,
                      met_dir, f):
    '''
    metilene_qval.0.05.out:

    | chr   | start     | stop      | q-value    | m.m.d.    | #CpGs
    | ----- | --------- | --------- | ---------- | --------- | ------
    | chr1  | 52842924  | 52843589  | 7.9501e-09 | -0.123670 | 6
    | chr1  | 100538915 | 100539953 | 0.01409    | -0.080123 | 8

    | p (MWU) | p (2D KS)
    --------- | -----
    | 0.21505 | 0.33872
    | 0.26316 | 0.34329

    '''
    table_name = os.path.join(
        OUTPUT_PATH, PROJECT_title, DRUGS_title, gender, met_dir,
        'metilene_qval.0.05.out')
    try:
        DF_met_tabl = pd.read_csv(table_name, sep='\t', header=None)
    except Exception as e:
        # TODO write in report that there is no table to read from....
        print(f'Exception thrown:\n{e}')
        print(f'table not available:\n{table_name}')
        print('returning')
        f.write('no columns to parse from: ')
        f.write(table_name)
        return

    # difficulties to access int label 3 for sort_values via by= parameter,
    # rename those labels:

    # Int64Index([0, 1, 2, 3, 4, 5, 6, 7], dtype='int64')
    DF_met_tabl.columns = ['chr', 'start', 'stop', 'q-value', 'm_m_d', 'CpGs',
                           'p(MWU)', 'p(2DKS)']
    DF_met_tabl.sort_values(by='q-value', inplace=True)
    DF_met_tabl['DMR length'] = DF_met_tabl['stop'] - DF_met_tabl['start']
    DF_met_tabl.reset_index(inplace=True, drop=True)
    # make table shorter:
    if len(DF_met_tabl) > 20:
        max_ind = 20
        more_dmr_than_max = True
    else:
        max_ind = len(DF_met_tabl)
        more_dmr_than_max = False
    # max_ind = len(DF_met_tabl)
    # more_dmr_than_max = False

    f.write('  \n\n')
    if more_dmr_than_max:
        f.write('### The first ')
        f.write(str(max_ind))
        f.write(' of ')
        f.write(str(len(DF_met_tabl)))
        f.write(' DMR ranges found')
    else:
        f.write('### All ')
        f.write(str(len(DF_met_tabl)))
        f.write(' DMR ranges found')
    f.write(', sorted after q-value, for ')
    f.write(gender)
    f.write(' genders ')
    f.write('  \n\n')

    # ##### begin table header
    f.write(' Chr | DMR | q-value | ')
    f.write('mean methylation diff | #CpGs \n')
    f.write('-|-|-|-|-  \n')
    # ##### end table header
    # ##### begin table content
    for i in range(0, max_ind):
        row_as_dict = DF_met_tabl.loc[i, :].to_dict()
        f.write(row_as_dict['chr'])
        f.write(' | ')
        # f.write(str(row_as_dict['start']))
        # f.write(' | ')
        dmr = str(row_as_dict['start']) + '-' + str(row_as_dict['stop'])
        f.write(dmr)
        f.write(' | ')
        f.write(str(row_as_dict['q-value']))
        f.write(' | ')
        f.write(str(row_as_dict['m_m_d']))
        f.write(' | ')
        f.write(str(row_as_dict['CpGs']))
        f.write('\n')
    # # ##### end table content
    f.write('  \n\n')
    f.write(r'\pagebreak')
    f.write('  \n\n')


def create_overview_table(OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID,
                          DRUGS_title, PROJECT_title, gender):
    '''
    create the overview on what cases were invoked actually, therefore the
    header intersect_header.tsv is sufficient
    to create the age distribution plots, filter the used cases out of the
    meta_info.dat in the single project folders
    '''
    projects = []

    for project in PROJECT_DRUG_UUID:
        projects.append(project)
    header_path = os.path.join(
        OUTPUT_PATH, PROJECT_title, DRUGS_title, gender,
        'intersect_header.tsv')
    # header = pd.read_csv(os.path.join(OUTPUT_PATH, PROJECT_title,
    #                         # DRUGS_title, 'intersect_header.tsv'), sep='\t')
    with open(header_path, 'r') as f:
        header = f.readline()
    # leave chromosome, Start, End, REF, out:
    header = header.strip().split('\t')[4:]
    header_dict = {
        'case_id': [], 'vital_status': [], 'gender': [], 'project': [],
        'drugnames': []}
    for head in header:
        temp_split = head.split(';')
        header_dict['vital_status'].append(temp_split[0])
        header_dict['case_id'].append(temp_split[1])
        header_dict['drugnames'].append(temp_split[2])
        header_dict['gender'].append(temp_split[3])
        header_dict['project'].append(temp_split[4])
    invoked_cases_DF = pd.DataFrame(header_dict)
    # # since the gender is missing, we can get that information right away out
    # # of the single project meta_info_druglist_merged_drugs_combined.tsv and
    # # merge on case_id
    # since gender is not missing anymore, just count the amount of the present
    # genders:
    # # but for the age at diagnosis we need it anyway:
    #
    # ####################
    first = True
    if isinstance(PROJECT_DRUG_UUID, dict):
        for project in PROJECT_DRUG_UUID:
            meta_path = os.path.join(
                OUTPUT_PATH, project,
                'meta_info_druglist_merged_drugs_combined.tsv')
            if first:
                try:
                    meta_info_DF = pd.read_csv(meta_path, sep='\t')
                    first = False
                except FileNotFoundError:
                    print('File not found: {}'.format(meta_path))
                    print('REPORT can not be created')
                    os._exit(0)
            else:
                meta_info_DF = pd.concat(
                    [meta_info_DF, pd.read_csv(meta_path, sep='\t')])
    else:
        meta_path = os.path.join(
            OUTPUT_PATH, PROJECT_title,
            'meta_info_druglist_merged_drugs_combined.tsv')
        try:
            meta_info_DF = pd.read_csv(meta_path, sep='\t')
        except FileNotFoundError:
            print('File not found: {}'.format(meta_path))
            print('REPORT can not be created')
            os._exit(0)
    # now join
    meta_info_DF = meta_info_DF.rename({'case_ids': 'case_id'}, axis=1)
    # 'cases.0.demographic.gender': 'gender',
    # 'cases.0.demographic.vital_status': 'vital_status',
    # 'cases.0.demographic.year_of_birth': 'year_of_birth',
    # 'cases.0.demographic.year_of_death': 'year_of_death',
    # 'cases.0.diagnoses.0.age_at_diagnosis': 'age_at_diagnosis',
    # 'cases.0.diagnoses.0.days_to_last_follow_up':
    # 'days_to_last_follow_up'}, axis=1)

    # (Pdb) invoked_cases_DF.shape # (236, 3)
    # (Pdb) meta_info_DF.shape # (470, 35)
    # meta_info_DF['age_at_diagnosis'] = meta_info_DF['age_at_diagnosis'] / 365
    invoked_cases_DF = invoked_cases_DF.set_index('case_id')
    meta_info_DF = meta_info_DF.set_index('case_id')
    # meta_info_DF contains also follow up year_of_birth and much more, we need
    # here just the age_at_diagnosis for the barplot (joining the whole
    # dataframe would clash with the same colnamings from invoked_cases_DF)

    complete_DF = invoked_cases_DF.join(meta_info_DF['age_at_diagnosis'])
    return(complete_DF)


def gene_region_table(OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID,
                      DRUGS_title, PROJECT_title, search_flag, met_dir,
                      gender):
    '''
    regions_DF: using the already saved gtf annotation file
    cand_table_list: holds the candidates, depending on the search flag (either
    candidates with, *cand*-> project specifix, or *CAND*-> across projects)
    '''
    temp_DF = pd.read_csv(os.path.join(OUTPUT_PATH,
                                       'gencode.v36.annotation.gtf.gz'),
                          sep='\t', skiprows=5, header=None)
    temp_DF = temp_DF[temp_DF[2] == 'gene']
    # filter out Genesymbol (gene_name) and description (gene_type) of the
    # table
    temp_DF['ENSG'] = temp_DF[8].apply(
        lambda x: re.search(r'ENSG\d*', x.split(';')[0]).group())
    temp_DF['Description'] = temp_DF[8].apply(
        lambda x: x.split(';')[1].replace('gene_type "', '').strip('"'))
    temp_DF['Genesymbol'] = temp_DF[8].apply(
        lambda x: x.split(';')[3].replace('gene_name "', '').strip('"'))
    regions_DF = temp_DF.loc[:, [3, 4, 'ENSG', 0, 'Genesymbol', 'Description']]
    regions_DF.rename({3: 'Start', 4: 'End', 0: 'chr'}, axis=1, inplace=True)
    cand_table_list = glob.glob(
        os.path.join(
            OUTPUT_PATH, PROJECT_title, DRUGS_title, gender, met_dir,
            'boxplot_range', 'boxplot_means_range_*' + search_flag
            + '_linegraph.tsv'))
    range_list = []
    for i in range(0, len(cand_table_list)):
        temp_list_of_str = os.path.basename(
            cand_table_list[i]).split('_')[3].split('-')
        temp_list_int = [int(x) for x in temp_list_of_str]
        cand_DF = pd.read_csv(cand_table_list[i], sep='\t')
        chr_tle = cand_DF.loc[0, 'chr']
        temp_list_int.append(chr_tle)
        range_list.append(temp_list_int)
    regions_not_found = []
    first = True
    regions_found_DF = pd.DataFrame()
    for range_pos in range(0, len(range_list)):
        # we check whether the start position of our found DMR lies
        # inbetween the start and and of the supplied regions table:
        # or if the end position of DMR lies inbetween Start and end
        # check just if chr also matches
        # check the subset of the df, which holds the right chromosome:
        regions_subs = regions_DF[regions_DF['chr'] == range_list[
            range_pos][2]]
        # also initialise the list, where regions are saved without annotations
        # found to them
        bool_DF = (
            (regions_subs['Start'] <= range_list[range_pos][1]) & (
                regions_subs['End'] >= range_list[range_pos][1])) | ((
                    regions_subs['Start'] <= range_list[range_pos][0]) & (
                        regions_subs['End'] >= range_list[range_pos][0]))
        # print("\nrange_list[range_pos]:\n", range_list[range_pos])
        # since we used a subset , the DF length is changed, us the index,
        # where the value is True: bool_DF[bool_DF].index
        if bool_DF.nunique() == 2:
            if first:
                regions_found_DF = regions_DF.loc[bool_DF[bool_DF].index, :]
                regions_found_DF = regions_found_DF.assign(DMR=str(
                    range_list[range_pos][0]) + ' - ' + str(
                        range_list[range_pos][1]))
                first = False
            else:
                temp_DF = regions_DF.loc[bool_DF[bool_DF].index, :]
                temp_DF = temp_DF.assign(DMR=str(
                    range_list[range_pos][0]) + ' - ' + str(
                        range_list[range_pos][1]))
                regions_found_DF = pd.concat(
                    [regions_found_DF, temp_DF])
        else:
            # to include also the DMR's which cannot be annotated, create here
            # a list of those ranges:
            # subset of the
            regions_not_found.append(range_list[range_pos])

    if regions_found_DF.empty:
        return(regions_found_DF, regions_not_found)
    regions_found_DF.drop_duplicates('Genesymbol', inplace=True)
    regions_found_DF.rename({'Genesymbol': 'Code', 'Description': 'gene_type'},
                            axis=1, inplace=True)
    regions_found_DF.reset_index(drop=True, inplace=True)

    return(regions_found_DF, regions_not_found)


def print_plots(OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID,
                DRUGS_title, PROJECT_title, search_flag, met_dir, gender, f,
                CAND_gender_dict):
    if search_flag == 'CAND':
        rate = 'Best'
    else:
        rate = 'Good'

        # # cand_gender_dict = {}
        # # for gender in gender_list:
        #     # glob_path = glob.glob(
        #         # os.path.join(
        #             # OUTPUT_PATH,
        #             # PROJECT_title,
        #             # DRUGS_title,
        #             # gender, met_dir, 'boxplot_range', 'boxplot*cand.pdf'))
        #     # CAND_found = 0
        #     # for path in glob_path:
        #         # CAND_found += 1
        #     # cand_gender_dict.update({gender: str(CAND_found)})

    f.write('## {} candidate ranges, __{}__ gender:  \n\n'.format(
        rate, gender))
    if search_flag == 'CAND':
        f.write('* at every position here, a higher methylation median'
                + ' of the mean-beta-values in alive cases occur than in'
                + ' the dead cases of **all** projects, or'
                + ' vice versa:  \n\n')
    else:
        f.write('* at every position here, a higher methylation median'
                + ' of the mean-beta-values in alive cases occur than in'
                + ' the dead cases **per** project, or'
                + ' vice versa:  \n\n')
        f.write(' * in contrast to the best candidates, the restriction')
        f.write(' that the methylation median lies above (or below) at ')
        f.write(' every postion within the DMR is here unburdened to a ')
        f.write(' single project.  \n\n')

    f.write('* {} DMR\'s found by metilene:  \n\n'.format(
        CAND_gender_dict[gender]))
    # max_plots = 0
    glob_path = glob.glob(
        os.path.join(
            OUTPUT_PATH,
            PROJECT_title,
            DRUGS_title,
            gender, met_dir, 'boxplot_range', 'boxplot*{}.pdf'.format(
                search_flag)))
    if len(glob_path) == 0:
        f.write(
            'No candidate ranges found for the {} cases \n\n'.format(
                gender))

    # ############## gene region table:
    f.write('  \n\n')
    regions_found_DF, regions_not_found = gene_region_table(
        OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID,
        DRUGS_title, PROJECT_title, search_flag, met_dir, gender)
    # (Pdb) regions_not_found
    # [[74515181, 74516537, 'chr2'], [50338027, 50338575, 'chr3']]
    # (Pdb) regions_found_DF
    #         chr     Start       End         gene_type             ENSG   Code
    # 56753  chr22  23838628  23840628   protein_coding  ENSG00000099958  DERL3
    # Index(['chr', 'Start', 'End', 'whatever', 'strand', 'DMR', 'gene_type',
    #        'ENSG', 'Code'],
    # (Pdb) regions_found_DF.iloc[:,3:6]
    #     whatever strand                  DMR
    # 56753         0      -  23838304 - 23839478
    # gene_region_table returns
    f.write('  \n\n')
    f.write('### Out of the __')
    f.write(CAND_gender_dict[gender])
    f.write('__ DMR´s in {} cases '.format(gender))
    f.write(str(len(regions_found_DF)))
    f.write(' genes within that range could be')
    f.write(' annotated')
    f.write('  \n\n')
    # ############## gene region table:
    # table with ENSG description | chr | gene symbol | gene range |
    # DMR
    # ##### begin table header
    f.write(' ENSG description | chr | gene symbol| ENSG-Range')
    f.write(' | DMR \n')
    f.write('-|-|-|-|-  \n')
    # ##### end table header
    # ##### begin table content
    for i in regions_found_DF.index:
        if str(regions_found_DF.loc[i, 'ENSG']) == 'nan':
            ENSG = '-'
            ENSG_range = '-'
            ' - '.join([str(x) for x in regions_found_DF.loc[
                i, ['Start', 'End']].values.tolist()])
        else:
            ENSG = 'ENSG-' + str(regions_found_DF.loc[i, 'ENSG'])[4:]
            ENSG_range = ' - '.join([str(x) for x in regions_found_DF.loc[
                i, ['Start', 'End']].values.tolist()])
        chr_tle = str(regions_found_DF.loc[i, 'chr'])
        gene_symbol = str(regions_found_DF.loc[i, 'Code'])
        gene_type = str(regions_found_DF.loc[i, 'gene_type']).strip()
        if gene_symbol == 'nan':
            gene_symbol = '-'
        DMR = str(regions_found_DF.loc[i, 'DMR'])
        if DMR == 'nan':
            DMR = '-'
        if gene_type == 'nan':
            gene_type = '-'
        f.write(ENSG + '|' + chr_tle + '|' + gene_symbol + '|' +
                ENSG_range + '|' + DMR + '  \n')
    for i in regions_not_found:
        ENSG = '-'
        ENSG_range = '-'
        chr_tle = str(i[2])
        gene_symbol = '-'
        ENSG_range = '-'
        DMR = (str(i[0]) + ' - ' + str(i[1]))
        f.write(ENSG + '|' + chr_tle + '|' + gene_symbol + '|' +
                ENSG_range + '|' + DMR + '  \n')

    # # ##### end table content
    f.write('  \n\n')
    f.write(r'\pagebreak')
    f.write('  \n\n')

    # the output of the plots shall be sorted by the p-value derived by
    # metilene for the regions, this is hold in the respective table:
    # ############ BEGIN belonging q_value and CpGs for specific lineplot:
    region_qvalue_dict = {}
    for path in glob_path:
        temp_DF_linegraph = pd.read_csv(path.replace(
            '.pdf', '_linegraph.tsv'), sep='\t')
        q_value = temp_DF_linegraph.loc[0, 'q_value']
        CpGs = int(temp_DF_linegraph.loc[0, 'CpGs'])
        region_qvalue_dict.update({path: [q_value, CpGs]})
    path_qval_DF = pd.DataFrame.from_dict(
        region_qvalue_dict,
        orient='index', columns=['q_value', 'CpGs'])
    path_qval_DF.index.name = 'path'
    path_qval_DF.sort_values('q_value', inplace=True)
    # ############ END belonging q_value and CpGs for specific lineplot:

    # ############ BEGIN  plotting boxplots of DMR and lineplots of DMR:
    # if len(glob_path) != 0:
    #     f.write('### plots of the {} candidates:  \n'.format(rate))
    #     f.write('out of:  \n**{}**\
    #             \n\n'.format(
    #                 os.path.join(met_dir, 'boxplot_range\n')))
    #     f.write('sorted after q_value, maximum 10 plots:  \n\n')
    # count = 0
    # for path in path_qval_DF.index:
    #     if count == 10:
    #         break
    #     f.write('* __q_value = {}; CpGs = {}:__  \n\n'.format(
    #         path_qval_DF.loc[path, 'q_value'], str(
    #             path_qval_DF.loc[path, 'CpGs'])))
    #     f.write('![DMR](' + path + ') **{}**   \n'.format(
    #         os.path.basename(path)))
    #     line_plot_path = path.replace('.pdf', '_linegraph.pdf')
    #     f.write('![DMR](' + line_plot_path + ') **{}**   \n'.format(
    #         os.path.basename(line_plot_path)))
    #     f.write('  \n\n')
    #     f.write(r'\pagebreak')
    #     f.write('  \n\n')
    #     count += 1
    # ############ END  plotting boxplots of DMR and lineplots of DMR:

    # ############ BEGIN plotting kaplan meyer of a choosen postion of DMR:
    # ## DONE include here the respective complements of drugs not invoked
    # in our analyses: -> therefore the
    # # lifeline_plots.create_lifeline_plots() must be applied to the
    # # additional lifeline tables created by fct
    # # lifeline_plots.prepare_lifeline_plots_complement()
    # # -> those fct. are called in the main, before report creation, here the
    # respective tables are already present, get out which drug combination
    # is used in the complement:
    if len(glob_path) != 0:
        f.write('### plots of the {} candidates:  \n'.format(rate))
        f.write('out of:  \n**{}**\
                \n\n'.format(
            os.path.join(met_dir, 'boxplot_range\n')))
        f.write('sorted after q_value, maximum 10 plots:  \n\n')
    count = 0
    for path in path_qval_DF.index:
        if count == 10:
            break
        f.write('* __q_value = {}; CpGs = {}:__  \n\n'.format(
            path_qval_DF.loc[path, 'q_value'], str(
                path_qval_DF.loc[path, 'CpGs'])))
        f.write('![DMR](' + path + ') **{}**   \n'.format(
            os.path.basename(path)))
        line_plot_path = path.replace('.pdf', '_linegraph.pdf')
        f.write('![DMR](' + line_plot_path + ') **{}**   \n'.format(
            os.path.basename(line_plot_path)))
        f.write('  \n\n')
        f.write(r'\pagebreak')
        f.write('  \n\n')
        # # directly after those 2 plots set the lifeline plot in place to keep
        # all together:
        base_path = os.path.basename(path)
        pattern = re.compile(r'\d*-\d*')
        DMR = pattern.search(base_path).group(0)
        # with the DMR and the chromosome the right plot can be found at
        # lifeline_out/:
        glob_pattern = '*' + DMR + '_' + '*pdf'
        life_out_path = os.path.join(os.path.split(os.path.dirname(path))[0],
                                     'lifeline_out')
        try:
            kapl_plot = glob.glob(os.path.join(life_out_path, glob_pattern))[0]
        except Exception as e:
            print(f'exception: {e}')
            breakpoint()

        # get the position out of the pdf filename:
        position_temp = re.findall(r'[0-9]*_[^0-9]*pdf', kapl_plot)[0]
        position = re.search(r'\d*', position_temp).group(0)

        f.write('### Kaplan Meier Plot at position {} out of'.format(position))
        f.write('the DMR {}:  \n'.format(DMR))
        f.write('from directory:  \n**{}** \n\n'.format(
            os.path.join(met_dir, 'lifeline_out')))
        # f.write('  \n\n')
        f.write('![kapl_plot](' + kapl_plot + ') **{}**   \n'.format(
            os.path.basename(kapl_plot)))
        f.write('  \n\n')
        # f.write('  \n\n')
        # drug_list, nr_cases = return_lifeline_plot_drug_list(kapl_plot)
        # f.write(f'- __The Kaplan Meier plot includes {str(nr_cases)} cases')
        # f.write(' and invokes therapies with:__\n')
        # for drug in drug_list:
        #     f.write(f'  - {drug}\n')
        kapl_plot = kapl_plot.replace('.pdf', '_complement.pdf')
        drug_list, nr_cases = return_lifeline_plot_drug_list(kapl_plot)
        f.write(f'### The __complement__ plot includes {str(nr_cases)} cases')
        # f.write(' and invokes therapies with:__\n')
        # f.write(' and invokes therapies with:__\n')
        # for drug in drug_list:
        #     f.write(f'  - {drug}\n')
        f.write('  \n\n')
        f.write('![kapl_plot](' + kapl_plot + ') **{}**   \n'.format(
            os.path.basename(kapl_plot)))
        # f.write('  \n\n')
        f.write('  \n\n')
        f.write(r'\pagebreak')
        f.write('  \n\n')
        count += 1
    # ############ END plotting kaplan meyer of a choosen postion of DMR:


def create_pdf_from_DF(OUTPUT_PATH, SCRIPT_PATH,
                       DRUGS_title, PROJECT_title, DF, tex_file_name, name=''):

    logger = set_logger.set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title)
    tex_file_path = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                 tex_file_name)
    DF.name = name
    try:
        if isinstance(DF, pd.DataFrame):
            DF.to_latex(tex_file_path, index=False)
        else:
            DF.to_latex(tex_file_path)
        log_file = os.path.join(PROJECT_title, DRUGS_title, tex_file_name)
        logger.info('REPORT_17:\t{}'.format(log_file))
    except FileNotFoundError:
        print('Could not save {}'.format(tex_file_path))

    tex_template_path = os.path.join(
        OUTPUT_PATH, PROJECT_title, DRUGS_title, 'template_' + tex_file_name)
    tex_template_name = 'template_' + tex_file_name

    # now we need to write the table into the tex_template_path
    # content of the tex table needs to be flanked by latex commands, read the
    # table and write the content within those directives
    with open(tex_template_path, 'w') as f2, open(tex_file_path, 'r') as f3:
        f2.write(r'\documentclass[12pt, a4paper]{article}' + '\n')
        f2.write(r'\usepackage{booktabs}' + '\n')
        f2.write(r'\begin{document}' + '\n')
        f2.write(r'\pagenumbering{gobble}' + '\n')
        for line in f3:
            f2.write(line)
        f2.write(r'\end{document}' + '\n')

    # # the total amount shall be inserted as last value, therefore open file
    # # again, count lines and insert total at len()-3
    # with open(tex_template_path, 'r') as f2:
    #     contents = f2.readlines()

    # contents.insert(
    #     len(contents)-3, ' & & & & total: ' + str(total_cases
    #                                               ) + r'\\  ' + '\n')
    # with open(tex_template_path, 'w') as f2:
    #     contents = ''.join(contents)
    #     f2.write(contents)

    os.chdir(os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title))
    system_call = ['pdflatex', tex_template_name]
    try:
        subprocess.check_call(system_call)
        log_file = os.path.join(
            PROJECT_title,
            DRUGS_title, tex_template_name)
        logger.info('REPORT_17:\t{}'.format(log_file))
        log_file = os.path.join(
            PROJECT_title,
            DRUGS_title, tex_template_name.replace('.tex', '.pdf'))
        logger.info('REPORT_17:\t{}'.format(log_file))
    except subprocess.CalledProcessError:
        print('\nfailed to create the latex_template.pdf', end='')
        print(', is pdflatex available?')
        os._exit(0)

    # delete the auxiliary files
    os.remove(tex_template_name.replace('.tex', '.aux'))
    os.remove(tex_template_name.replace('.tex', '.log'))

    # now crop the margins around the created pdf:
    pdf_name = tex_template_name.replace('.tex', '.pdf')
    pdfcrop = os.path.join(SCRIPT_PATH, 'pdfcrop')
    system_call = [pdfcrop, pdf_name, pdf_name]
    try:
        subprocess.check_call(system_call)
    except subprocess.CalledProcessError:
        print('failed to create the ' + pdf_name + ', is pdfcrop in the \
              SCRIPT_PATH ?')
        os._exit(0)

    # if the table is too big, the first page is empty and the second holds the
    # cropped table, just take the last page of the document to be the
    # latex_template.pdf, which then should be included in the REPORT.pdf
    pdf = PdfFileReader(pdf_name)
    if pdf.getNumPages() > 1:
        pdf_writer = PdfFileWriter()
        pdf_writer.addPage(pdf.getPage(pdf.getNumPages() - 1))
        with open(pdf_name, 'wb') as out:
            pdf_writer.write(out)
    return(os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title, pdf_name))


def populate_can_gen_dict(OUTPUT_PATH,
                          PROJECT_title, DRUGS_title, gender_list, met_dir):
    '''
    cand_hash:
    {'CAND': {'both': '1', 'male': '12'}, 'cand': {'both': '4', 'male': '12'}}
    CAND_gender_dict:
    {'both': '4', 'male': '12'}
    '''

    cand_hash = {}
    for search_flag in ['CAND', 'cand']:
        CAND_gender_dict = {}
        for gender in gender_list:
            glob_path = glob.glob(
                os.path.join(
                    OUTPUT_PATH,
                    PROJECT_title,
                    DRUGS_title,
                    gender,
                    met_dir,
                    'boxplot_range',
                    'boxplot*{}.pdf'.format(search_flag)))
            CAND_found = 0
            for path in glob_path:
                CAND_found += 1
                CAND_gender_dict.update({gender: str(CAND_found)})
        cand_hash.update({search_flag: CAND_gender_dict})
    return (cand_hash, CAND_gender_dict)


def return_lifeline_plot_drug_list(life_line_plot_pdf):
    table_path = life_line_plot_pdf.replace('pdf', 'tsv')
    kapl_DF = pd.read_csv(table_path, sep='\t')
    nr_cases = len(kapl_DF)
    kapl_DF['drug'] = kapl_DF['drug'].replace(np.nan, 'Not Available')
    drug_list = list(set(kapl_DF['drug'].tolist()))
    return (sorted(drug_list), nr_cases)
