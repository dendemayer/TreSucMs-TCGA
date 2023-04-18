import subprocess       # to call pandoc
import os
import glob
from create_matrix_new import set_logger
import re
import pandas as pd
import seaborn as sns
from PyPDF2 import PdfFileReader, PdfFileWriter


# fct 13 in main
def create_report_pdf(OUTPUT_PATH, DRUGS, SCRIPT_PATH, PROJECT_DRUG_UUID,
                      threshold, cutoff, DRUGS_title):
    """
    :param: OUTPUT_PATH: path for DESeq2 pipeline outputs
    :type: OUTPUT_PATH: str
    :param: DRUGS: applied drug(s)
    :type: DRUGS: list of str
    :param: SCRIPT_PATH: path to the DESeq2_pipeline repo
    :type: SCRIPT_PATH: str
    :param: PROJECT_DRUG_UUID: hash holding a project to the UUID of the\
    belonging drugtable
    :type: PROJECT_DRUG_UUID: dict
    :param: threshold: parameter for the lifeline plots helping for the\
    classification of expression data
    :type: threshold: int

    creating a report file with all visual outputs created in an analysis,
    saved at OUTPUT_PATH/PROJECT_title/DRUGS_title/REPORT.pdf


    .. _function_13:
    .. code-block:: bash

        # for executing this step via terminal, issue -f 13 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 13 -t 0
        # if you want to include several thresholds you had created in prior
        # with the KaplanMeier function (-f 11 or within -A), apply them here
        # as well, for
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 13 -t 0 -t 50 -t
        100

    """
    # DRUGS_title = '_'.join(sorted(DRUGS))
    if isinstance(PROJECT_DRUG_UUID, dict):  # a multi project as dict
        projects = []
        for project in PROJECT_DRUG_UUID:
            projects.append(project)
        PROJECT_title = '_'.join(sorted(map(str.upper, projects)))
    else:
        PROJECT_title = PROJECT_DRUG_UUID
    pdf_file = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                            'REPORT.pdf')
    md_file = os.path.join(
        OUTPUT_PATH, PROJECT_title, DRUGS_title, 'REPORT.md')
    threshold = sorted(list(threshold))
    with open(md_file, 'w') as f:
        f.write('# Report file for the DESeq2_pipeline  \n\n')
        f.write('## Query:  \n\n')
        f.write('- PROJECT:\t')
        if isinstance(PROJECT_DRUG_UUID, dict):  # a multi project as dict
            f.write('__' + '; '.join(projects) + '__  \n')
        else:
            f.write('__' + PROJECT_title + '__  \n')
        f.write('- DRUGS:\t')
        f.write('__' + '; '.join(DRUGS) + '__  \n')
        # don't print threshold when it's not passed, and therefore is 0 and
        # therefore does nothing:
        if not len(threshold) == 1 and threshold[0] == 0:
            f.write('- threshold:\t')
            f.write('__' + str(threshold) + '__')
            f.write('  \n')
        f.write('- OUTPUT_PATH:\t')
        f.write('__')
        f.write(os.path.basename(OUTPUT_PATH))
        f.write('__  \n')
        # f.write('  \n')
        f.write('- SCRIPT_PATH:\t')
        f.write('__')
        f.write(os.path.basename(SCRIPT_PATH))
        f.write('__  \n')
        name_copied = os.path.join(
            'Snakes', 'snakemake_config_'
            + PROJECT_title + '_' + DRUGS_title + '.yaml')
        f.write(r'- your configfile within the SCRIPT\_PATH:')
        f.write('\n\n   ')
        f.write(name_copied)
        f.write('   \n\n')
        if cutoff > 0.0:
            f.write('- cutoff Parameter:\t')
            f.write('__')
            # no need for decimal point if whole number
            if cutoff % 1 == 0:
                f.write(str(int(cutoff)))
                cut_str = str(int(cutoff))
            else:
                f.write(str(cutoff))
                cut_str = str(cutoff)
            f.write('__  \n')
        f.write('   \n\n')
        # #################START Overview of applied cases: ###################

        # create_case_count gives a list of 3 values:
        # case_count_list[0] -> bool, if a deseq run was actually performed
        # case_count_list[1] -> total cases
        # case_count_list[2] -> DF to make bins for age at diagnosis
        # case_count_list[3] -> DF_info_counts frequency counts for all invoked
        # cases
        # create the DF:
        # take care that just non complement cases are included! DONE
        case_count_list = create_case_count(OUTPUT_PATH, PROJECT_title,
                                            DRUGS_title, SCRIPT_PATH)
        if case_count_list[0]:
            alive_count = case_count_list[2].loc[
                :, 'vital_status'].value_counts()['alive']
            dead_count = case_count_list[2].loc[
                :, 'vital_status'].value_counts()['dead']
            # Index(['vital_status', 'gender', 'PROJECT', 'drugnames',
            # 'age_at_diagnosis'], dtype='object')
            # create the pdf table, which is inserted here in the report:
            tex_file_name = 'latex_all_invoked_cases.tex'
            pdf_name = create_pdf_from_DF(
                OUTPUT_PATH, SCRIPT_PATH, DRUGS_title,
                PROJECT_title, case_count_list[3], tex_file_name, '# of cases')
            f.write('## Overview of your applied cases to the analysis  \n\n')

            # latex_output = os.path.join(OUTPUT_PATH, PROJECT_title,
            #                             DRUGS_title, 'latex_template.pdf')
            f.write(
                '![factors of your applied cases](' + pdf_name + ')  \n\n')
            f.write('__total cases: ' + str(case_count_list[1]) + '__  \n')
            f.write('__alive cases: ' + str(alive_count) + '__  \n')
            f.write('__dead cases: ' + str(dead_count) + '__  \n\n')
            gender_freq = len(
                case_count_list[2].loc[:, 'gender'].value_counts())
            if gender_freq > 1:
                female_count = case_count_list[2]['gender'].value_counts()[
                    'female']
                male_count = case_count_list[2][
                    'gender'].value_counts()['male']
                f.write('__female cases: ' + str(female_count) + '__  \n')
                f.write('__male cases: ' + str(male_count) + '__  \n\n')
        else:
            f.write('## Overview of your applied cases to the analysis  \n\n')
            f.write('### not enough cases available for a proper DESeq2 run'
                    + '  \n\n')

            # ### start of cutoff added cases:
        if cutoff > 0:
            path = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                'cutoff_cases_add_' + str(cutoff) + '.tsv')
            try:
                cutoff_add_DF = pd.read_csv(path, sep='\t')
                # ### the same frequency table like the case_count_list[1], but
                # just the cutoff dependent added cases:
                DF_freq_cutoff_add = cutoff_add_DF.T.loc[
                    ['vital_status', 'gender', 'PROJECT', 'drugnames'],
                    :].T.value_counts()
                tex_file_name = 'latex_cases_cutoff_add.tex'
                # create the pdf and include it in the report:
                pdf_name = create_pdf_from_DF(OUTPUT_PATH, SCRIPT_PATH,
                                              DRUGS_title, PROJECT_title,
                                              DF_freq_cutoff_add,
                                              tex_file_name, '# of cases')

                f.write(r'\pagebreak')
                f.write('\n')
                f.write('## Cases converted  due to the cutoff of '
                        + cut_str + ' years   \n\n')
                f.write('While having a vital status of dead, those cases ')
                f.write('outlived the cutoff limit and are ')
                f.write('therefore taken as alive cases into account to your ')
                f.write('analyses.  \n\n')
                f.write(
                    '![factors of added cases]('
                    + pdf_name + ')  \n\n')

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
                f.write('## Survivaltimes of the added cases (cutoff = {}\
                        ):   \n\n'.format(cut_str))
                f.write(
                    '![survival_times of added cases]('
                    + pdf_name + ')  \n\n')

            except FileNotFoundError:
                print('not able to find {}'.format(path))

        # #################END Overview of applied cases: ###################

        # ###################START barplot age at diagnose ####################
        # in
        # OUTPUT_PATH/PROJECT_title/Lifeline_plots_general/lifelines_table.tsv
        # we have every case of the project with age_at_diagnosis, restrict it
        # to the cases (drugquery) which is used here:

        if case_count_list[0]:
            DF_info_table = case_count_list[2]
            plot = sns.displot(
                data=DF_info_table,
                x='age_at_diagnosis', col='gender', kde=True)
            plot_path = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                     'displot.pdf')
            plot.savefig(plot_path)

            f.write(r'\pagebreak')
            f.write('  \n\n')
            f.write('## Distribution of age at diagnosis, all cases:  \n\n')
            f.write('![dis_plot](' + plot_path + ') \
                    {} in \n\n **{}**'.format(
                        os.path.basename(plot_path),
                        os.path.split(plot_path)[0].replace(OUTPUT_PATH +
                                                            os.path.sep, '')))
            f.write('\n\n')
        # ###################END barplot age at diagnose ######################

        # #############DESeq2 results#################################
        f.write(r'\pagebreak')
        f.write('  \n\n')
        f.write('## DESeq2 results, located in:  \n\n')
        f.write(
            '### {}/...  \n\n'.format(
                os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title)))
        f.write('* everything in:  \n')
        f.write('    * .../DESeq2_out_DRUG_**combi/** combines both, male and\
                female  \n')
        f.write('    * .../DESeq2_out_DRUG_combi_**female/** is restricted \
                to female cases  \n')
        f.write('    * .../DESeq2_out_DRUG_combi_**male/** is restricted\
                to male cases \n\n')

        for point_cloud in sorted(glob.glob(os.path.join(
            OUTPUT_PATH, PROJECT_title, DRUGS_title,
                'DESeq2_out_DRUG_combi*/DESeq2_results.pdf'))):
            dir_name = os.path.dirname(point_cloud)
            dir_name = os.path.split(dir_name)[1]
            f.write('![point cloud](' + point_cloud + ') \
                    {} in **{}**'.format(
                        os.path.basename(point_cloud), dir_name))
            f.write('\n\n')

        # #############gene tables with symbols################################
        f.write(r'\pagebreak')
        f.write('## Table with gene symbols:  \n\n')
        f.write('### The following tables contain\
                the most differentiate expressed genes between the groups of \
                alive and dead:  \n\n')
        f.write(r'* Genes in tables with \*__DECREASE__\* within the ')
        f.write('filenames are overall more expressed among the __dead__ \
                cases  \n\n')
        f.write(r'* Genes in tables with \*__INCREASE__\* within the ')
        f.write('filenames are overall more expressed among the\
                __alive__ cases  \n\n')
        f.write('* The complete table with exhaustive variables for the DESeq2\
                results like p_value and baseMean, as well as the belonging \
                gene infos like full gene name or the ensembl transcript\
                description can be found at the respective referred file\
                names.\n\n')
        # f.write('__TODO / Question__: how do we actually know if e.g. STRIT1\
        #     is stronger expressed in dead, the log2foldchange could also\
        #     emerge because of a lag of expression in alive cases...  \n\n')
        f.write(r'\pagebreak')
        f.write('  \n\n')

        for path in sorted(
            glob.glob(
                os.path.join(
                    OUTPUT_PATH,
                    PROJECT_title, DRUGS_title, 'DESeq2_out_DRUG_combi*'))):
            table_par_dir = os.path.basename(path)
            # f.write("\ntable in:  \n{}  \n\n".format(table_par_dir))
            for table in sorted(
                glob.glob(
                    os.path.join(
                        path, 'results_*_10_lgfch_ENSG_and_Gene_info.tsv'))):
                table_name = os.path.basename(table)
                ENSG_GENE_DF = pd.read_csv(table, sep='\t')
                ENSG_GENE_DF = ENSG_GENE_DF.sort_values(
                    by=['log2FoldChange'], ascending=False, ignore_index=True)
                # # print(ENSG_GENE_DF.loc[:, [
                #     # 'Unnamed: 0', 'symbol', 'log2FoldChange',
                #     'type_of_gene',
                #     # 'genomic_pos.chr']])

                f.write('ENSG description|gene symbol|log2FoldChange|padj|\
                        type of gene\n')
                f.write('|-|-|-|-|-|  \n')
                # padj is sometimes with e- notation, round that correctly:
                for i in range(0, len(ENSG_GENE_DF)):
                    gene_symbol = str(ENSG_GENE_DF.loc[i, 'symbol'])
                    if gene_symbol == 'nan':
                        gene_symbol = '-'
                    gene_type = str(ENSG_GENE_DF.loc[i, 'type_of_gene'])
                    if gene_type == 'nan':
                        gene_type = '-'
                    # shorten the scientific e notation
                    ENSG_name = str(ENSG_GENE_DF.loc[i, 'Unnamed: 0']
                                    ).replace('ENSG', 'ENSG-')
                    if re.search('e', str(ENSG_GENE_DF.loc[i, 'padj'])):
                        f.write('|' + ENSG_name + '|' + gene_symbol + '|' +
                                str(
                                    round(
                                        ENSG_GENE_DF.loc[
                                            i, 'log2FoldChange'], ndigits=2))
                                + '|' + '{:0.3e}'.format(
                                    ENSG_GENE_DF.loc[i, 'padj']) + '|' +
                                gene_type + '|'
                                + '  \n')
                    else:
                        f.write('|' + ENSG_name + '|' + gene_symbol + '|' +
                                str(round(
                                    ENSG_GENE_DF.loc[
                                        i, 'log2FoldChange'],
                                    ndigits=2)) + '|' +
                                str(
                                    round(
                                        ENSG_GENE_DF.loc[
                                            i, 'padj'], ndigits=4)) + '|' +
                                gene_type +
                                '|' + '  \n')
                f.write('Table: ' + '__' + table_name + '__' + '\n' +
                        'in ' + table_par_dir + '  \n\n')
            f.write(r'\pagebreak')
            f.write('  \n\n')

        # #############Heatmaps##########################################
        f.write(r'\pagebreak')
        f.write('## Heatmaps, located in:  \n\n')
        f.write(
            '### {}/...  \n\n'.format(
                os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title)))
        f.write('* everything in:  \n')
        f.write('    * .../DESeq2_out_DRUG_**combi/** combines both, male and\
                female  \n')
        f.write('    * .../DESeq2_out_DRUG_combi_**female/** is restricted to\
                female cases  \n')
        f.write('    * .../DESeq2_out_DRUG_combi_**male/** is restricted to\
                male cases  \n\n')

        f.write('### Explanation of the file naming:  \n\n\n')

        f.write('#### example filename:  \n\n\n')
        f.write('* __DESeq2_heatmap_log2fINCREASE_norm_60.pdf__  \n\n')

        f.write('#### composition of the filenames:  \n\n')
        f.write(r'* **\*\_log2fINCREASE\_\*** :')
        f.write('\n')
        f.write('    * genes with highest log2fold_change,\
                increased sorted \n')
        f.write('    * those are genes with higher expression in\
                alive cases:  \n')
        f.write('        * INCREASE -> higher expression in alive  \n')

        f.write(r'* **\*\_log2fDECREASE\_\*** :')
        f.write('  \n')
        f.write('    * genes with highest\
                log2fold_change, decreased sorted \n')
        f.write('    * those are genes with higher expression\
                in dead cases:  \n')
        f.write('        * DECREASE -> higher expression in dead  \n\n')

        # increase: alive up, dead down
        # decrease, alive down, dead up

        f.write(r'* **\*\_norm\_\***:')
        f.write('\n')
        f.write('    * normalized counts, extracted from the DESeq2 dataset\
                with counts() function, e.g. counts(dds,\
                                                    normalized=TRUE)\n\n')

        # f.write(r'* **\*\_nt\_\***:')
        # f.write('  \n')
        # f.write('    * norm_transform: shifted logarithm\
        #         transformation: ntd <-\
        #         normTransform(dds)  \n\n')

        # f.write(r'* **\*\_raw_counts\_\***:')
        # f.write('\n')
        # f.write('    * raw_counts: the count tables in original form\
        #         loaded from TCGA, without any transformation performed\
        #         by DESeq2  \n\n')

        f.write(r'* **\*\_60\***:')
        f.write('  \n')
        f.write('    * showing 60 genes in the plot, sorted in increasing or\
                decreasing order after log2foldchange  \n\n')

        for heatmap in sorted(glob.glob(os.path.join(
            OUTPUT_PATH, PROJECT_title, DRUGS_title,
                'DESeq2_out_DRUG_combi*/DESeq2_heatmap_log*norm*.pdf'))):
            dir_name = os.path.dirname(heatmap)
            dir_name = os.path.split(dir_name)[1]

            f.write('![heatmap](' + heatmap + ') \
                    {} in  \n**{}**'.format(
                        os.path.basename(heatmap), dir_name))
            f.write('\n\n')

        f.write(r'\pagebreak')
        f.write('\n')
        # f.write('* ntd: **shifted logarithm transformation**.')
        # f.write('    * The shifted logarithm\
        #         has elevated standard deviation in the lower count range ')
        # for heatmap_ntd in sorted(glob.glob(os.path.join(
        #     OUTPUT_PATH, PROJECT_title, DRUGS_title,
        #         'DESeq2_out_DRUG_combi*/DESeq2_heatmap_ntd*.pdf'))):
        #     dir_name = os.path.dirname(heatmap)
        #     dir_name = os.path.split(dir_name)[1]

        #     f.write('![heatmap](' + heatmap + ') \
        #             {} in **{}**'.format(os.path.basename(heatmap),
        #             dir_name))
        #     f.write('\n\n')

        # f.write('* vsd:  **Variance stabilizing transformation**  \n')
        # f.write('    * variance stabilized data, the standard deviation is\
        #         roughly constant along the whole dynamic range  \n')

        # for heatmap_ntd in sorted(glob.glob(os.path.join(
        #     OUTPUT_PATH, PROJECT_title, DRUGS_title,
        #         'DESeq2_out_DRUG_combi*/DESeq2_heatmap_vsd*.pdf'))):
        #     dir_name = os.path.dirname(heatmap)
        #     dir_name = os.path.split(dir_name)[1]

        #     f.write('![heatmap](' + heatmap + ') \
        #             {} in **{}**'.format(os.path.basename(heatmap),
        #             dir_name))
        #     f.write('\n\n')

        # ###############PCA plots##########################################
        f.write('## PCA plots:   \n\n')
        for pca in sorted(glob.glob(os.path.join(
            OUTPUT_PATH, PROJECT_title, DRUGS_title,
                'DESeq2_out_DRUG_combi*/DESeq2_pca_groups.pdf'))):
            dir_name = os.path.dirname(pca)
            dir_name = os.path.split(dir_name)[1]
            f.write('![pca plot](' + pca + ') \
                    {} in\n\n**{}**'.format(os.path.basename(pca), dir_name))
            f.write('\n\n')

        # ##############Lifeline plots########################################

        f.write('## Kaplan Meier plots for: {}  \n\n'.format(PROJECT_title))
        f.write('### every case of your project, without\
                restriction of your drug selection  \n\n')
        # leave the cumulative density out, its just mirrored to the other plot
        for life_line in sorted(glob.glob(os.path.join(
            OUTPUT_PATH, PROJECT_title,
                'Lifeline_plots_general/lifelines_survival*.pdf'))):
            f.write('![kaplan meier plot](' + life_line + ') \
                    {}'.format(os.path.basename(life_line)))
            f.write('  \n\n')
        # # life_line_general = os.path.join(
        #     # OUTPUT_PATH,
        #     # PROJECT_title,
        #     # 'Lifeline_plots_general/lifelines_survival_every_case.pdf')
        # # f.write('![kaplan meier plot](' + life_line_general + ') \
        #         # {}'.format(os.path.basename(life_line_general)))
        # # f.write('\n\n')

        # # first loop: threshold
        # #   second loop: male, female or both
        # #       third loop: ICREASE or DECREASE

        # first loop: threshold
        for threshold_iter in threshold:
            f.write('## Kaplan Meier plots for: {} and {} with\
                    threshold {}\n\n'.format(
                        PROJECT_title, DRUGS_title, threshold_iter))
            for gender_iter in sorted(glob.glob(
                    os.path.join(
                        OUTPUT_PATH,
                        PROJECT_title,
                        DRUGS_title, 'DRUG_combi*Lifeline_norm_counts'))):

                #   second loop: male, female or both
                gender_iter = os.path.split(gender_iter)[1]

                lifeline_num = len(
                    glob.glob(os.path.join(
                        OUTPUT_PATH, PROJECT_title, DRUGS_title,
                        gender_iter, 'threshold_' + str(threshold_iter,),
                        'lifelines_*INCREASE*.pdf')))
                if gender_iter == 'DRUG_combi_Lifeline_norm_counts':
                    f.write('### males and females combined: first {}\
                            '.format(lifeline_num))
                    f.write(' genes with the highest expression')
                    f.write(r'in __alive__\ cases:')
                    f.write('  \n\n ')
                elif re.search('female', gender_iter):
                    f.write('### female cases only: first {}'.format(
                        lifeline_num))
                    f.write(' genes with the highest expression in')
                    f.write(r' __alive__\ cases:')
                    f.write('  \n\n ')
                else:
                    f.write('### male cases only: first {} genes with the\
                            highest'.format(lifeline_num))
                    f.write(r' expression in __alive__\ cases:')
                    f.write('  \n\n ')

                # third loop: ICREASE or DECREASE
                for INC_iter in sorted(glob.glob(
                        os.path.join(
                            OUTPUT_PATH,
                            PROJECT_title,
                            DRUGS_title, gender_iter, 'threshold_' +
                            str(threshold_iter,),
                            'lifelines_*INCREASE*.pdf'))):
                    life_file = os.path.dirname(INC_iter)  # dir to filename
                    life_file = os.path.split(life_file)  # (head, tail) of
                    # path
                    life_file = os.path.split(life_file[0])  # head an tail out
                    # of head
                    life_file = life_file[1]  # tail
                    f.write('![kaplan meier plot](' + INC_iter + ') {}\
                            in\n\n**{}** (threshold: {})'.format(
                                os.path.basename(INC_iter), life_file,
                                threshold_iter))
                    f.write('\n\n')

                lifeline_num = len(
                    glob.glob(os.path.join(
                        OUTPUT_PATH, PROJECT_title, DRUGS_title,
                        gender_iter, 'threshold_' + str(threshold_iter,),
                        'lifelines_*DECREASE*.pdf')))
                if gender_iter == 'DRUG_combi_Lifeline_norm_counts':
                    f.write('### males and females combined: {}'.format(
                        lifeline_num))
                    f.write(' genes with the highest expression ')
                    f.write(r'in __dead__\ cases')
                    f.write(':  \n\n')
                elif re.search('female', gender_iter):
                    f.write('### female cases only: {}'.format(lifeline_num))
                    f.write('genes with the highest expression in')
                    f.write(r'__dead__\ cases:')
                    f.write('  \n\n')
                else:
                    f.write('### male cases only: {}'.format(lifeline_num))
                    f.write(' genes with the highest expression')
                    f.write(r'in __dead__\ cases:')
                    f.write('  \n\n ')
                # third loop: ICREASE or DECREASE
                for INC_iter in sorted(glob.glob(
                        os.path.join(
                            OUTPUT_PATH,
                            PROJECT_title,
                            DRUGS_title, gender_iter, 'threshold_' +
                            str(threshold_iter,),
                            'lifelines_*DECREASE*.pdf'))):
                    life_file = os.path.dirname(INC_iter)  # dir to filename
                    life_file = os.path.split(life_file)  # (head, tail) of
                    # path
                    life_file = os.path.split(life_file[0])  # head an tail out
                    # of head
                    life_file = life_file[1]  # tail
                    f.write('![kaplan meier plot](' + INC_iter + ') {}\
                            in\n\n**{}** (threshold: {})'.format(
                                os.path.basename(INC_iter), life_file,
                                threshold_iter))
                    f.write('\n\n')

        # ##############Drug frequencies############################
        # ## bar plot if len PROJECTS is ==1 ,else heatmap
        f.write('## Drugfrequencies:   \n\n')
        f.write('### Every drug/drugcombination available according to your\
                choice of projects, NOT restricted to your choice of\
                therapeutics :   \n\n')
        f.write('* to keep the plot well-arranged,\
                rows with a maximum value of 1 are filtered out   \n\n')
        if len(PROJECT_DRUG_UUID) == 1 or isinstance(PROJECT_DRUG_UUID, str):
            drug_frequency = os.path.join(
                OUTPUT_PATH, PROJECT_title,
                'DF_3t_both_with_DRUG_combi_frequency.pdf')
        else:
            drug_frequency = os.path.join(
                OUTPUT_PATH, PROJECT_title,
                'DRUG_combi_frequency_heatmap.pdf')
        f.write('![frequency](' + drug_frequency + ') \
                {} for {}'.format(os.path.basename(drug_frequency),
                                  PROJECT_title))
        f.write('\n\n')

        f.close()

    # log the md file:
    logger = set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title)
    log_file = md_file.replace(OUTPUT_PATH + os.path.sep, '')
    # change into the SCRIPT_PATH/iftex s.t. pandoc can find the iftex.sty, if
    # needed
    os.chdir(os.path.join(SCRIPT_PATH, os.path.pardir, 'resources/iftex'))
    pandoc_sequence = ['pandoc', md_file, '-f', 'markdown-implicit_figures',
                       '-t', 'latex', '-o', pdf_file]
    print('calling pandoc with:\n{}'.format(pandoc_sequence))
    # md_file and pdf_file are absolute anyway
    # in case the pdf file exists already, delete it in prior:
    if os.path.isfile(pdf_file):
        subprocess.check_call(['rm', pdf_file])
    try:
        subprocess.check_call(pandoc_sequence)
        # log the pdf
        log_file = pdf_file.replace(OUTPUT_PATH + os.path.sep, '')
        logger.info('REPORT_13:\t{}'.format(log_file))
    except subprocess.CalledProcessError:
        print('failed to create the report pdf, have you latex installed?')
        return


def create_case_count(OUTPUT_PATH, PROJECT_title, DRUGS_title, SCRIPT_PATH):
    '''
    creates a tuple of for:
    # case_count_list[0] -> bool, if a deseq run was actually performed
    # case_count_list[1] -> total cases
    # case_count_list[2] -> DF to make bins for age at diagnosis
    # case_count_list[3] -> DF to make bins for age at diagnosis
    unless no deseq analysis could be performed, then case_count_list[0] is
    false and tuple length is 2
    ## create a short overview table of the cases invoked in the
    analysis:
    all cases included in the deseq analysis can be found in the:
    /OUTPUT_PATH/PROJECT_title/DRUGS_title/DRUG_combi_*_INFO.tsv
    depending on presence of male or female cases, collect all available
    tables and unify on case_id col:

    set_index(keys, drop=True, append=False, inplace=False,
    verify_integrity=False) pd.DataFrame().to_latex
    !!! don't count double!! if female and male is present, the aggregation
    of both must be omitted -> this is done with dropping duplicate case_id
    '''
    info_tables_pre = glob.glob(os.path.join(OUTPUT_PATH, PROJECT_title,
                                             DRUGS_title,
                                             'DRUG_combi_*_INFO.tsv'))
    info_tables = []
    for table in info_tables_pre:
        if not re.search('complement', (os.path.split(table)[1])):
            info_tables.append(table)
    # if no tables available, leave:
    if len(info_tables) == 0:
        return(False, 0)
    # first = True
    DF_list = []
    for table in info_tables:
        temp_DF = pd.read_csv(table, sep='\t')
        temp_DF = temp_DF.T
        temp_DF = temp_DF.set_index(4)
        temp_DF.columns = temp_DF.loc['case_id', :]
        temp_DF.columns.name = ''
        temp_DF = temp_DF.drop(labels='case_id')
        temp_DF.index.name = 'case_id'
        DF_list.append(temp_DF)
    DF_info_table = pd.concat(DF_list)
    # in case duplicate case_id arise, make them column, drop duplicates, set
    # back to index:
    DF_info_table = DF_info_table.reset_index().drop_duplicates(
        'case_id').set_index('case_id')
    DF_info_counts = DF_info_table.value_counts()
    total_cases = DF_info_table.value_counts().sum()
    DF_info_counts.name = '# of cases'
    print("\nDF_info_counts:\n", DF_info_counts)
    # convert the table to tex, the md format is useless,
    # compile tex to pdf with pdflatex, then include the pdf here in the
    # report:
    # ###########
    # to create the age_at_diagnosis barplot, merge on case_id with
    # OUTPUT_PATH/PROJECT_title/Lifeline_plots_general/lifelines_table.tsv
    try:
        all_cases_DF = pd.read_csv(os.path.join(OUTPUT_PATH, PROJECT_title,
                                                'Lifeline_plots_general',
                                                'lifelines_table.tsv'),
                                   sep='\t', index_col='case_id')
    except FileNotFoundError:
        print('Not enough cases available')
        return(False, False, False)
    DF_info_table = pd.merge(DF_info_table, all_cases_DF['age_at_diagnosis'],
                             left_index=True, right_index=True)
    # case_count_list[0] -> bool, if a deseq run was actually performed
    # case_count_list[1] -> total cases
    # case_count_list[2] -> DF to make bins for age at diagnosis
    # case_count_list[3] -> DF to make bins for age at diagnosis
    return(True, total_cases, DF_info_table, DF_info_counts)


def create_pdf_from_DF(OUTPUT_PATH, SCRIPT_PATH,
                       DRUGS_title, PROJECT_title, DF, tex_file_name, name=''):

    logger = set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title)
    tex_file_path = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title,
                                 tex_file_name)
    DF.name = name
    try:
        if isinstance(DF, pd.DataFrame):
            DF.to_latex(tex_file_path, index=False)
        else:
            DF.to_latex(tex_file_path)
        log_file = os.path.join(PROJECT_title, DRUGS_title, tex_file_name)
        logger.info('REPORT_13:\t{}'.format(log_file))
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

    work_path = os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title)
    os.chdir(work_path)
    # use the pdflatex in the SCRIPT_PATH/../bins dir
    # temp_script = os.path.split(SCRIPT_PATH)[0]
    # bin_path = os.path.join(temp_script, 'bins')
    # pdf_latex_bin = os.path.join(bin_path, 'pdflatex')
    # for now use the systems pdflatex in /usr/bin/pdflatex
    system_call = ['pdflatex', tex_template_name]
    try:
        subprocess.check_call(system_call)
        log_file = os.path.join(
            PROJECT_title,
            DRUGS_title, tex_template_name)
        logger.info('REPORT_13:\t{}'.format(log_file))
        log_file = os.path.join(
            PROJECT_title,
            DRUGS_title, tex_template_name.replace('.tex', '.pdf'))
        logger.info('REPORT_13:\t{}'.format(log_file))
    except subprocess.CalledProcessError:
        print('\nfailed to create the latex_template.pdf', end='')
        print(', is pdflatex available?')
        os._exit(0)

    # delete the auxiliary files
    os.remove(tex_template_name.replace('.tex', '.aux'))
    os.remove(tex_template_name.replace('.tex', '.log'))

    # now crop the margins around the created pdf:
    pdf_name = tex_template_name.replace('.tex', '.pdf')
    pdfcrop = os.path.join(SCRIPT_PATH, os.path.pardir, 'bins/pdfcrop')
    pdfcrop_path = os.path.join(SCRIPT_PATH, os.path.pardir, 'bins')
    system_call = [pdfcrop, pdf_name, pdf_name]
    try:
        subprocess.check_call(system_call)
    except subprocess.CalledProcessError:
        print(f'failed to create the {pdf_name} , is ', end='')
        print('pdfcrop in {}'.format(pdfcrop_path))
        os._exit(0)

    # if the table is too big, the first page is empty and the second holds the
    # cropped table, just take the last page of the document to be the
    # latex_template.pdf, which then should be included in the REPORT.pdf
    pdf = PdfFileReader(pdf_name)
    if pdf.getNumPages() > 1:
        pdf_writer = PdfFileWriter()
        pdf_writer.addPage(pdf.getPage(pdf.getNumPages()-1))
        with open(pdf_name, 'wb') as out:
            pdf_writer.write(out)
    return(os.path.join(OUTPUT_PATH, PROJECT_title, DRUGS_title, pdf_name))
