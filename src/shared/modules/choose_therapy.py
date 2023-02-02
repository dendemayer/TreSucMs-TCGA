import pandas as pd
import os
import subprocess       # to clear the console after choose steps


def Choose_project():
    '''
    interactively requesting the Projects that shall be applied to the deseq
    approach
    '''
    # print("Choose your project out of the list:")
    all_projects_list = {'TCGA-CESC': 'Cervical Squamous Cell Carcinoma and' +
                         ' Endocervical Adenocarcinoma',
                         'TCGA-HNSC': 'Head and Neck Squamous Cell Carcinoma',
                         'TCGA-LUSC': 'Lung Squamous Cell Carcinoma',
                         'TCGA-ESCA': 'Esophageal Carcinoma',
                         'TCGA-BRCA': 'Breast Invasive Carcinoma',
                         'TCGA-GBM': 'Glioblastoma Multiforme',
                         'TCGA-OV': 'Ovarian Serous Cystadenocarcinaoma',
                         'TCGA-LUAD': 'Lung Adenocarcinoma',
                         'TCGA-UCEC': 'Uterine Corpus Endometrial Carinoma',
                         'TCGA-KIRC': 'kindney renal clear cell carcinoma',
                         'TCGA-LGG': 'brain lower grade glioma',
                         'TCGA-THCA': 'thyroid carcinoma',
                         'TCGA-PRAD': 'prostate adenocarcinoma',
                         'TCGA-SKCM': 'skin cutaneous melanoma',
                         'TCGA-COAD': 'colon adenocarcinoma',
                         'TCGA-STAD': 'stomach adenocarcinoma',
                         'TCGA-BLCA': 'bladder urothelial carcinoma',
                         'TCGA-LIHC': 'liver hepatocellular carcinoma',
                         'TCGA-KIRP': 'kidney renal papillary cell carcinoma',
                         'TCGA-SARC': 'sarcoma',
                         'TCGA-PAAD': 'pancreatic adenocarcinoma',
                         'TCGA-PCPG': 'pheochromocytoma and paraganglioma',
                         'TCGA-READ': 'rectum adenocarcinoma',
                         'TCGA-TGCT': 'testicular germcelltumors',
                         'TCGA-THYM': 'thymoma',
                         'TCGA-KICH': 'kidney chromophobe',
                         'TCGA-ACC': 'adrenochordical carcinoma',
                         'TCGA-MESO': 'mesothelioma',
                         'TCGA-UVM': 'uveal melanoma',
                         'TCGA-DLBC': 'lymphoid neoplasm diffuse large b-cell'
                         + ' lymphoma',
                         'TCGA-UCS': 'uterine carcinoma',
                         'TCGA-CHOL': 'cholangiocarcinoma'}

    project_list = []
    print('\nwhich projects' +
          ' do you want to include in your analysis:\n')
    for i in range(0, len(all_projects_list.keys())):
        print('{:>2}:\t{:<20}{:<70}'.format(
            i, [*all_projects_list][i], all_projects_list[
                [*all_projects_list][i]]))

    while True:
        inp = input(
            'enter your choices one by one,' +
            ' if you are finnished, just type "Enter": ')
        if inp == "":
            break
        else:
            if int(inp) not in range(0, len(all_projects_list)):
                print("please choose within {}".format(
                    range(0, len(all_projects_list)-1)))
                continue
            project_list.append([*all_projects_list][int(inp)])
    subprocess.check_call('clear')
    # make a set, in case a project was selected twice:
    project_list = set(project_list)
    project_list = list(project_list)
    print('\nyou choose:')
    print('PROJECTS:\t', project_list)
    return(project_list)


def Choose_drugs(SCRIPT_PATH, PROJECTS):
    '''
    :param: SCRIPT_PATH: path to the DESeq2_pipeline repo
    :type: SCRIPT_PATH: str
    :param: PROJECTS: list of projects chosen
    :type: PROJECTS: list of str

    interactively requesting the drugs which shall be applied to the deseq
    approach
    '''
    DF_drug_table = pd.read_csv(os.path.join(
        SCRIPT_PATH, os.path.pardir, 'resources/DRUG_combi_frequency_all.tsv'),
        sep='\t', index_col=0)
    # DF_drug_table = DF_drug_table.loc[:, PROJECTS].dropna(how='any')
    # by dropping any na, just drugs are provided which were applied in every
    # project, leverage that by throwing out drugs which were applied once
    # maximal:
    DF_drug_table = DF_drug_table.loc[:, PROJECTS].dropna(how='all')
    filt = DF_drug_table.max(axis=1) > 1
    DF_drug_table = DF_drug_table[filt]
    DF_drug_table.fillna(0, inplace=True)
    DF_drug_table = DF_drug_table.astype(int)
    '''
                            TCGA-HNSC   TCGA-LUSC  ...
    cisplatin               65.0        1.0  ...
    carboplatin,paclitaxel  26.0        14.0  ...
    '''
    drug_hash = {}
    '''
    {0: 'cisplatin', 1: 'carboplatin,paclitaxel',...
    '''
    for i in range(0, len(DF_drug_table.index)):
        if DF_drug_table.index[
                i] == '[not available]' or DF_drug_table.index[
                    i] == '[unknown]':
            continue
        drug_hash.update({i: DF_drug_table.index[i]})

    project_frequency_hash = {}
    '''
    {0: 'TCGA-CESC: 104.0 TCGA-LUSC: 1.0', 1: 'TCGA-CESC: 5.0 TCGA-LUSC: 14.0',
    '''
    for i in drug_hash:
        first_proj = True
        for project in PROJECTS:
            if first_proj:
                project_frequency_hash.update(
                    {i: [project + ':', str(
                        DF_drug_table.loc[drug_hash[i], project])]})
                first_proj = False
            else:
                project_frequency_hash[i].extend(
                    [project + ':', str(DF_drug_table.loc[drug_hash[i],
                                                          project])])
    for i in project_frequency_hash:
        project_frequency_hash[i] = ' '.join(project_frequency_hash[i])

    print('\nwhich therapy approach' +
          ' do you want to include in your analysis:\n')
    for i in drug_hash:
        print('{:>2}: {:<41}{:<50}'.format(
            i, drug_hash[i], project_frequency_hash[i]))

    key_list = []
    '''
    this list holds the (number) keys for drug_hash and project_frequency_hash
    '''
    print('')
    while True:
        inp = input(
            'enter your choices one by one,' +
            ' if you are finnished, just type "Enter": ')
        if inp == "":
            break
        else:
            if int(inp) not in range(0, len(drug_hash)):
                print("please choose within {}".format(
                    range(0, len(drug_hash) - 1)))
                continue
            key_list.append(inp)
    subprocess.check_call('clear')
    drug_list = []
    for i in key_list:
        # print('Nr.: {} \t {}'.format(i, drug_hash[int(i)]))
        # create a list of the drugcombinations which shall be applied to the
        # deseq run:
        drug_list.append(drug_hash[int(i)])

    return(drug_list)


def Choose_path_and_option(OUTPUT_PATH, PROJECTS, DRUGS, function,
                           SCRIPT_PATH, analyse_data, download_data,
                           threshold):
    '''
    :param: OUTPUT_PATH: path for DESeq2 pipeline outputs
    :type: OUTPUT_PATH: str
    :param: PROJECTS: list of projects chosen
    :type: PROJECTS: list of str
    :param: DRUGS: applied drug(s)
    :type: DRUGS: list of str
    :param: function: applied functions
    :type: function: int
    :param: SCRIPT_PATH: path to the DESeq2_pipeline repo
    :type: SCRIPT_PATH: str

    interactively requesting whether download steps, analysis or both should be
    performed
    '''
    while True:
        print('\nyou did not choose a analyse or download option,' +
              ' you can choose now out of:\n')
        print('0:\t download count data for {}'.format(str(PROJECTS)))
        print('1:\t perform DESeq2 run on {} with {}'.format(
            str(PROJECTS), str(DRUGS)) +
            '\n\tcount data must be already be downloaded in your output' +
            ' path')
        print('2:\t perform download and DESeq2 run')
        input_choice_analyse = input()
        if input_choice_analyse == '0':
            function = (1, 2, 3, 4, 5, 6)
            download_data = True
            while True:
                print('\npress enter if the default path:\n' +
                      '\n{}\n\nshall be set as OUTPUT_PATH path'.format(
                          OUTPUT_PATH) +
                      '\nor type in the OUTPUT_PATH:')
                input_choice = input()
                if input_choice == '':
                    break
                else:
                    OUTPUT_PATH = input_choice
                    print(
                        'the OUTPUT_PATH is set to: \
                        \n{}'.format(OUTPUT_PATH))
                    break
            break
        if input_choice_analyse == '1':
            function = (7, 8, 9, 10, 11, 12, 13, 14, 15)
            analyse_data = True
            while True:
                print('\npress enter if the default path:\n' +
                      '\n{}\n\nshall be set as OUTPUT_PATH path'.format(
                          OUTPUT_PATH) +
                      '\nor type in the OUTPUT_PATH:')
                input_choice = input()
                if input_choice == '':
                    break
                else:
                    OUTPUT_PATH = input_choice
                    print(
                        'the OUTPUT_PATH is set to \n{}'.format(OUTPUT_PATH))
                    break
            break
        if input_choice_analyse == '2':
            function = (100,)
            analyse_data = True
            download_data = True
            while True:
                print('press enter if the default path:\n' +
                      ' \n{}\n\nshall be set as OUTPUT_PATH path'.format(
                          OUTPUT_PATH) +
                      '\nor type in the OUTPUT_PATH')
                input_choice = input()
                if input_choice == '':
                    break
                else:
                    OUTPUT_PATH = input_choice
                    print(
                        'the OUTPUT_PATH is set to \n{}'.format(OUTPUT_PATH))
                    break
            break
        else:
            print('please provide your choice with the numbers' +
                  ' 0 or 1 or 2')
            continue
    first = True
    print(
        'do you want to keep the threshold with {}'.format(str(threshold))
        + ', then press enter, if not\n')
    while True:
        input_choice = input(
            'enter your choices one by one,' +
            ' if you are finnished, just type "Enter": ')
        if input_choice == '':
            break
        else:
            if first:
                threshold = (int(input_choice),)
            else:
                # threshold_old = threshold
                threshold = tuple(set(threshold + (int(input_choice),)))
            first = False
            continue

    while True:
        subprocess.check_call('clear')
        print('your analysis will start with:\n')
        print('PROJECTS:\t', PROJECTS)
        print('DRUGS:\t\t', DRUGS)
        print('OUTPUT_PATH:\t', OUTPUT_PATH)
        print('threshold\t', threshold)
        # print('SCRIPT_PATH:\t', SCRIPT_PATH)
        if input_choice_analyse == '0':
            print('\n- download count data for the selected projects -\n')
        elif input_choice_analyse == '1':
            print('\n- perform just DESeq2 analysis, the count data is' +
                  ' already downloaded -\n')
        else:
            print('\n- perform download of count data and DESeq2 run -\n')
        print('press Enter to start the analysis')
        input_choice = input()
        if input_choice == '':
            break
        else:
            continue
    return((OUTPUT_PATH, function, analyse_data, download_data, threshold))
