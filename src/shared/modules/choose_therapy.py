import pandas as pd
import os
import subprocess       # to clear the console after choose steps


def Choose_project():
    '''
    Choose_project.

    interactively requesting the Projects that shall be applied to the approach
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
            ' when you are done, simply press "Enter": ')
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
    project_list = sorted(map(str.upper, project_list))
    print('\nyou choose:')
    print('PROJECTS:\t', project_list)
    return(project_list)


def Choose_drugs(PROJECTS):
    """Choose_drugs.

    :param str SCRIPT_PATH: path to the DESeq2_pipeline repo
    :param list of str PROJECTS: list of projects chosen

    interactively requesting the drugs which shall be applied to the deseq
    approach
    """
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
    DESeq2_treatment_heatmap_table.tsv
    '''
    '''
                            TCGA-HNSC   TCGA-LUSC  ...
    cisplatin               65.0        1.0  ...
    carboplatin,paclitaxel  26.0        14.0  ...
    '''
    drug_hash = {}
    '''
    {0: 'cisplatin', 1: 'carboplatin,paclitaxel',...
    '''
    table_path = os.path.join(os.path.dirname(__file__), os.path.pardir, 'resources', 'DESeq2_treatment_heatmap_table.tsv')
    DF_drug_table = pd.read_table(table_path, index_col=0).loc[:, PROJECTS].dropna(how='all')
    filt = DF_drug_table.max(axis=1) > 1
    DF_drug_table = DF_drug_table[filt]
    DF_drug_table.fillna(0, inplace=True)
    DF_drug_table = DF_drug_table.astype(int)
    # the drugtable can be limited to the projects chosen and dropping complete
    # nan rows, and convert remaining nans to 0:
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
            ' when you are done, simply press "Enter": ')
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

    drug_list = sorted(map(str.lower, drug_list))
    return(drug_list)


def update_parameters(parameter, parameter_str):
    """update_parameters.

    :param str or tuple or list of strs parameter: depending on the datatype, the parameter which shall be assign is requested in interactive mode
    :param str parameter: default OUTPUT_PATH can be confirmed or changed in interactive mode
    :param int parameter: cores can be changed in interactive mode
    :param str parameter_str: description of the given parameter

    """
    # setting OUTPUT_PATH:
    if isinstance(parameter, str):
        print(f'do you want to keep the default OUTPUT_PATH of:\n{parameter}')
        print('if so, press ENTER, if not, enter your custom output path:')
        output_input = input()
        if output_input == '':
            return  parameter
        else:
            return  output_input
    # setting cores:
    if isinstance(parameter, int):
        print(f'do you want to keep the default number of cores invoked of 1?')
        print('if so, press ENTER, if not, enter the number of cores:')
        new_cores = input()
        if new_cores == '':
            return parameter
        else:
            return int(new_cores)
    # setting pipelines:
    if isinstance(parameter, list):
        if parameter_str == 'pipelines':
            while True:
                print('which pipeline do you want to include into your analysis')
                print('press ENTER if DESeq2 and metilene (default) or\n1 for DESeq2 or \n2 for metilene')
                pipelines = input()
                if pipelines == '':
                    return parameter
                if pipelines == "1":
                    return ['DESeq2']
                if pipelines == "2":
                    return ['metilene']
                else:
                    print('just press ENTER or apply\n1 for DESeq2 or \n2 for metilene')
                    continue
    # setting cutoff:
    if isinstance(parameter, tuple):
        if parameter_str == 'cutoff':
            print('do you want to add one or multiple cutoffs?')
            print('it is recommend to choose cutoff values between 5 and 10 years')
            print('if not, just press ENTER, if so enter the coutoffs one by one:')
            while True:
                cutoff = input()
                if cutoff == '':
                    return parameter
                # sanity check, the intput value must be casteable into into
                # float:
                try:
                    float(cutoff)
                except ValueError:
                    print('the cutoff value must be number, starting over')
                    print('do you want to add one or multiple cutoffs?')
                    print('it is recommend to choose cutoff values between 5 and 10 years')
                    print('if not, just press ENTER, if so enter the coutoffs one by one:')
                    continue
                parameter = parameter + (float(cutoff), )
    # setting threshold:
    if isinstance(parameter, tuple):
        if parameter_str == 'threshold':
            print('do you want to add one or multiple thresholds?')
            print('it is recommend to choose threshold values which do not exceed a value of 50')
            print('if not, just press ENTER, if so enter the thresholds one by one:')
            while True:
                threshold = input()
                if threshold == '':
                    return parameter
                # sanity check, the intput value must be casteable into into
                # float:
                try:
                    float(threshold)
                except ValueError:
                    print('the threshold value must be number, starting over')
                    print('do you want to add one or multiple thresholds?')
                    print('it is recommend to choose threshold values between 5 and 10 years')
                    print('if not, just press ENTER, if so enter the coutoffs one by one:')
                    continue
                parameter = parameter + (float(threshold), )
