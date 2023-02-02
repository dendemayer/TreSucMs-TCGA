import os
import logging
import glob
import shutil
import re


def set_logger(OUTPUT_PATH, PROJECT_title, DRUGS_title):
    '''
    :param: OUTPUT_PATH: path for metilene pipeline outputs
    :type: OUTPUT_PATH: str
    :param: PROJECT_title: merged project title out of multiple projects
    :type: PROJECT_title: str
    :param: DRUGS_title: merged drug title out of multiple drugs
    :type: DRUGS_title: str

    every paths and options are set, configure here the logfiles, with which
    the snakemake config files are going to be created we create 2 loggers, in
    case just one project is applied the logs are written in
    PROJECT/DRUGS_title/test_log.log in case multi project is applied, the logs
    are written in PROJECT_title/DRUGS_title/test_log.log with that, it is
    clear which config file shall be createt out of the logfiles present in one
    outputpath (the drugs path must therefore be created from the first fct, to
    write the log file also, the dir of the logfile must be logged, s.t.
    snakemake knows where the input file for the final snakemake configuratioin
    file is located
    '''

    # in that dir)
    #   setting a level
    #   adding a filehandler
    #   adding a formatter to the filehandler
    # logging.basicConfig(
    # filename=os.path.join(
    # OUTPUT_PATH,
    # PROJECT_title,
    # DRUGS_title,
    # 'logfile.tsv'),
    # level=logging.INFO,
    # format='%(asctime)s\t%(levelname)s\t%(message)s')

    logger = logging.getLogger('logger')
    if logger.hasHandlers():
        logger.removeHandler(logger.handlers[0])
    formatter = logging.Formatter(
        '%(asctime)s\t%(levelname)s\t%(name)s\t%(message)s\t')
    handler = logging.FileHandler(
        os.path.join(
            OUTPUT_PATH, PROJECT_title, DRUGS_title, 'test_log.log'))
    handler.setFormatter(formatter)
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    return logger
    # if(logger_single.hasHandlers()):
    # logger_single.removeHandler(logger_single.handlers[0])
    # handler = logging.FileHandler(os.path.join(
    # OUTPUT_PATH, PROJECT, DRUGS_title,'test_log.log'))
    # handler.setFormatter(formatter)
    # logger_single.addHandler(handler)

    # # logger_multi = logging.getLogger('multi_logger')
    # # logger_multi.setLevel(logging.INFO)


def snake_meta(PROJECT, FILE_TYPE, OUTPUT_PATH, DRUGS_title, SCRIPT_PATH,
               PROJECT_title):
    '''
    copy the respective
    SCRIPT_PATH/Snakes/meta_infos/PROJECT_title/DRUGS_title/meta_info.dat
    into the actual OUTPUT_PATH/PROJECT/ path
    '''
    # its ok to glob the DRUGS_title_cutoff value, the meta_info.dat ist not
    # affected by cutoff
    temp_dir = PROJECT_title + '_' + DRUGS_title
    meta_source = glob.glob(os.path.join(SCRIPT_PATH, 'Snakes', 'meta_infos',
                                         temp_dir,
                                         PROJECT + '_meta_info.dat'))[0]
    dest = os.path.join(OUTPUT_PATH, PROJECT, 'meta_info.dat')
    shutil.copy2(meta_source, dest)


def apply_list(value):
    # print(drug)
    if re.search("palixtaxel", value):
        return "paclitaxel"
    if re.search("premetrexed", value):
        return "pemetrexed"
    if re.search("ironotecan", value):
        return "irinotecan"
    if (re.search("5-flurouracil", value) or
        re.search("5-fu", value) or re.search("5fu", value) or
        re.search("5 fu", value) or
        re.search("fluorouracil", value) or
        re.search("fluorouracil", value) or
            re.search("fluoruracil", value)):
        return "5-fluorouracil"
    if re.search("cisplatinum", value):
        return "cisplatin"
    if re.search("vinorelbin", value):
        return "vinorelbine"

    # ##harmonize synonyms:
    if re.search("cisplatin-xrt", value):
        return "cisplatin"
    if re.search("carboplatinum", value):
        return "carboplatin"
    if re.search("taxotere", value):
        return "docetaxel"
    if re.search("taxol", value) or re.search("abraxane", value):
        return "paclitaxel"
    if re.search("navelbine", value):
        return "vinorelbine"
    if re.search("cetuximab", value):
        return "erbitux"
    if re.search("erlotinib", value):
        return "tarceva"
    if re.search("panitumumab", value):
        return "vectibix"
    if re.search("paraplatin", value):
        return "carboplatin"
    if re.search("vepesid", value):
        return "etoposide"
    if re.search("gemzar", value):
        return "gemcitabine"
    if re.search("ironotecan", value):
        return "irinotecan"
    if re.search("alimta", value):
        return "pemetrexed"
    # for BRCA merge arimidex and anastrozole:
    if re.search("anastrozole", value):
        return "arimidex"
    else:
        return value
