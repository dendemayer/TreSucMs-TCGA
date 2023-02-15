import pandas as pd
import os
import requests
import json
import time
import subprocess
"""
searching through the manifest file to get the UUID of the requested file:
"""

manifest_file = snakemake.input[0]
aux_out = snakemake.output[0]
aux_file = os.path.split(aux_out)[1]
DF_mani = pd.read_table(manifest_file)
# UUID_list = [(DF_mani.set_index('filename').loc[aux_file, 'id'])]
UUID = DF_mani.set_index('filename').loc[aux_file, 'id']
if isinstance(UUID, pd.Series):
    UUID_list = DF_mani.set_index('filename').loc[aux_file, 'id'].to_list()
else:
    UUID_list = [(DF_mani.set_index('filename').loc[aux_file, 'id'])]

suffix = 1
for UUID in UUID_list:
    while True:
        # it happens that we have 2 different UUID for the exact same filename
        # (?!?) -> save both, the first name is then saved as expected and
        data_endpt = "https://api.gdc.cancer.gov/data/"
        params = {"ids": UUID}
        try:
            response = requests.post(
                data_endpt, data=json.dumps(
                    params), headers={"Content-Type": "application/json"})
        #  *** TypeError: Object of type Series is not JSON serializable
        except Exception as e:
            print(f'Exception {e} with file {aux_out}, retrying')
            continue
            # breakpoint()
            # open(aux_out, 'a').close()
            # os._exit(0)

        # The file name can be found in the header within the
        # Content-Disposition key.
        try:
            response_head_cd = response.headers["Content-Disposition"]
            break
        except KeyError:
            iterations = 0
            print(f'problems while loading {aux_file}, retrying')
            continue
            # while iterations < 5:
            #     print(f'problems while loading {aux_file}')
            #     print('TCGA GDC dataportal not reachable at the moment\
            #         check the website: https://portal.gdc.cancer.gov/\
            #         for potentially network maintenance news and try again later')
            #     print(f'(trying {5-iterations} times again...)')
            #     print('waiting 5 Minutes before new attempt')
            #     continue
                # time.sleep(300)
                # try:
                #     try:
                #         response = requests.post(
                #             data_endpt, data=json.dumps(
                #                 params), headers={"Content-Type": "application/json"})
                #     #  *** TypeError: Object of type Series is not JSON serializable
                #     except TypeError:
                #         print(f'cannot find aux file {aux_out}, skipping:')
                #         breakpoint()
                #         open(aux_out, 'a').close()
                #         os._exit(0)
                #     response_head_cd = response.headers["Content-Disposition"]
                # except KeyError:
                #     iterations += 1
            os._exit(0)

    # if the connection breaks while writing the file, just try again until its
    # working again
    while True:
        with open(aux_out, "wb") as output_file:
            try:
                output_file.write(response.content)
                break
            except Exception as e:
                print('Exception catched: {e},\
                      trying again to write file {aux_out}')
                print('wait 5 min:')
                time.sleep(300)
                continue
    # now the md5sum can be checked:
    md5sum = DF_mani.set_index('id').loc[UUID, 'md5']
    md5sum_file = subprocess.check_output( ['md5sum', aux_out]).decode('utf-8').split(' ')[0]
    if md5sum_file == md5sum:
        print(f'md5sum of file {aux_out} with:\n{md5sum}\nis verified')
    else:
        print(f'cannot confirm md5sum of file {aux_out}!!!')
        breakpoint()

    aux_out = aux_out + '_' + str(suffix)
    suffix += 1


# unpack all auxiliary files and put them into PROJECT/aux_files/ dir:
# file_name holds the just downloaded tar.gz, composed out of date and
# other numerical flags, ex:gdc_download_20220524_114836.301313.tar.gz
# this filename is generic and changes with every download, do not log it!
# instead log the content from this archive.
# after extracting, the generic .gz file can be deleted
########
# os.makedirs(os.path.join(OUTPUT_PATH, PROJECT, 'aux_files'), exist_ok=True)
# tar = tarfile.open(file_name, 'r:gz')
# tar.extractall(os.path.join(OUTPUT_PATH, PROJECT, 'aux_files'))
# tar.close()
# list_aux_tables = glob.glob(os.path.join(OUTPUT_PATH, PROJECT,
#                                            # 'aux_files', '*', '*'))
# now mv every table one dir above and delete the remaining empty dir:
# # for aux_table in list_aux_tables:
#     # src = aux_table
#     # dest = os.path.join(
#         # os.path.split(
#             # aux_table)[0], os.path.pardir, os.path.split(aux_table)[1])
#     # os.replace(src, dest)
# compl_ls = os.listdir(os.path.join(OUTPUT_PATH, PROJECT, 'aux_files'))
# # for element in compl_ls:
#     # dir_full_path = os.path.join(OUTPUT_PATH, PROJECT, 'aux_files',
#                                     # element)
#     # if os.path.isdir(dir_full_path):
#         # os.removedirs(dir_full_path)
# # log every file within aux_files:
# logger = set_logger.set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
# # for aux_file in os.listdir(os.path.join(
#         # OUTPUT_PATH, PROJECT, 'aux_files')):
#     # logger.info('download_clinical_tables_2\t{}'.format(
#         # os.path.join(PROJECT, 'aux_files', aux_file)))
# # os.remove(file_name)