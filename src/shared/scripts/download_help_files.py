import requests
import os
import subprocess
# get the belonging md5sum out of the config file
out_file = snakemake.wildcards[1]
out_dest = snakemake.output[0]
server_file_name = snakemake.config[out_file][0]
md5sum = snakemake.config[out_file][1]
url = snakemake.config[out_file][2]
in_dest = os.path.join(url, server_file_name)
# manifest_out =
# dest_manifest = os.path.join(OUTPUT_PATH, api_manifest)
data = requests.get(in_dest)
# logger = set_logger.set_logger(OUTPUT_PATH, PROJECT, DRUGS_title)
with open(out_dest, 'wb') as f:
    try:
        f.write(data.content)
    except BaseException:
        print(f'problems while downloading the {out_file}')
        print(f'is {url} online?')
        print('exiting the program')
        os._exit(0)

md5sum_file = subprocess.check_output(
    ['md5sum', out_dest]).decode('utf-8').split(' ')[0]
if md5sum_file == md5sum:
    print(f'md5sum of\n{md5sum}\nis verified')
