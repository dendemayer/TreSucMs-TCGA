import os
import difflib

# ###  scenarios:
# pipeline runs first time (config run):
# - in par dir of tcga_data the all_data_files.txt is saved
# pipeline runs for multiple time:
# - it can be checked if the newly downloaded files are the same as the
# previous, with help of already saved all_data_files.txt
# pipilen runs with snake -> the all_data_files.txt is not present, but it
# shall be fixed which files shall be downloaded (not just those which are
# specified in config, but also those which are additionally downloaded ->
# - from the Snakemakefile, create the all_data_files.txt within the pardir
# than it is automatically checked that all files match as if it would be a
# second run

OUTPUT_PATH = '/scr/dings/PEVO/NEW_downloads_2/metilene_config_run_CESC'
PROJECT_title = 'TCGA-CESC'
data_dir = PROJECT_title + '_data_files'
data_path = os.path.join(OUTPUT_PATH, PROJECT_title, data_dir)
file_list = os.listdir(data_path)
file_list.sort()
all_files_path = os.path.join(data_path, os.path.pardir, 'all_data_files.txt')

# if all_data_files.txt exists allready,
if os.path.isfile(all_files_path):
    print('\nyou already downloaded the data files for this project,' +
          ' checking for consistency:\n')
    # compare the created files list with the list already saved:
    with open(all_files_path, 'r') as f_1:
        file_list_old = [line.strip() for line in f_1.readlines()]
        file_list_old.sort()
    # now newly loaded data is hold in list file_list and previous loaded data
    # is hold in file_list_old, both should be identical:
    diff_list = difflib.unified_diff(file_list, file_list_old,
                                     fromfile='new_data',
                                     tofile='old_data', lineterm='', n=0)
    # diff list is a generator object, if diffs are found, the unpacked object
    # is longer than 0
    if len([*diff_list]) == 0:
        print('no differences found in previously loaded data, everythings' +
              ' fine')
    else:
        print('\ninconsistencies detected:\n')
        for line in difflib.unified_diff(
                                        file_list, file_list_old,
                                        fromfile='new_data',
                                        tofile='old_data', lineterm='', n=0):
            print(line)
        print('\nif you started from a snake_config.yaml, manually check ' +
              'what went wrong with the downloaded data and start again...' +
              ' exiting now.\n')
        os.sys.exit(0)

# if all_data_files.txt is not present allready, write it:
else:
    with open(all_files_path, 'w') as f:
        f.writelines('\n'.join(file_list))
