import os, shutil, re

for folder_name, sub_folders, file_names in os.walk('./plots-5.0_30'):
    reg = re.compile(r'\(\d+\)')
    for file_name in file_names:
        file_name_new = reg.search(file_name).group()
        file_name_new = file_name_new.lstrip('(').rstrip(')')
        shutil.move('./plots-5.0_30/' + file_name, './plots-5.0_30/' + file_name_new.zfill(5) + '.png')