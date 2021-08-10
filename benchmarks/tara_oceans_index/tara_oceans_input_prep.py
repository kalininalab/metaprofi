"""
For preparing Tara Ocenas dataset input_data.txt file for MetaProFi
"""
import os

base = "/data/tara_oceans"
metaprofi_input_file = "/data/tara_oceans/input_data.txt"
files_list = []

# Gather all the sub-dir list
folder_list = [dir for dir in os.listdir(base) if os.path.isdir(f"{base}/{dir}")]
for dir_ in folder_list:
    sub_l = [dir_]
    for file in os.listdir(f"{base}/{dir_}"):
        if file.endswith(".gz"):
            sub_l.append(f"{base}/{dir_}/{file}")
    files_list.append(sub_l)

# MetaProFi input file construction
with open(metaprofi_input_file, "w") as outf:
    for f in files_list:
        try:
            outf.write(f"{f[0]}: {f[1]}; {f[2]}\n")
        except:
            outf.write(f"{f[0]}: {f[1]}\n")
