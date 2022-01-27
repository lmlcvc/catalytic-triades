import os
import configparser

# Ukloniti svaki atom koji nije:
# * OG u SER
# * SG u CYS
# * CG u HIS/ASP/GLU
# * OD1 ili OD2

def clean_files():
    config = configparser.ConfigParser()
    config.read('config.ini')
    config = config['default']

    dirpath = config['location']
    destpath = config['transformed_location']
    keypath = config['location'] + config['keywords']

    if not os.path.isdir(destpath):
        os.makedirs(destpath)

    with open(keypath, 'r') as keys:
        key_lines = [s.strip() for s in keys.readlines()]

    for filename in os.listdir(dirpath):
        if filename.endswith(".pdb"):
            file_fullname = os.path.join(dirpath, filename)
            destname = os.path.join(destpath, filename)

            with open(file_fullname, 'r') as file:
                atoms = file.readlines()
                with open(destname, 'w') as dest:
                    for key in key_lines:
                        text = ""
                        for atom in atoms:
                            if key in atom[11:20].strip().replace(" ", "") \
                                    and 'ASN' not in atom[11:20].strip().replace(" ", ""):
                                text += atom
                        if len(text) > 0:
                            dest.write(text)
