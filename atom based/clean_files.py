import os
import configparser
from Bio.PDB.PDBParser import PDBParser

parser = PDBParser(PERMISSIVE=1, QUIET=True)

"""
NUC – BASE: any O/S to N
BASE – ACID: any O to N
ACID – NUC: previously measured O/S to previously measured O
"""


def clean_files():
    """
    Removes each atom that doesn't fit triad point criteria.
    Stores NUC, ACID, BASE points to respective files.

    :return: None
    """

    config = configparser.ConfigParser()
    config.read(os.path.join(os.pardir, 'config.ini'))
    config = config['atom']

    directory = config['location']
    output_directory = config['transformed_location']

    if not os.path.isdir(output_directory):
        print('here')
        os.makedirs(output_directory)
        print(os.path.isdir(output_directory))

    for filename in os.listdir(directory):
        if filename.endswith(".pdb"):
            print(os.path.join(directory, filename))
            file_absolute_path = os.path.join(directory, filename)

            nuc_path = os.path.join(output_directory, 'nuc_' + filename)
            base_path = os.path.join(output_directory, 'base_' + filename)
            acid_path = os.path.join(output_directory, 'acid_' + filename)

            with open(file_absolute_path, 'r') as file:
                atoms = file.readlines()

                with open(nuc_path, 'w+') as nuc, \
                        open(acid_path, 'w+') as acid, \
                        open(base_path, 'w+') as base:
                    nuc_list, acid_list, base_list = [], [], []

                    for atom in atoms:
                        if len(atom) > 77:  # only look at rows with proper atom formatting

                            # add NUCs (any O/S) to NUC list
                            if atom[13] == "O" \
                                    or atom[13] == "S":
                                nuc_list.append(atom)

                            # add ACIDs (any O) to ACID list
                            if atom[13] == "O":
                                acid_list.append(atom)

                            # add BASEs (any N) to BASE list
                            if atom[13] == "N":
                                base_list.append(atom)

                    # write atom lists to NUC/ACID/BASE files
                    nuc.write(''.join(nuc_list))
                    acid.write(''.join(acid_list))
                    base.write(''.join(base_list))
