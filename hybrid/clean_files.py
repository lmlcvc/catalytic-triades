import os
import configparser
from Bio.PDB.PDBParser import PDBParser

parser = PDBParser(PERMISSIVE=1, QUIET=True)

"""
NUC – BASE: OG, SG
BASE – ACID: OD1, OD2, O, OE1, OE2
ACID – NUC: ND1, ND2
"""

NUC_CLASSIFICATION = ["OG", "SG"]
ACID_CLASSIFICATION = ["OD1", "OD2", "O", "OE1", "OE2"]
BASE_CLASSIFICATION = ["ND1", "NE2"]


def clean_files():
    """
    Removes each atom that doesn't fit triad point criteria.
    Stores NUC, ACID, BASE points to respective files.

    :return: None
    """

    config = configparser.ConfigParser()
    config.read(os.path.join(os.pardir, 'config.ini'))
    config = config['hybrid']

    directory = config['location']
    output_directory = config['transformed_location']

    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    for filename in os.listdir(directory):
        if filename.endswith(".pdb"):
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
                        atom_classification = atom[13:16].strip().replace(" ", "")

                        # add NUCs (OG/SG) to NUC list
                        if atom_classification in NUC_CLASSIFICATION:
                            nuc_list.append(atom)

                        # add ACIDs (OD1, OD2, O, OE1, OE2) to ACID list
                        if atom_classification in ACID_CLASSIFICATION:
                            acid_list.append(atom)

                        # add BASEs (ND1, NE2) to BASE list
                        if atom_classification in BASE_CLASSIFICATION:
                            base_list.append(atom)

                    # write atom lists to NUC/ACID/BASE files
                    nuc.write(''.join(nuc_list))
                    acid.write(''.join(acid_list))
                    base.write(''.join(base_list))
