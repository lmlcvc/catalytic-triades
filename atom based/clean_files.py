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
    config = config['default']

    directory = config['location']
    output_directory = config['transformed_location_atom']
    print(directory)
    print(output_directory)

    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    nuc_classifications = ['O', 'OD1', 'OD2']

    for filename in os.listdir(directory):
        if filename.endswith(".pdb"):
            print(os.path.join(directory, filename))
            structure = parser.get_structure(filename, os.path.join(directory, filename))

            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            # atoms_tmp.append((atom, residue))
                            print(atom.get_id())

            file_fullname = os.path.join(directory, filename)
            nuc_path = os.path.join(output_directory, 'nuc_' + filename)
            base_path = os.path.join(output_directory, 'base_' + filename)
            acid_path = os.path.join(output_directory, 'acid_' + filename)

            """with open(file_fullname, 'r') as file:
                atoms = file.readlines()

                with open(nuc_path, 'w') as nuc, \
                        open(acid_path, 'w') as acid, \
                        open(base_path, 'w') as base:
                    nuc_list, acid_list, base_list = [], [], []

                    for atom in atoms:
                        atom_classificaiton = atom[11:16].strip().replace(" ", "")

                        # add NUCs (OG SER / SG CYS) to NUC list
                        if "OGSER" in atom_classificaiton \
                                or "SGCYS" in atom_classificaiton:
                            nuc_list.append(atom)

                        # add ACIDs (CG HIS / CG ASP / CG GLU) to ACID list
                        if ("OD1" in atom_classificaiton
                            or "OD2" in atom_classificaiton) \
                                and "ASN" not in atom_classificaiton:
                            acid_list.append(atom)

                        # add BASEs (CG HIS / CG ASP / CG GLU) to BASE list
                        if "CGHIS" in atom_classificaiton \
                                or "CGASP" in atom_classificaiton \
                                or "CGGLU" in atom_classificaiton:
                            base_list.append(atom)

                    # write atom lists to NUC/ACID/BASE files
                    nuc.write(''.join(nuc_list))
                    acid.write(''.join(acid_list))
                    base.write(''.join(base_list))"""
