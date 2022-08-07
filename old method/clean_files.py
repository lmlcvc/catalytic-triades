import os
import configparser

config = configparser.ConfigParser()
config.read(os.path.join(os.pardir, 'config.ini'))
config = config['default']

dirpath = config['location']
destpath = config['transformed_location']
keypath = config['location'] + config['keywords']


def clean_files():
    """
    Removes each atom that doesn't fit triad point criteria.
    Stores NUC, ACID, BASE points to respective files.

    :return: None
    """

    if not os.path.isdir(destpath):
        os.makedirs(destpath)

    for filename in os.listdir(dirpath):
        if filename.endswith(".pdb"):
            file_fullname = os.path.join(dirpath, filename)
            nuc_path = os.path.join(destpath, 'nuc_' + filename)
            base_path = os.path.join(destpath, 'base_' + filename)
            acid_path = os.path.join(destpath, 'acid_' + filename)

            with open(file_fullname, 'r') as file:
                atoms = file.readlines()

                with open(nuc_path, 'w') as nuc, \
                        open(acid_path, 'w') as acid, \
                        open(base_path, 'w') as base:
                    nuc_list, acid_list, base_list = [], [], []

                    for atom in atoms:
                        atom_classificaiton = atom[11:20].strip().replace(" ", "")

                        # add NUCs (OG SER / SG CYS) to NUC list
                        if "OGSER" in atom_classificaiton \
                                or "SGCYS" in atom_classificaiton:
                            nuc_list.append(atom)

                        # add ACIDs (OD1, OD2) to ACID list
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
                    base.write(''.join(base_list))
