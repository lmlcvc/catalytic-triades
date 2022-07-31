import logging
import math
import os

import pandas as pd
from Bio.PDB.PDBParser import PDBParser

parser = PDBParser(PERMISSIVE=1, QUIET=True)


def create_folder(output_directory):
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
    return


def write_file(output_directory, header, protein, text, descriptor=''):
    file = open(output_directory + "/" + descriptor + str(protein) + '.csv', 'w')
    file.write(header)
    file.write(text)
    file.close()
    return


def find_angle(u, v, w):
    """
    :param u: surrounding point 1 coordinates
    :param v: angle point coordinates
    :param w: surrounding point 2 coordinates
    :return: angle in v [°]
    """

    a = u - v
    b = v - w
    c = w - u

    try:
        return round(math.degrees(
            math.acos((a * a + b * b - c * c) / (2 * a * b))), 2)
    except ValueError as e:
        logging.warning(e, u, v, w, a, b, c)
        return -1


def store_atoms(directory):
    """
    Store NUC, ACID and BASE atoms to respective lists, depending on name of file they're read from.
    """

    nuc_atoms, acid_atoms, base_atoms = {}, {}, {}

    for filename in os.listdir(directory):
        if filename.endswith('.pdb'):
            atoms_tmp = []  # stores file's atoms

            # write NUC atoms to NUC list
            if filename.startswith('nuc'):
                protein = filename[4:8]
                structure = parser.get_structure(protein, os.path.join(directory, filename))

                for model in structure:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                atoms_tmp.append((atom, residue))

                # store newly created list to NUC dict by protein
                nuc_atoms[protein] = atoms_tmp

            # write ACID atoms to ACID list
            if filename.startswith('acid'):
                protein = filename[5:9]
                structure = parser.get_structure(protein, os.path.join(directory, filename))

                for model in structure:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                atoms_tmp.append((atom, residue))

                # store newly created list to ACID dict by protein
                acid_atoms[protein] = atoms_tmp

            # write BASE atoms to BASE list
            if filename.startswith('base'):
                protein = filename[5:9]
                structure = parser.get_structure(protein, os.path.join(directory, filename))

                for model in structure:
                    for chain in model:
                        for residue in chain:
                            for atom in residue:
                                atoms_tmp.append((atom, residue))

                # store newly created list to BASE dict by protein
                base_atoms[protein] = atoms_tmp

    return nuc_atoms, acid_atoms, base_atoms


def store_triads_protein(directory, similar=True):
    protein_df = {}

    for filename in os.listdir(directory):
        if filename.startswith('similar') and similar is False:
            continue

        protein = filename.replace('.csv', '')
        triads_df = pd.read_csv(os.path.join(directory, filename), header=0)

        protein_df[protein] = triads_df

    return protein_df


def read_triads_df(directory, similar=True):
    triads_df_list = []
    header = ['NUC', 'ACID', 'BASE', 'Dist_Nuc_Acid', 'Dist_Acid_Base']  # TODO: change depending on columns nr

    for filename in os.listdir(directory):
        if filename.startswith('similar') and similar is False:
            continue

        df = pd.read_csv(os.path.join(directory, filename), header=None)
        df.columns = header
        triads_df_list.append(df)

    return pd.concat(triads_df_list, axis=0, ignore_index=True)


def get_triad_ranges_old(triads_all_df):
    """
    :param triads_all_df: Dataframe of all triads in all proteins
    :return: Dict of min and max values for old method range implementation
    """

    return {"Dist_Nuc_Acid_min": triads_all_df["Dist_Nuc_Acid"].min(),
            "Dist_Nuc_Acid_max": triads_all_df["Dist_Nuc_Acid"].max(),
            "Dist_Acid_Base_min": triads_all_df["Dist_Acid_Base"].min(),
            "Dist_Acid_Base_max": triads_all_df["Dist_Acid_Base"].max(),
            "Dist_Base_Nuc_min": triads_all_df["Dist_Base_Nuc"].min(),
            "Dist_Base_Nuc_max": triads_all_df["Dist_Base_Nuc"].max(),
            "Angle_Nuc_min": triads_all_df["Angle_Nuc"].min(), "Angle_Nuc_max": triads_all_df["Angle_Nuc"].max(),
            "Angle_Acid_min": triads_all_df["Angle_Acid"].min(), "Angle_Acid_max": triads_all_df["Angle_Acid"].max()}
