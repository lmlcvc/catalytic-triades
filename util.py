import logging
import math
import os
from collections import defaultdict

import numpy as np
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


def read_triads_dict(directory, similar=True):
    triads_df_dict = {}
    triads_df_dict = defaultdict(lambda: pd.DataFrame(), triads_df_dict)

    header = ['NUC', 'ACID', 'BASE', 'Dist_Nuc_Acid', 'Dist_Acid_Base']  # TODO: change depending on columns nr

    for filename in os.listdir(directory):
        if filename.startswith('similar') and similar is False:
            continue

        protein = filename.replace("similar_", "").replace(".csv", "")

        df = pd.read_csv(os.path.join(directory, filename), header=None)
        df.columns = header
        triads_df_dict[protein] = pd.concat([triads_df_dict[protein], df], axis=0, ignore_index=True)

    return triads_df_dict


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


def get_iteration_info(population, algo_type, iteration, header, output_directory):
    pop = []
    for i in population:
        pop.append([g for g in i.x])
        pop[-1].append(i.f)

    population_df = pd.DataFrame(pop, columns=header)

    best = population_df[population_df.fitness == population_df.fitness.max()].iloc[-1]  # best result (just the first)
    best_df = pd.DataFrame(best, index=header).T

    print(population_df)
    print(best_df)
    population_df.to_csv(os.path.join(output_directory, algo_type + iteration), header=True, index=False)

    return best_df


def levenshtein(a, b):
    """
    Calculates the Levenshtein distance between a and b.
    a: list
    b: list
    :return: Levenshtein distance of a and b
    """
    n, m = len(a), len(b)
    if n > m:
        # Make sure n <= m, to use O(min(n,m)) space
        a, b = b, a
        n, m = m, n

    current = range(n + 1)
    for i in range(1, m + 1):
        previous, current = current, [i] + [0] * n
        for j in range(1, n + 1):
            add, delete = previous[j] + 1, current[j - 1] + 1
            change = previous[j - 1]
            if a[j - 1] != b[i - 1]:
                change = change + 1
            current[j] = min(add, delete, change)

    return current[n]


def list_similarities(x):
    """
    Calculate Levenhstein distance of all 2-element combinations from list.
    List similarity is equal to their average
    x: list
    """

    similarities = []

    for i in range(0, len(x)):
        for j in range((i + 1), len(x)):
            similarities.append(levenshtein(x[i], x[j]))

    return np.round(sum(similarities) / len(similarities), 2)


def store_fitness_trends(fitness, filepath):
    fitness_df = pd.DataFrame(fitness, columns=["fitness"])
    fitness_df.to_csv(filepath, header=True, index=False)
