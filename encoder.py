import configparser
import math
import os

import pandas as pd

import util

config = configparser.ConfigParser()
config.read(os.path.join(os.pardir, 'config.ini'))
config = config['default']

encoded_directory = config['encoded_location']
encoded_minmax_directory = config['encoded_minmax_location']


def encode_atoms(triangles, method="old"):
    """
    Encode atoms as arrays representing the population in GA application.
    Write encoded items to files.

    :param triangles: List of triangles from CSV
    :param method: Name of method in which atoms are used, defines their encoding type
    """

    # TODO: urednost koda, OO?

    if method == "old":
        for triangle in triangles:
            genes = []

            # Nuc atom gene
            if triangle["Nuc_name"] == "OG":
                genes.append(0)
            elif triangle["Nuc_name"] == "SG":
                genes.append(1)

            # Acid atom gene
            if triangle["Acid_name"] == "OD1":
                genes.append(0)
            elif triangle["Acid_name"] == "OD2":
                genes.append(1)

            # Base atom gene
            if triangle["Base_aa"] == "HIS":
                genes.append(0)
            elif triangle["Base_aa"] == "ASP":
                genes.append(1)
            elif triangle["Base_aa"] == "GLU":
                genes.append(2)

            # - - - - - - - - - -

            # Nuc-Acid distance
            genes.append((1 - (8.57 - triangle["Dist_Nuc_Acid"]) / 2.84) * 10)

            # Acid-Base distance
            genes.append((1 - (4.62 - triangle["Dist_Acid_Base"]) / 1.96) * 10)

            # - - - - - - - - - -

            # Nuc angle
            genes.append((1 - (42.22 - triangle["Angle_Nuc"]) / 30.76) * 10)

            # Acid angle
            genes.append((1 - (65.85 - triangle["Angle_Acid"]) / 55.68) * 10)

            individual = ''.join(genes)


def old_encode(triads, all_distances=False, angles=False, distance_categories=20, angle_categories=20):
    """
    Encode triads to csv files of integer rows, taking min and max accepted values as categor reference points.

    :param triads: Dict of triads by protein
    :param all_distances: Boolean - Whether all 3 distances should be included in encoding
    :param angles: Whether angles should be included in encoding
    :param distance_categories: Number of distance categories in min-max range
    :param angle_categories: Number of angle categories in min-max range
    :return:
    """

    # make df of all triads
    triads_all_df = pd.concat(triads.values(), ignore_index=True)

    # get min & max value for all
    ranges = util.get_triad_ranges_old(triads_all_df)

    # TODO: prepisati ovo dolje da pišu te min i max vrijednosti

    util.create_folder(encoded_directory)

    for protein in triads.keys():
        population = []

        for index, triad in triads[protein].iterrows():
            genes = []

            # Nuc atom gene
            if triad["Nuc_name"] == "OG":
                genes.append(0)
            elif triad["Nuc_name"] == "SG":
                genes.append(1)

            # Acid atom gene
            if triad["Acid_name"] == "OD1":
                genes.append(0)
            elif triad["Acid_name"] == "OD2":
                genes.append(1)

            # Base atom gene
            if triad["Base_aa"] == "HIS":
                genes.append(0)
            elif triad["Base_aa"] == "ASP":
                genes.append(1)
            elif triad["Base_aa"] == "GLU":
                genes.append(2)

            # - - - - - - - - - -

            # Nuc-Acid distance
            nuc_acid_range = ranges["Dist_Nuc_Acid_max"] - ranges["Dist_Nuc_Acid_min"]
            genes.append(
                math.floor((1 - (ranges["Dist_Nuc_Acid_max"] - triad["Dist_Nuc_Acid"])
                            / nuc_acid_range) * distance_categories))

            # Acid-Base distance
            acid_base_range = ranges["Dist_Acid_Base_max"] - ranges["Dist_Acid_Base_min"]
            genes.append(math.floor((1 - (ranges["Dist_Acid_Base_max"] - triad["Dist_Acid_Base"])
                                     / acid_base_range) * distance_categories))

            # Base-Nuc distance (optional)
            if all_distances:
                base_nuc_range = ranges["Dist_Base_Nuc_max"] - ranges["Dist_Base_Nuc_min"]
                genes.append(math.floor((1 - (ranges["Dist_Base_Nuc_max"] - triad["Dist_Base_Nuc"]) /
                                         base_nuc_range) * distance_categories))

            # - - - - - - - - - -

            if angles:
                # Nuc angle
                nuc_angle_range = ranges["Angle_Nuc_max"] - ranges["Angle_Nuc_min"]
                genes.append((1 - (ranges["Angle_Nuc_max"] - triad["Angle_Nuc"]) /
                              nuc_angle_range) * angle_categories)

                # Acid angle
                acid_angle_range = ranges["Angle_Acid_max"] - ranges["Angle_Acid_min"]
                genes.append((1 - (ranges["Angle_Acid_max"] - triad["Angle_Acid"]) /
                              acid_angle_range) * angle_categories)

            genes.append('\n')
            genes = [str(i) for i in genes]
            individual = "".join([",".join(genes[:-1]), genes[-1]])
            population.append(individual)

            util.write_file(encoded_directory, '', protein, ''.join(population))


def old_encode_minmax(triads, all_distances=False, angles=False, distance_categories=20, angle_categories=20):
    """
    Encode triads to csv files of integer rows, taking min and max accepted values as categor reference points.

    :param triads: Dict of triads by protein
    :param all_distances: Boolean - Whether all 3 distances should be included in encoding
    :param angles: Whether angles should be included in encoding
    :param distance_categories: Number of distance categories in min-max range
    :param angle_categories: Number of angle categories in min-max range
    :return:
    """

    util.create_folder(encoded_minmax_directory)

    for protein in triads.keys():
        population = []

        for index, triad in triads[protein].iterrows():
            genes = []

            # Nuc atom gene
            if triad["Nuc_name"] == "OG":
                genes.append(0)
            elif triad["Nuc_name"] == "SG":
                genes.append(1)

            # Acid atom gene
            if triad["Acid_name"] == "OD1":
                genes.append(0)
            elif triad["Acid_name"] == "OD2":
                genes.append(1)

            # Base atom gene
            if triad["Base_aa"] == "HIS":
                genes.append(0)
            elif triad["Base_aa"] == "ASP":
                genes.append(1)
            elif triad["Base_aa"] == "GLU":
                genes.append(2)

            # - - - - - - - - - -

            # Nuc-Acid distance
            genes.append(math.floor((1 - (8.57 - triad["Dist_Nuc_Acid"]) / 2.84) * distance_categories))

            # Acid-Base distance
            genes.append(math.floor((1 - (4.62 - triad["Dist_Acid_Base"]) / 1.96) * distance_categories))

            # Base-Nuc distance (optional)
            if all_distances:
                genes.append(math.floor((1 - (4.86 - triad["Dist_Base_Nuc"]) / 2.12) * distance_categories))

            # - - - - - - - - - -

            if angles:
                # Nuc angle
                genes.append((1 - (42.22 - triad["Angle_Nuc"]) / 30.76) * angle_categories)

                # Acid angle
                genes.append((1 - (65.85 - triad["Angle_Acid"]) / 55.68) * angle_categories)

            genes.append('\n')
            genes = [str(i) for i in genes]
            individual = "".join([",".join(genes[:-1]), genes[-1]])
            population.append(individual)

            util.write_file(encoded_minmax_directory, '', protein, ''.join(population))
            # TODO: custom headers depending on parameters


""" Types of individual string encoding, depending on the problem question. """
ENCODING_TYPES = ['DISTANCES']


class Encoder:

    def __init__(self, **kwargs):
        pass

    def encode(self, individual):
        # TODO: return occurrences of this pattern throughout the population
        pass
