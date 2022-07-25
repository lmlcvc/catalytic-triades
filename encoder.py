import configparser
import math
import os

import util

config = configparser.ConfigParser()
config.read(os.path.join(os.pardir, 'config.ini'))
config = config['default']

encoded_directory = config['encoded_location']


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

            util.write_file(encoded_directory, '', protein, ''.join(population))
            # TODO: custom headers depending on parameters


""" Types of individual string encoding, depending on the problem question. """
ENCODING_TYPES = ['DISTANCES']


class Encoder:

    def __init__(self, **kwargs):
        pass

    def encode(self, individual):
        # TODO: return occurrences of this pattern throughout the population
        pass
