"""Read transformed files, make an Atom object list for each protein, find triades in each protein."""
import math
import os
import configparser
import logging
from Bio.PDB.PDBParser import PDBParser

parser = PDBParser(PERMISSIVE=1, QUIET=True)

config = configparser.ConfigParser()
config.read('config.ini')
config = config['default']

directory = config['transformed_location']
output_directory = config['output_location']

old_directory = r'C:\Users\Marina\Desktop\Lana\triades\uoop_projekt\transformed_old'

"""
Udaljenosti:
STDEVx3: 1.42 Nuc-Acid, 1.06 Nuc-Base, 0.98 Base-Acid
AVG: 7.15- Nuc-Acid, 4.86 - Nuc-Base, 3.64 - Base-Acid

Kutevi:
STDEVx3: 15.38 Nuc,  41.27Base, 27.84 Acid
AVG: 26.84 - Nuc, 115.14 - Base, 38.01 - Acid"""

NUC_ACID_MIN = 7.15 - 1.42
NUC_ACID_MAX = 7.15 + 1.42
NUC_BASE_MIN = 4.86 - 1.06
NUC_BASE_MAX = 4.86 + 1.06
ACID_BASE_MIN = 3.64 - 0.98
ACID_BASE_MAX = 3.64 + 0.98

ANGLE_NUC_MIN = 26.84 - 15.38
ANGLE_NUC_MAX = 26.84 + 15.38
ANGLE_BASE_MIN = 115.15 - 41.27
ANGLE_BASE_MAX = 115.15 + 41.27
ANGLE_ACID_MIN = 38.01 - 27.84
ANGLE_ACID_MAX = 38.01 + 27.84

HEADER = 'Nuc_name,Nuc_posX,Nuc_posY,Nuc_posZ,Nuc_aa,Nuc_aaID,Acid_name,Acid_posX,Acid_posY,Acid_posZ,' \
         'Acid_aa,Acid_aaID,Base_name,Base_posX,Base_posY,Base_posZ,Base_aa,Base_aaID,Dist_Nuc_Acid,' \
         'Dist_Acid_Base,Dist_Base_Nuc,Angle_Nuc,Angle_Acid,Angle_Base\n'


def store_atoms():
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


def find_angle(u, v, w):
    a = u - v
    b = v - w
    c = w - u

    try:
        return round(math.degrees(
            math.acos((a * a + b * b - c * c) / (2 * a * b))), 2)
    except ValueError as e:
        logging.warning(e, u, v, w, a, b, c)
        return -1


def parse_triangle_descriptors(nuc, base, acid, nuc_angle, acid_angle, base_angle):
    nuc_coord = nuc[0].get_coord()
    base_coord = base[0].get_coord()
    acid_coord = acid[0].get_coord()

    return (str(nuc[0].get_name()) + ',' + str(nuc_coord[0]) + ',' + str(nuc_coord[1]) + ','
            + str(nuc_coord[2]) + ',' + str(nuc[1].get_resname()) + ',' + str(nuc[1].get_id()[1]) + ','
            + str(acid[0].get_name()) + ',' + str(acid_coord[0]) + ',' + str(acid_coord[1]) + ','
            + str(acid_coord[2]) + ',' + str(acid[1].get_resname()) + ',' + str(acid[1].get_id()[1]) + ','
            + str(base[0].get_name()) + ',' + str(base_coord[0]) + ',' + str(base_coord[1]) + ','
            + str(base_coord[2]) + ',' + str(base[1].get_resname()) + ',' + str(base[1].get_id()[1]) + ','
            + str(nuc[0] - acid[0]) + ',' + str(acid[0] - base[0]) + ',' + str(base[0] - nuc[0]) + ','
            + str(nuc_angle) + ',' + str(acid_angle) + ',' + str(base_angle) + '\n')


def create_folder():
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
    return


def write_file(protein, text, descriptor):
    file = open(output_directory + "/" + descriptor + str(protein) + '.csv', 'w')
    file.write(HEADER)
    file.write(text)
    file.close()
    return


def store_triads(nuc_atoms, acid_atoms, base_atoms):
    for protein in nuc_atoms.keys():

        nucs = nuc_atoms[protein]
        acids = acid_atoms[protein]
        bases = base_atoms[protein]

        text_list = []
        similar_text_list = []

        for nuc in nucs:
            for base in bases:
                for acid in acids:

                    # if triad candidates found, test angles
                    nuc_angle = find_angle(acid[0], nuc[0], base[0])
                    acid_angle = find_angle(nuc[0], acid[0], base[0])
                    base_angle = find_angle(nuc[0], base[0], acid[0])

                    if ((ANGLE_NUC_MIN <= nuc_angle <= ANGLE_NUC_MAX)
                            and (ANGLE_ACID_MIN <= acid_angle <= ANGLE_ACID_MAX)
                            and (ANGLE_BASE_MIN <= base_angle <= ANGLE_BASE_MAX)):

                        # if angle fits criteria, test distances
                        if ((NUC_ACID_MIN <= nuc[0] - acid[0] <= NUC_ACID_MAX)
                                and (NUC_BASE_MIN <= nuc[0] - base[0] <= NUC_BASE_MAX)
                                and (ACID_BASE_MIN <= acid[0] - base[0] <= ACID_BASE_MAX)):
                            # if distances fit, mark as a triad found
                            text_list.append(
                                parse_triangle_descriptors(nuc, base, acid, nuc_angle, acid_angle,
                                                           base_angle))

                        # if distances don't fit but angles do, mark as a similar triangle
                        else:
                            similar_text_list.append(parse_triangle_descriptors(nuc, base, acid,
                                                                                nuc_angle,
                                                                                acid_angle,
                                                                                base_angle))

        write_file(protein, ''.join(text_list), '')
        write_file(protein, ''.join(similar_text_list), 'similar_')

        """atoms = protein_atoms[protein]
        text_list = []
        similar_text_list = []
        for nuc in atoms:
            if (nuc[0].get_name() == 'OG' and nuc[1].get_resname() == 'SER') or (
                    nuc[0].get_name() == 'SG' and nuc[1].get_resname() == 'CYS'):
                for base in atoms:
                    if base[0].get_name() == 'CG' and (
                            base[1].get_resname() == 'HIS' or base[1].get_resname() == 'ASP'
                            or base[1].get_resname() == 'GLU'):
                        for acid in atoms:
                            if acid[0].get_name() == 'OD1' or acid[0].get_name() == 'OD2':

                                # if triad candidates found, test angles
                                nuc_angle = find_angle(acid[0], nuc[0], base[0])
                                acid_angle = find_angle(nuc[0], acid[0], base[0])
                                base_angle = find_angle(nuc[0], base[0], acid[0])

                                if ((ANGLE_NUC_MIN <= nuc_angle <= ANGLE_NUC_MAX)
                                        and (ANGLE_ACID_MIN <= acid_angle <= ANGLE_ACID_MAX)
                                        and (ANGLE_BASE_MIN <= base_angle <= ANGLE_BASE_MAX)):

                                    # if angle fits criteria, test distances
                                    if ((NUC_ACID_MIN <= nuc[0] - acid[0] <= NUC_ACID_MAX)
                                            and (NUC_BASE_MIN <= nuc[0] - base[0] <= NUC_BASE_MAX)
                                            and (ACID_BASE_MIN <= acid[0] - base[0] <= ACID_BASE_MAX)):
                                        # if distances fit, mark as a triad found
                                        text_list.append(
                                            parse_triangle_descriptors(nuc, base, acid, nuc_angle, acid_angle,
                                                                       base_angle))

                                    # if distances don't fit but angles do, mark as a similar triangle
                                    else:
                                        similar_text_list.append(parse_triangle_descriptors(nuc, base, acid, nuc_angle,
                                                                                            acid_angle,
                                                                                            base_angle))
        write_file(protein, ''.join(text_list), '')
        write_file(protein, ''.join(similar_text_list), 'similar_')"""


def find_triades():
    nuc_atoms, acid_atoms, base_atoms = store_atoms()

    create_folder()
    store_triads(nuc_atoms, acid_atoms, base_atoms)
