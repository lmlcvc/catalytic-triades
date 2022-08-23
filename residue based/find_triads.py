"""Read transformed files, make an Atom object list for each protein, find triads in each protein."""
import os
import configparser
import util
from Bio.PDB.PDBParser import PDBParser

parser = PDBParser(PERMISSIVE=1, QUIET=True)

config = configparser.ConfigParser()
config.read(os.path.join(os.pardir, 'config.ini'))
config = config['residue']

directory = config['transformed_location']
output_directory = config['output_location']

ACID_BASE_MAX = 8
BASE_NUC_MAX = 8
NUC_ACID_MAX = 11

ANGLE_ACID_MAX = 55
ANGLE_BASE_MAX = 130
ANGLE_NUC_MAX = 40

HEADER = 'Nuc_name,Nuc_posX,Nuc_posY,Nuc_posZ,Nuc_aa,Nuc_aaID,' \
         'Acid_name,Acid_posX,Acid_posY,Acid_posZ,Acid_aa,Acid_aaID,' \
         'Base_name,Base_posX,Base_posY,Base_posZ,Base_aa,Base_aaID,' \
         'Dist_Nuc_Acid,Dist_Acid_Base,Dist_Base_Nuc,Angle_Nuc,Angle_Acid,Angle_Base\n'


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
                    nuc_angle = util.find_angle(acid[0], nuc[0], base[0])
                    acid_angle = util.find_angle(nuc[0], acid[0], base[0])
                    base_angle = util.find_angle(nuc[0], base[0], acid[0])

                    if ((nuc_angle <= ANGLE_NUC_MAX)
                            and (acid_angle <= ANGLE_ACID_MAX)
                            and (base_angle <= ANGLE_BASE_MAX)):

                        # if angle fits criteria, test distances
                        if ((nuc[0] - acid[0] <= NUC_ACID_MAX)
                                and (nuc[0] - base[0] <= BASE_NUC_MAX)
                                and (acid[0] - base[0] <= ACID_BASE_MAX)):
                            # if distances fit, mark as a triad found
                            text_list.append(
                                parse_triangle_descriptors(nuc, base, acid, nuc_angle, acid_angle,
                                                           base_angle))

                        # if distances don't fit but angles do, mark as a similar triangle
                        else:
                            similar_text_list.append(parse_triangle_descriptors(nuc, base, acid, nuc_angle,
                                                                                acid_angle,
                                                                                base_angle))

        util.write_file(output_directory, HEADER, protein, ''.join(text_list), '')
        util.write_file(output_directory, HEADER, protein, ''.join(similar_text_list), 'similar_')


def find_triads():
    nuc_atoms, acid_atoms, base_atoms = util.store_atoms(directory)

    util.create_folder(output_directory)
    store_triads(nuc_atoms, acid_atoms, base_atoms)
