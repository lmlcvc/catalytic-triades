"""Read transformed files, make an Atom object list for each protein, find triades in each protein."""
import configparser
import util
from Bio.PDB.PDBParser import PDBParser

parser = PDBParser(PERMISSIVE=1, QUIET=True)

config = configparser.ConfigParser()
config.read('config.ini')
config = config['default']

directory = config['transformed_location']
output_directory = config['output_location']

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


# TODO: lesser triads found when separating nuc/acid/base than before?
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

        util.write_file(output_directory, HEADER, protein, ''.join(text_list), '')
        util.write_file(output_directory, HEADER, protein, ''.join(similar_text_list), 'similar_')

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


def find_triads():
    nuc_atoms, acid_atoms, base_atoms = util.store_atoms(directory)

    util.create_folder(output_directory)
    store_triads(nuc_atoms, acid_atoms, base_atoms)
