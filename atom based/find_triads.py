"""Read transformed files, make an Atom object list for each protein, find triads in each protein."""
import os
import configparser
import util
from Bio.PDB.PDBParser import PDBParser

parser = PDBParser(PERMISSIVE=1, QUIET=True)

config = configparser.ConfigParser()
config.read(os.path.join(os.pardir, 'config.ini'))
config = config['atom']

directory = config['transformed_location']
output_directory = config['output_location']

ACID_BASEN1_MAX = 4
BASEN1_BASEN2_MAX = 2.5
BASEN2_NUC_MAX = 4.5
NUC_ACID_MAX = 9

ANGLE_ACID_MAX = 50
ANGLE_N1_MAX = 180
ANGLE_N2_MAX = 180
ANGLE_NUC_MAX = 33

HEADER = 'Nuc_name,Nuc_posX,Nuc_posY,Nuc_posZ,Nuc_aa,Nuc_aaID,' \
         'Acid_name,Acid_posX,Acid_posY,Acid_posZ,Acid_aa,Acid_aaID,' \
         'BaseN1_name,BaseN1_posX,BaseN1_posY,BaseN1_posZ,BaseN1_aa,BaseN1_aaID,' \
         'BaseN2_name,BaseN2_posX,BaseN2_posY,BaseN2_posZ,BaseN2_aa,BaseN2_aaID,' \
         'Dist_Acid_N2,Dist_N2_Nuc,Dist_Nuc_Acid,Dist_Acid_N1,Dist_N1_N2,' \
         'Angle_Nuc,Angle_Acid,Angle_N1,Angle_N2\n'


def parse_triangle_descriptors(nuc, n1, n2, acid, nuc_angle, n1_angle, n2_angle, acid_angle):
    nuc_coord = nuc[0].get_coord()
    n1_coord = n1[0].get_coord()
    n2_coord = n2[0].get_coord()
    acid_coord = acid[0].get_coord()

    return (str(nuc[0].get_name()) + ',' + str(nuc_coord[0]) + ',' + str(nuc_coord[1]) + ','
            + str(nuc_coord[2]) + ',' + str(nuc[1].get_resname()) + ',' + str(nuc[1].get_id()[1]) + ','
            + str(acid[0].get_name()) + ',' + str(acid_coord[0]) + ',' + str(acid_coord[1]) + ','
            + str(acid_coord[2]) + ',' + str(acid[1].get_resname()) + ',' + str(acid[1].get_id()[1]) + ','
            + str(n1[0].get_name()) + ',' + str(n1_coord[0]) + ',' + str(n1_coord[1]) + ','
            + str(n1_coord[2]) + ',' + str(n1[1].get_resname()) + ',' + str(n1[1].get_id()[1]) + ','
            + str(n2[0].get_name()) + ',' + str(n2_coord[0]) + ',' + str(n2_coord[1]) + ','
            + str(n2_coord[2]) + ',' + str(n2[1].get_resname()) + ',' + str(n2[1].get_id()[1]) + ','
            + str(acid[0] - n2[0]) + ',' + str(n2[0] - nuc[0]) + ',' + str(nuc[0] - acid[0]) + ','
            + str(acid[0] - n1[0]) + ',' + str(n1[0] - n2[0]) + ','
            + str(nuc_angle) + ',' + str(acid_angle) + ',' + str(n1_angle) + ',' + str(n2_angle) + '\n')


def find_candidates(nuc_atoms, acid_atoms, base_atoms):
    """
    Find and store triad candidates by each protein.

    :param nuc_atoms: nuc candidate atoms in each protein
    :param acid_atoms: acid candidate atoms in each protein
    :param base_atoms: base (N1/N2) candidate atoms in each protein
    :return: dict of {protein : triad candidates in that protein}
    """

    protein_candidates = {}

    for protein in nuc_atoms.keys():

        nucs = nuc_atoms[protein]
        acids = acid_atoms[protein]
        bases = base_atoms[protein]

        acid_n1 = []
        acid_n1_n2 = []
        acid_n1_n2_nuc = []

        # find ACID-N1 pairs that fit distance criteria
        for acid in acids:
            for base in bases:
                if acid[0] - base[0] < ACID_BASEN1_MAX:
                    tmp_dict = {"acid": acid, "n1": base}
                    acid_n1.append(tmp_dict)

        # find (ACID-N1)-N2 sets that fit distance criteria
        # N2 must not be the same atom as N1
        for combination in acid_n1:
            for base in bases:
                if (combination["n1"] != base) and (combination["n1"][0] - base[0] < BASEN1_BASEN2_MAX):
                    tmp_dict = {"acid": combination["acid"], "n1": combination["n1"], "n2": base}
                    acid_n1_n2.append(tmp_dict)

        # find (ADIC-N1-N2)-NUC sets that fit distance criteria
        # NUC must not be the same atom as ACID
        for combination in acid_n1_n2:
            for nuc in nucs:
                if (combination["acid"] != nuc) and (combination["n2"][0] - nuc[0] < BASEN2_NUC_MAX):
                    tmp_dict = {"acid": combination["acid"], "n1": combination["n1"], "n2": combination["n2"],
                                "nuc": nuc}
                    acid_n1_n2_nuc.append(tmp_dict)

        protein_candidates[protein] = acid_n1_n2_nuc

    return protein_candidates


def store_triads(protein_candidates):
    """
    Filter out triad candidates that do not fit angle criteria.
    Write filtered and formatted triad list to csv file.

    :param protein_candidates: list of all triad candidates by protein
    """

    for protein in protein_candidates.keys():
        text_list = []

        for candidate in protein_candidates[protein]:
            nuc = candidate["nuc"]
            n1 = candidate["n1"]
            n2 = candidate["n2"]
            acid = candidate["acid"]

            nuc_angle = util.find_angle(n2[0], nuc[0], acid[0])
            n1_angle = util.find_angle(acid[0], n1[0], n2[0])
            n2_angle = util.find_angle(acid[0], n2[0], nuc[0])
            acid_angle = util.find_angle(n1[0], acid[0], nuc[0])

            if (acid_angle < ANGLE_ACID_MAX) \
                    and (n1_angle < ANGLE_N1_MAX) \
                    and (n2_angle < ANGLE_N2_MAX) \
                    and (nuc_angle < ANGLE_NUC_MAX):
                text_list.append(
                    parse_triangle_descriptors(nuc, n1, n2, acid, nuc_angle, n1_angle, n2_angle, acid_angle))

        util.write_file(output_directory, HEADER, protein, ''.join(text_list))


def find_triads():
    nuc_atoms, acid_atoms, base_atoms = util.store_atoms(directory)

    util.create_folder(output_directory)
    protein_candidates = find_candidates(nuc_atoms, acid_atoms, base_atoms)
    store_triads(protein_candidates)
