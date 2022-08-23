"""Read transformed files, make an Atom object list for each protein, find triads in each protein."""
import os
import configparser
import util
from Bio.PDB.PDBParser import PDBParser

parser = PDBParser(PERMISSIVE=1, QUIET=True)

config = configparser.ConfigParser()
config.read(os.path.join(os.pardir, 'config.ini'))
config = config['hybrid']

directory = config['transformed_location']
output_directory = config['output_location']

ACID_ND1_MAX = 4
ND1_NE2_MAX = 2.5
NE2_NUC_MAX = 4.5
NUC_ACID_MAX = 9

ANGLE_ACID_MAX = 50
ANGLE_N1_MAX = 180
ANGLE_N2_MAX = 180
ANGLE_NUC_MAX = 33

HEADER = 'Nuc_name,Nuc_posX,Nuc_posY,Nuc_posZ,Nuc_aa,Nuc_aaID,' \
         'Acid_name,Acid_posX,Acid_posY,Acid_posZ,Acid_aa,Acid_aaID,' \
         'ND1_name,ND1_posX,ND1_posY,ND1_posZ,ND1_aa,ND1_aaID,' \
         'NE2_name,NE2_posX,NE2_posY,NE2_posZ,NE2_aa,NE2_aaID,' \
         'Dist_Acid_NE2,Dist_NE2_Nuc,Dist_Nuc_Acid,Dist_Acid_ND1,Dist_ND1_NE2,' \
         'Angle_Nuc,Angle_Acid\n'


def parse_triangle_descriptors(nuc, nd1, ne2, acid, nuc_angle, acid_angle):
    nuc_coord = nuc[0].get_coord()
    nd1_coord = nd1[0].get_coord()
    ne2_coord = ne2[0].get_coord()
    acid_coord = acid[0].get_coord()

    return (str(nuc[0].get_name()) + ',' + str(nuc_coord[0]) + ',' + str(nuc_coord[1]) + ','
            + str(nuc_coord[2]) + ',' + str(nuc[1].get_resname()) + ',' + str(nuc[1].get_id()[1]) + ','
            + str(acid[0].get_name()) + ',' + str(acid_coord[0]) + ',' + str(acid_coord[1]) + ','
            + str(acid_coord[2]) + ',' + str(acid[1].get_resname()) + ',' + str(acid[1].get_id()[1]) + ','
            + str(nd1[0].get_name()) + ',' + str(nd1_coord[0]) + ',' + str(nd1_coord[1]) + ','
            + str(nd1_coord[2]) + ',' + str(nd1[1].get_resname()) + ',' + str(nd1[1].get_id()[1]) + ','
            + str(ne2[0].get_name()) + ',' + str(ne2_coord[0]) + ',' + str(ne2_coord[1]) + ','
            + str(ne2_coord[2]) + ',' + str(ne2[1].get_resname()) + ',' + str(ne2[1].get_id()[1]) + ','
            + str(acid[0] - ne2[0]) + ',' + str(ne2[0] - nuc[0]) + ',' + str(nuc[0] - acid[0]) + ','
            + str(acid[0] - nd1[0]) + ',' + str(nd1[0] - ne2[0]) + ','
            + str(nuc_angle) + ',' + str(acid_angle) + '\n')


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

        acid_nd1 = []
        acid_nd1_ne2 = []
        acid_nd1_ne2_nuc = []

        # find ACID-N1 pairs that fit distance criteria
        for acid in acids:
            for base in bases:
                if base[0].get_name() == "ND1" \
                        and acid[0] - base[0] < ACID_ND1_MAX:
                    tmp_dict = {"acid": acid, "ND1": base}
                    acid_nd1.append(tmp_dict)

        # find (ACID-ND1)-NE2 sets that fit distance criteria
        # N2 must not be the same atom as N1
        for combination in acid_nd1:
            for base in bases:
                if base[0].get_name() == "NE2" \
                        and combination["ND1"][0] - base[0] < ND1_NE2_MAX:
                    tmp_dict = {"acid": combination["acid"], "ND1": combination["ND1"], "NE2": base}
                    acid_nd1_ne2.append(tmp_dict)

        # find (ADIC-ND1-NE2)-NUC sets that fit distance criteria
        for combination in acid_nd1_ne2:
            for nuc in nucs:
                if combination["acid"] != nuc \
                        and combination["NE2"][0] - nuc[0] < NE2_NUC_MAX:
                    tmp_dict = {"acid": combination["acid"], "ND1": combination["ND1"], "NE2": combination["NE2"],
                                "nuc": nuc}
                    acid_nd1_ne2_nuc.append(tmp_dict)

        protein_candidates[protein] = acid_nd1_ne2_nuc

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
            nd1 = candidate["ND1"]
            ne2 = candidate["NE2"]
            acid = candidate["acid"]

            nuc_angle = util.find_angle(ne2[0], nuc[0], acid[0])
            acid_angle = util.find_angle(nd1[0], acid[0], nuc[0])

            if (acid_angle < ANGLE_ACID_MAX) \
                    and (nuc_angle < ANGLE_NUC_MAX):
                text_list.append(
                    parse_triangle_descriptors(nuc, nd1, ne2, acid, nuc_angle, acid_angle))

        util.write_file(output_directory, HEADER, protein, ''.join(text_list))


def find_triads():
    nuc_atoms, acid_atoms, base_atoms = util.store_atoms(directory)

    util.create_folder(output_directory)
    protein_candidates = find_candidates(nuc_atoms, acid_atoms, base_atoms)
    store_triads(protein_candidates)
