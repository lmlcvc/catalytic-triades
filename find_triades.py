"""Read transformed files, make an Atom object list for each protein, find triades in each protein."""
import math
import os
import configparser
import numpy as np
import vg
from Bio.PDB.PDBParser import PDBParser

# QUIET=True zbog PDBConstructionWarninga koji nam dolaze jer on vidi greške u strukturi peptida
# budući da se bavimo samo pojedinim atomima, te su nam greške nevažne?
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


def store_atoms(protein_atoms_tmp):
    for filename in os.listdir(directory):
        if filename.endswith('.pdb'):
            protein = filename[0:4]
            structure = parser.get_structure(protein, os.path.join(directory, filename))

            atoms_tmp = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            atoms_tmp.append((atom, residue))
            protein_atoms_tmp[protein] = atoms_tmp


def find_angle(u, v, w):
    a = u - v
    b = v - w
    c = w - u

    return math.degrees(math.acos((a*a + b*b - c*c) / (2 * a * b)))


def parse_triangle_descriptors(nuc, base, acid):
    return ('Nuc: ' + str(nuc[0].get_name()) + ' ' + str(nuc[0].get_coord())
    + ' ' + str(nuc[1].get_resname()) + ' ' + str(nuc[1].get_id()) + '\n'
    + 'Acid: ' + str(acid[0].get_name()) + ' ' + str(acid[0].get_coord())
    + ' ' + str(acid[1].get_resname()) + ' ' + str(acid[1].get_id()) + '\n'
    + 'Base: ' + str(base[0].get_name()) + ' ' + str(base[0].get_coord())
    + ' ' + str(base[1].get_resname()) + ' ' + str(base[1].get_id()) + '\n'
    + 'Nuc - Acid: ' + str(nuc[0] - acid[0]) + '\n'
    + 'Acid - Base: ' + str(acid[0] - base[0]) + '\n'
    + 'Base - Nuc: ' + str(base[0] - nuc[0]) + '\n'
    + 'Nuc angle: ' + str(find_angle(acid[0], nuc[0], base[0])) + '° \n'
    + 'Acid angle: ' + str(find_angle(nuc[0], acid[0], base[0])) + '° \n'
    + 'Base angle: ' + str(find_angle(nuc[0], base[0], acid[0])) + '° \n\n')

def store_triads():
    for protein in protein_atoms.keys():
        atoms = protein_atoms[protein]
        text = ""
        for nuc in atoms:
            if (nuc[0].get_name() == 'OG' and nuc[1].get_resname() == 'SER') or (
                    nuc[0].get_name() == 'SG' and nuc[1].get_resname() == 'CYS'):
                for base in atoms:
                    if base[0].get_name() == 'CG' and (
                            base[1].get_resname() == 'HIS' or base[1].get_resname() == 'ASP'
                            or base[1].get_resname() == 'GLU'):
                        for acid in atoms:
                            if acid[0].get_name() == 'OD1' or acid[0].get_name() == 'OD2':
                                if ((NUC_ACID_MIN <= nuc[0] - acid[0] <= NUC_ACID_MAX)
                                        and (NUC_BASE_MIN <= nuc[0] - base[0] <= NUC_BASE_MAX)
                                        and (ACID_BASE_MIN <= acid[0] - base[0] <= ACID_BASE_MAX)):
                                    text += parse_triangle_descriptors(nuc, base, acid)
                                # triad = [nuc, base, acid]
                                # triads.append((protein, triad))
        file = open(output_directory + "/" + str(protein) + ".txt", "w")
        file.write(text)
        file.close()
        # print(text)
        # break


if __name__ == "__main__":
    protein_atoms = {}
    store_atoms(protein_atoms)

    triads = []
    store_triads()

    print(triads)