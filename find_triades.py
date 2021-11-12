"""Read transformed files, make an Atom object list for each protein, find triades in each protein."""

import os
import configparser
from Bio.PDB.PDBParser import PDBParser

# QUIET=True zbog PDBConstructionWarninga koji nam dolaze jer on vidi greške u strukturi peptida
# budući da se bavimo samo pojedinim atomima, te su nam greške nevažne?
parser = PDBParser(PERMISSIVE=1, QUIET=True)

config = configparser.ConfigParser()
config.read('config.ini')
config = config['default']

directory = config['transformed_location']
output_directory = config['output_location']

NUC_ACID_MIN = 6.68
NUC_ACID_MAX = 7.43
NUC_BASE_MIN = 4.4
NUC_BASE_MAX = 4.83
ACID_BASE_MIN = 2.6
ACID_BASE_MAX = 2.8


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


def store_triads():
    for protein in protein_atoms.keys():
        atoms = protein_atoms[protein]
        for nuc in atoms:
            if (nuc[0].get_name() == 'OG' and nuc[1].get_resname() == 'SER') or (
                    nuc[0].get_name() == 'SG' and nuc[1].get_resname() == 'CYS'):
                for base in atoms:
                    if base[0].get_name() == 'CG' and (
                            base[1].get_resname() == 'HIS' or base[1].get_resname() == 'ASP'
                            or base[1].get_resname() == 'GLU'):
                        for acid in atoms:
                            if acid[0].get_name() == 'OD1' or acid[0].get_name() == 'OD2':
                                # print(nuc[0] - acid[0])
                                # print(nuc[0].get_coord() - acid[0].get_coord())
                                # print(nuc[0].get_coord())
                                # print(acid[0].get_coord())
                                # print(nuc[0].get_name())
                                # print(acid[0].get_name())
                                # break
                                if ((NUC_ACID_MIN <= nuc[0] - acid[0] <= NUC_ACID_MAX)
                                        and (NUC_BASE_MIN <= nuc[0] - base[0] <= NUC_BASE_MAX)
                                        and (ACID_BASE_MIN <= acid[0] - base[0] <= ACID_BASE_MAX)):
                                    triad = [nuc, base, acid]
                                    triads.append((protein, triad))


if __name__ == "__main__":
    protein_atoms = {}
    store_atoms(protein_atoms)

    triads = []
    store_triads()

    print(triads)

    # distance = atom1 - atom2
