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

if __name__ == "__main__":
    protein_atoms = {}

    for filename in os.listdir(directory):
        if filename.endswith('.pdb'):
            protein = filename[0:4]
            structure = parser.get_structure(protein, os.path.join(directory, filename))

            atoms_tmp = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            atoms_tmp.append(atom)
            protein_atoms[protein] = atoms_tmp
