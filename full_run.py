import time

from pymatgen.io.cif import CifParser
from pymatgen.symmetry.groups import PointGroup, SpaceGroup
from FindingGenerators.finding_generators import AtomPermutator, GeneratorFinder
from ase import Atoms, Atom
from ase.geometry import cellpar_to_cell
from ase.io import write
import numpy as np
import math
import json

def parse_inp(rescale: tuple[int, int, int]):
    parser = CifParser("rocksalt_alnit.cif")
    structure = parser.parse_structures(primitive=False)[0]*rescale
    a, b, c = structure.lattice.abc
    alpha, beta, gamma = structure.lattice.angles
    lattice_matrix = structure.lattice.matrix
    return structure, lattice_matrix

def run_permutation(structure):

    structure_dict = {}

    for i, site in enumerate(structure):
        site = site.as_dict()
        structure_dict[i] = {}

        frac_coord = np.array(site["abc"])

        for j in range(3):
            frac_coord[j] -= math.floor(frac_coord[j])
            if np.allclose(frac_coord[j], 1) or np.allclose(frac_coord[j], 0):
                frac_coord[j] = 0
            else:
                frac_coord[j] = round(frac_coord[j], 7)

        structure_dict[i] = {"element": site["species"][0]["element"], "frac_coord": frac_coord}

    structure_dict[i] = {"element": site["species"][0]["element"], "frac_coord": site["abc"]}

    ap = AtomPermutator(structure_dict)
    ap.dope_cell(2, "Sc")
    ap.run_permutation("/Users/oliver/Documents/programming/FindingGenerators/all_configurations.json")

def run_reduction():
    sg = SpaceGroup("Fm-3m")

    rotations = [op.rotation_matrix for op in sg.symmetry_ops]
    translations = [op.translation_vector for op in sg.symmetry_ops]

    i = 0
    transformations_dict = {}
    for rot, trans in zip(rotations, translations):
        #print("Rotation:")
        #print(rot)
        #print("Translation")
        #print(trans)
        #print(".................")

        transformations_dict[i] = {"trans": trans, "rot": rot}

        i += 1

    gf = GeneratorFinder("/Users/oliver/Documents/programming/FindingGenerators/all_configurations.json", transformations_dict)
    gf.start_reduction("/Users/oliver/Documents/programming/FindingGenerators/ired_configurations.json")

def rebuild_to_cif(lattice_matrix):
    with open("ired_configurations.json") as f:
        configurations = json.load(f)

    ase = []
    for config_key in configurations:
        symbols = []
        frac_coords = []
        for atom in configurations[config_key]:
            symbol, frac_coord = atom["element"], atom["frac_coord"]
            symbols.append(symbol)
            frac_coords.append(frac_coord)

        ase_atoms = Atoms(symbols=symbols, cell=lattice_matrix, scaled_positions=frac_coords, pbc=True)
        ase.append(ase_atoms)
        write(f"alnit_111_doped/4/{config_key}_4.xyz", ase_atoms)

if __name__ == "__main__":
    structures, lattice_matrix = parse_inp((3, 3, 3))
    run_permutation(structures)
    #run_reduction()
    #rebuild_to_cif(lattice_matrix)