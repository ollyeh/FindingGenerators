import finding_generators.finding_generators_binding as fg
import typing
from ase import Atoms
from pymatgen.io.cif import CifParser
import numpy as np
from numpy.typing import NDArray
import math
from pymatgen.symmetry.groups import PointGroup, SpaceGroup

class AtomPermutator:
    def __init__(self, input: str | Atoms, rescale_factors: tuple[int, int, int] = (1,1,1)):
        self.structure_dict: dict[int, dict[str, str | NDArray[float]]] = {}
        self.lattice_matrix: NDArray[tuple[np.float64]]
        if isinstance(input, str):
            self.prepare_structure_from_str(input, rescale_factors)
        elif isinstance(input, Atoms):
            raise NotImplementedError("Implementation is coming soon. Stick to using a cif file.")
        else:
            raise TypeError(f"Got {type(input)} as input.")

        self.atom_permutator = fg.AtomPermutator(self.structure_dict)


    def prepare_structure_from_str(self, cif_path: str, rescale_factors: tuple[int, int, int]) -> None:
        structure = CifParser(cif_path).parse_structures(primitive=False)[0]*rescale_factors
        self.lattice_matrix = structure.lattice.matrix
        for i, site in enumerate(structure):
            site = site.as_dict()
            self.structure_dict[i] = {}

            frac_coord = np.array(site["abc"])

            for j in range(3):
                frac_coord[j] -= math.floor(frac_coord[j])
                if np.allclose(frac_coord[j], 1) or np.allclose(frac_coord[j], 0):
                    frac_coord[j] = 0
                else:
                    frac_coord[j] = round(frac_coord[j], 7)

            self.structure_dict[i] = {"element": site["species"][0]["element"], "frac_coord": frac_coord}

    def dope_cell(self, n_swaps: int, element: str) -> None:
        self.atom_permutator.dope_cell(n_swaps, element)

    def run_permutation(self, all_configurations_path: str):
        self.atom_permutator.run_permutation(all_configurations_path)
        # after this we do not need the permutator anymore -> delete it by overwriting
        self.atom_permutator = None

class GeneratorFinder:
    def __init__(self, all_configurations_path: str, space_group: str):
        if not all_configurations_path.endswith(".json"):
            if set(all_configurations_path.split(".")[1:]) != set(".json"):
                raise TypeError(f"Need json object.")

        transformations_dict = {}
        sg = SpaceGroup(space_group)
        for i, op in enumerate(sg.symmetry_ops):
            transformations_dict[i] = {"trans": op.translation_vector, "rot": op.rotation_matrix}

        self.generator_finder = fg.GeneratorFinder(all_configurations_path, transformations_dict)

    def start_reduction(self, ired_configurations_path: str):
        self.generator_finder.start_reduction(ired_configurations_path)
        # after this we do not need the finder anymore -> delete it by overwriting
        self.generator_finder = None
