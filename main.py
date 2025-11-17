import time

from pymatgen.io.cif import CifParser
from pymatgen.analysis.local_env import VoronoiNN, get_neighbors_of_site_with_index
from pymatgen.symmetry.groups import PointGroup, SpaceGroup
from package import finding_generators as fg
import numpy as np
import math
import json


voronoi = VoronoiNN(compute_adj_neighbors=False)

parser = CifParser("zeug.cif")

structures = parser.parse_structures(primitive=False)

sg = SpaceGroup("Fm-3m")
# x -> x*rotation + translation

rotations = [op.rotation_matrix for op in sg.symmetry_ops]
translations = [op.translation_vector for op in sg.symmetry_ops]

structure_dict = {}
neighbor_dict = {}

with open("temp.json", "r") as f:
    neighbor_dict = json.load(f)
    f.close()

neighbor_dict = {int(key): val for (key, val) in zip(neighbor_dict.keys(), neighbor_dict.values())}

for i, site in enumerate(rescaled_structure := structures[0]*(5,5,5)):
    site = site.as_dict()
    structure_dict[i] = {}

    #neighbor_dict[i] = [int(i) for i in np.unique([int(n["site_index"]) for n in voronoi.get_nn_info(rescaled_structure, i)])]

    frac_coord = np.array(site["abc"])

    for j in range(3):
        frac_coord[j] -= math.floor(frac_coord[j])
        if np.allclose(frac_coord[j], 1) or np.allclose(frac_coord[j], 0):
            frac_coord[j] = 0
        else:
            frac_coord[j] = round(frac_coord[j], 7)

    structure_dict[i] = {"element": site["species"][0]["element"], "frac_coord": frac_coord}

#with open("temp.json", "w") as f:
#    json.dump(neighbor_dict, f)
#    f.close()

structure_dict[i] = {"element": site["species"][0]["element"], "frac_coord": site["abc"]}

print(neighbor_dict)

framework = fg.Framework(structure_dict, neighbor_dict, crystal_system="cubic")
framework.dope_cell(100, "Sc")
#framework.display_cell()
#framework.display_cell()
framework.run_simulation(True)
#framework.dope_cell(4, "Sc")
