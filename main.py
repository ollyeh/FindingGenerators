from pymatgen.io.cif import CifParser
from package import finding_generators as fg
import numpy as np
import math
import time

parser = CifParser("rocksalt_alnit.cif")

structures = parser.parse_structures(primitive=False)
structure_dict = {}

for i, site in enumerate(structures[0]*(1,1,1)):
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

#structure_dict[i] = {"element": site["species"][0]["element"], "frac_coord": site["abc"]}

framework = fg.Framework(structure_dict, crystal_system="cubic")
#framework.dope_cell(4, "Sc")
framework.display_cell()
framework.dope_cell(4, "Sc")
