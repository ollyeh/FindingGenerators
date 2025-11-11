from pymatgen.io.cif import CifParser
from package import finding_generators as fg

parser = CifParser("rocksalt_alnit.cif")

structures = parser.parse_structures(primitive=False)
structure_dict = {}
for i, site in enumerate(structures[0]*(1,1,1)):
    site = site.as_dict()
    structure_dict[i] = {}
    structure_dict[i] = {"element": site["species"][0]["element"], "frac_coord": site["abc"]}

framework = fg.Framework(structure_dict, crystal_system="cubic")
framework.dope_cell(4, "Sc")
framework.display_cell()