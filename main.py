from pymatgen.io.cif import CifParser
import package.fg
parser = CifParser("rocksalt_alnit.cif")

structures = parser.parse_structures(primitive=False)
atoms_dict = {}
for i, site in enumerate(structures[0]*(1,1,1)):
    site = site.as_dict()
    print(site)
    atoms_dict[i] = {}
    atoms_dict[i] = {"element": site["species"][0]["element"], "frac_coord": site["abc"]}
