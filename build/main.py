import mod

from pymatgen.core import Structure

if __name__ == "__main__":
    structure = Structure.from_file("AlN.cif")

    structure_dict = {}
    for i, site in enumerate(structure.sites):
        structure_dict[i] = {}
        symbol = site.species_string  # chemical symbol
        coords = site.coords  # Cartesian coordinates in Å
        structure_dict[i]["symbol"] = symbol
        structure_dict[i]["coords"] = coords

    swapper = mod.Swapper(structure_dict)
