from pymatgen.symmetry.groups import PointGroup, SpaceGroup
from package import finding_generators as fg

#sg = SpaceGroup("P-1")
sg = SpaceGroup("Fm-3m")
# x -> x*rotation + translation

rotations = [op.rotation_matrix for op in sg.symmetry_ops]
translations = [op.translation_vector for op in sg.symmetry_ops]

i = 0
transformations_dict = {}
for rot, trans in zip(rotations, translations):
    print("Rotation:")
    print(rot)
    print("Translation")
    print(trans)
    print(".................")

    transformations_dict[i] = {"trans": trans, "rot": rot}

    i += 1


fg = fg.GeneratorFinder("/Users/oliver/Documents/programming/FindingGenerators/configurations.json", transformations_dict)
fg.start_search()