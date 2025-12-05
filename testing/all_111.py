from finding_generators import AtomPermutator, GeneratorFinder
import json

if __name__ == "__main__":
    for n in range(1, 8):
        ap = AtomPermutator("../rocksalt_alnit.cif", (1,1,1))
        ap.dope_cell(n, "Sc")
        ap.run_permutation("all_configurations.json")
        gf = GeneratorFinder("all_configurations.json", "Fm-3m")
        gf.start_reduction(f"ired_configurations_{n}.json")
    for name in [f"ired_configurations_{n}.json" for n in range(1,8)]:
        print(name)
        with open(name) as f:
            js = json.load(f)
            parse_dict()
            f.close()



