from finding_generators import AtomPermutator, GeneratorFinder

if __name__ == "__main__":
    ap = AtomPermutator("../rocksalt_alnit.cif", (2,2,2))
    ap.dope_cell(5, "Mo")
    ap.run_permutation("all_configurations.json")
    gf = GeneratorFinder("all_configurations.json", "Fm-3m")
    gf.start_reduction("ired_configurations.json")
