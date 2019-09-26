from deriver.api import Deriver


my_deriver = Deriver()

seeds = [
    "CC(NC)CC1=CC=C(OCO2)C2=C1",
    "CCCCCc1cc(c2c(c1)OC([C@H]3[C@H]2C=C(CC3)C)(C)C)O",
    "CN(CCC1=CNC2=C1C=CC=C2)C",
    "CN(C)CCc1c[nH]c2cccc(O)c12",
    "O(c1cc(cc(OC)c1OC)CCN)C",
    "CCN(CC)C(=O)[C@H]1CN([C@@H]2Cc3c[nH]c4c3c(ccc4)C2=C1)C",
    "COc1cc(CCN)c(OC)cc1Br",
    "O=C(OC)[C@H]2[C@@]3(CC[C@H]4C(=O)O[C@H](c1ccoc1)C[C@@]4([C@H]3C(=O)[C@@H](OC(=O)C)C2)C)C"
]

my_deriver.set_seeds(seeds)  # tell Deriver the starting molecules to build from
my_deriver.enable_and_expand_filter()  # enable molecule filtering
my_deriver.manual_filter_set("MW", 100, 500)  # set the molecular weight limits to be 100-500, rather than the default
my_deriver.set_filter_molecules(["CCCC"])  # if this molecule is generated don't keep it
my_deriver.set_must_have_patterns(["[#6]-[#6](=[#8])-[#8]"])  # only include molecules with this motif

my_deriver.derive_selfies(n_children=10_000)

print(f"Produced {len(my_deriver.data.all_good_selfies_children)} "
      f"good children:\n{my_deriver.data.all_good_selfies_children}")

# now you can change the seeds and generate new children, for a second generation for example:
my_deriver.set_seeds(my_deriver.data.all_good_selfies_children[0:5])
# disable the pattern requirement
my_deriver.set_must_have_patterns(None)
# make then filter 1000 children
my_deriver.derive_selfies(n_children=1_000)
print(f"\nGeneration 2:\n"
      f"Produced {len(my_deriver.data.all_good_selfies_children)} "
      f"good children:\n{my_deriver.data.all_good_selfies_children}")