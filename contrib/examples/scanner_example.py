from deriver.api import Deriver
import random
random.seed(44)
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
results = my_deriver.scan_selfies()
print(results)

