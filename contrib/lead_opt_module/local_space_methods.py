from deriver.api import Deriver
import random
random.seed(42)

test_smiles = ["O=C(C)Oc1ccccc1C(=O)O"]

# test crem
crem_db = "crem/replacements02_sc2.5.db"  # obtained from the crem github repo, see readme.txt
myderiver = Deriver()
myderiver.set_seeds(test_smiles)
myderiver.set_crem_source_db(crem_db=crem_db)
myderiver.enable_and_expand_filter()
myderiver.derive_local_space()

print(myderiver.data.all_good_local_children)
