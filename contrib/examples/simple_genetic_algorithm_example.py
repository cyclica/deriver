from deriver.api import Deriver
from deriver.child_filter import get_filter_values
from rdkit import Chem

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
my_deriver.manual_filter_set("MW", 100, 1000)  # set the molecular weight limits to be 100-1000, rather than the default


def get_fitness(smile):
    mol = Chem.MolFromSmiles(smile)
    values = get_filter_values(mol)
    fitness = abs(500 - values["MW"]) + abs(3 - values["logP"]) + abs(2 - values["num_rings"])
    return fitness

top_mols_per_generation = []
for generation in range(10):  # go for ten generations
    # get the first generation
    good_children = my_deriver.derive_selfies(n_children=10_000)
    scored_children = [(child, get_fitness(child)) for child in good_children]
    scored_children.sort(key=lambda tup: tup[1])  # sort by the score, ascending
    best_children = scored_children[:25]  # the 25 best children
    top_mols_per_generation.append(best_children)
    print(f"Finished generation {generation}! Best children were: {best_children}")
    best_smiles = [smile for smile, score in best_children]
    # now set the seeds for generation i+1
    my_deriver.set_seeds(best_smiles)

print("Finished!")
print(top_mols_per_generation)


"""
Finished generation 9! 
Best children were: [
('C[C]CC1C=COOO[C@H]1OC(=O)CCC1[C](C)[CH]C(=O)[C@@H](OC(C)=O)C[C@H]1C(=O)Br', 0.07576447999998948), 
('C[C]C1CC=COO[C@H]1OC(=O)CCC1[C](C)[CH]OC(=O)[C@@H](OC(C)=O)C[C@H]1C(=O)Br', 0.07576447999998948), 
('C[C]CC1C=COOO[C@H]1OC(=O)CCC1[C](C)[CH]C(=O)[C@@H](OC(C)=O)C[C@H]1C(=O)Br', 0.07576447999998948), 
('[C]CC1C=COOO[C@H]1OC(=O)CO[C]C1C[C@H](C(=O)Br)C[C@H](OC(C)=O)C(C)[C]1C', 0.09066447999998939), 
('[C]CC1C(=CO)OO[C@H]1OC(=O)CO[C]C1C[C@H](C(=O)Br)C[C@H](OC(C)=O)C(C)[C]1C', 0.0916244799999868), 
('C[C]N1CC(=CBr)O[C@H]1OC(=O)CC1O[C](C)[CH]C(=O)[C@@H](OC(=O)O)C[C@H]1C(=N)CC', 0.13326785999995394), 
('C[C]C1CC2=C(Br)[C@H](C(=N)ON)C[C@H](CC(=O)O)C(=O)[CH][C](C)CCOC(=O)O[C@@H]1O2', 0.13328785999983772), 
('C[C]N1CC2=C(Br)[C@H](C(=N)OC)C[C@H](OC(=O)O)C(=O)[CH][C](C)CCCC(=O)O[C@@H]1O2', 0.13486785999984097), 
('C[C]C1CC(=CBr)O[C@H]1OC(=O)C1CC[C](C)[CH]C(=O)[C@@H](OC(=O)O)C[C@H]1C(=N)ON', 0.134887859999953), 
('C=C(OC)[C@@H]1C[C@H](OC(C)=N)C(=O)[CH][C@@]12OCC[CH]C(=O)C[C@H](ONBr)C[CH][C](C)O2', 0.15098336799991863), 
('C=C(C)O[C@@H]1C(=C)[CH]C(O)(O)OC2C[CH]C(=O)O[C@H](COBr)C[CH][C]2CC[C]OC1C', 0.16242998800001018), 
('C[C]N1CC(=CBr)O[C@H]1OC(=O)CCC1[C](C)[CH]C(=O)[C@@H](OC(=O)O)C[C@H]1C(=N)OC', 0.16808785999989517), 
('C[C]N1CC(=CBr)O[C@H]1OC(=O)CCC1[C](C)[CH]C(=O)[C@@H](OC(=O)O)C[C@H]1C(=N)OC', 0.16808785999989517), 
('C=C1[CH][C](C)C(CCC(=O)O[C@@H]2OC(=CBr)NN2[C]C)[C@H](C(=O)OC)C[C@@H]1OC(=O)O', 0.16955785999995232), 
('COC(=O)[C@@H]1C[C@H](OC(C)=O)C(=O)[CH][C@@]1(C)CC1C[CH]C(=O)O[C@H](OCBr)C[CH][C]1C', 0.17251998799995416), 
('COC(=O)[C@@H]1C[C@H](OC(C)=O)C(=O)[CH][C@H]1CCCC1[CH]C(=O)O[C@H](OCBr)C[CH][C]1C', 0.17251998799995416), 
('COC(=O)[C@@H]1C[C@H](OC(C)=O)C(=O)[CH][C@H]1CC1CC[CH]C(=O)O[C@H](OCBr)C[CH][C]1C', 0.17251998799995416), 
('COC(=O)[C@@H]1C[C@H](CC(C)=O)C(=O)[CH][C@H]1CCC1O[CH]C(=O)O[C@H](OCBr)C[CH][C]1C', 0.17251998799995416), 
('COC(=O)[C@@H]1C[C@H](OOC=O)C(=O)[CH][C@H]1CCC1C[CH]C(=O)O[C@H](CCBr)C[CH][C]1C', 0.17251998799995416), 
('COC(=O)[C@@H]1C[C@H](OC(C)=O)C(=O)[CH][C@H]1CCCC[CH]C(=O)O[C@H]1O[CH][C](C)C1CBr', 0.17251998799995416), 
('CCC(=O)[C@@H]1C[C@H](OC(C)=O)C(=O)[CH][C@H]1CC1CO[CH]C(=O)O[C@H](CCBr)O[CH][C]1C', 0.17251998800001012), 
('COC(=O)[C@@H]1C[C@H](OC(C)=O)C(=O)[CH][C@H]1C1CCC[CH]C(=O)O[C@H](CCBr)O[CH][C]1C', 0.172519988000011), 
('COC(=O)[C@@H]1C[C@H](OC(C)=O)C(=O)[CH][C@H]1CCC1C[CH]C(=O)O[C@H](CCBr)O[CH][C]1C', 0.172519988000011), 
('COC(=O)[C@@H]1C[C@H](OC(C)=O)C(=O)[CH][C@H]1CCC1C[CH]C(=O)O[C@H](CCBr)O[CH][C]1C', 0.172519988000011), 
('COC(=O)[C@@H]1C[C@H](OC(C)=O)C(=O)[CH][C@H]1CCC1C[CH]C(=O)O[C@H](CCBr)O[CH][C]1C', 0.172519988000011)
]
"""