"""
Modified from code by Emilie S. Henault and Jan H. Jensen 2019
"""
import random
from .jensen_selfies_crossover import mol2string, string2mol
from selfies import get_semantic_robust_alphabet
from rdkit import Chem, rdBase
rdBase.DisableLog('rdApp.error')
symbols = [symbol.strip("[]") for symbol in get_semantic_robust_alphabet()]


def mutate(mol):
    child = mol2string(mol)
    for _ in range(50):
        mutated_gene = random.randint(0, len(child) - 1)
        random_symbol_number = random.randint(0, len(symbols)-1)
        new_child = list(child)
        new_child[mutated_gene] = symbols[random_symbol_number]
        new_child_mol = string2mol(new_child)
        return new_child_mol

    return mol
