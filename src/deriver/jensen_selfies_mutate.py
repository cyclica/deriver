"""
Written by Emilie S. Henault and Jan H. Jensen 2019, copied 04/2020
"""
import random

from .jensen_crossover import mol_OK
from .jensen_selfies_crossover import mol2string, string2mol

from rdkit import Chem, rdBase
rdBase.DisableLog('rdApp.error')


def get_symbols():
    symbols = ['C', 'Branch1_2', 'epsilon', 'Branch1_3', '=C', 'O', '#N', '=O', 'N', 'Ring1',
               'Branch1_1', 'F', '=N', '#C', 'C@@H', 'S', 'Branch2_2', 'Ring2', 'Branch2_3',
               'Branch2_1', 'Cl', 'O-', 'C@H', 'NH+', 'C@', 'Br', '/C', '/O', 'NH3+', '=S', 'NH2+',
               'C@@', '=N+', '=NH+', 'N+', '\\C', '\\O', '/N', '/S', '\\S', 'S@', '\\O-', 'N-', '/NH+',
               'S@@', '=NH2+', '/O-', 'S-', '/S-', 'I', '\\N', '\\Cl', '=P', '/F', '/C@H', '=OH+',
                '\\S-', '=S@@', '/C@@H', 'P', '=S@', '\\C@@H', '/S@', '/Cl', '=N-', '/N+', 'NH-',
                '\\C@H', 'P@@H', 'P@@', '\\N-', 'Expl\\Ring1', '=P@@', '=PH2', '#N+', '\\NH+', 'P@',
                'P+', '\\N+', 'Expl/Ring1', 'S+', '=O+', '/N-', 'CH2-', '=P@', '=SH+', 'CH-', '/Br',
                '/C@@', '\\Br', '/C@', '/O+', '\\F', '=S+', 'PH+', '\\NH2+', 'PH', '/NH-', '\\S@', 'S@@+',
                '/NH2+', '\\I']

    return symbols


def mutate(mol, mutation_rate):
    if random.random() > mutation_rate:
        return mol
    Chem.Kekulize(mol, clearAromaticFlags=True)
    child = mol2string(mol)
    symbols = get_symbols()
    for i in range(50):
        mutated_gene = random.randint(0, len(child) - 1)
        random_symbol_number = random.randint(0, len(symbols)-1)
        new_child = list(child)
        new_child[mutated_gene] = symbols[random_symbol_number]
        new_child_mol = string2mol(new_child)
        if mol_OK(new_child_mol):
            return new_child_mol

    return mol
