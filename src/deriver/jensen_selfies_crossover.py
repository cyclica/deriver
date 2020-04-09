"""
Written by Emilie S. Henault and Jan H. Jensen 2019, copied 04/2020
"""

from rdkit import Chem
import random
from selfies import encoder, decoder
from .jensen_crossover import mol_OK

from rdkit import rdBase
rdBase.DisableLog('rdApp.error')


def cut_point(parent):
    m = random.randint(0, len(parent) - 1)
    return m


def mol2string(mol):
    Chem.Kekulize(mol, clearAromaticFlags=True)
    smiles = Chem.MolToSmiles(mol, canonical=False)
    return encoder(smiles).split('][')


def string2mol(string):
    string = ']['.join(string)
    try:
        smiles = decoder(string,PrintErrorMessage=False)
    except:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol
    except:
        return None


def crossover(parent_a_mol, parent_b_mol):
    parent_a, parent_b = mol2string(parent_a_mol), mol2string(parent_b_mol)

    for _ in range(50):
        cut_point_a = cut_point(parent_a)
        cut_point_b = cut_point(parent_b)
        a1 = parent_a[0:cut_point_a]
        b2 = parent_b[cut_point_b:len(parent_b)]
        child_string = a1 + b2
        child_mol = string2mol(child_string)
        if mol_OK(child_mol):
            return child_mol

    return None
