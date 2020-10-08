"""
Written by Jan H. Jensen 2018, copied 04/2020
"""
from rdkit.Chem import AllChem

import random
import numpy as np
from .jensen_crossover import ring_OK, mol_OK

from rdkit import Chem, rdBase
rdBase.DisableLog('rdApp.error')


def delete_atom():
    choices = ['[*:1]~[D1]>>[*:1]', '[*:1]~[D2]~[*:2]>>[*:1]-[*:2]',
               '[*:1]~[D3](~[*;!H0:2])~[*:3]>>[*:1]-[*:2]-[*:3]',
               '[*:1]~[D4](~[*;!H0:2])(~[*;!H0:3])~[*:4]>>[*:1]-[*:2]-[*:3]-[*:4]',
               '[*:1]~[D4](~[*;!H0;!H1:2])(~[*:3])~[*:4]>>[*:1]-[*:2](-[*:3])-[*:4]']
    p = [0.25, 0.25, 0.25, 0.1875, 0.0625]
  
    return np.random.choice(choices, p=p)


def append_atom():
    choices = [['single',['C','N','O','F','S','Cl','Br'],7*[1.0/7.0]],
               ['double',['C','N','O'],3*[1.0/3.0]],
               ['triple',['C','N'],2*[1.0/2.0]] ]
    p_BO = [0.60, 0.35, 0.05]

    index = np.random.choice(list(range(3)), p=p_BO)

    BO, atom_list, p = choices[index]
    new_atom = np.random.choice(atom_list, p=p)

    if BO == 'single':
      rxn_smarts = '[*;!H0:1]>>[*:1]X'.replace('X','-'+new_atom)
    if BO == 'double':
      rxn_smarts = '[*;!H0;!H1:1]>>[*:1]X'.replace('X','='+new_atom)
    if BO == 'triple':
      rxn_smarts = '[*;H3:1]>>[*:1]X'.replace('X','#'+new_atom)

    return rxn_smarts


def insert_atom():
    choices = [['single',['C','N','O','S'],4*[1.0/4.0]],
               ['double',['C','N'],2*[1.0/2.0]],
               ['triple',['C'],[1.0]] ]
    p_BO = [0.60,0.35,0.05]

    index = np.random.choice(list(range(3)), p=p_BO)

    BO, atom_list, p = choices[index]
    new_atom = np.random.choice(atom_list, p=p)

    if BO == 'single':
      rxn_smarts = '[*:1]~[*:2]>>[*:1]X[*:2]'.replace('X',new_atom)
    if BO == 'double':
      rxn_smarts = '[*;!H0:1]~[*:2]>>[*:1]=X-[*:2]'.replace('X',new_atom)
    if BO == 'triple':
      rxn_smarts = '[*;!R;!H1;!H0:1]~[*:2]>>[*:1]#X-[*:2]'.replace('X',new_atom)

    return rxn_smarts


def change_bond_order():
    choices = ['[*:1]!-[*:2]>>[*:1]-[*:2]','[*;!H0:1]-[*;!H0:2]>>[*:1]=[*:2]',
               '[*:1]#[*:2]>>[*:1]=[*:2]','[*;!R;!H1;!H0:1]~[*:2]>>[*:1]#[*:2]']
    p = [0.45,0.45,0.05,0.05]

    return np.random.choice(choices, p=p)


def delete_cyclic_bond():
    return '[*:1]@[*:2]>>([*:1].[*:2])'


def add_ring():
    choices = ['[*;!r;!H0:1]~[*;!r:2]~[*;!r;!H0:3]>>[*:1]1~[*:2]~[*:3]1',
               '[*;!r;!H0:1]~[*!r:2]~[*!r:3]~[*;!r;!H0:4]>>[*:1]1~[*:2]~[*:3]~[*:4]1',
               '[*;!r;!H0:1]~[*!r:2]~[*:3]~[*:4]~[*;!r;!H0:5]>>[*:1]1~[*:2]~[*:3]~[*:4]~[*:5]1',
               '[*;!r;!H0:1]~[*!r:2]~[*:3]~[*:4]~[*!r:5]~[*;!r;!H0:6]>>[*:1]1~[*:2]~[*:3]~[*:4]~[*:5]~[*:6]1']
    p = [0.05,0.05,0.45,0.45]

    return np.random.choice(choices, p=p)


def change_atom(mol):
    choices = ['#6','#7','#8','#9','#16','#17','#35']
    random.shuffle(choices)
    X = ""
    for choice in choices:
        if choice in Chem.MolToSmarts(mol):
            X = choice
            break
    random.shuffle(choices)
    Y = ""
    for choice in choices:
        if choice != X:
            Y = choice
            break

    if Y == "" or X == "":
        return None
    return '[X:1]>>[Y:1]'.replace('X',X).replace('Y',Y)


def reactor(mol):

    ia = lambda x: insert_atom()
    cbo = lambda x: change_bond_order()
    dcb = lambda x: delete_cyclic_bond()
    ar = lambda x: add_ring()
    da = lambda x: delete_atom()
    ca = lambda x: change_atom(x)
    aa = lambda x: append_atom()

    rxn_func_list = [ia, cbo, dcb, ar, da, ca, aa]
    random.shuffle(rxn_func_list)

    for rxn_func in rxn_func_list:
        rxn_smarts = []
        for _ in range(10):
            rxn_smarts.append(rxn_func(mol))
        rxn_smarts = list(set(rxn_smarts))
        for smarts in rxn_smarts:
            if smarts is None:
                continue
            rxn = AllChem.ReactionFromSmarts(smarts)
            yield rxn


def mutate(mol):
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=True), sanitize=True)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    my_reactor = reactor(mol)
    for rxn in my_reactor:
        new_mol_trial = rxn.RunReactants((mol,))
        for m in new_mol_trial:
            m = m[0]
            if mol_OK(m) and ring_OK(m):
                return m

    return None
