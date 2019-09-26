from rdkit import Chem
from selfies import encoder, decoder
import re
import random
import logging as logger


ALPHABET = ['Branch1_1', 'Branch1_2', 'Branch1_3', 'Ring1', 'Branch2_1', 'Branch2_2', 'Branch2_3',
            'Ring2', 'Branch3_1', 'Branch3_2', 'Branch3_3', 'O', '=O', 'N', '=N', 'C', '=C',
            '#C', 'S', '=S', 'P', 'F', 'C@Hexpl', 'C@@Hexpl', 'C@expl', 'C@@expl', 'H']


def selfies_substitution(*,
                         parent_smiles: str,
                         n_children: int = 100,
                         mut_rate: float = 0.03,
                         mut_min: int = 1,
                         mut_max: int = 2
                         ):

    """
    This function takes a parent molecule, and generates a number of children derived from it, via converting
    the parent into SELFIES format and substituting symbols for other symbols, based on primitive chemical rules.
    :param parent_smiles: The smiles of the parent molecule.
    :param n_children: How many children to produce.
    :param mut_rate: How frequent should mutations be? 0.0 to 1.0
    :param mut_min: What is the min number of mutations to allow in a given child, relative to the parent?
    :param mut_max: Same as above but the max.
    :return:
    """

    # get the parent mol set up properly with defined aromaticity
    parent_mol = Chem.MolFromSmiles(parent_smiles, sanitize=True)
    Chem.rdmolops.Kekulize(parent_mol)
    parent_smiles = Chem.MolToSmiles(parent_mol, isomericSmiles=True, kekuleSmiles=True)
    logger.info(f"Generating children from: {parent_smiles}")

    children = []  # finished children
    spawns = 0  # counter for children produced
    parent_selfies = encoder(parent_smiles)
    symbols = re.findall(r"[^[]*\[([^]]*)\]", parent_selfies)   # get the SELFIES symbols into a list

    # follow these rules for which types of substitution should be allowed
    allowed_subs = {
        "C": ["N", "O", "H"],
        "=C": ["=N", "N", "O", "=O", "S"],
        "F": ["Cl", "Br", "I", "O", "N", "C", "H"],
        "Cl": ["F", "Br", "I", "O", "N", "C", "H"],
        "Br": ["Cl", "F", "I", "O", "N", "C", "H"],
        "I": ["Cl", "Br", "F", "O", "N", "C", "H"],
        "H": ["C", "O", "N", "S", "=C", "=O", "=S"],
        "O": ["S", "N", "=O", "=N", "C", "=C", "Cl", "F", "Br", "I", "H"],
        "=O": ["=S", "=N", "=C", "O"],
        "N": ["O", "C", "H"],
        "=N": ["=O", "O", "S", "=C", "C"],
        "#N": ["#C"],
        "S": ["O", "N", "C", "=O", "=N", "H"]
    }

    while spawns < n_children:  # try to produce the correct number of children
        muts = 0
        mutations = []  # which parts of the SELFIES to remove
        mut_symbols = symbols.copy()  # don't manipulate the original
        mut_positions = list(range(len(symbols)))  # need the index

        while muts < mut_min:
            random.shuffle(mut_positions)  # shuffle the order so that mutations will be random
            for pos in mut_positions:  # try to mutate
                if random.random() <= mut_rate:
                    # ignore special symbols, leave them alone
                    if not (symbols[pos] == "epsilon" or "Branch" in symbols[pos] or "Ring" in symbols[pos]):
                        if symbols[pos] in allowed_subs:  # if we have rules for how to handle this symbol
                            mutations.append(pos)  # record which symbol this is
                            muts += 1  # record the intention to mutate
                            if muts == mut_max:
                                # when we're done, stop looking
                                break

        # for each planned mutation, actually execute it, on the non-shuffled original SELFIES list
        for index in sorted(mutations, reverse=True):
            mut_symbols[index] = random.choice(allowed_subs[mut_symbols[index]])

        # convert the new child into a SELFIES (rather than a list)
        child_symbols = [f"[{symb}]" for symb in mut_symbols]
        child = "".join(child_symbols)
        child_smiles = decoder(child)  # get the smiles

        # test that smiles is valid
        try:
            # same as parent, have to have explicit aromaticity
            child_mol = Chem.MolFromSmiles(child_smiles, sanitize=True)
            Chem.rdmolops.Kekulize(child_mol)
            child_smiles = Chem.MolToSmiles(child_mol, isomericSmiles=True, kekuleSmiles=True)
            assert child_mol  # if MolToSmiles fails, it will be a None
            if child_smiles == parent_smiles:  # ignore this child if it's the same as the parent
                continue
        except Exception:  # pylint: disable=broad-except
            logger.warning(f"Produced improper SELFIES. Ignoring and trying again. Details below:")
            logger.warning(f"Child SELFIES: {child}")
            logger.warning(f"Parent SELFIES:{parent_selfies}")
            logger.warning(f"Child SMILES: {child_smiles}")
            logger.warning(f"Parent SMILES: {parent_smiles}")
            continue
        # Every good child deserves fudge
        children.append(child_smiles)
        spawns += 1  # update our counter

    return children


def random_selfies_generation(*, n_symbols: int = 50):
    # should be a generator
    yield NotImplementedError


def selfies_insertion(*,
                         parent_smiles: str,
                         n_children: int = 100,
                         mut_rate: float = 0.03,
                         mut_min: int = 1,
                         mut_max: int = 2
                         ):

    """
    This function takes a parent molecule, and generates a number of children derived from it, via converting
    the parent into SELFIES format and inserting other symbols.
    :param parent_smiles: The smiles of the parent molecule.
    :param n_children: How many children to produce.
    :param mut_rate: How frequent should mutations be? 0.0 to 1.0
    :param mut_min: What is the min number of mutations to allow in a given child, relative to the parent?
    :param mut_max: Same as above but the max.
    :return:
    """

    # get the parent mol set up properly with defined aromaticity
    parent_mol = Chem.MolFromSmiles(parent_smiles, sanitize=True)
    Chem.rdmolops.Kekulize(parent_mol)
    parent_smiles = Chem.MolToSmiles(parent_mol, isomericSmiles=True, kekuleSmiles=True)
    logger.info(f"Generating children from: {parent_smiles}")

    children = []  # finished children
    spawns = 0  # counter for children produced
    parent_selfies = encoder(parent_smiles)
    symbols = re.findall(r"[^[]*\[([^]]*)\]", parent_selfies)   # get the SELFIES symbols into a list


    while spawns < n_children:  # try to produce the correct number of children
        muts = 0
        mutations = []  # which parts of the SELFIES to remove
        mut_symbols = symbols.copy()  # don't manipulate the original
        mut_positions = list(range(len(symbols)))  # need the index

        while muts < mut_min:
            random.shuffle(mut_positions)  # shuffle the order so that mutations will be random
            for pos in mut_positions:  # try to mutate
                if random.random() <= mut_rate:
                    mutations.append(pos)  # record which symbol this is
                    muts += 1  # record the intention to mutate
                    if muts == mut_max:
                        # when we're done, stop looking
                        break

        # for each planned mutation, actually execute it, on the non-shuffled original SELFIES list
        for index in sorted(mutations, reverse=True):
            mut_symbols.insert(index, random.choice(ALPHABET))

        # convert the new child into a SELFIES (rather than a list)
        child_symbols = [f"[{symb}]" for symb in mut_symbols]
        child = "".join(child_symbols)
        child_smiles = decoder(child)  # get the smiles

        # test that smiles is valid
        try:
            # same as parent, have to have explicit aromaticity
            child_mol = Chem.MolFromSmiles(child_smiles, sanitize=True)
            Chem.rdmolops.Kekulize(child_mol)
            child_smiles = Chem.MolToSmiles(child_mol, isomericSmiles=True, kekuleSmiles=True)
            assert child_mol  # if MolToSmiles fails, it will be a None
            if child_smiles == parent_smiles:  # ignore this child if it's the same as the parent
                continue
        except Exception:  # pylint: disable=broad-except
            logger.warning(f"Produced improper SELFIES. Ignoring and trying again. Details below:")
            logger.warning(f"Child SELFIES: {child}")
            logger.warning(f"Parent SELFIES:{parent_selfies}")
            logger.warning(f"Child SMILES: {child_smiles}")
            logger.warning(f"Parent SMILES: {parent_smiles}")
            continue
        # Every good child deserves fudge
        children.append(child_smiles)
        spawns += 1  # update our counter

    return children


def selfies_deletion(*,
                         parent_smiles: str,
                         n_children: int = 100,
                         mut_rate: float = 0.03,
                         mut_min: int = 1,
                         mut_max: int = 2
                         ):

    """
    This function takes a parent molecule, and generates a number of children derived from it, via converting
    the parent into SELFIES format and deleting symbols.
    :param parent_smiles: The smiles of the parent molecule.
    :param n_children: How many children to produce.
    :param mut_rate: How frequent should mutations be? 0.0 to 1.0
    :param mut_min: What is the min number of mutations to allow in a given child, relative to the parent?
    :param mut_max: Same as above but the max.
    :return:
    """

    # get the parent mol set up properly with defined aromaticity
    parent_mol = Chem.MolFromSmiles(parent_smiles, sanitize=True)
    Chem.rdmolops.Kekulize(parent_mol)
    parent_smiles = Chem.MolToSmiles(parent_mol, isomericSmiles=True, kekuleSmiles=True)
    logger.info(f"Generating children from: {parent_smiles}")

    children = []  # finished children
    spawns = 0  # counter for children produced
    parent_selfies = encoder(parent_smiles)
    symbols = re.findall(r"[^[]*\[([^]]*)\]", parent_selfies)   # get the SELFIES symbols into a list


    while spawns < n_children:  # try to produce the correct number of children
        muts = 0
        mutations = []  # which parts of the SELFIES to remove
        mut_symbols = symbols.copy()  # don't manipulate the original
        mut_positions = list(range(len(symbols)))  # need the index

        while muts < mut_min:
            random.shuffle(mut_positions)  # shuffle the order so that mutations will be random
            for pos in mut_positions:  # try to mutate
                if random.random() <= mut_rate:
                    mutations.append(pos)  # record which symbol this is
                    muts += 1  # record the intention to mutate
                    if muts == mut_max:
                        # when we're done, stop looking
                        break

        # for each planned mutation, actually execute it, on the non-shuffled original SELFIES list
        for index in sorted(mutations, reverse=True):
            del mut_symbols[index]

        # convert the new child into a SELFIES (rather than a list)
        child_symbols = [f"[{symb}]" for symb in mut_symbols]
        child = "".join(child_symbols)
        child_smiles = decoder(child)  # get the smiles

        # test that smiles is valid
        try:
            # same as parent, have to have explicit aromaticity
            child_mol = Chem.MolFromSmiles(child_smiles, sanitize=True)
            Chem.rdmolops.Kekulize(child_mol)
            child_smiles = Chem.MolToSmiles(child_mol, isomericSmiles=True, kekuleSmiles=True)
            assert child_mol  # if MolToSmiles fails, it will be a None
            if child_smiles == parent_smiles:  # ignore this child if it's the same as the parent
                continue
        except Exception:  # pylint: disable=broad-except
            logger.warning(f"Produced improper SELFIES. Ignoring and trying again. Details below:")
            logger.warning(f"Child SELFIES: {child}")
            logger.warning(f"Parent SELFIES:{parent_selfies}")
            logger.warning(f"Child SMILES: {child_smiles}")
            logger.warning(f"Parent SMILES: {parent_smiles}")
            continue
        # Every good child deserves fudge
        children.append(child_smiles)
        spawns += 1  # update our counter

    return children
