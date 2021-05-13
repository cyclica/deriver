from rdkit import Chem
from selfies import encoder, decoder, get_semantic_robust_alphabet
import re
import random
from .config import logger

ALPHABET = [symbol.strip("[]") for symbol in get_semantic_robust_alphabet()]

# follow these rules for which types of substitution should be allowed
ALLOWED_SUBS = {
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
    "S": ["O", "N", "C", "=O", "=N", "H"],
    "=S": ["=O", "=N", "=C", "O"],
    "#C": ["#N"]
}


def selfies_substitution(*,
                         parent_smiles: str,
                         n_children: int = 100,
                         mut_rate: float = 0.03,
                         mut_min: int = 1,
                         mut_max: int = 2,
                         max_trials: int = 100
                         ):

    """
    This function takes a parent molecule, and generates a number of children derived from it, via converting
    the parent into SELFIES format and substituting symbols for other symbols, based on primitive chemical rules.
    :param parent_smiles: The smiles of the parent molecule.
    :param n_children: How many children to produce.
    :param mut_rate: How frequent should mutations be? 0.0 to 1.0
    :param mut_min: What is the min number of mutations to allow in a given child, relative to the parent?
    :param mut_max: Same as above but the max.
    :param max_trials: number of attempts to create valid mutations before moving on, useful for pathological selfies
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

        t = 0
        while (muts < mut_min) and (t <= max_trials):
            random.shuffle(mut_positions)  # shuffle the order so that mutations will be random
            for pos in mut_positions:  # try to mutate
                if pos not in mutations:
                    if random.random() <= mut_rate:
                        # ignore special symbols, leave them alone
                        if not (symbols[pos] == "epsilon" or "Branch" in symbols[pos] or "Ring" in symbols[pos]):
                            if symbols[pos] in ALLOWED_SUBS:  # if we have rules for how to handle this symbol
                                mutations.append(pos)  # record which symbol this is
                                muts += 1  # record the intention to mutate
                                if muts == mut_max:
                                    # when we're done, stop looking
                                    break
            t += 1
        if t > max_trials:
            logger.warning(f"Failed to produce any selfies after {max_trials} trials. Returning empty list.")
            return list()

        # for each planned mutation, actually execute it, on the non-shuffled original SELFIES list
        for index in sorted(mutations, reverse=True):
            mut_symbols[index] = random.choice(ALLOWED_SUBS[mut_symbols[index]])

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


def selfies_insertion(*,
                         parent_smiles: str,
                         n_children: int = 100,
                         mut_rate: float = 0.03,
                         mut_min: int = 1,
                         mut_max: int = 2,
                         max_trials: int = 100
                         ):

    """
    This function takes a parent molecule, and generates a number of children derived from it, via converting
    the parent into SELFIES format and inserting other symbols.
    :param parent_smiles: The smiles of the parent molecule.
    :param n_children: How many children to produce.
    :param mut_rate: How frequent should mutations be? 0.0 to 1.0
    :param mut_min: What is the min number of mutations to allow in a given child, relative to the parent?
    :param mut_max: Same as above but the max.
    :param max_trials: number of attempts to create valid mutations before moving on, useful for pathological selfies
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

        t = 0
        while (muts < mut_min) and (t <= max_trials):
            random.shuffle(mut_positions)  # shuffle the order so that mutations will be random
            for pos in mut_positions:  # try to mutate
                if pos not in mutations:
                    if random.random() <= mut_rate:
                        mutations.append(pos)  # record which symbol this is
                        muts += 1  # record the intention to mutate
                        if muts == mut_max:
                            # when we're done, stop looking
                            break

            t += 1
        if t > max_trials:
            logger.warning(f"Failed to produce any selfies after {max_trials} trials. Returning empty list.")
            return list()

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
                         mut_max: int = 2,
                         max_trials: int = 100
                         ):

    """
    This function takes a parent molecule, and generates a number of children derived from it, via converting
    the parent into SELFIES format and deleting symbols.
    :param parent_smiles: The smiles of the parent molecule.
    :param n_children: How many children to produce.
    :param mut_rate: How frequent should mutations be? 0.0 to 1.0
    :param mut_min: What is the min number of mutations to allow in a given child, relative to the parent?
    :param mut_max: Same as above but the max.
    :param max_trials: number of attempts to create valid mutations before moving on, useful for pathological selfies
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

        t = 0
        while (muts < mut_min) and (t <= max_trials):
            random.shuffle(mut_positions)  # shuffle the order so that mutations will be random
            for pos in mut_positions:  # try to mutate
                if pos not in mutations:
                    if random.random() <= mut_rate:
                        mutations.append(pos)  # record which symbol this is
                        muts += 1  # record the intention to mutate
                        if muts == mut_max:
                            # when we're done, stop looking
                            break

            t += 1
        if t > max_trials:
            logger.warning(f"Failed to produce any selfies after {max_trials} trials. Returning empty list.")
            return list()

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


def random_selfies_generator(*, n_symbols: int = 100):

    rand_gen_alphabet = (list(ALLOWED_SUBS.keys()) * 2) + (["C"] * 50) + ALPHABET

    while True:
        symbols = [random.choice(rand_gen_alphabet) for i in range(n_symbols)]
        # print(symbols)
        # convert the new child into a SELFIES (rather than a list)
        child_symbols = [f"[{symb}]" for symb in symbols]
        child = "".join(child_symbols)
        child_smiles = decoder(child)  # get the smiles
        if len(child_smiles) < 10:
            continue
        elif Chem.MolFromSmiles(child_smiles) is None:
            continue
        else:
            yield child_smiles


def selfies_scanner(*, parent_smiles: str, safe_mode: bool = False):
    # get the parent mol set up properly with defined aromaticity
    parent_mol = Chem.MolFromSmiles(parent_smiles, sanitize=True)
    Chem.rdmolops.Kekulize(parent_mol)
    parent_smiles = Chem.MolToSmiles(parent_mol, isomericSmiles=True, kekuleSmiles=True)
    if safe_mode:
        if parent_mol.GetNumAtoms() > 100:
            logger.warning(f"NOT making children from: {parent_smiles}, safe mode on, too many atoms.")
            return []
    logger.info(f"Generating children from: {parent_smiles}")

    children = []  # finished children
    spawns = 0  # counter for children produced
    parent_selfies = encoder(parent_smiles)
    symbols = re.findall(r"[^[]*\[([^]]*)\]", parent_selfies)  # get the SELFIES symbols into a list

    for i, symb in enumerate(symbols):
        if not (symb == "epsilon" or "Branch" in symb or "Ring" in symb):
            if symb in ALLOWED_SUBS:  # if we have rules for how to handle this symbol
                for replacement in ALLOWED_SUBS[symb]:
                    mut_symbols = symbols.copy()  # don't manipulate the original
                    mut_symbols[i] = replacement
                    child_symbols = [f"[{symb}]" for symb in mut_symbols]
                    child = "".join(child_symbols)
                    child_smiles = decoder(child)  # get the smiles
                    # test that smiles is valid
                    try:
                        # same as parent, have to have explicit aromaticity
                        child_mol = Chem.MolFromSmiles(child_smiles, sanitize=True)
                        if child_mol is None:  # if MolToSmiles fails, it will be a None
                            continue
                        Chem.rdmolops.Kekulize(child_mol)
                        child_smiles = Chem.MolToSmiles(child_mol, isomericSmiles=True, kekuleSmiles=True)
                        if child_smiles == parent_smiles:  # ignore this child if it's the same as the parent
                            continue
                    except Exception as e:  # pylint: disable=broad-except
                        print(e)
                        # logger.warning(f"Produced improper SELFIES. Ignoring and trying again. Details below:")
                        # logger.warning(f"Child SELFIES: {child}")
                        # logger.warning(f"Parent SELFIES:{parent_selfies}")
                        # logger.warning(f"Child SMILES: {child_smiles}")
                        # logger.warning(f"Parent SMILES: {parent_smiles}")
                        continue
                    # Every good child deserves fudge
                    children.append(child_smiles)
                    spawns += 1  # update our counter

    return children
