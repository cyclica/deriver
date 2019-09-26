from .config import pains_smarts
from rdkit import Chem
from rdkit.Chem import BRICS as Brics
from rdkit.Chem import Crippen as crip
from rdkit.Chem import Lipinski as lip
from rdkit.Chem import rdMolDescriptors as desc
from rdkit.Chem import rdmolops
from typing import List


def get_atom_props(mol):
    """
    this function takes a mol and calculates the number of carbons, charges, and the size of the largest ring
    """

    num_carbons = 0
    num_charges = 0
    max_ring_size = 0

    # loop over all the atoms
    for at in mol.GetAtoms():
        if at.GetSymbol() == "C":
            num_carbons += 1
        if at.GetFormalCharge() != 0:
            num_charges += abs(at.GetFormalCharge())

        # I suspect that there will never be a ring larger than 100, but this is technically an arbitrary cutoff
        # this really could be optimized, but it's not rate-limiting enough to bother right now
        for size in range(3, 101):
            if at.IsInRingSize(size):
                # which is larger, lazy method
                max_ring_size = max((size, max_ring_size))

    return num_carbons, num_charges, max_ring_size


def get_filter_values(mol):
    """
    calculate the values, for a given molecule, that are used to filter
    return as a dictionary
    """

    assert isinstance(mol, Chem.Mol)

    values = {}
    values["MW"] = desc.CalcExactMolWt(mol)
    values["logP"] = crip.MolLogP(mol)
    values["HBA"] = lip.NumHAcceptors(mol)
    values["HBD"] = lip.NumHDonors(mol)
    values["tPSA"] = desc.CalcTPSA(mol)
    values["rot_bonds"] = lip.NumRotatableBonds(mol)
    values["rigid_bonds"] = mol.GetNumBonds() - values["rot_bonds"]  # assume mutual exclusion
    values["num_rings"] = lip.RingCount(mol)
    values["num_hetero_atoms"] = lip.NumHeteroatoms(mol)
    values["charge"] = rdmolops.GetFormalCharge(mol)  # trusting this charge calculation method
    values["num_carbons"], values["num_charges"], values["max_ring_size"] = get_atom_props(mol)
    try:
        values["hc_ratio"] = float(values["num_hetero_atoms"]) / float(values["num_carbons"])
    except ZeroDivisionError:
        values["hc_ratio"] = 100000000  # if there are zero carbons
    values["fc"] = len(list(Brics.FindBRICSBonds(mol)))  # how many BRICS bonds, related to complexity
    values["is_good"] = True  # default to true, but not yet observed
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]  # get all the atoms, and make the list unique (only types)
    atoms = set(atoms)
    atoms = list(atoms)
    values["atoms"] = atoms
    values["num_chiral_centers"] = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    values["rejections"] = []  # empty list to store the reasons for rejection

    return values


def apply_filter(final_limits, children, must_have_patterns: List[str] = None):

    """
    Apply the filters to a set of molecules.
    :param final_limits: The drug-like parameters to use as bounds for the filters.
    :param children: The Chem.Mol objects to actually filter.
    :param must_have_patterns: An optional list of SMARTS patterns to ensure at least one of which must be present.
    :return: A dictionary with top-level keys being SMILES strings, and next level being filter criteria.
    """

    assert isinstance(children, list)
    assert isinstance(children[0], Chem.Mol)
    assert isinstance(final_limits, dict)

    # compute the values to check for all the children
    child_values = {}
    for child in children:
        smile = Chem.MolToSmiles(child, isomericSmiles=True)
        # avoid filtering the same smile multiple times
        if smile in child_values:
            pass
        else:
            child_values[smile] = get_filter_values(child)

    # actually filter the children here
    # record each reason why a child violates the filter rules
    for smile in child_values:
        # ranges
        if child_values[smile]["MW"] <= final_limits["MW"][0]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("MW")
        if child_values[smile]["MW"] >= final_limits["MW"][1]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("MW")
        if child_values[smile]["num_carbons"] <= final_limits["num_carbons"][0]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("num_carbons")
        if child_values[smile]["num_carbons"] >= final_limits["num_carbons"][1]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("num_carbons")
        if child_values[smile]["num_hetero_atoms"] <= final_limits["num_hetero_atoms"][0]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("num_hetero_atoms")
        if child_values[smile]["num_hetero_atoms"] >= final_limits["num_hetero_atoms"][1]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("num_hetero_atoms")
        if child_values[smile]["hc_ratio"] <= final_limits["hc_ratio"][0]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("hc_ratio")
        if child_values[smile]["hc_ratio"] >= final_limits["hc_ratio"][1]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("hc_ratio")
        if child_values[smile]["charge"] <= final_limits["charge"][0]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("charge")
        if child_values[smile]["charge"] >= final_limits["charge"][1]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("charge")
        if child_values[smile]["logP"] <= final_limits["logP"][0]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("logP")
        if child_values[smile]["logP"] >= final_limits["logP"][1]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("logP")

        # upper limits
        if child_values[smile]["HBA"] >= final_limits["HBA"]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("HBA")
        if child_values[smile]["HBD"] >= final_limits["HBD"]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("HBD")
        if child_values[smile]["tPSA"] >= final_limits["tPSA"]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("tPSA")
        if child_values[smile]["rot_bonds"] >= final_limits["rot_bonds"]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("rot_bonds")
        if child_values[smile]["rigid_bonds"] >= final_limits["rigid_bonds"]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("rigid_bonds")
        if child_values[smile]["num_rings"] >= final_limits["num_rings"]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("num_rings")
        if child_values[smile]["max_ring_size"] >= final_limits["max_ring_size"]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("max_ring_size")
        if child_values[smile]["num_charges"] >= final_limits["num_charges"]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("num_charges")
        if child_values[smile]["num_chiral_centers"] >= final_limits["num_chiral_centers"]:
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("num_chiral_centers")

        # PAINS
        if filter_pains(smile):
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("PAINS")

        # filter smarts patterns
        if must_have_patterns:
            if not filter_smarts(smile=smile, patterns=must_have_patterns):
                child_values[smile]["is_good"] = False
                child_values[smile]["rejections"].append("missing_pattern")

    return child_values


pains_patts = None  # scope this here, embarrasingly. should be improved by a proper representation for the filter


def filter_pains(smile):
    """
    function for filtering based on pains criteria
    :param smile:
    :return:
    """
    # read in the PAINS patterns
    global pains_patts  # pylint: disable=global-statement
    if pains_patts is None:
        pains_patts = [Chem.MolFromSmarts(x, mergeHs=True) for x in pains_smarts]

    # if there is a PAINS match, stop searching and return True
    mol = Chem.MolFromSmiles(smile)
    for patt in pains_patts:
        if mol.HasSubstructMatch(patt):
            return True
    return False


def filter_smarts(*, smile: str, patterns: List[str]):

    """
    Check a molecule for the presence of any of a set of SMARTS patterns.
    :param smile: SMILES string of the molecule to be checked
    :param patterns: list of SMARTS patterns to look for
    :return: True if molecule DOES have AT LEAST ONE of the patterns
    """

    # get mols of the search patterns
    filter_mols = [Chem.MolFromSmarts(x, mergeHs=True) for x in patterns]
    # if there is a match, stop searching and return True
    test_mol = Chem.MolFromSmiles(smile)
    for patt in filter_mols:
        if test_mol.HasSubstructMatch(patt):
            return True
    # else, no match
    return False
