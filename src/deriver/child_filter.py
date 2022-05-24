from rdkit import Chem
from rdkit.Chem import BRICS as Brics
from rdkit.Chem import Crippen as crip
from rdkit.Chem import Lipinski as lip
from rdkit.Chem import rdMolDescriptors as desc
from rdkit.Chem import rdmolops
from rdkit.Chem import rdqueries
from typing import List
import pandas as pd
import os
import copy

alerts = pd.read_csv(
    os.path.join(os.path.dirname(__file__), "data/alert_collection.csv")
)
pains = alerts.loc[alerts["rule_set_name"] == "PAINS", "smarts"]
glaxo = alerts.loc[alerts["rule_set_name"] == "Glaxo", "smarts"]
always_filter = {"PAINS": pains, "glaxo": glaxo}


def get_filter_values(mol):
    """
    calculate the values, for a given molecule, that are used to filter
    return as a dictionary
    """

    assert isinstance(mol, Chem.Mol)

    atoms = mol.GetAtoms()

    values = {}
    values["MW"] = desc.CalcExactMolWt(mol)
    values["logP"] = crip.MolLogP(mol)
    values["HBA"] = lip.NumHAcceptors(mol)
    values["HBD"] = lip.NumHDonors(mol)
    values["tPSA"] = desc.CalcTPSA(mol)
    values["rot_bonds"] = lip.NumRotatableBonds(mol)
    # assume mutual exclusion
    values["rigid_bonds"] = mol.GetNumBonds() - values["rot_bonds"]
    values["num_rings"] = lip.RingCount(mol)
    values["max_ring_size"] = max(
        (len(r) for r in mol.GetRingInfo().AtomRings()), default=0
    )
    values["num_carbons"] = len(
        mol.GetAtomsMatchingQuery(rdqueries.AtomNumEqualsQueryAtom(6))
    )
    values["num_heavy_atoms"] = desc.CalcNumHeavyAtoms(mol)
    values["num_hetero_atoms"] = lip.NumHeteroatoms(mol)
    values["charge"] = rdmolops.GetFormalCharge(mol)
    values["num_charges"] = sum((abs(atom.GetFormalCharge()) for atom in atoms))
    try:
        values["hc_ratio"] = float(values["num_hetero_atoms"]) / float(
            values["num_carbons"]
        )
    except ZeroDivisionError:
        values["hc_ratio"] = 100000000  # if there are zero carbons
    # how many BRICS bonds, related to complexity
    values["fc"] = len(list(Brics.FindBRICSBonds(mol)))
    # get all the atoms, and make the list unique (only types)
    values["atoms"] = list({atom.GetSymbol() for atom in atoms})
    values["num_chiral_centers"] = len(
        Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    )
    values["fSP3"] = desc.CalcFractionCSP3(mol)
    # default to true, but not yet observed
    values["is_good"] = True
    values["rejections"] = []  # empty list to store the reasons for rejection

    return values


def apply_filter(
    final_limits,
    children,
    must_have_patterns: List[str] = None,
    must_not_have_patterns: List[str] = None,
):

    """
    Apply the filters to a set of molecules.
    :param final_limits: The drug-like parameters to use as bounds for the filters.
    :param children: The Chem.Mol objects to actually filter.
    :param must_have_patterns: An optional list of SMARTS patterns to ensure at least one of which must be present.
    :param must_not_have_patterns: An optional list of SMARTS patterns where none can be present.
    :return: A dictionary with top-level keys being SMILES strings, and next level being filter criteria.
    """

    assert isinstance(children, list)
    assert isinstance(children[0], Chem.Mol)
    assert isinstance(final_limits, dict) or (final_limits is None)
    if not must_not_have_patterns:
        must_not_have_patterns = []

    # compute the values to check for all the children
    child_values = {}
    for child in children:
        # TODO: this is wasteful when it's already precalculated?
        smile = Chem.MolToSmiles(child, isomericSmiles=True)
        # avoid filtering the same smile multiple times
        if smile in child_values:
            pass
        else:
            if final_limits is None:
                child_values[smile] = {"is_good": True}
            else:
                child_values[smile] = get_filter_values(child)
    if final_limits is None:
        return child_values

    # actually filter the children here
    # record each reason why a child violates the filter rules
    for smile in child_values:
        # ranges
        if (
            child_values[smile]["MW"] <= final_limits["MW"][0]
            or child_values[smile]["MW"] >= final_limits["MW"][1]
        ):
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("MW")
        if (
            child_values[smile]["num_heavy_atoms"] <= final_limits["num_heavy_atoms"][0]
            or child_values[smile]["num_heavy_atoms"]
            >= final_limits["num_heavy_atoms"][1]
        ):
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("num_heavy_atoms")
        if (
            child_values[smile]["num_carbons"] <= final_limits["num_carbons"][0]
            or child_values[smile]["num_carbons"] >= final_limits["num_carbons"][1]
        ):
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("num_carbons")
        if (
            child_values[smile]["num_hetero_atoms"]
            <= final_limits["num_hetero_atoms"][0]
        ) or (
            child_values[smile]["num_hetero_atoms"]
            >= final_limits["num_hetero_atoms"][1]
        ):
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("num_hetero_atoms")
        if (
            child_values[smile]["hc_ratio"] <= final_limits["hc_ratio"][0]
            or child_values[smile]["hc_ratio"] >= final_limits["hc_ratio"][1]
        ):
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("hc_ratio")
        if (
            child_values[smile]["charge"] <= final_limits["charge"][0]
            or child_values[smile]["charge"] >= final_limits["charge"][1]
        ):
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("charge")
        if (
            child_values[smile]["logP"] <= final_limits["logP"][0]
            or child_values[smile]["logP"] >= final_limits["logP"][1]
        ):
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("logP")
        if (
            child_values[smile]["fSP3"] <= final_limits["fSP3"][0]
            or child_values[smile]["fSP3"] >= final_limits["fSP3"][1]
        ):
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("fSP3")

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
        if (
            child_values[smile]["num_chiral_centers"]
            >= final_limits["num_chiral_centers"]
        ):
            child_values[smile]["is_good"] = False
            child_values[smile]["rejections"].append("num_chiral_centers")

        # exclude smarts
        smarts_dict = copy.deepcopy(always_filter)
        smarts_dict.update({"user_excluded": must_not_have_patterns})
        for name in smarts_dict.keys():
            if filter_smarts(smile=smile, patterns=smarts_dict[name]):
                child_values[smile]["is_good"] = False
                child_values[smile]["rejections"].append(name)

        # require smarts
        if must_have_patterns:
            if not filter_smarts(smile=smile, patterns=must_have_patterns):
                child_values[smile]["is_good"] = False
                child_values[smile]["rejections"].append("missing_pattern")

    return child_values


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
