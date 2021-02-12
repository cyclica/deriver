""""Fragment submodule, for fragmenting and calculating fragmentation coefficients."""

import time
from argparse import ArgumentParser
from collections import deque
from playhouse.sqlite_ext import SqliteExtDatabase
from .lib_read import lib_read
from peewee import chunked
from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import BRICS
from .config import logger
from copy import deepcopy
import re
from .clean_frag_db import clean

# set up rdkit logger to ignore WARNINGS (it was so irritating (and possibly responsible for a slight slowdown))
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)

# these atoms are allowed in the fragments, other atoms are not drug-like enough for our purposes
# TODO: could be made into an argument/parameter to be user-provided
allowed_atoms = [
    "C", "N", "S", "O", "P",
    "F", "Br", "Cl", "Na", "K",
    "Mg", "Ca", "B", "Si", "I",
    "Fe", "H", "*"
]


# count atoms (a thread-safe counter that can provide unique IDs) (was more important when multithreaded)
class Counter(object):
    def __init__(self, counter=0):
        self.val = counter

    def increment(self, n=1):
        val = deepcopy(self.val)
        self.val = val + n
        return val


ID_COUNTER = Counter()


class RecomposerMol(object):
    """
    This is a graph-based description of a molecule, which allows for better handling of stereochemistry when
    breaking and reattaching bonds. Otherwise, stereochemistry gets pretty junked up. This also permits fragments
    to be built 'up' from the smallest pieces, producing more fragments, faster, than using the builtin rdkit
    methods. This function based off draft code developed by Harsh Patel @harsh6646
    """
    atomCounter = None
    pseudoIndex = {}  # where are the pseudoatoms if any
    atoms = {}  # same but for atoms
    fragments = {}  # here are all the fragments, keys are FC values
    init_mol = None  # rdkit mol object of this molecule
    edges = []  # basically where are all the bonds
    fc0_index = {}  # given a specific pseudoatom, which fc0 fragment has it
    frag_cache = []  # previously seen fragments (no need to build up from them twice) (can't be a set, dicts inside)
    brics_legend = {}  # which BRICS pseudoatom type is each pseudoatom, needed to rebuild at the end

    def __eq__(self, item):
        if isinstance(item, RecomposerMol):
            return self.pseudoIndex == item.pseudoIndex
        else:
            return self.pseudoIndex == item

    def __hash__(self):
        return self.pseudoIndex

    # track atoms
    class Atom():
        index = None
        atom_num = None
        isotope = None
        chiralTag = None
        aromaticTag = None
        bonds = None
        explicitHs = None

        def __init__(self, *, index:int, atom_num:int, isotope:int,
                     chiralTag:Chem.rdchem.ChiralType, explicitHs:int, bonds:[],
                     aromaticTag: bool):
            self.index = index
            self.atom_num = atom_num
            self.isotope = isotope
            self.chiralTag = chiralTag
            self.explicitHs = explicitHs
            self.aromaticTag = aromaticTag
            self.bonds = bonds

    # track bonds
    class Bond():
        bondType = None
        beginAtomIdx = None
        endAtomIdx = None
        # may want to add bond dir and stereo if needed

        def __init__(self, bondType, beginAtomIdx, endAtomIdx):
            self.bondType = bondType
            self.beginAtomIdx = beginAtomIdx
            self.endAtomIdx = endAtomIdx

    @staticmethod
    def fromMol(mol):
        """
        Build a RecomposerMol object from an rdkit Mol object
        :param mol: rdkit Mol object
        :return: a RecomposerMol object
        """
        newMol = RecomposerMol()
        atoms = {}
        pseudoIndex = {}
        # add info for each atom in the mol into the atom dict
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            bonds = atom.GetBonds()
            aromaticTag = atom.GetIsAromatic()
            newBonds = [newMol.Bond(bond.GetBondType(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in bonds]
            newAtom = newMol.Atom(index=idx, atom_num=atom.GetAtomicNum(),
                                         isotope=atom.GetIsotope(), chiralTag=atom.GetChiralTag(),
                                         explicitHs=atom.GetNumExplicitHs(), bonds=newBonds,
                                         aromaticTag=aromaticTag)
            atoms[idx] = newAtom
            if atom.GetAtomicNum() == 0:  # if it's a pseudoatom
                pseudoIndex[atom.GetIsotope()] = idx
        newMol.init_mol = mol
        newMol.atoms = atoms  # setting these like this prevents reference errors
        newMol.atomCounter = Counter(mol.GetNumAtoms())
        newMol.pseudoIndex = pseudoIndex  # likewise here
        newMol.frag_cache = [{}]  # the parent mol always has an empty pseudoindex
        return newMol

    def addPseudoAtom(self, isotope):
        new_bonds = []
        # get a counter to use as an id
        idx = self.atomCounter.increment()
        new_atom = self.Atom(index=idx, atom_num=0, isotope=isotope,
                            chiralTag=Chem.rdchem.ChiralType.CHI_UNSPECIFIED, explicitHs=0,
                            bonds=new_bonds, aromaticTag=False)
        self.atoms[idx] = new_atom
        self.pseudoIndex[isotope] = idx  # track in order to replace at the end
        return idx

    def breakAndReplace(self, beginIdx, endIdx, beginValue, endValue):
        """
        Replace a bond with two bonds to pseudoatoms. Edit in place.
        :param beginIdx: start atom
        :param endIdx:  end atom
        :param beginValue:  pseudoatom type for start atom
        :param endValue: pseudoatom type for end atom
        :return:
        """
        # just edit in place
        begin_prime = self.addPseudoAtom(beginValue)
        end_prime = self.addPseudoAtom(endValue)

        # replace first bond to bond with pseudoatom
        for i, bond in enumerate(self.atoms[beginIdx].bonds):
            # what follows is a bunch of checks to detect how the rdkit object bond and atom order matches the desired
            if (bond.beginAtomIdx, bond.endAtomIdx) == (beginIdx, endIdx):
                self.atoms[beginIdx].bonds[i].endAtomIdx = begin_prime
                new_bonds = [self.Bond(bondType=bond.bondType, beginAtomIdx=beginIdx, endAtomIdx=begin_prime)]
                self.atoms[begin_prime].bonds = new_bonds
                break
            elif (bond.beginAtomIdx, bond.endAtomIdx) == (endIdx, beginIdx):
                self.atoms[beginIdx].bonds[i].beginAtomIdx = begin_prime
                new_bonds = [self.Bond(bondType=bond.bondType, beginAtomIdx=begin_prime, endAtomIdx=beginIdx)]
                self.atoms[begin_prime].bonds = new_bonds
                break
        # replace second bond with bond with pseudoatom
        for i, bond in enumerate(self.atoms[endIdx].bonds):
            if (bond.beginAtomIdx, bond.endAtomIdx) == (beginIdx, endIdx):
                self.atoms[endIdx].bonds[i].beginAtomIdx = end_prime
                new_bonds = [self.Bond(bondType=bond.bondType, beginAtomIdx=end_prime, endAtomIdx=endIdx)]
                self.atoms[end_prime].bonds = new_bonds
                break
            elif (bond.beginAtomIdx, bond.endAtomIdx) == (endIdx, beginIdx):
                self.atoms[endIdx].bonds[i].endAtomIdx = end_prime
                new_bonds = [self.Bond(bondType=bond.bondType, beginAtomIdx=endIdx, endAtomIdx=end_prime)]
                self.atoms[end_prime].bonds = new_bonds
                break

    def get_fc0_frags(self):
        """
        This function breaks every bond and enumerates all the FC0 fragments. This is a needed first step, since these
        are used to build up all the more complex fragments.
        :return:
        """
        fc0_frags = []
        brics = BRICS.FindBRICSBonds(self.init_mol)
        brics_counter = Counter(counter=100)  # counter starts at 100 since it's not going to be in use as a p.atom
        edges = []
        for bond, brics_type in brics:
            i, j = bond
            x, y = brics_type
            pseudo_i = brics_counter.increment()
            pseudo_j = brics_counter.increment()  # this is always done in pairs so that eg. 100 and 101 go together
            self.brics_legend[pseudo_i] = int(x)  # remember to store the BRICS typing, we need it later
            self.brics_legend[pseudo_j] = int(y)
            edges.append((pseudo_i, pseudo_j))  # keep a record of this bond
            self.breakAndReplace(i, j, pseudo_i, pseudo_j)
        self.edges = edges

        # this while loop finds parts of a graph that are no longer connected, meaning FC0 fragments in this case
        untouched = set(self.atoms.keys())  # keep a list of atoms we haven't investigated yet
        while untouched:  # if there are still atoms we haven't processed
            currSet = set()
            currSet.add(untouched.pop())
            newMol = RecomposerMol()
            newAtoms = {}
            newPseudoAtoms = {}
            new_brics_legend = {}
            new_brics_legend.update(self.brics_legend)
            while currSet:
                # perform graph search
                currAtom = self.atoms[currSet.pop()]
                for bond in currAtom.bonds:
                    idx = None
                    if currAtom.index == bond.endAtomIdx:
                        idx = bond.beginAtomIdx
                    elif currAtom.index == bond.beginAtomIdx:
                        idx = bond.endAtomIdx
                    # add connected atoms to curr mol atom set
                    if idx in untouched:
                        currSet.add(idx)
                        untouched.remove(idx)
                    else:
                        continue
                newAtoms[currAtom.index] = self._copyAtom(currAtom)  # watch out for reference errors!
                if currAtom.atom_num == 0:
                    newPseudoAtoms[currAtom.isotope] = currAtom.index
            newMol.atoms = newAtoms
            newMol.pseudoIndex = newPseudoAtoms
            newMol.brics_legend = new_brics_legend
            fc0_frags.append(newMol)
            self.frag_cache.append(newMol)
        # now need to index the fragments to make it easier to bond them
        fc0_index = {}
        for i, frag in enumerate(fc0_frags):
            for pseudo_atom in frag.pseudoIndex.keys():
                fc0_index[pseudo_atom] = i  # each FC0 fragment should be accessible by pseudoatom
        self.fc0_index = fc0_index
        self.fragments[0] = fc0_frags

    def react(self, *, frag1, frag2, p_atom1, p_atom2):
        """
        basically the opposite of breaking the bonds. Remove the pseudoatoms and put a bond there.
        :param frag1: RecomposerMol
        :param frag2: RecomposerMol
        :param p_atom1: the pseudoatom on fragment 1
        :param p_atom2: the pseudoatom on fragment 2
        :return: a new RecomposerMol
        """

        # merge the two mols
        new_mol = RecomposerMol()
        newAtoms = {}
        new_pi = {}
        new_bl = {}
        new_bl.update(frag1.brics_legend)  # the brics legend is so useful, keep it
        for atom_idx in frag1.atoms:
            newAtoms[atom_idx] = frag1._copyAtom(frag1.atoms[atom_idx])
        for atom_idx in frag2.atoms:
            newAtoms[atom_idx] = frag2._copyAtom(frag2.atoms[atom_idx])

        new_mol.atoms = newAtoms
        pi1 = deepcopy(frag1.pseudoIndex)
        new_pi.update(pi1)
        pi2 = deepcopy(frag2.pseudoIndex)
        new_pi.update(pi2)

        # new mol is ready, now bond it and remove the pseudoatoms
        # get the real atoms that are to be attached
        p_a_idx = new_pi[p_atom1]
        p_b_idx = new_pi[p_atom2]
        a_idx = None
        b_idx = None

        # need to check each end of both bonds to figure out which is the right atom to replace
        for bond in new_mol.atoms[p_a_idx].bonds:  # should only ever be one bond for any pseudoatom
            if bond.beginAtomIdx == p_a_idx:
                a_idx = bond.endAtomIdx
            else:
                a_idx = bond.beginAtomIdx

        for bond in new_mol.atoms[p_b_idx].bonds:  # should only ever be one bond for any pseudoatom
            if bond.beginAtomIdx == p_b_idx:
                b_idx = bond.endAtomIdx
            else:
                b_idx = bond.beginAtomIdx

        # now make the new bonds
        # a
        for i, bond in enumerate(new_mol.atoms[a_idx].bonds):
            if (bond.beginAtomIdx, bond.endAtomIdx) == (a_idx, p_a_idx):
                new_mol.atoms[a_idx].bonds[i].endAtomIdx = b_idx
                break
            elif (bond.beginAtomIdx, bond.endAtomIdx) == (p_a_idx, a_idx):
                new_mol.atoms[a_idx].bonds[i].beginAtomIdx = b_idx
                break
        # b
        for i, bond in enumerate(new_mol.atoms[b_idx].bonds):
            if (bond.beginAtomIdx, bond.endAtomIdx) == (b_idx, p_b_idx):
                new_mol.atoms[b_idx].bonds[i].endAtomIdx = a_idx
                break
            elif (bond.beginAtomIdx, bond.endAtomIdx) == (p_b_idx, b_idx):
                new_mol.atoms[b_idx].bonds[i].beginAtomIdx = a_idx
                break

        # now remove the pseudoatoms
        del new_mol.atoms[p_a_idx]
        del new_mol.atoms[p_b_idx]
        del new_pi[p_atom1]
        del new_pi[p_atom2]
        new_mol.pseudoIndex = new_pi
        new_mol.brics_legend = new_bl

        return new_mol

    def get_fc1_frags(self):
        """
        This function gets (and sets) all the FC1 fragments. This is another special case, since you can just combine
        all the pairs of FC0 fragments that match, and that is definitively all the FC1 fragments.
        :return:
        """
        fc1_frags = []
        for p_atom1, p_atom2 in self.edges:  # for every bond (now broken)
            # get the FC0 fragments
            frag1 = self.fragments[0][self.fc0_index[p_atom1]]
            frag2 = self.fragments[0][self.fc0_index[p_atom2]]
            # combine them
            new_mol = self.react(frag1=frag1, frag2=frag2, p_atom1=p_atom1, p_atom2=p_atom2)
            if new_mol:
                if new_mol not in self.frag_cache:
                    fc1_frags.append(new_mol)
                    self.frag_cache.append(new_mol)
        self.fragments[1] = fc1_frags

    def _copyAtom(self, atom):
        newBonds = [self.Bond(bond.bondType, bond.beginAtomIdx, bond.endAtomIdx) for bond in atom.bonds]
        return self.Atom(index=atom.index, atom_num=atom.atom_num, isotope=atom.isotope,
                         chiralTag=atom.chiralTag, explicitHs=atom.explicitHs,
                         bonds=newBonds, aromaticTag=atom.aromaticTag)

    def get_next_fc_level_frags(self, *, fc):
        """
        Given an FC value for which the fragments already exist, compute the next level up of complexity
        :param fc:
        :return: nothing, edit in place
        """
        new_frags = []
        for old_frag in self.fragments[fc]:
            p_atoms = deepcopy(list(old_frag.pseudoIndex.keys()))
            for p_atom_label in p_atoms:
                if p_atom_label % 2 == 0:
                    other_p_atom = p_atom_label+1
                else:
                    other_p_atom = p_atom_label - 1
                frag2 = self.fragments[0][self.fc0_index[other_p_atom]]
                new_mol = self.react(frag1=old_frag, frag2=frag2, p_atom1=p_atom_label, p_atom2=other_p_atom)
                # it is possible to have multiple paths create the same fragment
                # detect here and ignore (this keeps the search space much smaller)
                if new_mol not in self.frag_cache:
                    new_frags.append(new_mol)
                    self.frag_cache.append(new_mol)
        self.fragments[fc+1] = new_frags

    def get_all_fragments(self, max_fc=7):
        """
        Given a target FC value, calculate all the fragments up to and including that FC. Format them in a way that is
        ready to be written to the database.
        :param max_fc:
        :return: frag_list, heritage_list, atoms_list, pseudo_atoms_list (lists of dictionaries matching models)
        """
        self.get_fc0_frags()  # always do this first
        max_fc = min(len(self.edges), max_fc)
        if len(self.edges) > 1:  # if there is more than one brics bond
            self.get_fc1_frags()
        if len(self.edges) > 2:  # if there are more than two brics bonds
            for i in range(1, max_fc):
                self.get_next_fc_level_frags(fc=i)

        # these will store the dictionaries the database code expects
        frag_list = []
        heritage_list = []
        atoms_list = []
        pseudo_atoms_list = []

        # using this function to avoid copying code
        def insert_mol_info(mol, fc, curr_id, parent_id):

            seen_ps = set()
            seen_atom = set()
            num_pseudo_atoms = 0
            num_unique_pseudo_atoms = 0
            has_only_good_atoms = True

            # sanity checking to see if this is a fragment we want to use for making molecules
            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                if symbol not in allowed_atoms:
                    has_only_good_atoms = False
                    break
                atomic_num = atom.GetAtomicNum()
                ps_num = atom.GetIsotope()

                # storing info about each atom
                if symbol not in seen_atom and atomic_num != 0:
                    atoms_list.append({
                        'frag': curr_id,
                        'atom': symbol
                    })
                    seen_atom.add(symbol)
                if atomic_num == 0:
                    num_pseudo_atoms += 1
                    pseudo_atoms_list.append({
                        'frag': curr_id,
                        'pseudo_atom': ps_num
                    })
                    if ps_num not in seen_ps:
                        num_unique_pseudo_atoms += 1
                        seen_ps.add(ps_num)

            # yeah we'll take it
            if has_only_good_atoms:
                frag_list.append({
                    'id': curr_id,
                    'smile': Chem.MolToSmiles(mol, isomericSmiles=True),
                    'frag_coeff': fc,
                    'num_pseudo_atoms': num_pseudo_atoms,
                    'num_unique_pseudo_atoms': num_unique_pseudo_atoms
                })
                heritage_list.append({
                    'frag': curr_id,
                    'parent': parent_id
                })

        # insert parent first, since the fragments refer to the parent_id
        parent_id = ID_COUNTER.increment()
        curr_id = parent_id
        insert_mol_info(self.init_mol, len(self.edges), curr_id, parent_id)  # parents have a frag_id == parent_id
        curr_id = ID_COUNTER.increment()
        # insert all sub fragments
        for fc in range(max_fc):
            frag_group = self.fragments[fc]
            for frag in frag_group:
                insert_mol_info(frag.toMol(), fc, curr_id, parent_id)
                curr_id = ID_COUNTER.increment()
        return frag_list, heritage_list, atoms_list, pseudo_atoms_list

    def toMol(self):
        """
        Turn a RecomposerMol object back into an RDKIT Chem.Mol object, with correct stereochemistry
        and pseudoatoms if relevant
        :return: Chem.Mol
        """
        emol = Chem.RWMol(Chem.MolFromSmiles(""))  # create an empty editable mol
        atomsIndices = sorted(list(self.atoms.keys()))
        idxMapping = {}
        # add all the atoms back
        for i, idx in enumerate(atomsIndices):
            atom = self.atoms[idx]
            newAtom = Chem.Atom(atom.atom_num)
            if atom.atom_num == 0:
                newAtom.SetIsotope(self.brics_legend[atom.isotope])
            else:
                newAtom.SetIsotope(atom.isotope)
            newAtom.SetNoImplicit(True)
            newAtom.SetNumExplicitHs(atom.explicitHs)
            newAtom.SetChiralTag(atom.chiralTag)
            newAtom.SetIsAromatic(atom.aromaticTag)
            emol.AddAtom(newAtom)
            idxMapping[idx] = i
        # add all the bonds back
        for idx in atomsIndices:
            for bond in self.atoms[idx].bonds:
                begin = bond.beginAtomIdx
                end = bond.endAtomIdx
                if idx == begin:
                    emol.AddBond(idxMapping[begin], idxMapping[end], bond.bondType)
        return emol.GetMol()


def libgen(mol_list, output_name):
    """
    function to generate a database format library of fragments from a mol, list of mol objects, .smi, or .sdf file
    :param mol_list: list of molecules, a single molecule, or a filename of molecules to read
    :type mol_list: str|Chem.Mol|[Chem.Mol]
    :param output_name: name of the database to use?
    :type output_name: str
    :return:
    """
    # if a file not a list then read into list
    if isinstance(mol_list, str) and mol_list.endswith(".smi"):
        mol_list = Chem.SmilesMolSupplier(mol_list, delimiter="\t", titleLine=False)
    elif isinstance(mol_list, str) and mol_list.endswith(".sdf"):
        mol_list = Chem.SDMolSupplier(mol_list)
    elif type(mol_list) == Chem.Mol:
        mol_list = [mol_list]
    elif type(mol_list) == list:
        assert type(mol_list[0]) == Chem.Mol
    else:
        raise Exception("Did you provide a list of mol objects? Input type error.")

    fragment_dict_deque = deque()
    heritage_dict_deque =deque()
    atoms_dict_deque = deque()
    pseudoatoms_dict_deque = deque()
    logger.info("Fragmenting:")
    n = len(mol_list)
    i = 0
    t0 = time.time()
    for mol in mol_list:
        re_mol = RecomposerMol.fromMol(mol=mol)
        frag_list, heritage_list, atoms_list, pseudo_atoms_list = re_mol.get_all_fragments(7)
        fragment_dict_deque.extend(frag_list)
        heritage_dict_deque.extend(heritage_list)
        atoms_dict_deque.extend(atoms_list)
        pseudoatoms_dict_deque.extend(pseudo_atoms_list)
        logger.info("DONE: %d/%d %.f" % (
            i, n, 1000 * (time.time() - t0) / (i + 1)
        ))
        i += 1

    logger.info("Done")
    logger.info("Saving:")

    # create the database for the output
    db = SqliteExtDatabase(output_name, pragmas={
        'cache_size': -1024 * 64,  # 64MB page-cache.
        'journal_mode': 'wal',  # Use WAL-mode (you should always use this!).
        'foreign_keys': 0,
        'wal_autocheckpoint': 10,
    })

    db.connect()
    # get the models
    Fragment, Heritage, PseudoAtoms, Atoms = lib_read(db)

    Fragment.create_table(safe=True)
    Heritage.create_table(safe=True)
    PseudoAtoms.create_table(safe=True)
    Atoms.create_table(safe=True)
    with db.atomic():
        if len(fragment_dict_deque) > 0:
            for ents in chunked(fragment_dict_deque, 200):
                query = Fragment.replace_many(ents)
                query.execute()
            for ents in chunked(heritage_dict_deque, 200):
                query = Heritage.replace_many(ents)
                query.execute()
            for ents in chunked(pseudoatoms_dict_deque, 200):
                query = PseudoAtoms.replace_many(ents)
                query.execute()
            for ents in chunked(atoms_dict_deque, 200):
                query = Atoms.replace_many(ents)
                query.execute()
    db.close()
    clean(output_name)

    return 1


def get_pseudo_atoms(smile):
    """
    function to get BRICS pseudoatoms from a smile, and return the number of pseudoatoms, the num unique, and the atoms
    :param smile:
    :return:
    """
    num_pseudo_atoms = 0
    num_unique_pseudo_atoms = 0
    pseudo_atoms = []
    pattern = re.compile(r'\[([0-9]+)\*\]')
    try:
        pseudo_atoms = pattern.findall(smile)
        num_pseudo_atoms = len(pseudo_atoms)
        num_unique_pseudo_atoms = len(set(pseudo_atoms))
    except Exception:  # pylint: disable=broad-except
        pass

    return num_pseudo_atoms, num_unique_pseudo_atoms, pseudo_atoms


if __name__ == "__main__":
    # one of the only CLIs here. Just for convenience. Calls libgen on an input file. Now automatically cleaned
    parser = ArgumentParser()
    parser.add_argument('-i',
                        dest="input",
                        type=str,
                        help='.smi or .sdf file of molecules to be fragmented.',
                        required=True)
    parser.add_argument("-o",
                        dest="output",
                        type=str,
                        help="Database location for new fragment library.",
                        required=True)

    options = parser.parse_args()
    runtime = libgen(options.input, options.output)
    logger.info(runtime)
