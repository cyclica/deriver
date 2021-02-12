import os

import deriver.fragment as fragment
from .brics_mating_rules import brics_mating_rules
from .child_filter import apply_filter
from .fragment_index import frag_index
from .lib_read import lib_read
from peewee import SqliteDatabase
from rdkit.Chem.BRICS import reverseReactions
from rdkit import Chem
import random
from .config import logger


def react(user_frag, mate_mols):
    """
    a reactor for brics fragments
    :param user_frag:
    :param mate_mols:
    :return:
    """

    for mate_frag in mate_mols:

        # list of reactions in the Brics system, try them all in a random order (to prevent systematic bias)
        random.shuffle(reverseReactions, random=random.random)
        for rxn in reverseReactions:
            # run the appropriate reactions
            ps = rxn.RunReactants((user_frag, mate_frag))
            # give up the reacted fragment
            if len(ps) > 0:
                # logger.debug(f"Mated {Chem.MolToSmiles(user_frag)} and "
                #              f"{Chem.MolToSmiles(mate_frag)}")
                for p in ps:
                    yield p[0]
            ps = rxn.RunReactants((mate_frag, user_frag))
            # give up the reacted fragment
            if len(ps) > 0:
                # logger.debug(f"Mated {Chem.MolToSmiles(user_frag)} and "
                #              f"{Chem.MolToSmiles(mate_frag)}")
                for p in ps:
                    yield p[0]

        # logger.error(f"No reactions took place between {Chem.MolToSmiles(user_frag)} and "
        #              f"{Chem.MolToSmiles(mate_frag)}")


def mate(db,
         user_smile,
         user_num_pseudo_atoms,
         missing_piece_fc,
         missing_piece_len,
         permissivity=1.0,
         max_num_children=25,
         final_limits=None,
         must_have_patterns=None,
         must_not_have_patterns=None
         ):
    """
    mate()

    'main' function for this module
    generates children using a library of processed fragments, and a specific user fragment
        db = location of the fragment library database
        user_frag = Fragment object from lib_read.py
        missing_piece_fc = the estimated FC of the pieces removed from the user parent, to create this user frag
        missing_piece_len = see above but for length
        permissivity = scaling factor that controls how deviant molecules are allowed to be:
                        1~parity to parent, 0~smallest possible additions, >1~allow adding larger/more complex fragments
        max_num_children = the number of children to return for this fragment
        max_pseudoatoms = not user facing, use caution. The starting value for how many pseudoatoms a mate frag is allowed

    :param db:
    :param user_smile:
    :param user_num_pseudo_atoms:
    :param missing_piece_fc:
    :param missing_piece_len:
    :param permissivity:
    :param max_num_children:
    :param max_pseudoatoms:
    :return:

    """

    # now get the list of pseudoatoms in the user fragment
    _, _, p_atoms = fragment.get_pseudo_atoms(user_smile)
    user_frag_mol = Chem.MolFromSmiles(user_smile)
    if user_frag_mol is None:
        logger.error(f"Could not process this fragment into a Mol object, skipping: {user_smile}")
        return [], {}

    # the final list of all the children
    finished_children = list()

    # a counter
    num_children = 0

    # for the set of pseudoatoms in this fragment, get the list of pseudoatoms that can react
    allowed_atoms = list()
    for atom in p_atoms:
        atom = int(atom)
        allowed_atoms.extend(brics_mating_rules[atom])
    allowed_atoms = set(allowed_atoms)

    # default settings for max_pseudoatoms assigned here
    # careful changing things, everything is more or less in a delicate balance
    max_pseudoatoms = round(missing_piece_fc * permissivity * user_num_pseudo_atoms) - 1.0

    # get the mate fragments (first pass)
    mate_mols = get_mate_fragments(db,
                                   missing_piece_fc=missing_piece_fc,
                                   missing_piece_len=missing_piece_len,
                                   permissivity=permissivity,
                                   max_pseudoatoms=max_pseudoatoms,
                                   allowed_atoms=allowed_atoms,
                                   max_num_children=max_num_children)

    # get the generator that will produce the children
    children = react(user_frag_mol, mate_mols)

    # catch the end of the generator
    done_mating = False

    # loop through generator until correct number of children generated, or until generator craps out
    while (not done_mating) and (num_children < max_num_children):
        try:
            # get the next child
            try:
                child = next(children)
            except StopIteration:
                done_mating = True
                break
            child_smile = Chem.MolToSmiles(child, isomericSmiles=True)
            child = Chem.MolFromSmiles(child_smile)
            _, _, p_atoms = fragment.get_pseudo_atoms(child_smile)
            num_p_atoms = len(p_atoms)

            # within loop values of parameters, to iteratively change as the loop continues
            # hopefully safe against infinite loops
            max_pseudoatoms_tmp = round(max_pseudoatoms / user_num_pseudo_atoms)
            permissivity_tmp = permissivity
            missing_piece_len_tmp = missing_piece_len

            # this part is where any children that still have pseudoatoms will be "completed"
            # basically, increasingly smaller and less complex pieces will be attached, until no pseudoatoms remain
            retries = 0
            while num_p_atoms > 0:
                # if retries > 5:
                # logger.warning(f"Continuing to mate {user_smile} child # {num_children}, "
                #                f"{child_smile}, filling in child, attempt # {retries}:")
                # decrease overall permissivity
                # for the set of pseudoatoms in this fragment, get the list of pseudoatoms that can react
                allowed_atoms_tmp = list()
                for atom in p_atoms:
                    atom = int(atom)
                    allowed_atoms_tmp.extend(brics_mating_rules[atom])
                allowed_atoms_tmp = set(allowed_atoms_tmp)
                max_pseudoatoms_tmp = max_pseudoatoms_tmp - 1
                permissivity_tmp = round(permissivity_tmp * 0.4, 2)
                mate_mols = get_mate_fragments(db,
                                               missing_piece_fc=missing_piece_fc,
                                               missing_piece_len=missing_piece_len_tmp,
                                               permissivity=permissivity_tmp,
                                               max_pseudoatoms=max_pseudoatoms_tmp,
                                               allowed_atoms=allowed_atoms_tmp,
                                               max_num_children=max_num_children)

                # react the original child with a new mate frag
                complete_children = react(child, mate_mols)

                # get the result as a child, in order to mask the original child
                # this allows convergence to a complete child to occur
                try:
                    child = next(complete_children)
                    child_smile = Chem.MolToSmiles(child, isomericSmiles=True)
                    child = Chem.MolFromSmiles(child_smile)
                except Exception as e:  # pylint: disable=broad-except
                    if retries > 5:
                        # logger.warning(f"Couldn't find any children for {child_smile}.\tTrying again.")
                        permissivity_tmp = round(permissivity_tmp * 2, 2)
                    retries += 1
                    if retries < 10:
                        continue
                    else:
                        logger.error(f"Couldn't find any children for {child_smile}. Aborting.")
                        break

                child_smile = Chem.MolToSmiles(child, isomericSmiles=True)

                # recalculate new  estimated missing piece length
                # basically, how much was added since the original fragment (when missing_piece_len was defined)?
                # remove that much from the estimate of missing_piece_len
                new_diff = len(child_smile) - len(user_smile)
                missing_piece_len_tmp = missing_piece_len - new_diff

                # loop condition, if there are still pseudoatoms, go again
                _, _, p_atoms = fragment.get_pseudo_atoms(child_smile)
                num_p_atoms = len(p_atoms)

            # need to check if the child is complete or not
            child_smile = Chem.MolToSmiles(child, isomericSmiles=True)
            _, _, p_atoms = fragment.get_pseudo_atoms(child_smile)
            num_p_atoms = len(p_atoms)
            if num_p_atoms > 0:
                num_children += 1
                continue

            # fix some chemical details of the generated child and store
            child.UpdatePropertyCache(strict=True)
            Chem.SanitizeMol(child)
            finished_children.append(child)
            num_children += 1
            # logger.debug(f"Finished with child {user_smile} child # {num_children}, {child_smile}.")

        except Exception as e:  # pylint: disable=broad-except
            # generator exhausted, or some other failure
            # (shouldn't matter much, just leave the loop, there will be others)
            logger.error(f"Couldn't mate {user_smile} for some reason.")
            logger.error(f"Error: {e}")
            done_mating = True
            # DO NOT RAISE AN ERROR HERE I SWEAR
            # LOOP HAS TO CONTINUE HERE, ERROR IS KNOWN AND I'M NOT WORRIED ABOUT IT

    logger.info(f"Finished Mating {user_smile}. Produced {len(finished_children)} children.")
    # pre-filter the children of this parent
    pre_filter_values = apply_filter(final_limits,
                                     finished_children,
                                     must_have_patterns=must_have_patterns,
                                     must_not_have_patterns=must_not_have_patterns)

    return finished_children, pre_filter_values


def get_mate_fragments(db,
                       missing_piece_fc,
                       missing_piece_len,
                       permissivity,
                       max_pseudoatoms,
                       allowed_atoms,
                       max_num_children):
    """
    for a given set of criteria, get compatible mate fragments
    db = location of the fragment library database
    missing_piece_fc = the estimated FC of the pieces removed from the user parent, to create this user frag
    missing_piece_len = see above but for length
    permissivity = scaling factor that controls how deviant molecules are allowed to be:
                    1~parity to parent, 0~smallest possible additions, >1~allow adding larger/more complex fragments
    max_num_children = the number of children to return for this fragment
    max_pseudoatoms = not user facing, use caution. The starting value for how many pseudoatoms a mate frag is allowed

    :param db:
    :param missing_piece_fc:
    :param missing_piece_len:
    :param permissivity:
    :param max_pseudoatoms:
    :param allowed_atoms:
    :param max_num_children:
    :return:
    """

    # sanity checks
    # ensure that all the input parameters make sense
    # minimum number of pseudoatoms is 1
    if max_pseudoatoms <= 0.0:
        max_pseudoatoms = 1.0

    # minimum smile length cutoff to use is 15 (chosen apropos of nothing)
    # important to have one of these here, as sometimes the estimated length can be negative for a variety of reasons
    if missing_piece_len < 20.0:
        missing_piece_len = 20.0

    # if the missing piece FC is zero, and we are not yet at the "base case" of missing_piece_len == 0
    # and we are not allowing any more pseudoatoms, then use the missing_piece_len as a rough
    # metric to allow more pseudoatoms, and increase the allowed FC
    # this will ensure that if we start with a very small piece and attach a small piece, then we are more likely
    # to attach multiple pieces in order to try to get back the rough size of the user parent
    if missing_piece_fc == 0.0 and missing_piece_len > 20.0 and max_pseudoatoms == 1.0:
        max_pseudoatoms = round(missing_piece_len / 10)
        missing_piece_fc = (max_pseudoatoms - 1.0)

    # load the fragment library
    if not os.path.exists(db):
        raise Exception("Database file not found: " + db)
    lib = SqliteDatabase(db)
    Fragment, _, _, _ = lib_read(lib)
    with lib:
        # define a length cutoff dynamically, with a minimum
        length_cutoff = ((missing_piece_len * 1.5) / (0.5 * max_pseudoatoms)) * permissivity
        if length_cutoff < 20.0:
            length_cutoff = 20.0

        # this is the part where we avoid the "ORDER BY RANDOM() LIMIT" query by using the FragmentIndex.
        # It saves quite a bit of time, at a moderate cost of memory.
        # returns an iterator rather than a list, so we only fetch as many fragments as we need.
        fx = frag_index(db)
        lim = max_num_children * 100
        fids = fx.mate_query(
            permissivity, max_pseudoatoms,
            missing_piece_fc, missing_piece_len,
            length_cutoff, allowed_atoms,
            limit=lim
        )

        def mate_mols():
            for fid in fids:
                my_mate = Fragment.select().where(Fragment.id == fid)[0]
                mol = Chem.MolFromSmiles(my_mate.smile)
                yield mol

        fast_iterator = mate_mols()
        # Don't close db when we return fast iterator, it will be closed when the iterator closes.
        return fast_iterator
