from .config import logger
import uuid
from .child_filter import get_filter_values, apply_filter
from .config import drug_like_params
from rdkit import Chem
from .selfies_methods import (
    selfies_substitution,
    selfies_deletion,
    selfies_insertion,
    random_selfies_generator,
    selfies_scanner,
)
from typing import List, Set
from collections import defaultdict
from .fragment import libgen
from peewee import SqliteDatabase
from .lib_read import lib_read
import re
from .fragment_index import frag_index
from .mate import mate
import os
import random
from .jensen_crossover import crossover as crossover_gb
from .jensen_mutate import mutate as mutate_gb
from .jensen_selfies_crossover import crossover as selfies_crossover_gb
from .jensen_selfies_mutate import mutate as selfies_mutate_gb
from crem.crem import grow_mol as crem_grow
from crem.crem import mutate_mol as crem_mutate
from copy import deepcopy


class Deriver(object):
    def __init__(self):

        logger.info("Initializing a new Deriver object!")
        self.data = self._Data()

    class _Data(object):

        """
        This object is meant to be a home for all the data and parameters used by deriver.
        It should be fine to read things from here, but do not set anything directly. Instead,
        use the methods that begin with "set_" in order to change these values. This
        prevents incompatible changes from being made silently.
        """

        def __init__(self):
            self.seed_smiles = None
            self.seed_mols = None
            self.filter_params = deepcopy(drug_like_params)
            self.filter = False
            self.child_db = None
            self.all_good_selfies_children = None
            self.all_good_scanner_children = None
            self.all_good_selfies_gb_children = None
            self.all_good_smiles_gb_children = None
            self.all_good_local_children = None
            self.filter_molecules = None
            self.must_have_patterns = None
            self.must_not_have_patterns = None
            self.heritage = defaultdict(list)
            self.track_heritage = True
            # BRICS specific
            self.seed_frags = None  # these are the fragments of the seed molecules
            self.fragment_source_db = None  # this is the location of the fragment DB
            self.seed_frag_db = None  # the is the DB where the seed_frags are stored and info about them
            self.all_good_brics_children = (
                None  # this is where the good (filtered) BRICS children are saved
            )
            self.crem_source_db = None  # the location of the crem fragment database used by the local space methods

    def set_seeds(self, seeds: list):

        """
        set the seeds that are used to generate new molecules
        :param seeds:
        :return:
        """

        if isinstance(seeds[0], Chem.rdchem.Mol):
            self.data.seed_smiles = []
            self.data.seed_mols = []
            for seed in seeds:
                smile = Chem.MolToSmiles(seed, isomericSmiles=True)
                self.data.seed_smiles.append(smile)
                self.data.seed_mols.append(seed)
        elif isinstance(seeds[0], str):
            self.data.seed_smiles = []
            self.data.seed_mols = []
            for smile in seeds:
                seed = Chem.MolFromSmiles(smile, sanitize=True)
                iso_smile = Chem.MolToSmiles(seed, isomericSmiles=True)
                self.data.seed_smiles.append(iso_smile)
                self.data.seed_mols.append(seed)
        else:
            logger.error(
                "Seeds must be provided as an iterable of Mol objects, or SMILES strings."
            )
            return 0

        return 1

    def set_must_have_patterns(self, must_have_patterns: List[str]):
        if isinstance(must_have_patterns, list):
            assert isinstance(must_have_patterns[0], str)
        elif must_have_patterns is None:
            pass
        else:
            raise TypeError(
                "must_have_patterns must be None or a list of SMARTS strings"
            )
        self.data.must_have_patterns = must_have_patterns

    def set_must_not_have_patterns(self, must_not_have_patterns: List[str]):
        if isinstance(must_not_have_patterns, list):
            assert isinstance(must_not_have_patterns[0], str)
        elif must_not_have_patterns is None:
            pass
        else:
            raise TypeError(
                "must_have_patterns must be None or a list of SMARTS strings"
            )
        self.data.must_not_have_patterns = must_not_have_patterns

    def set_filter_molecules(self, filter_molecules: Set[str]):
        assert isinstance(filter_molecules, set)
        assert isinstance(iter(filter_molecules).__next__(), str)
        self.data.filter_molecules = filter_molecules
        return 1

    def toggle_heritage_tracking(self):
        self.data.track_heritage = not self.data.track_heritage
        logger.info(f"Heritage tracking is now {self.data.track_heritage}")
        return 1

    def enable_and_expand_filter(self, seeds: list = None):

        """
        Turn on the filtering of produced molecules, and expand the ranges to include the provided
        seed molecules, or the stored seed molecules.
        :param seeds:
        :return:
        """

        if seeds is None:
            if self.data.seed_smiles is None:
                logger.info("Turning on filter using default filter values:")
                logger.info(self.data.filter_params)
                self.data.filter = True
                return 1
            else:
                logger.info(
                    f"Expanding filter ranges using current seed set: {self.data.seed_smiles}"
                )
                seeds = self.data.seed_mols
        for seed_mol in seeds:
            parent_values = get_filter_values(seed_mol)
            self._update_filter_params(parent_values)
        logger.info("Done! Turning on filter using new filter values:")
        logger.info(self.data.filter_params)
        self.data.filter = True
        return 1

    def manual_filter_set(self, var: str, val1, val2=None):

        """
        Use this function to change the filter values safely.
        :param var:
        :param val1:
        :param val2:
        :return:
        """

        if var not in self.data.filter_params:
            logger.error(f"{var} is not a valid filter parameter, try again.")
            raise Exception

        logger.info(
            f"Updating {var} from {self.data.filter_params[var]} to {val1}, {val2}"
        )

        if isinstance(self.data.filter_params[var], tuple):

            if val2 is None:
                logger.error(
                    f"{var} requires TWO values, lower and upper bounds, to be set."
                )
                raise Exception

            elif val2 <= val1:
                logger.error(
                    f"{var} requires TWO values, lower and upper bounds, to be set. "
                    f"The second must be LARGER than the first."
                )
                raise Exception

            else:
                # actually set the values!
                if isinstance(self.data.filter_params[var][0], int):
                    val1 = int(val1)
                    val2 = int(val2)
                self.data.filter_params[var] = (val1, val2)

        else:
            # here set these values
            if isinstance(self.data.filter_params[var], int):
                val1 = int(val1)
            self.data.filter_params[var] = val1

        logger.info(f"Done! Current filter parameters are: {self.data.filter_params}")

        return 1

    def _update_filter_params(self, parent_values):

        filter_params = {
            # ranges
            "MW": (
                min(parent_values["MW"], self.data.filter_params["MW"][0]),
                max(parent_values["MW"], self.data.filter_params["MW"][1]),
            ),
            "num_heavy_atoms": (
                min(
                    parent_values["num_heavy_atoms"],
                    self.data.filter_params["num_heavy_atoms"][0],
                ),
                max(
                    parent_values["num_heavy_atoms"],
                    self.data.filter_params["num_heavy_atoms"][1],
                ),
            ),
            "num_carbons": (
                min(
                    parent_values["num_carbons"],
                    self.data.filter_params["num_carbons"][0],
                ),
                max(
                    parent_values["num_carbons"],
                    self.data.filter_params["num_carbons"][1],
                ),
            ),
            "num_hetero_atoms": (
                min(
                    parent_values["num_hetero_atoms"],
                    self.data.filter_params["num_hetero_atoms"][0],
                ),
                max(
                    parent_values["num_hetero_atoms"],
                    self.data.filter_params["num_hetero_atoms"][1],
                ),
            ),
            "hc_ratio": (
                min(parent_values["hc_ratio"], self.data.filter_params["hc_ratio"][0]),
                max(parent_values["hc_ratio"], self.data.filter_params["hc_ratio"][1]),
            ),
            "charge": (
                min(parent_values["charge"], self.data.filter_params["charge"][0]),
                max(parent_values["charge"], self.data.filter_params["charge"][1]),
            ),
            "logP": (
                min(parent_values["logP"], self.data.filter_params["logP"][0]),
                max(parent_values["logP"], self.data.filter_params["logP"][1]),
            ),
            "fSP3": (
                min(parent_values["fSP3"], self.data.filter_params["fSP3"][0]),
                max(parent_values["fSP3"], self.data.filter_params["fSP3"][1]),
            ),
            # upper limits
            "HBA": max(parent_values["HBA"], self.data.filter_params["HBA"]),
            "HBD": max(parent_values["HBD"], self.data.filter_params["HBD"]),
            "tPSA": max(parent_values["tPSA"], self.data.filter_params["tPSA"]),
            "rot_bonds": max(
                parent_values["rot_bonds"], self.data.filter_params["rot_bonds"]
            ),
            "rigid_bonds": max(
                parent_values["rigid_bonds"], self.data.filter_params["rigid_bonds"]
            ),
            "num_rings": max(
                parent_values["num_rings"], self.data.filter_params["num_rings"]
            ),
            "max_ring_size": max(
                parent_values["max_ring_size"], self.data.filter_params["max_ring_size"]
            ),
            "num_charges": max(
                parent_values["num_charges"], self.data.filter_params["num_charges"]
            ),
            "num_chiral_centers": max(
                parent_values["num_chiral_centers"],
                self.data.filter_params["num_chiral_centers"],
            ),
        }

        self.data.filter_params = filter_params

    def derive_selfies(
        self,
        n_children: int = 100,
        mut_rate: float = 0.03,
        mut_min: int = 1,
        mut_max: int = 2,
    ):

        good_children = []
        all_filtered_children = {}
        self.data.all_good_selfies_children = []
        n_seeds = len(self.data.seed_smiles)
        if n_children < n_seeds:
            n_children_per_seed = 1
        else:
            n_children_per_seed = round(n_children / n_seeds) + 1

        logger.info(
            f"Mutating SELFIES to create {n_children_per_seed} children per seed."
        )

        if self.data.filter:
            filter_params = self.data.filter_params
        else:
            filter_params = None

        for seed in self.data.seed_smiles:
            children = selfies_substitution(
                parent_smiles=seed,
                n_children=round(n_children_per_seed * 0.7),  # three internal methods
                mut_rate=mut_rate,
                mut_min=mut_min,
                mut_max=mut_max,
            )
            if self.data.track_heritage:
                self.data.heritage[seed] += children
            child_mols = [
                Chem.MolFromSmiles(child, sanitize=True) for child in children
            ]

            children = selfies_insertion(
                parent_smiles=seed,
                n_children=round(n_children_per_seed * 0.15),  # three internal methods
                mut_rate=mut_rate,
                mut_min=mut_min,
                mut_max=mut_max,
            )
            if self.data.track_heritage:
                self.data.heritage[seed] += children
            child_mols += [
                Chem.MolFromSmiles(child, sanitize=True) for child in children
            ]

            children = selfies_deletion(
                parent_smiles=seed,
                n_children=round(n_children_per_seed * 0.15),  # three internal methods
                mut_rate=mut_rate,
                mut_min=mut_min,
                mut_max=mut_max,
            )
            if self.data.track_heritage:
                self.data.heritage[seed] += children
            child_mols += [
                Chem.MolFromSmiles(child, sanitize=True) for child in children
            ]

            # filter children
            filtered_children = apply_filter(
                filter_params,
                child_mols,
                self.data.must_have_patterns,
                self.data.must_not_have_patterns,
            )
            all_filtered_children.update(filtered_children)

            for child in filtered_children:
                if filtered_children[child]["is_good"]:
                    # check the cache
                    if self.data.filter_molecules:
                        if child not in self.data.filter_molecules:
                            good_children.append(child)
                        else:
                            logger.debug(f"skipping previously seen molecule: {child}")
                    else:
                        good_children.append(child)

        logger.info(f"Generated {len(good_children)} 'good' children.")
        self.data.all_good_selfies_children = good_children
        return good_children, all_filtered_children

    def random_selfies(self, *, n_symbols: int = 100, n_molecules: int = 100):

        good_children = []
        rand_selfies_gen = random_selfies_generator(n_symbols=n_symbols)
        self.data.all_good_selfies_children = []

        if self.data.filter:
            filter_params = self.data.filter_params
        else:
            filter_params = None

        while len(good_children) < n_molecules:
            child_mols = [
                Chem.MolFromSmiles(next(rand_selfies_gen))
                for i in range(n_molecules - len(good_children))
            ]
            # filter children
            filtered_children = apply_filter(
                filter_params,
                child_mols,
                self.data.must_have_patterns,
                self.data.must_not_have_patterns,
            )

            for child in filtered_children:
                if filtered_children[child]["is_good"]:
                    # check the cache
                    if self.data.filter_molecules:
                        if child not in self.data.filter_molecules:
                            good_children.append(child)
                        else:
                            logger.debug(f"skipping previously seen molecule: {child}")
                    else:
                        good_children.append(child)

        self.data.all_good_selfies_children = good_children
        return good_children

    def scan_selfies(self, safe_mode: bool = False):

        """
        Return all possible single substitution children for all the seeds.
        """
        if self.data.filter:
            filter_params = self.data.filter_params
        else:
            filter_params = None
        good_children = []
        self.data.all_good_scanner_children = []
        all_filtered_children = {}

        for seed in self.data.seed_smiles:
            children = selfies_scanner(parent_smiles=seed, safe_mode=safe_mode)
            if len(children) == 0:
                continue
            if self.data.track_heritage:
                self.data.heritage[seed] += children

            filtered_children = apply_filter(
                filter_params,
                [Chem.MolFromSmiles(child) for child in children],
                self.data.must_have_patterns,
                self.data.must_not_have_patterns,
            )
            all_filtered_children.update(filtered_children)

            for child in filtered_children:
                if filtered_children[child]["is_good"]:
                    # check the cache
                    if self.data.filter_molecules:
                        if child not in self.data.filter_molecules:
                            good_children.append(child)
                        else:
                            logger.debug(f"skipping previously seen molecule: {child}")
                    else:
                        good_children.append(child)

        self.data.all_good_scanner_children = good_children
        return good_children, all_filtered_children

    def derive_gb(self, n_children: int = 100, representation="selfies"):

        assert len(self.data.seed_smiles) > 0
        children = []
        good_children = []
        if representation == "selfies":
            self.data.all_good_selfies_gb_children = []
            crossover_fn = selfies_crossover_gb
            mutation_fn = selfies_mutate_gb
        elif representation == "smiles":
            self.data.all_good_smiles_gb_children = []
            crossover_fn = crossover_gb
            mutation_fn = mutate_gb
        else:
            raise ValueError(
                'Must specify derivation kind as one of "smiles" or "selfies"'
            )

        if self.data.filter:
            filter_params = self.data.filter_params
        else:
            filter_params = None

        # parent_a_smiles, parent_b_smiles = (None, None)  # not used. Delete?
        if len(self.data.seed_smiles) > 1:
            do_crossover = True
            new_child = None
        else:
            do_crossover = False
            new_child = self.data.seed_mols[0]

        for _ in range(n_children):
            if do_crossover:
                parent_a, parent_b = random.sample(self.data.seed_mols, 2)
                new_child = crossover_fn(parent_a, parent_b)
            if new_child is not None:
                mutated_child = mutation_fn(new_child)
                if mutated_child is None:
                    continue
            else:
                continue
            children.append(mutated_child)

        filtered_children = apply_filter(
            filter_params,
            children,
            self.data.must_have_patterns,
            self.data.must_not_have_patterns,
        )

        # bugfix for empty strings
        try:
            del filtered_children[""]
        except KeyError:
            pass

        for child in filtered_children:
            if filtered_children[child]["is_good"]:
                # check the cache
                if self.data.filter_molecules:
                    if child not in self.data.filter_molecules:
                        good_children.append(child)
                    else:
                        logger.debug(f"skipping previously seen molecule: {child}")
                else:
                    good_children.append(child)

        if representation == "smiles":
            self.data.all_good_smiles_gb_children = good_children
        else:
            self.data.all_good_selfies_gb_children = good_children
        logger.info(f"Generated {len(good_children)} 'good' children.")
        return good_children, filtered_children

    def set_fragment_source_db(self, frag_db):

        """
        set the location for the fragment database that is used to mate molecules
        :param frag_db:
        :return:
        """

        self.data.fragment_source_db = frag_db
        return 1

    def set_crem_source_db(self, crem_db: str):

        """
        set the location for the fragment database that is used by crem
        :param crem_db:
        :return:
        """
        self.data.crem_source_db = crem_db
        return 1

    def _process_seeds_for_brics(self):
        """
        This function parses the seed molecules and gets the BRICS fragments they make, then cleans them
        :return:
        """
        logger.info("Processing seeds to create scaffold fragments:")

        self.data.seed_frag_db = f"seed_frags_{uuid.uuid4()}.db"
        self.data.seed_frags = []

        # Databases are used in lieu of alternatives (like dataframes) in order to operate on larger datasets
        # where memory may be a bottleneck, and to ensure that only one source of truth exists for performing
        # this calculation.
        libgen(
            self.data.seed_mols, self.data.seed_frag_db
        )  # libgen is the basic command that generates frag libraries
        seed_frag_db = SqliteDatabase(
            self.data.seed_frag_db
        )  # load the db we just made
        Fragment, Heritage, _, _ = lib_read(
            seed_frag_db
        )  # we only care about fragment and heritage at this point
        seed_frag_db.connect()

        # get all the fragments from the user molecule
        user_frags = (
            Fragment.select()
            .join(Heritage, on=Heritage.frag)
            .where(
                Heritage.parent != Heritage.frag
            )  # only get fragments, not intact molecules
            .order_by(Fragment.frag_coeff.desc())
        )  # largest and most complex frags first

        # for every fragment from the user provided parent mol
        for user_frag in user_frags:

            # we want to ignore really small fragments, by counting atom symbols
            smaller_smile = re.sub(
                r"\[[0-9]+\*\]", "", user_frag.smile
            )  # ignore pseudoatoms
            smaller_smile = re.sub(
                r"[0-9]", "", smaller_smile
            )  # ignore numbers in general
            smaller_smile = re.sub(
                r"[\(\)=#@\-\]\[]+", "", smaller_smile
            )  # ignore a bunch of other symbols
            if len(smaller_smile) < 4:  # if there are less than four atoms
                logger.warning(f"Skipping user_frag {user_frag.smile} due to size.")
                continue

            # using this fragment and the whole parent molecule, estimate the "missing" FC and size
            try:
                parent = Heritage.get(Heritage.frag_id == user_frag.id).parent
            # todo: actual exception is deriver.lib_read.FragmentDoesNotExist, check if we can except just this case
            except Exception as e:  # pylint: disable=broad-except
                logger.warning(f"Encountered exception {e}")
                logger.warning(
                    "If this exception describes a missing parent in the Heritage table, this bug"
                    "is known and is being handled as intended."
                )
                continue

            # if parent is not None:
            missing_piece_fc = (
                parent.frag_coeff - user_frag.frag_coeff
            ) - 1.0  # -1.0 because two pieces combine
            missing_piece_len = len(parent.smile) - len(
                user_frag.smile
            )  # approximation
            # else:
            #       missing_piece_fc = 3.0  # approximation
            #       missing_piece_len = 40  # approximation

            # this is what we are going to keep
            seed_frag = (
                user_frag.smile,
                user_frag.num_pseudo_atoms,
                missing_piece_fc,
                missing_piece_len,
                parent.smile,
            )

            self.data.seed_frags.append(seed_frag)

        seed_frag_db.close()
        os.remove(
            self.data.seed_frag_db
        )  # we don't really care to keep the seed fragment database
        logger.info(
            f"Done! There are {len(self.data.seed_frags)} seed fragments ready to be mated."
        )
        # it is faster to keep these in memory rather than using the database

    def derive_brics(self, n_children: int = 100, permissivity: float = 1.0):
        """

        :param n_children: How many children do you want, in total. This is an approximation, not exact.
        :param permissivity: How unlike the parent molecules is the child allowed to be, higher is generally larger
        :return: (all_good_children [a list of smiles], all_filtered_children [a dict of values about the molecules])
        """
        # process the seeds
        self._process_seeds_for_brics()

        # get the "maximum number of children" per fragment
        n_seed_frags = len(self.data.seed_frags)
        if n_seed_frags == 0:
            logger.warning("No seed fragments! Cannot derive brics from these seeds.")
            return [], {}

        if self.data.filter:
            filter_params = self.data.filter_params
        else:
            filter_params = None

        if n_children < n_seed_frags:
            children_per_seed_frag = 1
        else:
            children_per_seed_frag = round(n_children / n_seed_frags) + 1

        logger.info(
            f"Creating/reading a fragment index for {self.data.fragment_source_db}"
        )
        # generate the frag index once, first, so it doesn't get generated in each pool process
        # this index serves to dramatically speed up queries
        frag_index(self.data.fragment_source_db)

        # again this is more of a guideline
        logger.info(
            f"Mating to create {children_per_seed_frag} children per seed frag."
        )

        all_filtered_children = (
            dict()
        )  # the filter returns a dictionary of calculated pk values and filter status

        for (
            seed_frag_smile,
            seed_frag_num_pa,
            missing_p_fc,
            missing_p_len,
            parent_smile,
        ) in self.data.seed_frags:
            try:
                res = mate(
                    self.data.fragment_source_db,
                    seed_frag_smile,
                    seed_frag_num_pa,
                    missing_p_fc,
                    missing_p_len,
                    permissivity,
                    children_per_seed_frag,
                    filter_params,
                    self.data.must_have_patterns,
                    self.data.must_not_have_patterns,
                )

                (
                    _,
                    filter_values,
                ) = res  # we only care about the filter dict, since it has everything
                all_filtered_children.update(filter_values)  # update our master dict
                if self.data.track_heritage:
                    self.data.heritage[parent_smile] += list(
                        filter_values.keys()
                    )  # this keeps track of heritage
            except IndexError as e:
                # This bug has never really been explored that much.
                logger.warning(
                    f"Error when trying to mate a molecule, ignoring this molecule. Error: {e}"
                )

        all_good_children = []
        self.data.all_good_brics_children = (
            []
        )  # every time you call `derive_brics` it deletes any old results
        for child in all_filtered_children:
            if all_filtered_children[child]["is_good"]:
                # check the cache of previously seen molecules (which we want to avoid reproducing)
                if self.data.filter_molecules:
                    if child not in self.data.filter_molecules:
                        all_good_children.append(child)
                    else:
                        logger.debug(f"skipping previously seen molecule: {child}")
                else:
                    # there is no provided list of molecules to skip
                    all_good_children.append(child)

        logger.info(
            f"Generated {len(self.data.all_good_brics_children)} 'good' children."
        )
        self.data.all_good_brics_children = all_good_children
        return all_good_children, all_filtered_children

    def derive_local_space(
        self, approx_children_per_seed: int = 1000, min_inc: int = -2, max_inc: int = 2
    ):

        if self.data.crem_source_db is None:
            raise AttributeError(
                "No crem source db. Please use `.set_crem_source_db()` to provide a source db. "
                "See readme for more information."
            )

        if self.data.filter:
            filter_params = self.data.filter_params
        else:
            filter_params = None

        children = []
        good_children = []
        # first we make the molecules by using grow to replace hydrogens, and mutate to do everything else
        for i, seed_mol in enumerate(self.data.seed_mols):
            logger.info(f"Growing children for {self.data.seed_smiles[i]}:")
            grown_children_smiles = crem_grow(
                seed_mol, self.data.crem_source_db, return_mol=False
            )
            grown_children_mols = [
                Chem.MolFromSmiles(smile, sanitize=True)
                for smile in grown_children_smiles
            ]
            children.extend(grown_children_mols)
            logger.info("Done!")

            logger.info(f"Mutating children for {self.data.seed_smiles[i]}:")
            mutate_children_smiles = crem_mutate(
                seed_mol,
                self.data.crem_source_db,
                return_mol=False,
                max_replacements=approx_children_per_seed,
                min_size=1,
                max_size=5,
                min_inc=min_inc,
                max_inc=max_inc,
            )
            mutate_children_mols = [
                Chem.MolFromSmiles(smile, sanitize=True)
                for smile in mutate_children_smiles
            ]
            children.extend(mutate_children_mols)
            logger.info("Done!")

        if self.data.filter:
            logger.info("Applying filters to local space children:")

        filtered_children = apply_filter(
            filter_params,
            children,
            self.data.must_have_patterns,
            self.data.must_not_have_patterns,
        )

        # bugfix for empty strings
        try:
            del filtered_children[""]
        except KeyError:
            pass

        logger.info("Done!")

        for child in filtered_children:
            if filtered_children[child]["is_good"]:
                # check the cache
                if self.data.filter_molecules:
                    if child not in self.data.filter_molecules:
                        good_children.append(child)
                    else:
                        logger.debug(f"skipping previously seen molecule: {child}")
                else:
                    good_children.append(child)

        self.data.all_good_local_children = good_children
        logger.info(
            f"Generated {len(self.data.all_good_local_children)} 'good' children."
        )

        return good_children, filtered_children
