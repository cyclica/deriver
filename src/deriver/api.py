import logging as logger
from .child_filter import get_filter_values, apply_filter
from .config import drug_like_params
from rdkit import Chem
from .selfies_methods import selfies_substitution, selfies_deletion, selfies_insertion, random_selfies_generator, \
    selfies_scanner
from typing import List, Set
from collections import defaultdict


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
            self.filter_params = drug_like_params
            self.filter = False
            self.child_db = None
            self.all_good_selfies_children = None
            self.all_good_scanner_children = None
            self.filter_molecules = None
            self.must_have_patterns = None
            self.heritage = defaultdict(list)

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
                self.data.seed_smiles.append(smile)
                self.data.seed_mols.append(seed)
        else:
            logger.error("Seeds must be provided as an iterable of Mol objects, or SMILES strings.")
            return 0

        return 1

    def set_must_have_patterns(self, must_have_patterns: List[str]):
        if isinstance(must_have_patterns, list):
            assert isinstance(must_have_patterns[0], str)
        elif must_have_patterns is None:
            pass
        else:
            raise TypeError("must_have_patterns must be None or a list of SMARTS strings")
        self.data.must_have_patterns = must_have_patterns

    def set_filter_molecules(self, filter_molecules: Set[str]):
        assert isinstance(filter_molecules, set)
        assert isinstance(iter(filter_molecules).__next__(), str)
        self.data.filter_molecules = filter_molecules
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
                logger.info(f"Expanding filter ranges using current seed set: {self.data.seed_smiles}")
                seeds = self.data.seed_mols
        for seed_mol in seeds:
            parent_values = get_filter_values(seed_mol)
            self._update_filter_params(parent_values)
        logger.info(f"Done! Turning on filter using new filter values:")
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

        logger.info(f"Updating {var} from {self.data.filter_params[var]} to {val1}, {val2}")

        if isinstance(self.data.filter_params[var], tuple):

            if val2 is None:
                logger.error(f"{var} requires TWO values, lower and upper bounds, to be set.")
                raise Exception

            elif val2 <= val1:
                logger.error(f"{var} requires TWO values, lower and upper bounds, to be set. "
                             f"The second must be LARGER than the first.")
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
            "MW": (min(parent_values["MW"], self.data.filter_params["MW"][0]),
                   max(parent_values["MW"], self.data.filter_params["MW"][1])),
            "num_carbons": (min(parent_values["num_carbons"], self.data.filter_params["num_carbons"][0]),
                            max(parent_values["num_carbons"], self.data.filter_params["num_carbons"][1])),
            "num_hetero_atoms": (min(parent_values["num_hetero_atoms"],
                                     self.data.filter_params["num_hetero_atoms"][0]),
                                 max(parent_values["num_hetero_atoms"],
                                     self.data.filter_params["num_hetero_atoms"][1])),
            "hc_ratio": (min(parent_values["hc_ratio"], self.data.filter_params["hc_ratio"][0]),
                         max(parent_values["hc_ratio"], self.data.filter_params["hc_ratio"][1])),
            "charge": (min(parent_values["charge"], self.data.filter_params["charge"][0]),
                       max(parent_values["charge"], self.data.filter_params["charge"][1])),
            "logP": (min(parent_values["logP"], self.data.filter_params["logP"][0]),
                     max(parent_values["logP"], self.data.filter_params["logP"][1])),

            # upper limits
            "HBA": max(parent_values["HBA"], self.data.filter_params["HBA"]),
            "HBD": max(parent_values["HBD"], self.data.filter_params["HBD"]),
            "tPSA": max(parent_values["tPSA"], self.data.filter_params["tPSA"]),
            "rot_bonds": max(parent_values["rot_bonds"], self.data.filter_params["rot_bonds"]),
            "rigid_bonds": max(parent_values["rigid_bonds"], self.data.filter_params["rigid_bonds"]),
            "num_rings": max(parent_values["num_rings"], self.data.filter_params["num_rings"]),
            "max_ring_size": max(parent_values["max_ring_size"], self.data.filter_params["max_ring_size"]),
            "num_charges": max(parent_values["num_charges"], self.data.filter_params["num_charges"]),
            "num_chiral_centers": max(parent_values["num_chiral_centers"],
                                      self.data.filter_params["num_chiral_centers"])

        }

        self.data.filter_params = filter_params

    def derive_selfies(self, n_children: int = 100, mut_rate: float = 0.03, mut_min: int = 1, mut_max: int = 2):

        good_children = []
        all_filtered_children = {}
        self.data.all_good_selfies_children = []
        n_seeds = len(self.data.seed_smiles)
        if n_children < n_seeds:
            n_children_per_seed = 1
        else:
            n_children_per_seed = round(n_children / n_seeds) + 1

        logger.info(f"Mutating SELFIES to create {n_children_per_seed} children per seed.")

        if self.data.filter:
            filter_params = self.data.filter_params
        else:
            logger.warning("Warning: No filter has been set, so all child molecules will be labeled"
                           " as 'good' regardless of quality. Please call Deriver.set_filter() first"
                           " in order to use a filter for drug-likeness.")
            filter_params = None

        for seed in self.data.seed_smiles:
            children = selfies_substitution(parent_smiles=seed,
                                                n_children=round(n_children_per_seed * 0.7),  # three internal methods
                                                mut_rate=mut_rate,
                                                mut_min=mut_min,
                                                mut_max=mut_max)
            self.data.heritage[seed] += children
            child_mols = [Chem.MolFromSmiles(child, sanitize=True) for child in children]

            children = selfies_insertion(parent_smiles=seed,
                                                n_children=round(n_children_per_seed * 0.15),  # three internal methods
                                                mut_rate=mut_rate,
                                                mut_min=mut_min,
                                                mut_max=mut_max)
            self.data.heritage[seed] += children
            child_mols += [Chem.MolFromSmiles(child, sanitize=True) for child in children]

            children = selfies_deletion(parent_smiles=seed,
                                                n_children=round(n_children_per_seed * 0.15),  # three internal methods
                                                mut_rate=mut_rate,
                                                mut_min=mut_min,
                                                mut_max=mut_max)
            self.data.heritage[seed] += children
            child_mols += [Chem.MolFromSmiles(child, sanitize=True) for child in children]

            # filter children
            filtered_children = apply_filter(filter_params, child_mols, self.data.must_have_patterns)
            all_filtered_children.update(filtered_children)

            for child in filtered_children:
                if filtered_children[child]["is_good"]:
                    # check the cache
                    if self.data.filter_molecules:
                        if child not in self.data.filter_molecules:
                            good_children.append(child)
                            self.data.all_good_selfies_children.append(child)
                        else:
                            logger.debug(f"skipping previously seen molecule: {child}")
                    else:
                        good_children.append(child)
                        self.data.all_good_selfies_children.append(child)

        logger.info(f"Generated {len(good_children)} 'good' children.")

        return good_children, all_filtered_children

    def random_selfies(self, *, n_symbols: int = 100, n_molecules: int = 100):

        good_children = []
        rand_selfies_gen = random_selfies_generator(n_symbols=n_symbols)
        self.data.all_good_selfies_children = []

        if self.data.filter:
            filter_params = self.data.filter_params
        else:
            logger.warning("Warning: No filter has been set, so all child molecules will be labeled"
                           " as 'good' regardless of quality. Please call Deriver.set_filter() first"
                           " in order to use a filter for drug-likeness.")
            filter_params = None

        while len(good_children) < n_molecules:
            child_mols = [Chem.MolFromSmiles(next(rand_selfies_gen)) for i in range(n_molecules - len(good_children))]
            # filter children
            filtered_children = apply_filter(filter_params, child_mols, self.data.must_have_patterns)

            for child in filtered_children:
                if filtered_children[child]["is_good"]:
                    # check the cache
                    if self.data.filter_molecules:
                        if child not in self.data.filter_molecules:
                            good_children.append(child)
                            self.data.all_good_selfies_children.append(child)
                        else:
                            logger.debug(f"skipping previously seen molecule: {child}")
                    else:
                        good_children.append(child)
                        self.data.all_good_selfies_children.append(child)

        return good_children

    def scan_selfies(self):

        """
        Return all possible single substitution children for all the seeds.
        """
        if self.data.filter:
            filter_params = self.data.filter_params
        else:
            logger.warning("Warning: No filter has been set, so all child molecules will be labeled"
                           " as 'good' regardless of quality. Please call Deriver.set_filter() first"
                           " in order to use a filter for drug-likeness.")
            filter_params = None
        good_children = []
        self.data.all_good_scanner_children = []
        all_filtered_children = {}

        for seed in self.data.seed_smiles:
            children = selfies_scanner(parent_smiles=seed)
            self.data.heritage[seed] += children

            filtered_children = apply_filter(filter_params,
                                             [Chem.MolFromSmiles(child) for child in children],
                                             self.data.must_have_patterns)
            all_filtered_children.update(filtered_children)

            for child in filtered_children:
                if filtered_children[child]["is_good"]:
                    # check the cache
                    if self.data.filter_molecules:
                        if child not in self.data.filter_molecules:
                            good_children.append(child)
                            self.data.all_good_scanner_children.append(child)
                        else:
                            logger.debug(f"skipping previously seen molecule: {child}")
                    else:
                        good_children.append(child)
                        self.data.all_good_scanner_children.append(child)

        return good_children, all_filtered_children
