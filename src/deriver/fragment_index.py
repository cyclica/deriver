# coding=utf-8
#
# Provide an index into the fragment library that finds fragments to mate with much faster than
# a database query. The problem with the database query is the "ORDER BY RANDOM() LIMIT" clause,
# which steadfastly resists query optimization.
#
# script by Andreas Windemuth

import os
import pickle
from typing import Set
import numpy as np
from peewee import SqliteDatabase
from .config import logger
from .util import lazy
from .fragment_caches import percent, CacheL2
from .fragment_queries import BackupQuery1, BackupQuery2, \
    BackupQuery3, BackupQuery4, MateQuery, FragArray, MateQueryParams
from .lib_read import lib_read

# This defines the backup query cascade. It should be possible to add and remove queries or rearrange them,
#   except after that comparison with the old mate query will no longer be possible after that
query_cascade = [MateQuery(), BackupQuery1(), BackupQuery2(), BackupQuery3(), BackupQuery4()]

# Singletons for fragment index
frag_index_cache = {}


def frag_index(db: str):
    global frag_index_cache  # pylint: disable=global-statement
    name = os.path.splitext(db)[0]
    if name not in frag_index_cache:
        fx = FragmentIndex.by_name(name)
        frag_index_cache[name] = fx
    else:
        fx = frag_index_cache[name]
    return fx


class FragmentIndex(object):
    """
    An index into the fragment base to make mate queries more efficient.
    """

    def __init__(self, frags: FragArray, checking: bool = False):
        self.id_map_value = {}
        self.fragments = frags
        self.id_map_calc = None
        if np.any(frags.id == 0):
            (ix,) = np.where(frags.id == 0)
            raise Exception("Missing fragments: %d " % len(ix) + str(ix))

        # set the cache up
        self.clear_cache()

    @classmethod
    def from_db(cls, db_name: str, limit: int = 0):
        """Read all the fragments from a fragment database"""
        lib = SqliteDatabase(db_name)
        Fragment, _, PseudoAtoms, _ = lib_read(lib)
        lib.connect()
        logger.info("Get count...")
        n_frags = len(Fragment.select())
        if limit > 0:
            n_frags = min(n_frags, limit)
        logger.info("  %d fragments. Loading...." % n_frags)
        fragments = FragArray.zeros(n_frags)
        idmap = {}
        n = 0
        for pa in PseudoAtoms.select().join(Fragment).iterator():
            frag = pa.frag
            ident = int(frag.id)
            k = int(pa.pseudo_atom)
            if ident in idmap:
                i = idmap[ident]
            else:
                i = len(idmap)
                if limit > 0 and i >= limit: break
                idmap[ident] = i
                fragments.id[i] = ident
                fragments.slen[i] = len(frag.smile)
                fragments.coeff[i] = int(frag.frag_coeff)
                fragments.npseudo[i] = int(frag.num_pseudo_atoms)
                if i % 100000 == 0:
                    logger.info(
                        "  %d " % n + percent(i, n_frags) + "  " + "%7d %7d %2d %8d  " % (ident, i, k, int(2 ** k)))
            fragments.pseudo[i] |= int(2 ** k)
            n += 1
        nf = len(idmap)
        flib = FragmentIndex(fragments.truncate(nf))
        flib.id_map_calc = idmap
        logger.info("Read %d/%d fragments, %d pseudo atoms." % (nf, n_frags, n))
        return flib

    @classmethod
    def by_name(cls, name: str, limit: int = 0):
        """Read all the fragments from a binary file, or create on from the fragment database"""
        if limit > 0:
            mod = "-%d" % limit
        else:
            mod = ""
        dbn = name + ".db"
        bfn = name + mod + ".bin"
        if os.path.exists(bfn):
            logger.info("Reading fragments from " + bfn)
            with open(bfn, 'rb') as f:
                fragments = pickle.load(f)
            return FragmentIndex(fragments)
        else:
            logger.info("Loading fragments from " + dbn)
            flib = FragmentIndex.from_db(dbn, limit)
            logger.info("Saving fragments to " + bfn)
            with open(bfn, 'wb') as f:
                pickle.dump(flib.fragments, f, protocol=pickle.HIGHEST_PROTOCOL)
            return flib

    @lazy
    def id_map(self):
        """A mapping to get the index of a fragment from its id"""
        if self.id_map_calc is None:
            self.id_map_value = {}
            for i, ident in enumerate(self.fragments.id):
                self.id_map_calc[ident] = i
        return self.id_map_calc

    def mate_query(self,
                   permissivity: float, max_pseudoatoms: int,
                   missing_piece_fc: int, missing_piece_len,
                   length_cutoff: int, allowed_atoms: Set[int],
                   limit: int
                   ):
        """
        Efficient mate query for deriver.
        :param permissivity:
        :param max_pseudoatoms:
        :param missing_piece_fc:
        :param missing_piece_len:
        :param length_cutoff:
        :param allowed_atoms:
        :param limit:
        :return:  np.array of fragment ids
        """

        # Compute the bit mask for the allowed pseudo atoms
        mask = 0
        for a in allowed_atoms:
            mask |= 2 ** a  # every pseudoatom type corresponds to a position in the bit vector

        params = MateQueryParams(
            permissivity, max_pseudoatoms,
            missing_piece_fc, missing_piece_len,
            length_cutoff, mask
        )

        def cascade(qs) -> np.ndarray:
            ixl: np.ndarray = None
            for k in range(len(qs)):
                ixl = qs[k].select(params)
                if len(ixl) > 0: break
            return ixl

        ix = cascade(self.queries)

        # Find random subset
        np.random.shuffle(ix)
        if limit > 0:
            return self.fragments.id[ix[:limit]]
        else:
            return self.fragments.id[ix]

    def clear_cache(self):
        """clear the cache"""
        self.queries = [CacheL2(self.fragments, q) for q in query_cascade]
