# coding=utf-8
#
# Alternate implementations of the cascading mate queries in mate.py, to go with the FragmentIndex..
#
# code by Andreas Windemuth

from typing import NamedTuple
import numpy as np
from .util import Vec


class MateQueryParams(NamedTuple):
    """Collection of all the parameters a query may depend on"""
    permissivity: float
    max_pseudoatoms: int
    missing_piece_fc: int
    missing_piece_len: int
    length_cutoff: int
    mask: int  # bit encoded version of allowed_atoms


class FragArray(NamedTuple):
    """An array of all the information about fragments needed to support mate queries.
       Note that this does not include SMILES, which need to still be read from the database.
       We could/should add smiles, but need to watch memory consumption.
    """
    id: Vec[int]  # fragment id
    slen: Vec[int]  # SMILES length
    coeff: Vec[int]  # frag_coeff (whatever it is...)
    npseudo: Vec[int]  # number of pseudo atoms
    pseudo: Vec[int]  # bitset of pseudo atoms

    @classmethod
    def zeros(cls, n):
        """Provide an empty array of a certain length"""

        def idim(): return np.zeros(n, dtype=int)

        return FragArray(id=idim(), slen=idim(), coeff=idim(), npseudo=idim(), pseudo=idim())

    def truncate(self, n):
        """Apply an index or index range.
           Can't call it __getitem__ as this conflicts with the tuple API and leads to infinite recursion
        """
        return FragArray(
            id=self.id[:n], slen=self.slen[:n], coeff=self.coeff[:n], npseudo=self.npseudo[:n], pseudo=self.pseudo[:n]
        )


# These are the cascading queries, in order of greatest stringency (except one seems superfluous?)

class MateQuery(object):
    def key(self, params: MateQueryParams) -> str:
        return "%d %d %d %d %d %d" % (
            int(params.max_pseudoatoms), int(params.max_pseudoatoms / 3.0),
            int(params.permissivity * params.missing_piece_fc), int(params.max_pseudoatoms / 3.0),
            params.length_cutoff, int(params.length_cutoff / 3.0)
        )

    def clause(self, fragments: FragArray, params: MateQueryParams):
        # logger.info("    npa: %d - %d" % (int(params.max_pseudoatoms / 3.0), int(params.max_pseudoatoms)))
        # logger.info("  coeff: %d - %d" % (int((params.permissivity / 3.0) * params.missing_piece_fc), int(params.permissivity * params.missing_piece_fc)))
        # logger.info("  smile: %d - %d" % (int(params.length_cutoff / 3.0), int(params.length_cutoff)))
        return (
                (fragments.npseudo <= int(params.max_pseudoatoms)) &
                (fragments.npseudo >= int(params.max_pseudoatoms / 3.0)) &
                (fragments.coeff <= int(params.permissivity * params.missing_piece_fc)) &
                (fragments.coeff >= int((params.permissivity / 3.0) * params.missing_piece_fc)) &
                (fragments.slen <= params.length_cutoff) &
                (fragments.slen >= int(params.length_cutoff / 3.0))
        )


class BackupQuery1(object):
    def key(self, params: MateQueryParams) -> str:
        return "%d %d %d" % (
            int(params.max_pseudoatoms),
            int(params.permissivity * params.missing_piece_fc),
            params.length_cutoff
        )

    def clause(self, fragments: FragArray, params: MateQueryParams):
        return (
                (fragments.npseudo <= int(params.max_pseudoatoms)) &
                (fragments.coeff <= int(params.permissivity * params.missing_piece_fc)) &
                (fragments.slen <= params.length_cutoff)
        )


def alt_length_cutoff(params: MateQueryParams):
    return max(20.0, (params.missing_piece_len * 1.5) / params.max_pseudoatoms)


class BackupQuery2(object):
    def key(self, params: MateQueryParams) -> str:
        return "%d %d %d" % (
            int(params.max_pseudoatoms),
            int(params.permissivity * params.missing_piece_fc),
            alt_length_cutoff(params)
        )

    def clause(self, fragments: FragArray, params: MateQueryParams):
        return (
                (fragments.npseudo <= int(params.max_pseudoatoms)) &
                (fragments.coeff <= int(params.permissivity * params.missing_piece_fc)) &
                (fragments.slen <= alt_length_cutoff(params))
        )


# This query appears to be more strict than the one above, and should therefore never be used ???

class BackupQuery3(object):
    def key(self, params: MateQueryParams) -> str:
        return "%d %d" % (int(params.max_pseudoatoms), alt_length_cutoff(params))

    def clause(self, fragments: FragArray, params: MateQueryParams):
        return (
                (fragments.npseudo <= int(params.max_pseudoatoms)) &
                (fragments.coeff == 0) &
                (fragments.slen <= alt_length_cutoff(params))
        )


class BackupQuery4(object):
    def key(self, params: MateQueryParams) -> str:
        return "single"  # clause does not depend on params

    def clause(self, fragments: FragArray, params: MateQueryParams):
        return (
                (fragments.npseudo == 1) & (fragments.coeff == 0)
        )
