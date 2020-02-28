# coding=utf-8
#
# Wrappers for query execution, with different levels of caching
#
# script by Andreas Windemuth

from typing import Callable, Any
import numpy as np
from .util import Vec, percent
from .fragment_queries import FragArray, MateQueryParams


class CacheL2(object):
    """Level 2 caching"""

    def __init__(self, fragments: FragArray, query):
        self.fragments = fragments
        self.query = query
        self.cache = Cache(self.l1_calc)

    def l1_calc(self, params):
        (ix,) = np.where(self.query.clause(self.fragments, params))
        pseudo = self.fragments.pseudo[ix]
        cache2 = Cache(self.l2_calc)
        return ix, pseudo, cache2

    def l2_calc(self, data):
        (pseudo, mask) = data
        (iy,) = np.where((pseudo & mask) != 0)
        return iy

    def select(self, params: MateQueryParams) -> Vec[int]:
        key1 = self.query.key(params)
        ix, pseudo, cache2 = self.cache.get(key1, params)
        key2 = params.mask
        iy = cache2.get(key2, (pseudo, params.mask))
        return ix[iy]

    def info(self):
        nh = 0
        nm = 0
        for _, _, c in self.cache.store.values():
            nh += c.n_hits
            nm += c.n_misses
        if nh + nm > 0:
            s = " %.2f%%" % (100.0 * nm / (nh + nm))
        else:
            s = ""
        return self.cache.info() + s


class Cache(object):
    """A basic general cache object. May put this into base/util later"""

    def __init__(self, calc: Callable[[Any], Any] = None):
        """
        Initializa a cache.
        :param calc: Function that will compute y from x
        """
        self.calc = calc
        self.store = {}
        self.n_hits = 0
        self.n_misses = 0
        self.size = 0

    def get(self, key, x=None):
        """Retrieve from cache if found, else compute and store"""
        if key in self.store:
            self.n_hits += 1
            y = self.store[key]
        else:
            self.n_misses += 1
            if self.calc is None:
                y = None
            else:
                y = self.calc(x)
                self.store[key] = y
        return y

    def set(self, key, y):
        self.store[key] = y

    def info(self):
        return percent(self.n_misses, self.n_misses + self.n_hits)
