from typing import TypeVar, Union, Sequence
import numpy as np
from rdkit import Chem

# Vector and matrix types
T = TypeVar('T')
Mat = Union[np.ndarray, Union[T]]
Vec = Union[np.ndarray, Union[T]]


def percent(k: int, n: int) -> str:
    """Format a fraction with percentage"""
    if n == 0:
        pc = ""
    else:
        pc = " (%.1f%%)" % (100.0 * k / n)
    return "%d/%d" % (k, n) + pc


# Decorators
class lazy(object):
    """
    meant to be used for lazy evaluation of an object attribute.
    property should represent non-mutable data, as it replaces itself.
    Usage:
        class Foo:
            @lazy
            def longish_computation(self):
                ... longish computation
                return result

        foo = Foo()
        print("Result: "+foo.longish_computation)  # This will take some time
        print("Result: "+foo.longish_computation)  # This won't
    """

    def __init__(self, fget):
        self.fget = fget
        self.func_name = fget.__name__

    def __get__(self, obj, cls):
        if obj is None:
            return None
        value = self.fget(obj)
        setattr(obj, self.func_name, value)
        return value
