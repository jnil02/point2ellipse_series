
"""
Caching utilities for function taking unsigned integer arguments.
"""

# External includes.
from collections.abc import Hashable


def cantor_pairing_two(k: int, l: int) -> int:
    """Cantor pairing for two non-negative integers."""
    assert k >= 0, "k must be >= 0"
    assert l >= 0, "l must be >= 0"
    return (k + l) * (k + l + 1) // 2 + l

def cantor_pairing(*args: int) -> int:
    """Chained Cantor pairing for an arbitrary number of non-negative integers."""
    if not args:
        raise ValueError("At least one integer required")
    assert all(isinstance(a, int) and a >= 0 for a in args), \
        "All arguments must be non-negative integers"
    key = args[0]
    for next_val in args[1:]:
        key = cantor_pairing_two(key, next_val)
    return key

def split_args(*args):
    """Split *args into (ints, others) where 'others' are hashable non-ints."""
    ints = []
    others = []
    for a in args:
        if isinstance(a, int):
            ints.append(a)
        elif isinstance(a, Hashable):
            others.append(a)
        else:
            raise TypeError(f"Argument {a!r} is not hashable")
    return ints, tuple(others)

class ints_cache:
    """Caches results of integer-arg + other hashable arg functions."""
    def __init__(self, f):
        self.f = f
        self.cache = {}

    def __call__(self, *args):
        ints, others = split_args(*args)
        if ints:
            int_key = cantor_pairing(*ints)
            key = (int_key,) + others
        else:
            key = others
        if key in self.cache:
            return self.cache[key]
        res = self.f(*args)
        self.cache[key] = res
        return res
