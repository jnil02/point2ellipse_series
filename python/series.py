
"""Series arithmetic classes.
"""

import abc
from collections.abc import Callable
import sympy as sp


class SeriesBase(abc.ABC):
    """Abstract base class for lazy series computations.

    This class and the typing could really be removed, but it has been
    introduced, at least partially, to support in C++ conversion.
    """
    @abc.abstractmethod
    def __getitem__(self, n: int) -> sp.core.Expr:
        """Compute the n:th series item.
        :param n: Series index.
        :return: n:th series item.
        """
        pass

    @abc.abstractmethod
    def __mul__(self, other: 'SeriesBase | sp.core.Expr') -> 'SeriesBase':
        """Multiply series with series or factor.
        :param other: Factor or series.
        :return: resulting series.
        """
        pass

    @abc.abstractmethod
    def __add__(self, other: 'SeriesBase') -> 'SeriesBase':
        """Add series to series.
        :param other: Other series.
        :return: resulting series.
        """
        pass


class Series(SeriesBase):
    """Class for representing a series and implement related series arithmetics.
    """

    def __init__(self, gen: Callable[[int], sp.core.Expr]):
        self.gen = gen  # Function taking an integer and returning coefficient.
        # Add 16 items in the cache to start with.
        self.items = [None] * 16

    def __getitem__(self, n: int) -> sp.core.Expr:
        # Double length of items cache while too short.
        while len(self.items) <= n:
            self.items = self.items + [None]*(len(self.items))
        # Compute and cache item if not cached.
        if self.items[n] is None:
            self.items[n] = self.gen(n)
        # Return cached item.
        return self.items[n]

    def __mul__(self, other: SeriesBase | sp.core.Expr) -> SeriesBase:
        if isinstance(other, Series):
            # Cauchy product.
            return Series(lambda n: sum([self[l] * other[n - l] for l in range(n + 1)]))
        elif isinstance(other, SeriesFactor):
            return Series(lambda n: self[n] * other[0])
        else:
            return Series(lambda n: self[n] * other)

    def __add__(self, other: SeriesBase) -> SeriesBase:
        return Series(lambda n: self[n] + other[n])


class SeriesFactor(SeriesBase):
    """Constant factor of a series as of The Series class.
    """

    def __init__(self, a: sp.core.Expr):
        self.a = a

    def __getitem__(self, n: int) -> sp.core.Expr:
        # A factor is also a series with only a constant term.
        if n == 0:
            return self.a
        return sp.Integer(0)

    def __mul__(self, other: SeriesBase | sp.core.Expr) -> SeriesBase:
        if isinstance(other, Series):
            return other * self.a
        elif isinstance(other, SeriesFactor):
            return SeriesFactor(self.a * other.a)
        else:
            return SeriesFactor(self.a * other)

    def __add__(self, other: SeriesBase) -> SeriesBase:
        # We could support this (see __getitem__) but it does not make sense.
        raise Exception("Cannot add to SeriesFactor.")

class SeriesEmpty(SeriesBase):
    """ Empty series object with the only supported operation being addition.
    """

    def __init__(self):
        pass

    def __getitem__(self, n: int):
        raise Exception("No series item available for empty series.")

    def __mul__(self, other: SeriesBase | sp.core.Expr) -> SeriesBase:
        raise Exception("Multiplication not defined for empty series.")

    def __add__(self, other: SeriesBase) -> SeriesBase:
        return other
