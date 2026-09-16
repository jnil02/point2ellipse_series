"""Generate test data for C++ consistency checks.

Outputs CSV files to test_data/ with columns: n,k,l,num,den
"""

import os
import csv

from dataclasses import dataclass
from typing import Callable, Optional, Sequence, Tuple

from coefficients import (d_phi, d_sin, d_cos, d_h, d_phi_evo2, c_phi_evo,
                          c_phi_pow_evo, c_sin_phi_evo, d_sin_phi_evo,
                          c_cos_phi_evo, d_cos_phi_evo, c_sin_phi_inv_evo,
                          c_N_evo, cp_evo_nkl, cp_evo_nkl2, c_h_evo2)

# Generate all indices up to and including this value.
M  = 5  # Max index
MP = 3  # Max index for power index. Lowever due to computational cost.

# Where to place the test data.
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'test_data')


def write_csv(filename, rows, header):
    """Write rows to a CSV file in TEST_DATA_DIR."""
    os.makedirs(TEST_DATA_DIR, exist_ok=True)
    path = os.path.join(TEST_DATA_DIR, filename)
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"Written {len(rows)} rows to {path}")


# (Coefficient) spec of the nested series.
@dataclass
class Spec:
    # Coefficient dependent on loop index variables.
    coef:  Callable
    # The nested for loops of the nested series.
    # Each loop is specified as (<name>, lo, hi) where lo, hi is a function
    # taking the proceeding loop variables as arguments, e.g. lambda k: k + 1.
    # The bounds of the first loop are just lambdas with no arguments returning
    # constants.
    loops: Sequence[Tuple]
    # Order in which the coefficient take the index variables, specified with a
    # list of the loop variable names. This cannot just be deferred to the
    # coefficient callable, e.g. rearranged via a lambda, since the order is
    # the same on the C-side and there we cannot rearrange them. None is
    # converted to the natural loop order.
    order: Optional[Tuple[str, ...]] = None   # call-arg & row order; default = loop order


def generate_csv(spec: Spec):
    names = [n for n, _, _ in spec.loops]
    order = spec.order or tuple(names)  # None means loop order.

    # Build all combinations of the nested for loops like:
    # for k in range(<constant>, <constant>):
    #   for l in range(lo(k), hi(k)):
    #       for n in range(lo(k,l), hi(k,l)):
    #           <coefficient>(k,l,n)
    # Each loop below corresponds to one for loop level.
    # Each loop below builds the combination of that level given the combination of the proceeding (outer) levels.
    # Example for loops (n, k, l):
    #   ()  ->  (n,)  ->  (n, k)  ->  (n, k, l)
    combos = [()]
    for _name, lo, hi in spec.loops:
        combos = [outer + (value,)
                  for outer in combos
                  for value in range(lo(*outer), hi(*outer))]

    # Rearrange loop variables contained in each combo to that of Spec.order.
    # Compute the corresponding coefficient and append to a list of
    # coefficients and indices.
    pos = {name: i for i, name in enumerate(names)}
    coefs = []
    for combo in combos:
        args = [combo[pos[n]] for n in order]
        c = spec.coef(*args)
        coefs.append(tuple(args) + (c.p, c.q))

    # Export to csv for C-side test consumption.
    write_csv(spec.coef.__name__ + '.csv', coefs, list(order) + ['num', 'den'])


SPECS = [
    Spec(d_phi, [('n', lambda: 0, lambda: M + 1),
                 ('k', lambda n: 1, lambda n: M + 1),
                 ('l', lambda n, k: max(n + 1, k), lambda n, k: n + k + 1)]),
    Spec(d_cos, [('n', lambda: 0, lambda: M + 1),
                 ('k', lambda n: 1, lambda n: M + 1),
                 ('l', lambda n, k: max(n, k), lambda n, k: n + k)]),
    Spec(d_sin, [('n', lambda: 0, lambda: M + 1),
                 ('k', lambda n: 1, lambda n: M + 1),
                 ('l', lambda n, k: max(n, k), lambda n, k: n + k + 1)]),
    Spec(d_h, [('n', lambda: 1, lambda: M + 1),
               ('k', lambda n: 0, lambda n: M + 1),
               ('l', lambda n, k: max(n, k + 1), lambda n, k: n + k + 1)]),

    Spec(d_phi_evo2,
         [('k', lambda: 0, lambda: M + 1),
          ('l', lambda k: 0, lambda k: (k-1)//2 + 1),
          ('n', lambda k, l: 0, lambda k, l: (k-1) // 2 + 1)]),
    Spec(c_phi_evo,
         [('l', lambda: 0, lambda: M + 1),
          ('k', lambda l: l + 1, lambda l: M + 1),
          ('n', lambda l, k: 1, lambda l, k: k + 1)],
         order=('k', 'l', 'n')),

    Spec(c_phi_pow_evo,
         [('i', lambda: 0, lambda: MP + 1),
          ('l', lambda i: 0, lambda i: M + 1),
          ('k', lambda i, l: l + i, lambda i, l: M + 1),
          ('n', lambda i, l, k: i, lambda i, l, k: k + 1)],
         order=('k', 'l', 'n', 'i')),
    Spec(c_sin_phi_evo,
         [('l', lambda: 0, lambda: M + 1),
          ('k', lambda l: l, lambda l: M + 1),
          ('n', lambda l, k: 0, lambda l, k: k + 1)],
         order=('k', 'l', 'n')),
    Spec(d_sin_phi_evo,
         [('k', lambda: 0, lambda: M + 1),
          ('l', lambda k: 0, lambda k: M + 1),
          ('n', lambda k, l: 1, lambda k, l: l // 2 + k + 1)]),
    Spec(c_cos_phi_evo,
         [('l', lambda: 0, lambda: M + 1),
          ('k', lambda l: l, lambda l: M + 1),
          ('n', lambda l, k: 1, lambda l, k: k + 1)],
         order=('k', 'l', 'n')),
    Spec(d_cos_phi_evo,
         [('k', lambda: 0, lambda: M + 1),
          ('l', lambda k: 0, lambda k: M + 1),
          ('n', lambda k, l: l % 2, lambda k, l: (l + 1) // 2 + k + 1)]),
    Spec(c_sin_phi_inv_evo,
         [('l', lambda: 0, lambda: M + 1),
          ('k', lambda l: l, lambda l: M + 1),
          ('n', lambda l, k: 0, lambda l, k: k + 1)],
         order=('k', 'l', 'n')),
    Spec(c_N_evo,
         [('l', lambda: 0, lambda: M + 1),
          ('k', lambda l: l, lambda l: M + 1),
          ('n', lambda l, k: 1, lambda l, k: k + 2)],
         order=('k', 'l', 'n')),
    Spec(cp_evo_nkl,
         [('l', lambda: 1, lambda: M + 1),
          ('k', lambda l: l, lambda l: M + 1),
          ('n', lambda l, k: 0, lambda l, k: k + 2)],
         order=('k', 'l', 'n')),
    Spec(cp_evo_nkl2,
         [('l', lambda: 1, lambda: M + 1),
          ('k', lambda l: l, lambda l: M + 1),
          ('n', lambda l, k: 1, lambda l, k: k + 2)],
         order=('k', 'l', 'n')),
    Spec(c_h_evo2,
         [('l', lambda: 1, lambda: M + 1),
          ('k', lambda l: l, lambda l: M + 1),
          ('n', lambda l, k: 1, lambda l, k: k + 2)],
         order=('k', 'l', 'n')),
]


if __name__ == '__main__':
    for spec in SPECS:
        generate_csv(spec)