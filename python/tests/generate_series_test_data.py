"""
Generates test data for series to ensure consistent definitions on the C++ side.
"""

import math
from symbols import sin_psi, rho_ae2, b_a
import fourier_series as fs
from generate_test_data import write_csv

A, B   = 1.0, 0.5
E2     = 1 - (B/A)**2
B_A    = B/A
ORDER  = 8
PSIS   = [44, 136, 224, 316]          # One per quadrant.

def evolute_rho(psi):                 # same formula as sweep_evo / plot
    ac = abs(A*math.cos(psi)); bs = abs(B*math.sin(psi))
    return (A**2 - B**2) / (ac**(2/3) + bs**(2/3))**1.5

# test points -> exact evo-series inputs
POINTS = []
for pdeg in PSIS:
    psi = math.radians(pdeg)
    rho = 0.5 * evolute_rho(psi)
    POINTS.append((math.sin(psi), rho/(A*E2), B_A))   # (sin_psi, rho_ae2, b_a)

# series under test: (filename, sympy callable returning expr in the 3 symbols)
SERIES = [
    ('phi_evo_sparse',      lambda: fs.phi_evo_sparse(ORDER, ORDER)),
    ('phi_evo_dense',       lambda: fs.phi_evo_dense(ORDER)),
    ('sin_phi_evo_sparse',  lambda: fs.sin_phi_evo_sparse(ORDER, ORDER)),
    ('sin_phi_evo_dense',   lambda: fs.sin_phi_evo_dense(ORDER)),
    ('cos_phi_evo_sparse',  lambda: fs.cos_phi_evo_sparse(ORDER, ORDER)),
    ('cos_phi_evo_dense',   lambda: fs.cos_phi_evo_dense(ORDER)),
    ('h_evo_sparse',        lambda: fs.h_evo_sparse(ORDER, ORDER)),
    ('h_evo_dense',         lambda: fs.h_evo_dense(ORDER)),
]

for name, make in SERIES:
    expr = make()
    rows = []
    for s, r, ba in POINTS:
        val = float(expr.subs({sin_psi: s, rho_ae2: r, b_a: ba}))
        rows.append((ORDER, repr(s), repr(r), repr(ba), repr(val)))
    write_csv(name + '_series.csv', rows,
              ['order', 'sin_psi', 'rho_ae2', 'b_a', 'value'])