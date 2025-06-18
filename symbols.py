
"""Symbols for variables in the symbolic point-to-ellipse series.
"""

# External include.
import sympy as sp

# varrho = a / rho, where 'a' is the ellipse semi-major axis and rho is the radial distance.
varrho = sp.symbols('varrho', Real=True, Positive=True)
# psi = atan(y,x), i.e. the polar angle.
psi = sp.symbols('psi', Real=True)
# sin(psi) and cos(psi).
sin_psi, cos_psi = sp.symbols('sin_psi, cos_psi', Real=True)
# Ellipse eccentricity squared.
e2 = sp.symbols('e2', Real=True, Positive=True)
