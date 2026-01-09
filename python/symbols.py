
"""Symbols for variables in the symbolic point-to-ellipse series.
"""

# External include.
import sympy as sp

# varrho = a / rho, where 'a' is the ellipse semi-major axis and rho is the radial distance.
# rho_ae2 = rho / (ae²)
varrho, rho_ae2 = sp.symbols('varrho, rho_ae2', Real=True, Positive=True)
# psi = atan(y,x), i.e. the polar angle.
psi = sp.symbols('psi', Real=True)
# sin(psi) and cos(psi).
sin_psi, cos_psi = sp.symbols('sin_psi, cos_psi', Real=True)
# Ellipse eccentricity squared e2=1-b²/a². b_a = b / a.
e2, b_a = sp.symbols('e2, b_a', Real=True, Positive=True)
