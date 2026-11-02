
"""Symbols for variables in the symbolic point-to-ellipse series.
"""

# External include.
import sympy as sp

# varrho = a / rho, where 'a' is the ellipse semi-major axis and rho is the radial distance.
# sigma = rho / (ae²)
rho, varrho, sigma = sp.symbols('rho, varrho, sigma', Real=True, Positive=True)
# psi = atan(y,x), i.e. the polar angle.
psi = sp.symbols('psi', Real=True)
# sin(psi) and cos(psi).
sin_psi, cos_psi = sp.symbols('sin_psi, cos_psi', Real=True)
# Ellipse eccentricity squared e2=1-b²/a². epsilon = b / a.
e2, epsilon = sp.symbols('e2, epsilon', Real=True, Positive=True)
