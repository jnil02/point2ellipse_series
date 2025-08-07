
"""Ellipse parameters and Ellipse to Cartesian coordinate transformation.
"""
import typing

# External include
from mpmath import mp

# Defining WGS84 constants.
mp_a = mp.mpf(6378137.0)  # Semi-major axis / Earth equatorial radius.
mp_f = mp.mpf(1.0) / mp.mpf(298.257223563)  # Flattening. f = (a - b) / a
# Derived WGS84 constants.
mp_b = mp_a - mp_f * mp_a  # Semi-minor axis / Earth polar radius.
mp_e2 = mp.mpf(1.0) - (mp_b * mp_b) / (mp_a * mp_a)  # First eccentricity squared.

def mp_ellipse_to_cartesian(phi : mp.mpf, h : mp.mpf) -> typing.Tuple[mp.mpf, mp.mpf]:
    """Multi-precision elliptical to Cartesian coordinate transformation.

    :param phi: Angle between semi-major axis and normal at foot point.
    :param h: Distance from foot point.
    :return: Cartesian coordinates (x,y).
    """
    sin_lat = mp.sin(phi)
    N = mp_a / mp.sqrt(mp.mpf(1.0) - mp_e2 * sin_lat * sin_lat)
    x = (N + h) * mp.cos(phi)
    y = ((mp.mpf(1.0) - mp_e2) * N + h) * sin_lat
    return x, y
