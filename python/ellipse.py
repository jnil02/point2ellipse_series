
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

def mp_ellipse_to_cartesian(phi: mp.mpf, h: mp.mpf) -> typing.Tuple[mp.mpf, mp.mpf]:
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

def mp_polar_to_cartesian(psi: mp.mpf, rho: mp.mpf) -> typing.Tuple[mp.mpf, mp.mpf]:
    """ Conversion from polar to cartesian coordinates.

    :param psi: Polar angle.
    :param rho: Polar distance.
    :return: Cartesian coordinate (x,y)
    """
    return rho * mp.cos(psi), rho * mp.sin(psi)

# Constants for transformation function below.
a2 = mp_a * mp_a
b2 = mp_b * mp_b
e4 = mp_e2 * mp_e2
b2am2 = b2 / a2
a3 = a2 * mp_a
a4 = a2 * a2
am2 = mp.mpf(1.) / a2
bam2 = mp_b / a2
b2am4 = b2 / a4
inv_6 = mp.mpf(1.) / mp.mpf(6.)
e2_2 = mp_e2 / mp.mpf(2)
e2bam3 = mp_e2 * mp_b / a3
M_PI_6 = mp.pi / mp.mpf(6.)
M_2_3 = mp.mpf(2.) / mp.mpf(3.)
bam1 = mp_b / mp_a  # = 1 - f
bem1 = mp_b / mp.sqrt(mp_e2)
e2bm1 = mp_e2 / mp_b
c = mp.sqrt(mp.mpf(2)) - mp.mpf(1.)

def mp_cartesian_to_ellipse(x, z):
    """Multiprecision Cartesian (ECEF) to geodetic (ellipse) coordinate transformation by Vermeille.

    For details, see:
    Vermeille, H, An analytical method to transform geocentric into geodetic
    coorindates. J. Geod. (2011) 85:105-117. DOI 10.1007/s00190-010-0419-x

    :param x: Cartesian x coordinate.
    :param z: Cartesian z coordinate.
    :return: Tuple of latitude and altitude.
    """
    t2 = x * x
    t = mp.sqrt(t2)
    p = t2 * am2
    q = z * z * b2am4
    r = (p + q - e4) * inv_6
    r38 = r * r * r * mp.mpf(8.)
    ev = r38 + e4 * p * q
    
    if (ev > 0):
        # Outside evolute.
        s = e2bam3 * mp.fabs(z) * t  # std::sqrt(e4*p*q)
        sqrt_ev = mp.sqrt(ev)
        c1 = sqrt_ev + s
        c2 = sqrt_ev - s
        u = r + mp.mpf(0.5) * mp.cbrt(c1 * c1) + mp.mpf(0.5) * mp.cbrt(c2 * c2)

        v = mp.sqrt(u * u + e4 * q)
        w = e2_2 * (u + v - q) / v
        k = (u + v) / (mp.sqrt(w * w + u + v) + w)
        D = k / (k + mp_e2) * t
        d = mp.sqrt(D * D + z * z)

        h = (k - b2am2) * d / k
        lat = 2. * mp.atan(z / (d + D))
        return lat, h
    elif (q != 0):
        # On or inside evolute and not on singular disc.
        s = e2bam3 * mp.fabs(z) * t  # std::sqrt(e4*p*q)
        up = M_2_3 * mp.atan(s / (mp.sqrt(-ev) + mp.sqrt(-r38)))
        u = mp.mpf(-4.) * r * mp.sin(up) * mp.cos(M_PI_6 + up)

        v = mp.sqrt(u * u + e4 * q)
        w = e2_2 * (u + v - q) / v
        k = (u + v) / (mp.sqrt(w * w + u + v) + w)
        D = k * t / (k + mp_e2)
        d = mp.sqrt(D * D + z * z)

        h = (k - b2am2) * d / k
        lat = 2. * mp.atan(z / (d + D))
        return lat, h
    else:
        # On the singular disc (including center of earth).
        # Values are taken to have positive latitude.
        h = -bem1 * mp.sqrt(mp_e2 - p)
        lat = mp.mpf(2.) * mp.atan(mp.sqrt(e4 - p) / (-e2bm1 * h + bam1 * mp.sqrt(p)))
        return lat, h
