"""Pytest tests for the far-field point-to-ellipse series evaluation.

Uses the WGS84 Earth ellipse with one test point per quadrant (lat
+-43.1 / +-136.9 deg, alt 10 000 m). All assertions compare a truncated
series (order MAX_ORDER) against the multi-precision reference value, so the
tolerance reflects truncation error, not floating-point noise.
"""

from mpmath import mp
mp.dps = 50

import pytest

from pytest_util import assert_close
import ellipse
from ellipse import mp_ellipse_to_cartesian, mp_e2, mp_a
from symbols import varrho, psi, sin_psi, cos_psi, e2
import expansions

# Order of the tested series. The convergence appear to have geometric convergence.
MAX_ORDER = 7
# Truncated-series tolerance: e2 ≈ 0.007, so e2^8 ≈ 6e-17; 1e-10 is generous.
TOL = mp.mpf("1e-10")

# One test point per quadrant. phi (the normal direction / geodetic latitude)
# is defined in (-180, 180] deg; |phi| > 90 places the point at x < 0 (the 2nd
# and 3rd quadrants), which exercises the full-circle Vermeille reference.
_PHI_DEG = ["43.1", "136.9", "-136.9", "-43.1"]   # Q1, Q2, Q3, Q4
_H   = mp.mpf("10000.")

# Reference values for the tests, one parametrization per quadrant.
@pytest.fixture(scope="module", params=_PHI_DEG, ids=[f"phi={d}" for d in _PHI_DEG])
def ref(request):
    phi = mp.mpf(request.param) / mp.mpf("180.") * mp.pi
    x, y = mp_ellipse_to_cartesian(phi, _H)
    rho      = mp.sqrt(x*x + y*y)
    psi_val  = mp.atan2(y, x)
    varrho_v = mp_a / rho
    return dict(phi=phi, h=_H, x=x, y=y, psi=psi_val, rho=rho, varrho=varrho_v)


def ev(expr, ref):
    """Substitute all numeric reference values into a symbolic series expression.

    Symbols absent from expr are silently ignored by sympy's subs.
    Returns an mpmath mpf preserving full precision.
    """
    result = (expr
              .subs(e2,      mp_e2)
              .subs(sin_psi, mp.sin(ref["psi"]))
              .subs(cos_psi, mp.cos(ref["psi"]))
              .subs(psi,     ref["psi"])
              .subs(varrho,  ref["varrho"]))
    return mp.mpf(result._mpf_)


# ---------------------------------------------------------------------------
# Ellipse <-> Cartesian round-trip
# ---------------------------------------------------------------------------

def test_cartesian_to_ellipse_roundtrip(ref):
    phi2, h2 = ellipse.mp_cartesian_to_ellipse(ref["x"], ref["y"])
    assert_close("roundtrip phi [rad]", ref["phi"], phi2,   mp.mpf("1e-40"))
    assert_close("roundtrip h [m]",     ref["h"],   h2,     mp.mpf("1e-40"))


# ---------------------------------------------------------------------------
# phi - psi
# ---------------------------------------------------------------------------

def test_phi_minus_psi_sin_pow(ref):
    expected = ref["phi"] - ref["psi"]
    result = ev(expansions.phi_in_sin_pow(MAX_ORDER, MAX_ORDER) * sin_psi * cos_psi, ref)
    assert_close("phi-psi  sin_pow", expected, result, TOL)


def test_phi_minus_psi_sin_mul(ref):
    expected = ref["phi"] - ref["psi"]
    result = ev(expansions.phi_in_sin_mul(MAX_ORDER, MAX_ORDER, MAX_ORDER), ref)
    assert_close("phi-psi  sin_mul", expected, result, TOL)


def test_phi_minus_psi_sin_pow2(ref):
    # phi_in_sin_pow2 absorbs cos(psi) as its positive square root sqrt(1-sin^2),
    # so it only represents the right half-plane (cos(psi) >= 0).
    if mp.cos(ref["psi"]) < 0:
        pytest.skip("phi_in_sin_pow2 is valid only for cos(psi) >= 0")
    expected = ref["phi"] - ref["psi"]
    result = ev(expansions.phi_in_sin_pow2(MAX_ORDER, MAX_ORDER), ref)
    # phi_in_sin_pow2 has worse convergence by design (see docstring); ~2e-7 at order 7.
    assert_close("phi-psi  sin_pow2", expected, result, mp.mpf("1e-5"))


# ---------------------------------------------------------------------------
# sin(phi) / sin(psi)
# ---------------------------------------------------------------------------

def test_sin_phi_over_sin_psi_sin_pow(ref):
    expected = mp.sin(ref["phi"]) / mp.sin(ref["psi"])
    result = ev(expansions.sin_phi_in_sin_pow(MAX_ORDER, MAX_ORDER), ref) + 1
    assert_close("sin(phi)/sin(psi)  sin_pow", expected, result, TOL)


def test_sin_phi_over_sin_psi_cos_mul(ref):
    expected = mp.sin(ref["phi"]) / mp.sin(ref["psi"])
    result = ev(expansions.sin_phi_in_cos_mul(MAX_ORDER, MAX_ORDER, MAX_ORDER), ref) + 1
    assert_close("sin(phi)/sin(psi)  cos_mul", expected, result, TOL)


def test_sin_phi_over_sin_psi_sin_pow2(ref):
    expected = mp.sin(ref["phi"]) / mp.sin(ref["psi"])
    result = ev(expansions.sin_phi_in_sin_pow2(MAX_ORDER, MAX_ORDER, MAX_ORDER), ref)
    assert_close("sin(phi)/sin(psi)  sin_pow2", expected, result, TOL)


# ---------------------------------------------------------------------------
# cos(phi) / cos(psi)
# ---------------------------------------------------------------------------

def test_cos_phi_over_cos_psi_sin_pow(ref):
    expected = mp.cos(ref["phi"]) / mp.cos(ref["psi"])
    result = ev(expansions.cos_phi_in_sin_pow(MAX_ORDER, MAX_ORDER), ref) + 1
    assert_close("cos(phi)/cos(psi)  sin_pow", expected, result, TOL)


def test_cos_phi_over_cos_psi_cos_mul(ref):
    expected = mp.cos(ref["phi"]) / mp.cos(ref["psi"])
    result = ev(expansions.cos_phi_in_cos_mul(MAX_ORDER, MAX_ORDER, MAX_ORDER), ref) + 1
    assert_close("cos(phi)/cos(psi)  cos_mul", expected, result, TOL)


def test_cos_phi_over_cos_psi_sin_pow2(ref):
    expected = mp.cos(ref["phi"]) / mp.cos(ref["psi"])
    result = ev(expansions.cos_phi_in_sin_pow2(MAX_ORDER, MAX_ORDER, MAX_ORDER), ref)
    assert_close("cos(phi)/cos(psi)  sin_pow2", expected, result, TOL)


# ---------------------------------------------------------------------------
# (h + a - rho) / a
# ---------------------------------------------------------------------------

def test_h_a_sin_pow(ref):
    expected = (ref["h"] + mp_a - ref["rho"]) / mp_a
    result = ev(expansions.h_in_sin_pow(MAX_ORDER, MAX_ORDER), ref)
    assert_close("(h+a-rho)/a  sin_pow", expected, result, TOL)


def test_h_a_cos_mul(ref):
    expected = (ref["h"] + mp_a - ref["rho"]) / mp_a
    result = ev(expansions.h_in_cos_mul(MAX_ORDER, MAX_ORDER, MAX_ORDER), ref)
    assert_close("(h+a-rho)/a  cos_mul", expected, result, TOL)


# ---------------------------------------------------------------------------
# h in metres (recovered from the h/a series)
# ---------------------------------------------------------------------------

def test_h_metres_sin_pow(ref):
    result = ev(expansions.h_in_sin_pow(MAX_ORDER, MAX_ORDER), ref) * mp_a + ref["rho"] - mp_a
    assert_close("h [m]  sin_pow", ref["h"], result, TOL)