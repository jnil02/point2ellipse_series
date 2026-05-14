"""Pytest tests for point-to-ellipse Fourier series evaluation.

Uses the WGS84 Earth ellipse with a fixed test point (lat 43.1°, alt 10 000 m).
All assertions compare a truncated series (order MAX_ORDER) against the
multi-precision reference value, so the tolerance reflects truncation error,
not floating-point noise.
"""

from mpmath import mp
mp.dps = 50

import pytest

import ellipse
from ellipse import mp_ellipse_to_cartesian, mp_e2, mp_a
from symbols import varrho, psi, sin_psi, cos_psi, e2
import fourier_series

# Order of the tested series. The convergence appear to have geometric convergence.
MAX_ORDER = 7
# Truncated-series tolerance: e2 ≈ 0.007, so e2^8 ≈ 6e-17; 1e-10 is generous.
TOL = mp.mpf("1e-10")

# Test point.
_PHI = mp.mpf("43.1") / mp.mpf("180.") * mp.pi
_H   = mp.mpf("10000.")

# Reference values for the tests.
@pytest.fixture(scope="module")
def ref():
    x, y = mp_ellipse_to_cartesian(_PHI, _H)
    rho      = mp.sqrt(x*x + y*y)
    psi_val  = mp.atan2(y, x)
    varrho_v = mp_a / rho
    return dict(phi=_PHI, h=_H, x=x, y=y, psi=psi_val, rho=rho, varrho=varrho_v)


def assert_close(title, expected, actual, tol):
    """Print a full-precision comparison and assert that actual is within tol of expected."""
    abs_err = abs(actual - expected)
    rel_err = abs_err / abs(expected) if expected != 0 else mp.inf
    print(f"\n{title}")
    print(f"  expected: {mp.nstr(expected, mp.dps)}")
    print(f"  actual:   {mp.nstr(actual,   mp.dps)}")
    print(f"  abs err:  {mp.nstr(abs_err, 3)}  (tol {mp.nstr(tol, 3)})")
    print(f"  rel err:  {mp.nstr(rel_err, 3)}")
    assert abs_err < tol


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
    result = ev(fourier_series.phi_in_sin_pow(MAX_ORDER, MAX_ORDER) * sin_psi * cos_psi, ref)
    assert_close("phi-psi  sin_pow", expected, result, TOL)


def test_phi_minus_psi_sin_mul(ref):
    expected = ref["phi"] - ref["psi"]
    result = ev(fourier_series.phi_in_sin_mul(MAX_ORDER, MAX_ORDER, MAX_ORDER), ref)
    assert_close("phi-psi  sin_mul", expected, result, TOL)


def test_phi_minus_psi_sin_pow2(ref):
    expected = ref["phi"] - ref["psi"]
    result = ev(fourier_series.phi_in_sin_pow2(MAX_ORDER, MAX_ORDER), ref)
    # phi_in_sin_pow2 has worse convergence by design (see docstring); ~2e-7 at order 7.
    assert_close("phi-psi  sin_pow2", expected, result, mp.mpf("1e-5"))


# ---------------------------------------------------------------------------
# sin(phi) / sin(psi)
# ---------------------------------------------------------------------------

def test_sin_phi_over_sin_psi_sin_pow(ref):
    expected = mp.sin(ref["phi"]) / mp.sin(ref["psi"])
    result = ev(fourier_series.sin_phi_in_sin_pow(MAX_ORDER, MAX_ORDER), ref) + 1
    assert_close("sin(phi)/sin(psi)  sin_pow", expected, result, TOL)


def test_sin_phi_over_sin_psi_cos_mul(ref):
    expected = mp.sin(ref["phi"]) / mp.sin(ref["psi"])
    result = ev(fourier_series.sin_phi_in_cos_mul(MAX_ORDER, MAX_ORDER, MAX_ORDER), ref) + 1
    assert_close("sin(phi)/sin(psi)  cos_mul", expected, result, TOL)


def test_sin_phi_over_sin_psi_sin_pow2(ref):
    expected = mp.sin(ref["phi"]) / mp.sin(ref["psi"])
    result = ev(fourier_series.sin_phi_in_sin_pow2(MAX_ORDER, MAX_ORDER, MAX_ORDER), ref)
    assert_close("sin(phi)/sin(psi)  sin_pow2", expected, result, TOL)


# ---------------------------------------------------------------------------
# cos(phi) / cos(psi)
# ---------------------------------------------------------------------------

def test_cos_phi_over_cos_psi_sin_pow(ref):
    expected = mp.cos(ref["phi"]) / mp.cos(ref["psi"])
    result = ev(fourier_series.cos_phi_in_sin_pow(MAX_ORDER, MAX_ORDER), ref) + 1
    assert_close("cos(phi)/cos(psi)  sin_pow", expected, result, TOL)


def test_cos_phi_over_cos_psi_cos_mul(ref):
    expected = mp.cos(ref["phi"]) / mp.cos(ref["psi"])
    result = ev(fourier_series.cos_phi_in_cos_mul(MAX_ORDER, MAX_ORDER, MAX_ORDER), ref) + 1
    assert_close("cos(phi)/cos(psi)  cos_mul", expected, result, TOL)


def test_cos_phi_over_cos_psi_sin_pow2(ref):
    expected = mp.cos(ref["phi"]) / mp.cos(ref["psi"])
    result = ev(fourier_series.cos_phi_in_sin_pow2(MAX_ORDER, MAX_ORDER, MAX_ORDER), ref)
    assert_close("cos(phi)/cos(psi)  sin_pow2", expected, result, TOL)


# ---------------------------------------------------------------------------
# (h + a - rho) / a
# ---------------------------------------------------------------------------

def test_h_a_sin_pow(ref):
    expected = (ref["h"] + mp_a - ref["rho"]) / mp_a
    result = ev(fourier_series.h_in_sin_pow(MAX_ORDER, MAX_ORDER), ref)
    assert_close("(h+a-rho)/a  sin_pow", expected, result, TOL)


def test_h_a_cos_mul(ref):
    expected = (ref["h"] + mp_a - ref["rho"]) / mp_a
    result = ev(fourier_series.h_in_cos_mul(MAX_ORDER, MAX_ORDER, MAX_ORDER), ref)
    assert_close("(h+a-rho)/a  cos_mul", expected, result, TOL)


# ---------------------------------------------------------------------------
# h in metres (recovered from the h/a series)
# ---------------------------------------------------------------------------

def test_h_metres_sin_pow(ref):
    result = ev(fourier_series.h_in_sin_pow(MAX_ORDER, MAX_ORDER), ref) * mp_a + ref["rho"] - mp_a
    assert_close("h [m]  sin_pow", ref["h"], result, TOL)