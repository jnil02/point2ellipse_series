"""Pytest tests for the inside-evolute point-to-ellipse Fourier series.

Uses the WGS84 Earth ellipse with a fixed test point inside the evolute
(polar angle 138°, radius 5000 m). All assertions compare a truncated
series (order MAX_ORDER_EVO) against the multi-precision reference value.
"""

from mpmath import mp
mp.dps = 50

import pytest

from pytest_util import assert_close
from ellipse import mp_polar_to_cartesian, mp_cartesian_to_ellipse, mp_e2, mp_a
from symbols import sin_psi, rho_ae2, b_a
import fourier_series

# Order of the tested series.
MAX_ORDER_EVO = 11

# convergence parameter rho_ae2 = rho/(a*e2) ≈ 0.117, so rho_ae2^12 ≈ 5e-12.
TOL_EVO = mp.mpf("1e-9")

# Test point: polar coords inside the ellipse evolute.
_PSI_EVO = mp.mpf("138.") / mp.mpf("180.") * mp.pi
_RHO_EVO = mp.mpf("5000.")


@pytest.fixture(scope="module")
def ref_evo():
    cx, cy = mp_polar_to_cartesian(_PSI_EVO, _RHO_EVO)
    phi, h = mp_cartesian_to_ellipse(cx, cy)
    sgn = mp.mpf("1") if _PSI_EVO > 0 else mp.mpf("-1")
    return dict(
        psi         = _PSI_EVO,
        rho         = _RHO_EVO,
        phi         = phi,
        h           = h,
        sgn         = sgn,
        abs_sin_psi = abs(mp.sin(_PSI_EVO)),
        abs_cos_psi = abs(mp.cos(_PSI_EVO)),
        rho_ae2_val = _RHO_EVO / (mp_a * mp_e2),
        b_a_val     = mp.sqrt(mp.mpf(1) - mp_e2),
    )


def ev_evo(expr, ref):
    """Substitute evo-series reference values into a symbolic expression.

    Uses abs(sin(psi)), rho/(a*e2), and b/a as the three substitution variables.
    Returns an mpmath mpf preserving full precision.
    """
    result = (expr
              .subs(sin_psi, ref["abs_sin_psi"])
              .subs(rho_ae2, ref["rho_ae2_val"])
              .subs(b_a,     ref["b_a_val"]))
    return mp.mpf(result._mpf_)


# ---------------------------------------------------------------------------
# (phi - sgn*pi/2) / (sgn * |cos(psi)|)
# ---------------------------------------------------------------------------

def test_phi_evo_dense(ref_evo):
    expected = ((ref_evo["phi"] - ref_evo["sgn"] * mp.pi / 2)
                / (ref_evo["sgn"] * ref_evo["abs_cos_psi"]))
    result = ev_evo(fourier_series.phi_evo_sin_pow_dense(MAX_ORDER_EVO, MAX_ORDER_EVO), ref_evo)
    assert_close("(phi - sgn*pi/2)/(sgn*|cos(psi)|)  evo_dense", expected, result, TOL_EVO)

def test_phi_evo_dense_m(ref_evo):
    expected = ((ref_evo["phi"] - ref_evo["sgn"] * mp.pi / 2)
                / (ref_evo["sgn"] * ref_evo["abs_cos_psi"]))
    result = ev_evo(fourier_series.phi_evo_sin_pow_dense_m(MAX_ORDER_EVO), ref_evo)
    assert_close("(phi - sgn*pi/2)/(sgn*|cos(psi)|)  evo_dense", expected, result, TOL_EVO)


# ---------------------------------------------------------------------------
# (sin(phi) - sgn) / sgn
# ---------------------------------------------------------------------------

def test_sin_phi_evo_dense(ref_evo):
    expected = (mp.sin(ref_evo["phi"]) - ref_evo["sgn"]) / ref_evo["sgn"]
    result = ev_evo(fourier_series.sin_phi_evo_dense(MAX_ORDER_EVO, MAX_ORDER_EVO), ref_evo)
    assert_close("(sin(phi) - sgn)/sgn  evo_dense", expected, result, TOL_EVO)


# ---------------------------------------------------------------------------
# cos(phi) / |cos(psi)|
# ---------------------------------------------------------------------------

def test_cos_phi_evo_dense(ref_evo):
    expected = mp.cos(ref_evo["phi"]) / ref_evo["abs_cos_psi"]
    result = ev_evo(fourier_series.cos_phi_evo_dense(MAX_ORDER_EVO, MAX_ORDER_EVO), ref_evo)
    assert_close("cos(phi)/|cos(psi)|  evo_dense", expected, result, TOL_EVO)


# ---------------------------------------------------------------------------
# h in metres
# ---------------------------------------------------------------------------

def test_h_evo(ref_evo):
    result = ev_evo(fourier_series.h_a_evo(MAX_ORDER_EVO, MAX_ORDER_EVO), ref_evo) * mp_a
    assert_close("h [m]  evo", ref_evo["h"], result, TOL_EVO * mp_a)


def test_h_evo_dense(ref_evo):
    result = ev_evo(fourier_series.h_a_evo_dense(MAX_ORDER_EVO, MAX_ORDER_EVO), ref_evo) * mp_a
    assert_close("h [m]  evo_dense", ref_evo["h"], result, TOL_EVO * mp_a)