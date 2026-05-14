from mpmath import mp


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