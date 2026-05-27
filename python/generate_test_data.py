"""Generate test data for C++ consistency checks.

Outputs CSV files to test_data/ with columns: n,k,l,num,den

All generators iterate over the full index range [0, MAX_INDEX] for all
indices. This means out-of-range indices are included, and Python's implicit
zeros for invalid indices are captured in the CSV. The C++ tests then verify
both the non-zero coefficients and the zeros.
"""

import os
import csv

from coefficients import (d_phi, d_sin, d_cos, d_h, d_phi_evo, c_phi_evo,
                          c_phi_pow_evo, c_sin_phi_evo, d_sin_phi_evo,
                          c_cos_phi_evo, d_cos_phi_evo,
                          c_sin_phi_inv_evo, a_mr, B_rt, C_mt, R, B_p, cp_evo_nkl,
                          ch_evo)

# Generate all indices up to and including this value.
MAX_INDEX = 5
MAX_INDEX_POW = 3

# Where to place the test data.
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'test_data')


def write_csv(filename, rows, header):
    """Write rows to a CSV file in TEST_DATA_DIR."""
    os.makedirs(TEST_DATA_DIR, exist_ok=True)
    path = os.path.join(TEST_DATA_DIR, filename)
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"Written {len(rows)} rows to {path}")


def generate_d_phi():
    """Generate test data for d_phi.

    Non-zero for k >= 1 and max(n+1, k) <= l <= n+k.
    Full range tested to verify zeros outside valid indices.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(1, MAX_INDEX + 1):
            for l in range(max(n+1,k), n+k + 1):
                c = d_phi(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_phi.csv', rows, ['n', 'k', 'l', 'num', 'den'])


def generate_d_sin():
    """Generate test data for d_sin.

    Non-zero for k >= 1 and max(n, k) <= l <= n+k.
    Full range tested to verify zeros outside valid indices.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(1, MAX_INDEX + 1):
            for l in range(max(n,k),n+k + 1):
                c = d_sin(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_sin.csv', rows, ['n', 'k', 'l', 'num', 'den'])


def generate_d_cos():
    """Generate test data for d_cos.

    Non-zero for k >= 1 and max(n, k) <= l <= n+k-1.
    Full range tested to verify zeros outside valid indices.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(1, MAX_INDEX + 1):
            for l in range(max(n,k), n+k):
                c = d_cos(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_cos.csv', rows, ['n', 'k', 'l', 'num', 'den'])


def generate_d_h():
    """Generate test data for d_h.

    Non-zero for max(n, k+1) <= l <= n+k.
    Full range tested to verify zeros outside valid indices.
    """
    rows = []
    for n in range(1, MAX_INDEX + 1):
        for k in range(MAX_INDEX + 1):
            for l in range(max(n,k+1), n+k + 1):
                c = d_h(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_h.csv', rows, ['n', 'k', 'l', 'num', 'den'])


def generate_d_phi_evo():
    """Generate test data for d_phi_evo.

    Non-zero for l <= n//2 + k.
    Full range tested to verify zeros outside valid indices.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(MAX_INDEX + 1):
            for l in range(n//2+k + 1):
                c = d_phi_evo(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_phi_evo.csv', rows, ['n', 'k', 'l', 'num', 'den'])


def generate_c_phi_evo():
    """Generate test data for c_phi_evo.

    Non-zero for k >= n+1, 1 <= l <= k, (k-n-1)%2==0, (l-k)%2==0.
    Full range tested to verify zeros outside valid indices.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(n+1, MAX_INDEX + 1):
            for l in range(1, k + 1):
                c = c_phi_evo(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('c_phi_evo.csv', rows, ['n', 'k', 'l', 'num', 'den'])


def generate_d_phi_pow_evo():
    """Generate test data for d_phi_pow_evo.

    Non-zero for (i+n-k)%2==0 and (l-k)%2==0.
    Uses a smaller MAX_INDEX due to computational cost.
    """
    rows = []
    for i in range(MAX_INDEX_POW + 1):
        for n in range(MAX_INDEX + 1):
            for k in range(n+i,MAX_INDEX + 1):
                for l in range(i, k + 1):
                    c = c_phi_pow_evo(n, k, l, i)
                    rows.append((n, k, l, i, c.p, c.q))
    write_csv('c_phi_pow_evo.csv', rows, ['n', 'k', 'l', 'i', 'num', 'den'])


def generate_c_sin_phi_evo():
    """Generate test data for c_sin_phi_evo (sparse, raw indices)."""
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(n, MAX_INDEX + 1):
            for l in range(k + 1):
                c = c_sin_phi_evo(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('c_sin_phi_evo.csv', rows, ['n', 'k', 'l', 'num', 'den'])


def generate_d_sin_phi_evo():
    """Generate test data for d_sin_phi_evo (dense, transformed indices).

    Valid range: n >= 0, k >= 0, 1 <= l <= n//2 + k.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(MAX_INDEX + 1):
            for l in range(1, n // 2 + k + 1):
                c = d_sin_phi_evo(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_sin_phi_evo.csv', rows, ['n', 'k', 'l', 'num', 'den'])


def generate_c_cos_phi_evo():
    """Generate test data for c_cos_phi_evo (sparse, raw indices)."""
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(n, MAX_INDEX + 1):
            for l in range(1, k + 1):
                c = c_cos_phi_evo(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('c_cos_phi_evo.csv', rows, ['n', 'k', 'l', 'num', 'den'])


def generate_d_cos_phi_evo():
    """Generate test data for d_cos_phi_evo (dense, transformed indices).

    Valid range: n >= 0, k >= 0, n%2 <= l <= (n+1)//2 + k.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(MAX_INDEX + 1):
            for l in range(n % 2, (n + 1) // 2 + k + 1):
                c = d_cos_phi_evo(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_cos_phi_evo.csv', rows, ['n', 'k', 'l', 'num', 'den'])


def generate_d_sin_phi_inv_evo():
    """Generate test data for d_sin_phi_inv_evo.

    Non-zero for (n-k)%2==0 and (n-l)%2==0.
    Full range tested to verify zeros outside valid indices.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(n, MAX_INDEX + 1):
            for l in range(k + 1):
                c = c_sin_phi_inv_evo(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('c_sin_phi_inv_evo.csv', rows, ['n', 'k', 'l', 'num', 'den'])


def generate_a_mr():
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(n + 1):
            c = a_mr(n,k)
            rows.append((n, k, c.p, c.q))
    write_csv('a_mr.csv', rows, ['n', 'k'])


def generate_B_rt():
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(n + 1):
            c = B_rt(n,k)
            rows.append((n, k, c.p, c.q))
    write_csv('B_rt.csv', rows, ['n', 'k'])


def generate_C_mt():
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(n + 1):
            c = C_mt(n,k)
            rows.append((n, k, c.p, c.q))
    write_csv('C_mt.csv', rows, ['n', 'k'])


def generate_R():
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(n, MAX_INDEX + 1):
            for l in range(k + 1):
                for i in range(l // 2 + 1):
                    c = R(n, k, l, i)
                    rows.append((n, k, l, i, c.p, c.q))
    write_csv('R.csv', rows, ['n', 'k', 'l', 'i', 'num', 'den'])


def generate_B_p():
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(n, MAX_INDEX + 1):
            for l in range(1, k+2):
                c = B_p(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('B_p.csv', rows, ['n', 'k', 'l', 'num', 'den'])


def generate_cp_evo_nkl():
    rows = []
    for n in range(1, MAX_INDEX + 1):
        for k in range(n, MAX_INDEX + 1):
            for l in range(k + 2):
                c = cp_evo_nkl(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('cp_evo_nkl.csv', rows, ['n', 'k', 'l', 'num', 'den'])


def generate_c_h_evo():
    rows = []
    for n in range(1, MAX_INDEX + 1):
        for k in range(n, MAX_INDEX + 1):
            for l in range(k + 2):
                c = ch_evo(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('c_h_evo.csv', rows, ['n', 'k', 'l', 'num', 'den'])


if __name__ == '__main__':
    generate_d_phi()
    generate_d_cos()
    generate_d_sin()
    generate_d_h()
    generate_d_phi_evo()
    generate_c_phi_evo()
    generate_d_phi_pow_evo()
    generate_c_sin_phi_evo()
    generate_d_sin_phi_evo()
    generate_c_cos_phi_evo()
    generate_d_cos_phi_evo()
    generate_d_sin_phi_inv_evo()
    generate_a_mr()
    generate_B_rt()
    generate_C_mt()
    generate_R()
    generate_B_p()
    generate_cp_evo_nkl()
    generate_c_h_evo()
