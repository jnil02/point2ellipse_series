"""Generate test data for C++ consistency checks.

Outputs CSV files to test_data/ with columns: n,k,l,num,den

All generators iterate over the full index range [0, MAX_INDEX] for all
indices. This means out-of-range indices are included, and Python's implicit
zeros for invalid indices are captured in the CSV. The C++ tests then verify
both the non-zero coefficients and the zeros.
"""

import os
import csv

from coefficients import d_phi, d_sin, d_cos, d_h, d_phi_evo, c_phi_evo

# Generate all indices up to and including this value.
MAX_INDEX = 5

# Where to place the test data.
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'test_data')


def write_csv(filename, rows):
    """Write rows to a CSV file in TEST_DATA_DIR."""
    os.makedirs(TEST_DATA_DIR, exist_ok=True)
    path = os.path.join(TEST_DATA_DIR, filename)
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['n', 'k', 'l', 'num', 'den'])
        writer.writerows(rows)
    print(f"Written {len(rows)} rows to {path}")


def generate_d_phi():
    """Generate test data for d_phi.

    Non-zero for k >= 1 and max(n+1, k) <= l <= n+k.
    Full range tested to verify zeros outside valid indices.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(MAX_INDEX + 1):
            for l in range(MAX_INDEX + 1):
                c = d_phi(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_phi.csv', rows)


def generate_d_sin():
    """Generate test data for d_sin.

    Non-zero for k >= 1 and max(n, k) <= l <= n+k.
    Full range tested to verify zeros outside valid indices.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(MAX_INDEX + 1):
            for l in range(MAX_INDEX + 1):
                c = d_sin(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_sin.csv', rows)


def generate_d_cos():
    """Generate test data for d_cos.

    Non-zero for k >= 1 and max(n, k) <= l <= n+k-1.
    Full range tested to verify zeros outside valid indices.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(MAX_INDEX + 1):
            for l in range(MAX_INDEX + 1):
                c = d_cos(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_cos.csv', rows)


def generate_d_h():
    """Generate test data for d_h.

    Non-zero for max(n, k+1) <= l <= n+k.
    Full range tested to verify zeros outside valid indices.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(MAX_INDEX + 1):
            for l in range(MAX_INDEX + 1):
                c = d_h(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_h.csv', rows)


def generate_d_phi_evo():
    """Generate test data for d_phi_evo.

    Non-zero for l <= n//2 + k.
    Full range tested to verify zeros outside valid indices.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(MAX_INDEX + 1):
            for l in range(MAX_INDEX + 1):
                c = d_phi_evo(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_phi_evo.csv', rows)


def generate_c_phi_evo():
    """Generate test data for c_phi_evo.

    Non-zero for k >= n+1, 1 <= l <= k, (k-n-1)%2==0, (l-k)%2==0.
    Full range tested to verify zeros outside valid indices.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(MAX_INDEX + 1):
            for l in range(MAX_INDEX + 1):
                c = c_phi_evo(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('c_phi_evo.csv', rows)


if __name__ == '__main__':
    generate_d_phi()
    generate_d_cos()
    generate_d_sin()
    generate_d_h()
    generate_d_phi_evo()
    generate_c_phi_evo()