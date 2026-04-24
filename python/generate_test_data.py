"""Generate test data for C++ consistency checks.

Outputs CSV files to test_data/ with columns: n,k,l,num,den
"""

import os
import sys
import csv

from coefficients import d_phi, d_sin, d_cos, d_h, d_phi_evo

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

    Valid indices: k >= 1 and max(n+1, k) <= l <= n+k.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(1, MAX_INDEX + 1):
            for l in range(max(n + 1, k), n + k + 1):
                c = d_phi(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_phi.csv', rows)


def generate_d_sin():
    """Generate test data for d_sin.

    Valid indices: k >= 1 and max(n, k) <= l <= n+k.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(1, MAX_INDEX + 1):
            for l in range(max(n, k), n + k + 1):
                c = d_sin(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_sin.csv', rows)


def generate_d_cos():
    """Generate test data for d_cos.

    Valid indices: k >= 1 and max(n, k) <= l <= n+k-1.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(1, MAX_INDEX + 1):
            for l in range(max(n, k), n + k):
                c = d_cos(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_cos.csv', rows)


def generate_d_h():
    """Generate test data for d_h.

    Valid indices: max(n, k+1) <= l <= n+k.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(MAX_INDEX + 1):
            for l in range(max(n, k + 1), n + k + 1):
                c = d_h(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_h.csv', rows)


def generate_d_phi_evo():
    """Generate test data for d_phi_evo.

    Valid indices: l <= n//2 + k.
    """
    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(MAX_INDEX + 1):
            for l in range(n // 2 + k + 1):
                c = d_phi_evo(n, k, l)
                rows.append((n, k, l, c.p, c.q))
    write_csv('d_phi_evo.csv', rows)

if __name__ == '__main__':
    generate_d_phi()
    generate_d_cos()
    generate_d_sin()
    generate_d_h()
    generate_d_phi_evo()
