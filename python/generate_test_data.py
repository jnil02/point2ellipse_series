"""Generate test data for C++ consistency checks.

Outputs test_data/d_phi_evo.csv with columns: n,k,l,num,den
Only valid indices (l <= n//2 + k) are included.
"""

import os
import sys
import csv

from coefficients import d_phi_evo

# Allow running from repo root or from python/.
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

# Generate all indices up to and including this value.
MAX_INDEX = 5

# Where to place the test data.
OUTPUT_PATH = os.path.join(os.path.dirname(__file__), '..', 'test_data', 'd_phi_evo.csv')

def main():
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)

    rows = []
    for n in range(MAX_INDEX + 1):
        for k in range(MAX_INDEX + 1):
            for l in range(n // 2 + k + 1):  # Only valid indices.
                c = d_phi_evo(n, k, l)  # c is a sympy Rational or Integer.
                rows.append((n, k, l, c.p, c.q)) # Indices, numerator and denominator.

    with open(OUTPUT_PATH, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['n', 'k', 'l', 'num', 'den'])
        writer.writerows(rows)

    print(f"Written {len(rows)} rows to {OUTPUT_PATH}")

if __name__ == '__main__':
    main()
