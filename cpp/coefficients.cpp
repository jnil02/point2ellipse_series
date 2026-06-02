
#include <cassert>
#include <map>

#include "series_substitution.hpp"
#include "cache.hpp"
#include "util.hpp"

#include "coefficients.hpp"

namespace point_to_ellipse_series {

using uint = unsigned int;

mpq_class d_phi(int n, int k, int l) {
	assert(n >= 0 && k >= 1 && l >= std::max(n + 1, k) && l <= n + k);
	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class d(0);
	for (int r = 0; r <= k - 1; ++r)
		for (int m = 0; m <= k - 1 - r; ++m)
			for (int q = 0; q <= r / 2; ++q)
				for (int p = 0; p <= m / 2; ++p) {
					mpz_class fact_q, fact_r2q;
					mpz_fac_ui(fact_q.get_mpz_t(), (unsigned long) q);
					mpz_fac_ui(fact_r2q.get_mpz_t(), (unsigned long) (r - 2 * q));
					mpz_class b_km1_mr, b_m1_2p1;
					mpz_bin_uiui(b_km1_mr.get_mpz_t(), (unsigned long) (k - 1), (unsigned long) (m + r));
					mpz_bin_uiui(b_m1_2p1.get_mpz_t(), (unsigned long) (m + 1), (unsigned long) (2 * p + 1));
					mpq_class outer = rf_half(k, r - q)
									  * mpq_class(powm1(n - l - k) * pow2(r - 2 * q));
					outer /= mpq_class(fact_q * fact_r2q * mpz_class(m + 1 + r));
					outer *= mpq_class(b_km1_mr * b_m1_2p1);
					outer.canonicalize();

					for (int t = l - n - k + m + r - p - q;
						 t <= std::min({l - k, r - q, l - n - k + m + r - q}); ++t) {
						d += outer
							 * int_bin(r - q, t)
							 * binomial_rational(mpq_class(k, 2) + (long) (r - q + l - k - t - 1),
												 (long) (l - k - t))
							 * int_bin(p, n - (m + r - q + l - k - t - p));
					}
				}

	mpq_class ret = d;
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

mpq_class d_phi2(int n, int k, int l) {
	assert(n >= 0 && k >= 1 && l >= k && l <= n + k);
	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class d(0);
	for (int r = 0; r <= k - 1; ++r)
		for (int m = 0; m <= k - 1 - r; ++m)
			for (int q = 0; q <= r / 2; ++q)
				for (int p = 0; p <= m / 2; ++p) {
					mpz_class fact_q, fact_r2q;
					mpz_fac_ui(fact_q.get_mpz_t(), (unsigned long) q);
					mpz_fac_ui(fact_r2q.get_mpz_t(), (unsigned long) (r - 2 * q));
					mpz_class b_km1_mr, b_m1_2p1;
					mpz_bin_uiui(b_km1_mr.get_mpz_t(), (unsigned long) (k - 1), (unsigned long) (m + r));
					mpz_bin_uiui(b_m1_2p1.get_mpz_t(), (unsigned long) (m + 1), (unsigned long) (2 * p + 1));
					mpq_class outer = rf_half(k, r - q)
									  * mpq_class(powm1(n - l - k) * pow2(r - 2 * q));
					outer /= mpq_class(fact_q * fact_r2q * mpz_class(m + 1 + r));
					outer *= mpq_class(b_km1_mr * b_m1_2p1);
					outer.canonicalize();

					for (int t = l - n - k + m + r - p - q;
						 t <= std::min(l - k, r - q); ++t) {
						d += outer
							 * int_bin(r - q, t)
							 * binomial_rational(mpq_class(k, 2) + (long) (r - q + l - k - t - 1),
												 (long) (l - k - t))
							 * binomial_rational(mpq_class(1, 2) + (long) p,
												 (long) (n - (m + r - q + l - k - t - p)));
					}
				}

	mpq_class ret = d;
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

mpq_class c_phi(int n, int k, int l) {
	assert(n >= 1 && k >= 1 && l >= std::max(n, k));
	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class d(0);
	for (int r = 0; r <= k - 1; ++r)
		for (int m = 0; m <= k - 1 - r; ++m)
			for (int p = 0; p <= m / 2; ++p)
				for (int q = 0; q <= r / 2; ++q) {
					// Precompute factors independent of t.
					mpz_class fact_q, fact_r2q;
					mpz_fac_ui(fact_q.get_mpz_t(), (unsigned long) q);
					mpz_fac_ui(fact_r2q.get_mpz_t(), (unsigned long) (r - 2 * q));
					mpz_class b_km1_mr, b_m1_2p1;
					mpz_bin_uiui(b_km1_mr.get_mpz_t(), (unsigned long) (k - 1), (unsigned long) (m + r));
					mpz_bin_uiui(b_m1_2p1.get_mpz_t(), (unsigned long) (m + 1), (unsigned long) (2 * p + 1));
					mpq_class rf = rf_half(k, r - q) * mpq_class(powm1(l - k));
					rf /= mpq_class(fact_q * fact_r2q * mpz_class(m + 1 + r));
					rf *= mpq_class(b_km1_mr * b_m1_2p1);
					rf.canonicalize();

					for (int t = 0; t <= std::min(l - k, r - q); ++t) {
						const int w = m + r - q + l - k - t - p;
						mpq_class h(0);
						for (int i = std::max(0, p - n + 1);
							 i <= std::min(p, p - n + 1 + w); ++i) {
							const int j = w + p - i + 1 - n;
							h += int_bin(2 * p + 1, i) * int_bin(1 + 2 * w, j) * mpq_class(powm1(j));
						}
						for (int i = p - w + n; i <= p; ++i) {
							const int j = w - p + i - n;
							h += int_bin(2 * p + 1, i) * int_bin(1 + 2 * w, j) * mpq_class(powm1(j));
						}
						for (int i = p - w - n; i <= p - n; ++i) {
							const int j = w - p + i + n;
							h -= int_bin(2 * p + 1, i) * int_bin(1 + 2 * w, j) * mpq_class(powm1(j));
						}
						if (h == mpq_class(0)) continue;

						mpz_class b_rq_t;
						mpz_bin_uiui(b_rq_t.get_mpz_t(), (unsigned long) (r - q), (unsigned long) t);
						mpz_class shift_val;
						mpz_ui_pow_ui(shift_val.get_mpz_t(), 2,
									  (unsigned long) (2 * (m + l - k - t) + r + 1));
						mpq_class outer = rf * mpq_class(b_rq_t) / mpq_class(shift_val);
						outer *= binomial_rational(mpq_class(k, 2) + (long) (r - q + l - k - t - 1),
												   (long) (l - k - t));
						outer.canonicalize();
						d += h * outer;
					}
				}

	mpq_class ret = d;
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

static E2Poly d_phi_pow_e2poly_se4(int n, int k, int i) {
	static UintsCache<E2Poly> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) i)) return *ret;

	static std::map<std::pair<int,int>, std::shared_ptr<TSeriesBase<LExpr>>> series_cache;
	auto key = std::make_pair(n, i);
	if (!series_cache.count(key))
		series_cache[key] = double_series_power_coeff_lexpr(n, i);

	LExpr lp = series_cache[key]->getItem(k);
	E2Poly result = lexpr_eval_e2poly(lp, [](int j, int m) {
		return a_nk_ser_lexpr(j, m, 1, d_phi);
	});

	return cache.insert(result, (uint) n, (uint) k, (uint) i);
}

mpq_class d_phi_pow_se4(int n, int k, int l, int i) {
	E2Poly ep = d_phi_pow_e2poly_se4(n, k, i);
	return (l < (int) ep.size()) ? ep[l] : mpq_class(0);
}

mpq_class d_phi_pow(int n, int k, int l, int i) {
	assert(n >= 0 && k >= i && l >= std::max(n + i, k) && l <= n + k && i >= 1);
	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l, (uint) i))
		return *ret;

	mpq_class ret = d_phi_pow_se4(n, k, l, i);  // ← LExpr/GMP pipeline (was: _se2)

	return cache.insert(ret, (uint) n, (uint) k, (uint) l, (uint) i);
}

mpq_class d_sin(int n, int k, int l) {
	assert(n >= 0 && k >= 1 && l >= std::max(n,k) && l <= n + k);
	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class d(0);
	const int i_max = std::min(k, 2 * n + 1);
	for (int i = 1; i <= i_max; ++i) {
		const int floor_i_2 = i / 2;
		const int ceil_i_2 = (i + 1) / 2;
		mpz_class fact;
		mpz_fac_ui(fact.get_mpz_t(), (unsigned long) i);
		for (int j = std::max(0, ceil_i_2 - l + n);
			 j <= std::min(ceil_i_2, std::min(n - floor_i_2, n+k-l-floor_i_2)); ++j) {
			mpz_class binom;
			mpz_bin_uiui(binom.get_mpz_t(), (unsigned long) ceil_i_2, (unsigned long) j);
			mpq_class coeff(powm1(floor_i_2 + j) * binom, fact);
			coeff.canonicalize();
			d += coeff * d_phi_pow(n - floor_i_2 - j, k, l, i);
		}
	}

	mpq_class ret = d;
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

static E2Poly d_sin_pow_e2poly_se4(int n, int k, int i) {
	static UintsCache<E2Poly> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) i)) return *ret;

	static std::map<std::pair<int,int>, std::shared_ptr<TSeriesBase<LExpr>>> series_cache;
	auto key = std::make_pair(n, i);
	if (!series_cache.count(key))
		series_cache[key] = double_series_power_coeff_lexpr(n, i);

	LExpr lp = series_cache[key]->getItem(k);
	E2Poly result = lexpr_eval_e2poly(lp, [](int j, int m) {
		return a_nk_ser_lexpr(j, m, 0, d_sin);
	});

	return cache.insert(result, (uint) n, (uint) k, (uint) i);
}

mpq_class d_sin_pow_se4(int n, int k, int l, int i) {
	E2Poly ep = d_sin_pow_e2poly_se4(n, k, i);
	return (l < (int) ep.size()) ? ep[l] : mpq_class(0);
}

mpq_class d_sin_pow(int n, int k, int l, int i) {
	assert(n >= 0 && k >= 0 && l >= 0 && l <= n + k && i >= 0);
	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l, (uint) i))
		return *ret;

	mpq_class ret = d_sin_pow_se4(n, k, l, i);  // ← LExpr/GMP pipeline (was: _se2)

	return cache.insert(ret, (uint) n, (uint) k, (uint) l, (uint) i);
}

mpq_class c_sin(int n, int k, int l) {
	assert(n >= 0 && k >= 1 && l >= std::max(n, k));
	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class ret = sin_pow_to_cos_mul(n, k, l, 0, 0, d_sin);
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

mpq_class d_cos(int n, int k, int l) {
	assert(n >= 0 && k >= 1 && l >= std::max(n, k) && l <= n + k - 1);
	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class d(0);
	const int i_max = std::min(k, 2 * n + 1);
	for (int i = 1; i <= i_max; ++i) {
		const int floor_i_2 = i / 2;
		const int ceil_i_2 = (i + 1) / 2;
		mpz_class fact;
		mpz_fac_ui(fact.get_mpz_t(), (unsigned long) i);
		for (int j = std::max(0, floor_i_2 - l + n);
			 j <= std::min(floor_i_2, std::min(n - ceil_i_2, n+k-l-ceil_i_2)); ++j) {
			mpz_class binom;
			mpz_bin_uiui(binom.get_mpz_t(), (unsigned long) floor_i_2, (unsigned long) j);
			mpq_class coeff(powm1(ceil_i_2 + j) * binom, fact);
			coeff.canonicalize();
			d += coeff * d_phi_pow(n - ceil_i_2 - j, k, l, i);
		}
	}

	mpq_class ret = d;
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

mpq_class c_cos(int n, int k, int l) {
	assert(n >= 0 && k >= 1 && l >= std::max(n, k));
	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class ret = sin_pow_to_cos_mul(n, k, l, 0, -1, d_cos);
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

mpq_class d_N_nkl(int n, int k, int l) {
	assert(n >= 1 && k >= 0 && l >= std::max(n, k + 1) && l <= n + k);
	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class d(0);
	for (int i = 1; i <= n; ++i) {
		mpq_class binom_half_i = binomial_rational(mpq_class(1, 2), (long) i)
								 * mpq_class(powm1(i));
		for (int j = 0; j <= std::min(2 * i, k); ++j) {
			mpz_class b2i_j;
			mpz_bin_uiui(b2i_j.get_mpz_t(), (unsigned long) (2 * i), (unsigned long) j);
			d += binom_half_i * mpq_class(b2i_j) * d_sin_pow(n - i, k, l - i, j);
		}
	}

	mpq_class ret = d;
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

mpq_class bp_nkl(int n, int k, int l) {
	assert(n >= 1 && k >= 1 && l >= std::max(n, k) && l <= n + k - 1);
	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class d(0);
	for (int i = 1; i <= std::min(n, k / 2); ++i) {
		mpz_class fact;
		mpz_fac_ui(fact.get_mpz_t(), (unsigned long)(2 * i));
		for (int j = std::max(0, n-l+i); j <= std::min(i, std::min(n - i, n+k-l-i)); ++j) {
			mpz_class binom;
			mpz_bin_uiui(binom.get_mpz_t(), (unsigned long) i, (unsigned long) j);
			mpq_class coeff(powm1(i + j) * binom, fact);
			coeff.canonicalize();
			d += coeff * d_phi_pow(n - i - j, k, l, 2 * i);
		}
	}

	mpq_class ret = d;
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

mpq_class d_h(int n, int k, int l) {
	assert(n >= 1 && k >= 0 && l >= std::max(n, k + 1) && l <= n+k);
	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class ret;
	if (k == 0) {
		ret = -d_N_nkl(n, k, l);
	} else {
		ret = bp_nkl(n, k + 1, l) - d_N_nkl(n, k, l);
	}
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

mpq_class c_h(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= std::max(n, k + 1));
	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class ret = sin_pow_to_cos_mul(n, k, l, 1, 0, d_h);
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

}  // namespace point_to_ellipse_series
