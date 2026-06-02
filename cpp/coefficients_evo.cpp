
#include <cassert>
#include <map>

#include "cache.hpp"
#include "util.hpp"
#include "series_substitution.hpp"

#include "coefficients_evo.hpp"
#include "stirling.hpp"

namespace point_to_ellipse_series {

using uint = unsigned int;

mpq_class d_phi_evo(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= 0 && l <= n / 2 + k);

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	const int s = n % 2;
	const int m = n + 1 + 2 * k;
	const int Q = s + 2 * l;  // Maximum value of q.

	// Every factor of the double sum is a non-negative integer, so the whole
	// accumulator c is an integer: rational arithmetic is deferred to the final
	// outer factor. Precompute the q-only factor T[q] = C(n-1+2k-q, Q-q), which
	// is independent of j (turns O(l^2) binomial calls into O(l)). The second
	// argument Q-q is always >= 0; the first is >= 0 except when Q-q == 0, where
	// C(top, 0) = 1, so a plain integer binomial covers every case.
	std::vector<mpz_class> T(Q + 1);
	for (int q = 0; q <= Q; ++q) {
		if (Q - q == 0)
			T[q] = 1;
		else
			mpz_bin_uiui(T[q].get_mpz_t(), (unsigned long) (n - 1 + 2 * k - q),
						 (unsigned long) (Q - q));
	}

	mpz_class c(0);
	// The inner b-sum reduces in closed form to b(j) = C(n/2+k-l, k-j) via the
	// alternating binomial identity sum_i (-1)^i C(j,i) C(x-i,k) = C(x-j,k-j).
	// It is nonzero only for max(0, l-n/2) <= j <= k, which bounds the loop.
	for (int j = std::max(0, l - n / 2); j <= std::min(l, k); ++j) {
		mpz_class b;
		mpz_bin_uiui(b.get_mpz_t(), (unsigned long) (n / 2 + k - l),
					 (unsigned long) (k - j));

		// Inner sum: sum_{q=2j}^{Q} 2^(q-2j) * C(q-j, j) * T[q]. Both 2^(q-2j)
		// and C(q-j, j) are advanced incrementally to avoid recomputation.
		mpz_class sum_(0);
		mpz_class shift(1);   // 2^(q-2j), starts at q = 2j.
		mpz_class bqj(1);     // C(q-j, j),  C(j, j) = 1 at q = 2j.
		for (int q = 2 * j, p = j; q <= Q; ++q, ++p) {
			sum_ += shift * bqj * T[q];
			shift *= 2;
			bqj *= (p + 1);        // C(p+1, j) = C(p, j) * (p+1) / (p+1-j).
			bqj /= (p + 1 - j);    // Exact integer division.
		}

		c += sum_ * b;
	}

	// Outer factor: (-1)^(n/2 + l + n + 1) * binomial(m/2, n/2 + k - l) / m
	mpq_class result =
			binomial_rational(mpq_class(m, 2), (long) (n / 2 + k - l))
			* mpq_class(mpz_class(powm1(n / 2 + l + n + 1)))
			/ mpq_class(mpz_class(m))
			* mpq_class(c);
	result.canonicalize();

	mpq_class ret = result;
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

// Old implementation
//mpq_class c_phi_evo(int n, int k, int l) {
//	assert(n >= 0 && k >= n + 1 && l >= 1 && l <= k);
//
//	// Parity constraints — coefficient is zero unless both hold.
//	if ((k - n - 1) % 2 != 0 || (l - k) % 2 != 0)
//		return {0};
//
//	static UintsCache<mpq_class> cache;
//	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
//		return *ret;
//
//	// Outer binomial factor: binomial(k/2, (k-l)/2) / k.
//	// (k-l)/2 >= 0 since (l-k)%2==0 and valid range has l <= k.
//	const long kl_2 = (long) ((k - l) / 2);
//	mpq_class outer = binomial_rational(mpq_class(k, 2), kl_2)
//					  / mpq_class(mpz_class(k));
//
//	mpq_class c(0);
//	for (int j = 0; j <= (l - 1) / 2; ++j) {
//
//		// Inner sum for b.
//		mpq_class b(0);
//		const int i_max = std::min(j, (n + 1 - l + 2 * j) / 2);
//		for (int i = 0; i <= i_max; ++i) {
//			const int ka = (n - l + 1 + 2 * j - 2 * i) / 2;
//			mpz_class bj_i;
//			mpz_bin_uiui(bj_i.get_mpz_t(), (unsigned long) j, (unsigned long) i);
//			b += mpq_class(mpz_class(powm1(ka)) * bj_i)
//				 * int_bin(j - i + (k - l) / 2, ka);
//		}
//
//		// Inner sum for s.
//		mpq_class s(0);
//		for (int q = 2 * j; q <= l - 1; ++q) {
//			mpz_class shift_val;
//			mpz_ui_pow_ui(shift_val.get_mpz_t(), 2, (unsigned long) (q - 2 * j));
//			mpz_class bqj_j;
//			mpz_bin_uiui(bqj_j.get_mpz_t(), (unsigned long) (q - j), (unsigned long) j);
//			s += mpq_class(shift_val * bqj_j) * int_bin(k - 2 - q, l - 1 - q);
//		}
//
//		c += outer * mpq_class(mpz_class(powm1(l - j))) * s * b;
//	}
//
//	mpq_class ret = c;
//	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
//}

mpq_class c_phi_evo(int n, int k, int l) {
	assert(n >= 0 && k >= n + 1 && l >= 1 && l <= k);

	// Parity constraints — coefficient is zero unless both hold.
	if ((k - n - 1) % 2 != 0 || (l - k) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	// Let a = (k-l)/2 and h = (n-l+1)/2 (both integers by parity).
	// Outer rational factor: binomial(k/2, a) / k.
	const int a = (k - l) / 2;
	const int h = (n - l + 1) / 2;
	mpq_class outer = binomial_rational(mpq_class(k, 2), (long) a)
					  / mpq_class(mpz_class(k));

	// Precompute q-only factor U[q] = C(k-2-q, l-1-q), independent of j.
	// l-1-q is always >= 0 for q in [0,l-1]. When q == l-1 the first arg may
	// be -1 (only when k == l); C(-1,0) = 1, handled explicitly.
	std::vector<mpz_class> U(l);
	for (int q = 0; q < l; ++q) {
		if (l - 1 - q == 0) U[q] = 1;
		else mpz_bin_uiui(U[q].get_mpz_t(), (unsigned long) (k - 2 - q),
						  (unsigned long) (l - 1 - q));
	}

	// b(j) collapses to (-1)^{h+j} * C(a, h+j) via the alternating binomial
	// identity sum_i (-1)^i C(j,i) C(x-i,r) = C(x-j,r-j), with x=j+a, r=a-h.
	// Combined with (-1)^{l-j}: total sign = (-1)^{l+h} = (-1)^{(n+l+1)/2},
	// constant in j (n+l+1 is always even since n and l have opposite parities).
	// C(a, h+j) is zero for h+j < 0 or h+j > a, giving tighter loop bounds.
	// All factors are integers; accumulate as mpz_class, apply outer at the end.
	const long sgn = powm1((long) (n + l + 1) / 2);
	const int j_min = std::max(0, -h);             // ensures h+j >= 0
	const int j_max = std::min((l - 1) / 2, a - h); // a-h = (k-n-1)/2; ensures h+j <= a

	mpz_class acc(0);
	for (int j = j_min; j <= j_max; ++j) {
		mpz_class bj;
		mpz_bin_uiui(bj.get_mpz_t(), (unsigned long) a, (unsigned long) (h + j));

		// s(j) = sum_{q=2j}^{l-1} 2^(q-2j) * C(q-j, j) * U[q], incremental.
		mpz_class s_sum(0), shift(1), bqj(1);
		for (int q = 2 * j, p = j; q < l; ++q, ++p) {
			s_sum += shift * bqj * U[q];
			shift *= 2;
			bqj *= (p + 1);        // C(p+1, j) = C(p, j) * (p+1) / (p+1-j).
			bqj /= (p + 1 - j);   // Exact integer division.
		}

		acc += bj * s_sum;
	}

	mpq_class result = outer * mpq_class(mpz_class(sgn)) * mpq_class(acc);
	result.canonicalize();

	mpq_class ret = result;
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

static E2Poly c_phi_pow_evo_e2poly_se4(int n, int k, int i) {
	static UintsCache<E2Poly> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) i)) return *ret;

	static std::map<std::pair<int,int>, std::shared_ptr<TSeriesBase<LExpr>>> series_cache;
	auto key = std::make_pair(n, i);
	if (!series_cache.count(key))
		series_cache[key] = double_series_power_coeff_lexpr(n, i);

	LExpr lp = series_cache[key]->getItem(k);
	E2Poly result = lexpr_eval_e2poly(lp, [](int j, int m) {
		return a_nk_C_lexpr(j, m, c_phi_evo);
	});

	return cache.insert(result, (uint) n, (uint) k, (uint) i);
}

mpq_class c_phi_pow_evo_se4(int n, int k, int l, int i) {
	E2Poly ep = c_phi_pow_evo_e2poly_se4(n, k, i);
	return (l < (int) ep.size()) ? ep[l] : mpq_class(0);
}

mpq_class c_phi_pow_evo(int n, int k, int l, int i) {
	assert(n >= 0 && k >= n + i && l >= i && l <= k && i >= 0);

	// Parity constraints from underlying c_phi_evo coefficients.
	if ((i + n - k) % 2 != 0 || (l - k) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l, (uint) i))
		return *ret;

	mpq_class ret = c_phi_pow_evo_se4(n, k, l, i);  // ← LExpr/GMP pipeline (was: _se2)

	return cache.insert(ret, (uint) n, (uint) k, (uint) l, (uint) i);
}

mpq_class c_sin_phi_evo(int n, int k, int l) {
	assert(n >= 0 && k >= n && l >= 0 && l <= k);

	// Parity constraints.
	if ((n - k) % 2 != 0 || (n - l) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class d(0);
	for (int i = 0; i <= l / 2; ++i) {
		mpz_class fact;
		mpz_fac_ui(fact.get_mpz_t(), (unsigned long)(2 * i));
		// j lower bound: ceil((n + 2*i - k) / 2) = (n + 2*i - k + 1) / 2
		const int j_min = std::max(0, (n + 2 * i - k + 1) / 2);
		const int j_max = n / 2;
		for (int j = j_min; j <= j_max; ++j) {
			mpz_class binom;
			mpz_bin_uiui(binom.get_mpz_t(), (unsigned long) i, (unsigned long) j);
			mpq_class coeff(mpz_class(powm1(i + j)) * binom, fact);
			coeff.canonicalize();
			d += coeff * c_phi_pow_evo(n - 2 * j, k, l, 2 * i);
		}
	}

	mpq_class ret = d;
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

mpq_class c_cos_phi_evo(int n, int k, int l) {
	assert(n >= 0 && k >= n && l >= 1 && l <= k);

	// Parity constraints.
	if ((n + 1 - k) % 2 != 0 || (l - k) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class d(0);
	for (int i = 0; i <= (l - 1) / 2; ++i) {
		mpz_class fact;
		mpz_fac_ui(fact.get_mpz_t(), (unsigned long)(2 * i + 1));
		// j lower bound: ceil((n + 2*i + 1 - k) / 2) = (n + 2*i + 1 - k + 1) / 2
		const int j_min = std::max(0, (n + 2 * i + 2 - k) / 2);
		const int j_max = std::min(i, n / 2);
		for (int j = j_min; j <= j_max; ++j) {
			mpz_class binom;
			mpz_bin_uiui(binom.get_mpz_t(), (unsigned long) i, (unsigned long) j);
			mpq_class coeff(mpz_class(powm1(i + 1 + j)) * binom, fact);
			coeff.canonicalize();
			d += coeff * c_phi_pow_evo(n - 2 * j, k, l, 2 * i + 1);
		}
	}

	mpq_class ret = d;
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

mpq_class c_sin_phi_inv_evo(int n, int k, int l) {
	assert(n >= 0 && k >= n && l >= 0 && l <= k);

	// Parity constraints — same as d_sin_phi_evo.
	if ((n - k) % 2 != 0 || (n - l) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class d(0);
	for (int i = 0; i <= l / 2; ++i) {
		mpz_class fact;
		mpz_fac_ui(fact.get_mpz_t(), (unsigned long)(2 * i));
		// j lower bound: ceil((n + 2*i - k) / 2) = (n + 2*i - k + 1) / 2
		const int j_min = std::max(0, (n + 2 * i - k + 1) / 2);
		const int j_max = n / 2;
		for (int j = j_min; j <= j_max; ++j) {
			mpz_class binom;
			mpz_bin_uiui(binom.get_mpz_t(), (unsigned long) i, (unsigned long) j);
			mpq_class coeff(E2(i) * mpz_class(powm1(j)) * binom, fact);
			coeff.canonicalize();
			d += coeff * c_phi_pow_evo(n - 2 * j, k, l, 2 * i);
		}
	}

	mpq_class ret = d;
	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

mpq_class a_mr(int m, int r) {
	assert(m >= 0 && r >= 0 && r <= m);

	if (m == 0 && r == 0)
		return {1, 1};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) m, (uint) r))
		return *ret;

	mpq_class a(0);

	for (int k = r; k <= m; ++k) {
		mpq_class b(0);

		for (int t = 1; t <= k; ++t) {
			mpz_class b2k_kt, t_2m;
			mpz_bin_uiui(b2k_kt.get_mpz_t(), (unsigned long) (2 * k), (unsigned long) (k - t));
			mpz_ui_pow_ui(t_2m.get_mpz_t(), (unsigned long) t, (unsigned long) (2 * m));
			b += mpq_class(mpz_class(powm1(k - t)) * b2k_kt * t_2m);
		}

		mpz_class s1 = stirling1_signed((uint) k, (uint) r);
		mpz_class fact_k, shift_kr;
		mpz_fac_ui(fact_k.get_mpz_t(), (unsigned long) k);
		mpz_ui_pow_ui(shift_kr.get_mpz_t(), 2, (unsigned long) (k - r));
		a += b * mpq_class(s1) / mpq_class(fact_k * shift_kr);
	}

	mpz_class fact_2m;
	mpz_fac_ui(fact_2m.get_mpz_t(), (unsigned long) (2 * m));
	mpq_class result = mpq_class(mpz_class(2 * powm1(m))) * a / mpq_class(fact_2m);
	result.canonicalize();

	mpq_class ret = result;
	return cache.insert(ret, (uint) m, (uint) r);
}

mpq_class B_rt(int r, int t) {
	assert(r >= 0 && t >= 0 && r >= t);

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) r, (uint) t))
		return *ret;

	mpq_class B(0);

	for (int k = t; k <= r; ++k) {
		mpz_class s2 = stirling2((uint) r, (uint) k);
		mpz_class bk_t;
		mpz_bin_uiui(bk_t.get_mpz_t(), (unsigned long) k, (unsigned long) t);
		B += mpq_class(s2 * mpz_class(powm1(k - t)) * bk_t) * rf_half(1, k);
	}

	mpq_class ret = B;
	return cache.insert(ret, (uint) r, (uint) t);
}

mpq_class C_mt(int m, int t) {
	assert(m >= 0 && t >= 0 && t <= m);

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) m, (uint) t))
		return *ret;

	mpq_class C(0);
	for (int r = t; r <= m; ++r)
		C += a_mr(m, r) * B_rt(r, t);

	mpq_class ret = C;
	return cache.insert(ret, (uint) m, (uint) t);
}

mpq_class R(int n, int k, int l, int i) {
	assert(n >= 0 && k >= n && l >= 0 && l <= k && i >= 0 && i <= l / 2);

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l, (uint) i))
		return *ret;

	mpq_class s(0);

	const int j_min = std::max(0, (n + 2 * i - k + 1) / 2);
	const int j_max = n / 2;

	for (int j = j_min; j <= j_max; ++j) {
		mpz_class binom;
		mpz_bin_uiui(binom.get_mpz_t(), (unsigned long) i, (unsigned long) j);
		s += mpq_class(mpz_class(powm1(j)) * binom)
			 * c_phi_pow_evo(n - 2 * j, k, l, 2 * i);
	}

	mpq_class ret = s;
	return cache.insert(ret, (uint) n, (uint) k, (uint) l, (uint) i);
}

mpq_class B_p(int n, int k, int p) {
	assert(n >= 0 && k >= n && p >= 0 && p <= k + 1);

	if ((n - k) % 2 != 0 || (n - p - 1) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) p))
		return *ret;

	mpq_class d(0);

	for (int l = 0; l <= k; ++l) {
		for (int i = 0; i <= l / 2; ++i) {
			const int t = l + 1 - p;

			if (t % 2 == 0 && t >= 0 && t <= 2 * i)
				d += C_mt(i, t / 2) * R(n, k, l, i);
		}
	}

	mpq_class ret = d;
	return cache.insert(ret, (uint) n, (uint) k, (uint) p);
}

mpq_class cp_evo_nkl(int n, int k, int l) {
	assert(n >= 1 && k >= n && l >= 0 && l <= k + 1);

	if ((n - k) % 2 != 0 || (n - l - 1) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class ret(0);

	if (l <= 1) {
		ret = c_sin_phi_inv_evo(n - 1, k - 1, l);
	} else if (2 <= l && l <= k - 1) {
		ret = c_sin_phi_inv_evo(n - 1, k - 1, l)
			  - c_sin_phi_inv_evo(n - 1, k - 1, l - 2);
	} else if (k <= l) {
		ret = -c_sin_phi_inv_evo(n - 1, k - 1, l - 2);
	}

	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

mpq_class c_h_evo(int n, int k, int l) {
	assert(n >= 0 && k >= n && l >= 0 && l <= k + 1);

	if ((n - k) % 2 != 0 || (n - l - 1) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	mpq_class ret(0);

	if (n == 0) {
		ret = -B_p(n, k, l);
	} else if (l == 0) {
		ret = cp_evo_nkl(n, k, l);
	} else {
		ret = cp_evo_nkl(n, k, l) - B_p(n, k, l);
	}

	return cache.insert(ret, (uint) n, (uint) k, (uint) l);
}

mpq_class d_h_evo(int n, int k, int l) {
	return c_h_evo(n, 2 * k + n, 2 * l + 1 - (n % 2));
}

mpq_class d_sin_phi_evo(int n, int k, int l) {
	return c_sin_phi_evo(n, n + 2 * k, 2 * l + (n % 2));
}

mpq_class d_cos_phi_evo(int n, int k, int l) {
	return c_cos_phi_evo(n, n + 1 + 2 * k, 2 * l + 1 - (n % 2));
}

}  // namespace point_to_ellipse_series
