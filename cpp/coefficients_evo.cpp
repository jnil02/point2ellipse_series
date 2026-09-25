
#include "detail/coefficients_evo_internal.hpp"
#include "detail/cache.hpp"
#include "detail/util.hpp"
#include "detail/series_substitution.hpp"
#include "detail/stirling.hpp"

#include <cassert>
#include <map>

namespace point_to_ellipse_series {

using uint = unsigned int;

mpq_class d_phi_evo(int k, int l, int n) {
	assert(l >= 0 && k >= 1 && n >= 0 && l <= (k - 1) / 2 && n <= (k - 1) / 2);

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) k, (uint) l, (uint) n))
		return *ret;

	const int s = k % 2;
	const int m_k = (k - 1) / 2;

	// The outer factor contains C(k/2, m_k-n). For a negative lower index
	// SymPy's binomial is zero, so the complete coefficient is zero.
	if (n > m_k) {
		mpq_class ret(0);
		return cache.insert(ret, (uint) k, (uint) l, (uint) n);
	}

	const int Q = 1 - s + 2 * n;  // Maximum value of q.

	// Precompute the q-only factor
	//
	//    T[q] = C(k-2-q, Q-q).
	//
	// Since n <= m_k, Q <= k-1. For q < Q the upper argument is therefore
	// non-negative. At q == Q the lower argument is zero, so T[Q] = 1 even
	// in the one case where the upper argument is -1.
	std::vector<mpz_class> T(Q + 1);
	for (int q = 0; q <= Q; ++q) {
		if (Q - q == 0)
			T[q] = 1;
		else
			mpz_bin_uiui(T[q].get_mpz_t(),
						 (unsigned long) (k - 2 - q),
						 (unsigned long) (Q - q));
	}

	mpz_class d(0);

	// The j-only factor is
	//
	//    C(m_k-n, l-n+j).
	//
	// Since m_k-n >= 0, it is nonzero only when
	//
	//    0 <= l-n+j <= m_k-n,
	//
	// giving the bounds below.
	const int j_min = std::max(0, n - l);
	const int j_max = std::min(n, m_k - l);

	for (int j = j_min; j <= j_max; ++j) {
		mpz_class bj;
		mpz_bin_uiui(bj.get_mpz_t(),
					 (unsigned long) (m_k - n),
					 (unsigned long) (l - n + j));

		// Inner sum:
		//
		//   sum_{q=2j}^{Q}
		//      2^(q-2j) * C(q-j, j) * T[q].
		//
		// Advance both q-dependent factors incrementally.
		mpz_class sum_(0);
		mpz_class shift(1);  // 2^(q-2j), starts at q = 2j.
		mpz_class bqj(1);    // C(q-j,j) = C(j,j) = 1 at q = 2j.

		for (int q = 2 * j, p = j; q <= Q; ++q, ++p) {
			sum_ += shift * bqj * T[q];

			shift *= 2;

			// C(p+1,j) = C(p,j) * (p+1) / (p+1-j).
			bqj *= (p + 1);
			bqj /= (p + 1 - j);  // Exact integer division.
		}

		d += sum_ * bj;
	}

	// Outer factor:
	//
	//   (-1)^(l+n+k) / k * C(k/2, m_k-n).
	//
	// Everything accumulated above is integral, so rational arithmetic is
	// deferred until here.
	mpq_class result =
			binomial_rational(mpq_class(k, 2), (long) (m_k - n))
			* mpq_class(powm1(l + n + k))
			/ mpq_class(mpz_class(k))
			* mpq_class(d);
	result.canonicalize();

	mpq_class ret = result;
	return cache.insert(ret, (uint) k, (uint) l, (uint) n);
}

mpq_class c_phi_evo(int k, int l, int n) {
	assert(l >= 0 && k >= l + 1 && n >= 1 && n <= k);

	// Parity constraints — coefficient is zero unless both hold.
	if ((k - l - 1) % 2 != 0 || (n - k) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) l, (uint) k, (uint) n))
		return *ret;

	// Let a = (k-l)/2 and h = (n-l+1)/2 (both integers by parity).
	// Outer rational factor: binomial(k/2, a) / k.
	const int a = (k - n) / 2;
	const int h = (l - n + 1) / 2;
	mpq_class outer = binomial_rational(mpq_class(k, 2), (long) a)
					  / mpq_class(mpz_class(k));

	// Precompute q-only factor U[q] = C(k-2-q, l-1-q), independent of j.
	// l-1-q is always >= 0 for q in [0,l-1]. When q == l-1 the first arg may
	// be -1 (only when k == l); C(-1,0) = 1, handled explicitly.
	std::vector<mpz_class> U(n);
	for (int q = 0; q < n; ++q) {
		if (n - 1 - q == 0) U[q] = 1;
		else mpz_bin_uiui(U[q].get_mpz_t(), (unsigned long) (k - 2 - q),
						  (unsigned long) (n - 1 - q));
	}

	// b(j) collapses to (-1)^{h+j} * C(a, h+j) via the alternating binomial
	// identity sum_i (-1)^i C(j,i) C(x-i,r) = C(x-j,r-j), with x=j+a, r=a-h.
	// Combined with (-1)^{l-j}: total sign = (-1)^{l+h} = (-1)^{(n+l+1)/2},
	// constant in j (n+l+1 is always even since n and l have opposite parities).
	// C(a, h+j) is zero for h+j < 0 or h+j > a, giving tighter loop bounds.
	// All factors are integers; accumulate as mpz_class, apply outer at the end.
	const mpz_class sgn = powm1((long) (l + n + 1) / 2);
	const int j_min = std::max(0, -h);             // ensures h+j >= 0
	const int j_max = std::min((n - 1) / 2, a - h); // a-h = (k-n-1)/2; ensures h+j <= a

	mpz_class acc(0);
	for (int j = j_min; j <= j_max; ++j) {
		mpz_class bj;
		mpz_bin_uiui(bj.get_mpz_t(), (unsigned long) a, (unsigned long) (h + j));

		// s(j) = sum_{q=2j}^{l-1} 2^(q-2j) * C(q-j, j) * U[q], incremental.
		mpz_class s_sum(0), shift(1), bqj(1);
		for (int q = 2 * j, p = j; q < n; ++q, ++p) {
			s_sum += shift * bqj * U[q];
			shift *= 2;
			bqj *= (p + 1);        // C(p+1, j) = C(p, j) * (p+1) / (p+1-j).
			bqj /= (p + 1 - j);   // Exact integer division.
		}

		acc += bj * s_sum;
	}

	mpq_class result = outer * mpq_class(sgn) * mpq_class(acc);
	result.canonicalize();

	mpq_class ret = result;
	return cache.insert(ret, (uint) l, (uint) k, (uint) n);
}

static E2Poly c_phi_pow_evo_e2poly_se4(int k, int l, int i) {
	static UintsCache<E2Poly> cache;
	if (auto *ret = cache.get((uint) k, (uint) l, (uint) i)) return *ret;

	static std::map<std::pair<int,int>, std::shared_ptr<TSeriesBase<LExpr>>> series_cache;
	auto key = std::make_pair(l, i);
	if (!series_cache.count(key))
		// Inside-evolute generator a_l starts at z^{l+1} implying a tightened
		// variant (Bell guard k>=i*(l+1)).
		series_cache[key] = double_series_power_coeff_evo_lexpr(l, i);

	LExpr lp = series_cache[key]->getItem(k);
	E2Poly result = lexpr_eval_e2poly(lp, [](int j, int m) {
		return a_nk_C_lexpr(j, m, c_phi_evo);
	});

	return cache.insert(result, (uint) k, (uint) l, (uint) i);
}

mpq_class c_phi_pow_evo_se4(int k, int l, int n, int i) {
	E2Poly ep = c_phi_pow_evo_e2poly_se4(k, l, i);
	return (n < (int) ep.size()) ? ep[n] : mpq_class(0);
}

mpq_class c_phi_pow_evo(int k, int l, int n, int i) {
	assert(l >= 0 && k >= l + i && n >= i && n <= k && i >= 0);

	// Parity constraints from underlying c_phi_evo coefficients.
	if ((i + l - k) % 2 != 0 || (n - k) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) l, (uint) k, (uint) n, (uint) i))
		return *ret;

	mpq_class ret = c_phi_pow_evo_se4(k, l, n, i);  // ← LExpr/GMP pipeline (was: _se2)

	return cache.insert(ret, (uint) l, (uint) k, (uint) n, (uint) i);
}

mpq_class c_sin_phi_evo(int k, int l, int n) {
	assert(l >= 0 && k >= l && n >= 0 && n <= k);

	// Parity constraints.
	if ((l - k) % 2 != 0 || (l - n) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) l, (uint) k, (uint) n))
		return *ret;

	mpq_class d(0);
	for (int i = 0; i <= n / 2; ++i) {
		mpz_class fact;
		mpz_fac_ui(fact.get_mpz_t(), (unsigned long)(2 * i));
		// j lower bound: ceil((n + 2*i - k) / 2) = (n + 2*i - k + 1) / 2
		const int j_min = std::max(0, (l + 2 * i - k + 1) / 2);
		const int j_max = l / 2;
		for (int j = j_min; j <= j_max; ++j) {
			mpz_class binom;
			mpz_bin_uiui(binom.get_mpz_t(), (unsigned long) i, (unsigned long) j);
			mpq_class coeff(powm1(i + j) * binom, fact);
			coeff.canonicalize();
			d += coeff * c_phi_pow_evo(k, l - 2 * j, n, 2 * i);
		}
	}

	mpq_class ret = d;
	return cache.insert(ret, (uint) l, (uint) k, (uint) n);
}

mpq_class c_cos_phi_evo(int k, int l, int n) {
	assert(l >= 0 && k >= l && n >= 1 && n <= k);

	// Parity constraints.
	if ((l + 1 - k) % 2 != 0 || (n - k) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) l, (uint) k, (uint) n))
		return *ret;

	mpq_class d(0);
	for (int i = 0; i <= (n - 1) / 2; ++i) {
		mpz_class fact;
		mpz_fac_ui(fact.get_mpz_t(), (unsigned long)(2 * i + 1));
		// j lower bound: ceil((n + 2*i + 1 - k) / 2) = (n + 2*i + 1 - k + 1) / 2
		const int j_min = std::max(0, (l + 2 * i + 2 - k) / 2);
		const int j_max = std::min(i, l / 2);
		for (int j = j_min; j <= j_max; ++j) {
			mpz_class binom;
			mpz_bin_uiui(binom.get_mpz_t(), (unsigned long) i, (unsigned long) j);
			mpq_class coeff(powm1(i + 1 + j) * binom, fact);
			coeff.canonicalize();
			d += coeff * c_phi_pow_evo(k, l - 2 * j, n, 2 * i + 1);
		}
	}

	mpq_class ret = d;
	return cache.insert(ret, (uint) l, (uint) k, (uint) n);
}

mpq_class c_sin_phi_inv_evo(int k, int l, int n) {
	assert(l >= 0 && k >= l && n >= 0 && n <= k);

	// Parity constraints — same as d_sin_phi_evo.
	if ((l - k) % 2 != 0 || (l - n) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) l, (uint) k, (uint) n))
		return *ret;

	mpq_class d(0);
	for (int i = 0; i <= n / 2; ++i) {
		mpz_class fact;
		mpz_fac_ui(fact.get_mpz_t(), (unsigned long)(2 * i));
		// j lower bound: ceil((n + 2*i - k) / 2) = (n + 2*i - k + 1) / 2
		const int j_min = std::max(0, (l + 2 * i - k + 1) / 2);
		const int j_max = l / 2;
		for (int j = j_min; j <= j_max; ++j) {
			mpz_class binom;
			mpz_bin_uiui(binom.get_mpz_t(), (unsigned long) i, (unsigned long) j);
			mpq_class coeff(E2(i) * powm1(j) * binom, fact);
			coeff.canonicalize();
			d += coeff * c_phi_pow_evo(k, l - 2 * j, n, 2 * i);
		}
	}

	mpq_class ret = d;
	return cache.insert(ret, (uint) l, (uint) k, (uint) n);
}

namespace detail {

mpq_class a_mr(int m, int r) {
	assert(r >= 0 && r <= m);

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
			b += mpq_class(powm1(k - t) * b2k_kt * t_2m);
		}

		mpz_class s1 = stirling1_signed((uint) k, (uint) r);
		mpz_class fact_k, shift_kr;
		mpz_fac_ui(fact_k.get_mpz_t(), (unsigned long) k);
		mpz_ui_pow_ui(shift_kr.get_mpz_t(), 2, (unsigned long) (k - r));
		a += b * mpq_class(s1) / mpq_class(fact_k * shift_kr);
	}

	mpz_class fact_2m;
	mpz_fac_ui(fact_2m.get_mpz_t(), (unsigned long) (2 * m));
	mpq_class result = mpq_class(2 * powm1(m)) * a / mpq_class(fact_2m);
	result.canonicalize();

	mpq_class ret = result;
	return cache.insert(ret, (uint) m, (uint) r);
}

mpq_class B_rt(int r, int t) {
	assert(t >= 0 && r >= t);

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) r, (uint) t))
		return *ret;

	mpq_class B(0);

	for (int k = t; k <= r; ++k) {
		mpz_class s2 = stirling2((uint) r, (uint) k);
		mpz_class bk_t;
		mpz_bin_uiui(bk_t.get_mpz_t(), (unsigned long) k, (unsigned long) t);
		B += mpq_class(s2 * powm1(k - t) * bk_t) * rf_half(1, k);
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

mpq_class R(int k, int l, int n, int i) {
	assert(l >= 0 && k >= l && n >= 0 && n <= k && i >= 0 && i <= n / 2);

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) k, (uint) l, (uint) n, (uint) i))
		return *ret;

	mpq_class s(0);

	const int j_min = std::max(0, (l + 2 * i - k + 1) / 2);
	const int j_max = l / 2;

	for (int j = j_min; j <= j_max; ++j) {
		mpz_class binom;
		mpz_bin_uiui(binom.get_mpz_t(), (unsigned long) i, (unsigned long) j);
		s += mpq_class(powm1(j) * binom)
			 * c_phi_pow_evo(k, l - 2 * j, n, 2 * i);
	}

	mpq_class ret = s;
	return cache.insert(ret, (uint) k, (uint) l, (uint) n, (uint) i);
}

} // namespace detail

mpq_class c_N_evo(int k, int l, int n) {
	assert(l >= 0 && k >= l && n >= 0 && n <= k + 1);

	if ((l - k) % 2 != 0 || (l - n - 1) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) l, (uint) k, (uint) n))
		return *ret;

	mpq_class d(0);

	for (int p = 0; p <= k; ++p) {
		for (int i = 0; i <= p / 2; ++i) {
			const int t = p + 1 - n;

			if (t % 2 == 0 && t >= 0 && t <= 2 * i)
				d += detail::C_mt(i, t / 2) * detail::R(k, l, p, i);
		}
	}

	mpq_class ret = d;
	return cache.insert(ret, (uint) l, (uint) k, (uint) n);
}

mpq_class cp_evo(int k, int l, int n) {
	assert(l >= 1 && k >= l && n >= 1 && n <= k + 1);

	if ((l - k) % 2 != 0 || (k + 1 - n) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) l, (uint) k, (uint) n))
		return *ret;

	mpq_class ret(0);

	if (n <= 2 && n <= k-1)
		// k = 0, l=0, n = 2
		ret = c_sin_phi_inv_evo(k - 1, l - 1, n);
	else if (3 <= n && n <= k-1)
		ret = c_sin_phi_inv_evo(k - 1, l - 1, n) - c_sin_phi_inv_evo(k - 1, l - 1, n - 2);
	else if (3 <= n && n <= k+1)
		ret = -c_sin_phi_inv_evo(k - 1, l - 1, n - 2);
	// else (ret = 0) will happen for e.g. 1,1,2.

	return cache.insert(ret, (uint) l, (uint) k, (uint) n);
}

mpq_class c_h_evo(int k, int l, int n) {
	assert(l >= 0 && k >= l && n >= 1 && n <= k + 1);

	if ((l - k) % 2 != 0 || (l - n - 1) % 2 != 0)
		return {0};

	static UintsCache<mpq_class> cache;
	if (auto *ret = cache.get((uint) l, (uint) k, (uint) n))
		return *ret;

	mpq_class ret(0);

	if (l == 0) {
		ret = -c_N_evo(k, l, n);
	} else {
		ret = cp_evo(k, l, n) - c_N_evo(k, l, n);
	}

	return cache.insert(ret, (uint) l, (uint) k, (uint) n);
}

mpq_class d_h_evo(int k, int l, int n) {
	assert(l >= 0 && k >= 2 && n >= 0 && l <= k / 2 && n <= k / 2);
	return c_h_evo(k, 2 * l + (k % 2), 2 * n + 1 + (k % 2));
}

mpq_class d_sin_phi_evo(int k, int l, int n) {
	assert(l >= 0 && k >= 2 && n >= 1 && l <= k / 2 && n <= k / 2);
	return c_sin_phi_evo(k, 2 * l + (k % 2), 2 * n + (k % 2));
}

mpq_class d_cos_phi_evo(int k, int l, int n) {
	assert(l >= 0 && k >= 1 && n >= (k - 1) % 2 && l <= (k - 1) / 2 && n <= k / 2);
	return c_cos_phi_evo(k, 2 * l + (k - 1) % 2, 2 * n + k % 2);
}

}  // namespace point_to_ellipse_series
