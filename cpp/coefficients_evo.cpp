
#include <symengine/expression.h>
#include <symengine/ntheory.h>
#include <symengine/rational.h>
#include <symengine/symbol.h>

#include "cache.hpp"
#include "util.hpp"
#include "series_substitution.hpp"

#include "coefficients_evo.hpp"
#include "stirling.hpp"

namespace point_to_ellipse_series {

using SymEngine::Expression;

using uint = unsigned int;

mpq_class d_phi_evo(int n, int k, int l) {
	assert(n >= 0 && k >= 0 && l >= 0 && l <= n / 2 + k);

	static auto cache = UintsCache<mpq_class>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	const int s = n % 2;
	const int m = n + 1 + 2 * k;

	mpq_class c(0);
	for (int j = 0; j <= l; ++j) {
		const int ka = n / 2 - l + j;

		// Inner sum for b.
		mpq_class b(0);
		const int i_max = std::min(j, ka);  // Empty range if ka < 0.
		for (int i = 0; i <= i_max; ++i) {
			mpz_class bj_i;
			mpz_bin_uiui(bj_i.get_mpz_t(), (unsigned long) j, (unsigned long) i);
			mpz_class bkaikk;
			mpz_bin_uiui(bkaikk.get_mpz_t(), (unsigned long) (ka - i + k),
						 (unsigned long) std::max(0, ka - i));
			b += mpq_class(mpz_class(powm1(i)) * bj_i * bkaikk);
		}

		// Inner sum for sum_.
		mpq_class sum_(0);
		for (int q = 2 * j; q <= s + 2 * l; ++q) {
			mpz_class shift_val;
			mpz_ui_pow_ui(shift_val.get_mpz_t(), 2, (unsigned long) (q - 2 * j));
			mpz_class bqj_j;
			mpz_bin_uiui(bqj_j.get_mpz_t(), (unsigned long) (q - j), (unsigned long) j);
			sum_ += mpq_class(shift_val * bqj_j) * int_bin(n - 1 + 2 * k - q, s + 2 * l - q);
		}

		c += sum_ * b;
	}

	// Outer factor: (-1)^(n/2 + l + n + 1) * binomial(m/2, n/2 + k - l) / m
	mpq_class result =
			binomial_rational(mpq_class(m, 2), (long) (n / 2 + k - l))
			* mpq_class(mpz_class(powm1(n / 2 + l + n + 1)))
			/ mpq_class(mpz_class(m))
			* c;
	result.canonicalize();

	mpq_class ret = result;
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

mpq_class c_phi_evo(int n, int k, int l) {
	assert(n >= 0 && k >= n + 1 && l >= 1 && l <= k);

	// Parity constraints — coefficient is zero unless both hold.
	if ((k - n - 1) % 2 != 0 || (l - k) % 2 != 0)
		return mpq_class(0);

	static auto cache = UintsCache<mpq_class>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l))
		return *ret;

	// Outer binomial factor: binomial(k/2, (k-l)/2) / k.
	// (k-l)/2 >= 0 since (l-k)%2==0 and valid range has l <= k.
	const long kl_2 = (long) ((k - l) / 2);
	mpq_class outer = binomial_rational(mpq_class(k, 2), kl_2)
					  / mpq_class(mpz_class(k));

	mpq_class c(0);
	for (int j = 0; j <= (l - 1) / 2; ++j) {

		// Inner sum for b.
		mpq_class b(0);
		const int i_max = std::min(j, (n + 1 - l + 2 * j) / 2);
		for (int i = 0; i <= i_max; ++i) {
			const int ka = (n - l + 1 + 2 * j - 2 * i) / 2;
			mpz_class bj_i;
			mpz_bin_uiui(bj_i.get_mpz_t(), (unsigned long) j, (unsigned long) i);
			b += mpq_class(mpz_class(powm1(ka)) * bj_i)
				 * int_bin(j - i + (k - l) / 2, ka);
		}

		// Inner sum for s.
		mpq_class s(0);
		for (int q = 2 * j; q <= l - 1; ++q) {
			mpz_class shift_val;
			mpz_ui_pow_ui(shift_val.get_mpz_t(), 2, (unsigned long) (q - 2 * j));
			mpz_class bqj_j;
			mpz_bin_uiui(bqj_j.get_mpz_t(), (unsigned long) (q - j), (unsigned long) j);
			s += mpq_class(shift_val * bqj_j) * int_bin(k - 2 - q, l - 1 - q);
		}

		c += outer * mpq_class(mpz_class(powm1(l - j))) * s * b;
	}

	mpq_class ret = c;
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

Expression d_phi_pow_evo_polynomial(int n, int k, int i) {
	assert(n >= 0 && k >= 0 && i >= 0);
	static auto cache = UintsCache<Expression>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) i))
		return *ret;

	// Series for b_{n,i}, take the k-th coefficient.
	Expression b_ni_k = double_series_power_coeff(n, i)->getItem(k);

	// Substitute a_{n,k} using a_nk_C with c_phi_evo.
	Expression d = a_nk_sub(b_ni_k, [](int n, int k) {
		return a_nk_C(n, k, c_phi_evo);
	});

	cache.insert(d, (uint) n, (uint) k, (uint) i);
	return d;
}

static mpq_class c_phi_pow_evo_se(int n, int k, int l, int i) {
	Expression poly = d_phi_pow_evo_polynomial(n, k, i);
	return expr_to_mpq(coeff_of(expand(poly).get_basic(), e2sym, l));
}

mpq_class c_phi_pow_evo(int n, int k, int l, int i) {
	assert(n >= 0 && k >= n + i && l >= i && l <= k && i >= 0);

	// Parity constraints from underlying c_phi_evo coefficients.
	if ((i + n - k) % 2 != 0 || (l - k) % 2 != 0)
		return mpq_class(0);

	static auto cache = UintsCache<mpq_class>();
	if (auto *ret = cache.get((uint) n, (uint) k, (uint) l, (uint) i))
		return *ret;

	mpq_class ret = c_phi_pow_evo_se(n, k, l, i);

	cache.insert(ret, (uint) n, (uint) k, (uint) l, (uint) i);
	return ret;
}

mpq_class c_sin_phi_evo(int n, int k, int l) {
	assert(n >= 0 && k >= n && l >= 0 && l <= k);

	// Parity constraints.
	if ((n - k) % 2 != 0 || (n - l) % 2 != 0)
		return mpq_class(0);

	static auto cache = UintsCache<mpq_class>();
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
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

mpq_class c_cos_phi_evo(int n, int k, int l) {
	assert(n >= 0 && k >= n && l >= 1 && l <= k);

	// Parity constraints.
	if ((n + 1 - k) % 2 != 0 || (l - k) % 2 != 0)
		return mpq_class(0);

	static auto cache = UintsCache<mpq_class>();
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
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

mpq_class c_sin_phi_inv_evo(int n, int k, int l) {
	assert(n >= 0 && k >= n && l >= 0 && l <= k);

	// Parity constraints — same as d_sin_phi_evo.
	if ((n - k) % 2 != 0 || (n - l) % 2 != 0)
		return mpq_class(0);

	static auto cache = UintsCache<mpq_class>();
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
	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

mpq_class a_mr(int m, int r) {
	assert(m >= 0 && r >= 0 && r <= m);

	if (m == 0 && r == 0)
		return {1, 1};

	static auto cache = UintsCache<mpq_class>();
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

		mpz_class s1 = se_integer_to_mpz(*stirling1_signed((uint) k, (uint) r));
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
	cache.insert(ret, (uint) m, (uint) r);
	return ret;
}

mpq_class B_rt(int r, int t) {
	assert(r >= 0 && t >= 0 && r >= t);

	static auto cache = UintsCache<mpq_class>();
	if (auto *ret = cache.get((uint) r, (uint) t))
		return *ret;

	mpq_class B(0);

	for (int k = t; k <= r; ++k) {
		mpz_class s2 = se_integer_to_mpz(*stirling2((uint) r, (uint) k));
		mpz_class bk_t;
		mpz_bin_uiui(bk_t.get_mpz_t(), (unsigned long) k, (unsigned long) t);
		B += mpq_class(s2 * mpz_class(powm1(k - t)) * bk_t) * rf_half(1, k);
	}

	mpq_class ret = B;
	cache.insert(ret, (uint) r, (uint) t);
	return ret;
}

mpq_class C_mt(int m, int t) {
	assert(m >= 0 && t >= 0 && t <= m);

	static auto cache = UintsCache<mpq_class>();
	if (auto *ret = cache.get((uint) m, (uint) t))
		return *ret;

	mpq_class C(0);
	for (int r = t; r <= m; ++r)
		C += a_mr(m, r) * B_rt(r, t);

	mpq_class ret = C;
	cache.insert(ret, (uint) m, (uint) t);
	return ret;
}

mpq_class R(int n, int k, int l, int i) {
	assert(n >= 0 && k >= n && l >= 0 && l <= k && i >= 0 && i <= l / 2);

	static auto cache = UintsCache<mpq_class>();
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
	cache.insert(ret, (uint) n, (uint) k, (uint) l, (uint) i);
	return ret;
}

mpq_class B_p(int n, int k, int p) {
	assert(n >= 0 && k >= n && p >= 0 && p <= k + 1);

	if ((n - k) % 2 != 0 || (n - p - 1) % 2 != 0)
		return mpq_class(0);

	static auto cache = UintsCache<mpq_class>();
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
	cache.insert(ret, (uint) n, (uint) k, (uint) p);
	return ret;
}

mpq_class cp_evo_nkl(int n, int k, int l) {
	assert(n >= 1 && k >= n && l >= 0 && l <= k + 1);

	if ((n - k) % 2 != 0 || (n - l - 1) % 2 != 0)
		return mpq_class(0);

	static auto cache = UintsCache<mpq_class>();
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

	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
}

mpq_class c_h_evo(int n, int k, int l) {
	assert(n >= 0 && k >= n && l >= 0 && l <= k + 1);

	if ((n - k) % 2 != 0 || (n - l - 1) % 2 != 0)
		return mpq_class(0);

	static auto cache = UintsCache<mpq_class>();
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

	cache.insert(ret, (uint) n, (uint) k, (uint) l);
	return ret;
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
