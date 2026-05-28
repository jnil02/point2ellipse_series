#pragma once

#include <unordered_map>
#include <gmpxx.h>

#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/pow.h>
#include <symengine/mp_class.h>

#include "coefficients.hpp"

namespace point_to_ellipse_series {

using SymEngine::Expression;
using SymEngine::Symbol;
using SymEngine::Add;
using SymEngine::Basic;
using SymEngine::RCP;
using SymEngine::Mul;
using SymEngine::Pow;
using SymEngine::Integer;

// Convert a SymEngine Integer to mpz_class via its decimal string representation.
inline mpz_class se_integer_to_mpz(const Integer& i) {
	mpz_class result;
	result.set_str(i.__str__(), 10);
	return result;
}

inline std::string mpz_to_str(const mpz_class& x) { return x.get_str(10); }

inline static Expression mpq_to_expr(const mpq_class& d)
{
	SymEngine::rational_class q(d.get_mpq_t());
	return Expression(SymEngine::Rational::from_mpq(std::move(q)));
}

inline static mpq_class expr_to_mpq(const Expression& d)
{
	RCP<const Basic> b = d.get_basic();

	if (is_a<SymEngine::Rational>(*b)) {
		auto r = rcp_static_cast<const SymEngine::Rational>(b);
		mpq_class q;
		mpq_set_num(q.get_mpq_t(), se_integer_to_mpz(*r->get_num()).get_mpz_t());
		mpq_set_den(q.get_mpq_t(), se_integer_to_mpz(*r->get_den()).get_mpz_t());
		return q;
	}
	else if (is_a<Integer>(*b)) {
		auto i = rcp_static_cast<const Integer>(b);
		return mpq_class(se_integer_to_mpz(*i));
	}
	else {
		throw std::runtime_error("Expression is not a Rational or Integer");
	}
}

/** Compute the n:th rising factorial of (k/2).
 *
 * @param k
 * @param n
 * @return
 */
inline mpq_class rf_half(unsigned long k, unsigned long n) {
	mpz_class r(1), s(k);
	for (unsigned long i = 0; i < n; ++i) { r *= s; s += 2; }
	mpz_class denom;
	mpz_ui_pow_ui(denom.get_mpz_t(), 2, n);
	mpq_class result(r, denom);
	result.canonicalize();
	return result;
}

/** (-1)^n
 *
 * See
 * https://stackoverflow.com/questions/29110752/what-is-the-correct-way-to-obtain-1n
 *
 * @param n
 * @return
 */
inline long powm1(long n) {
	return 1L - ((n & 1L) << 1);
}

/** Generalised integer binomial C(n, k) for arbitrary integer n and any k.
 *
 * Returns 0 for k < 0. For n < 0 uses C(n,k) = (-1)^k * C(k-n-1, k).
 *
 * @param n Integer (may be negative).
 * @param k Integer.
 * @return
 */
inline mpq_class int_bin(long n, long k) {
	if (k < 0) return mpq_class(0);
	if (n >= 0) {
		if (k > n) return mpq_class(0);
		mpz_class result;
		mpz_bin_uiui(result.get_mpz_t(), (unsigned long) n, (unsigned long) k);
		return mpq_class(result);
	}
	// n < 0: C(n,k) = (-1)^k * C(k-n-1, k)
	mpz_class result;
	mpz_bin_uiui(result.get_mpz_t(), (unsigned long) (k - n - 1), (unsigned long) k);
	return mpq_class(mpz_class(powm1(k)) * result);
}

/** Compute the binomial coefficient for rational n and integer k.
 *
 * Computes (n-k+1)*(n-k+2)*...*n / k!  Returns 0 for k < 0.
 *
 * @param n Rational number.
 * @param k Integer.
 * @return
 */
inline mpq_class binomial_rational(const mpq_class& n, long k) {
	if (k < 0) return mpq_class(0);
	if (k == 0) return mpq_class(1);
	mpq_class result(1);
	mpq_class d = n - mpq_class(k - 1);
	for (long i = 0; i < k; ++i) {
		result *= d;
		d += 1;
	}
	mpz_class fact;
	mpz_fac_ui(fact.get_mpz_t(), (unsigned long) k);
	result /= mpq_class(fact);
	result.canonicalize();
	return result;
}

/** Euler secant number E_{2n}.
 *
 * Mirrors Python E2(n). Uses a static cache for efficiency.
 * Defined by: E2(0)=1, E2(n) = sum_{j=0}^{n-1} (-1)^(n-j+1) * C(2n,2j) * E2(j)
 *
 * @param n Non-negative integer.
 * @return E_{2n} as mpz_class.
 */
inline mpz_class E2(int n) {
	static std::unordered_map<int, mpz_class> cache;
	auto it = cache.find(n);
	if (it != cache.end())
		return it->second;

	if (n == 0)
		return mpz_class(1);

	mpz_class result(0);
	for (int j = 0; j < n; ++j) {
		mpz_class binom;
		mpz_bin_uiui(binom.get_mpz_t(), (unsigned long)(2 * n), (unsigned long)(2 * j));
		result += mpz_class(powm1(n - j + 1)) * binom * E2(j);
	}

	cache[n] = result;
	return result;
}

// Extract the coefficient of `sym^exp` from an expanded expression.
static Expression
coeff_of(const Expression &expr, const RCP<const Symbol> &sym, int exp) {
	const RCP<const Basic> b = expr.get_basic();

	// Sum of terms: recurse over args and add results.
	if (is_a<Add>(*b)) {
		Expression acc(0);
		for (const auto &term: b->get_args()) {
			acc = acc + coeff_of(Expression(term), sym, exp);
		}
		return acc;
	}

	// For a single term, gather the power of `sym` and the remaining factor.
	long power = 0;
	Expression rest(1);

	if (is_a<Mul>(*b)) {
		for (const auto &f: b->get_args()) {
			if (is_a<Pow>(*f)) {
				auto pw = rcp_static_cast<const Pow>(f);
				const RCP<const Basic> base = pw->get_base();
				if (eq(*base, *sym)) {
					const auto &iexp = SymEngine::down_cast<const Integer &>(
							*pw->get_exp());
					power += iexp.as_int();
				} else {
					rest = rest * Expression(f);
				}
			} else if (is_a<Symbol>(*f) && eq(*f, *sym)) {
				power += 1;
			} else {
				rest = rest * Expression(f);
			}
		}
	} else if (is_a<Pow>(*b)) {
		auto pw = rcp_static_cast<const Pow>(b);
		const RCP<const Basic> base = pw->get_base();
		if (eq(*base, *sym)) {
			const auto &iexp = SymEngine::down_cast<const Integer &>(
					*pw->get_exp());
			power += iexp.as_int();
		} else {
			rest = Expression(b);
		}
	} else if (is_a<Symbol>(*b) && eq(*b, *sym)) {
		power = 1;
	} else {
		rest = Expression(b);
	}

	return (power == exp) ? rest : Expression(0);
}

static Expression
coeff_of2(const Expression &expr, const RCP<const Symbol> &sym, int exp) {
	const RCP<const Basic> b = expr.get_basic();

	if (is_a<Add>(*b)) {
		Expression acc(0);
		for (const auto &term: b->get_args()) {
			acc = acc + coeff_of2(Expression(term), sym, exp);
		}
		return acc;
	}

	long power = 0;
	Expression rest(1);

	if (is_a<Mul>(*b)) {
		for (const auto &f: b->get_args()) {
			if (is_a<Pow>(*f)) {
				auto pw = rcp_static_cast<const Pow>(f);
				const RCP<const Basic> base = pw->get_base();
				if (eq(*base, *sym)) {
					const auto &iexp = SymEngine::down_cast<const Integer &>(
							*pw->get_exp());
					power += iexp.as_int();
				} else {
					rest = rest * Expression(f);
				}
			} else if (is_a<Symbol>(*f) && eq(*f, *sym)) {
				power += 1;
			} else {
				rest = rest * Expression(f);
			}
		}
	} else if (is_a<Pow>(*b)) {
		auto pw = rcp_static_cast<const Pow>(b);
		const RCP<const Basic> base = pw->get_base();
		if (eq(*base, *sym)) {
			const auto &iexp = SymEngine::down_cast<const Integer &>(
					*pw->get_exp());
			power += iexp.as_int();
		} else {
			rest = Expression(b);
		}
	} else if (is_a<Symbol>(*b) && eq(*b, *sym)) {
		power = 1;
	} else {
		rest = Expression(b);
	}

	return (power == exp) ? rest : Expression(0);
}

} // namespace point_to_ellipse_series
