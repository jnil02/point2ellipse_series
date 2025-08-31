#pragma once

#include <cstdint>

using SymEngine::Expression;
using SymEngine::rational;
using SymEngine::Symbol;
using SymEngine::Add;
using SymEngine::Basic;
using SymEngine::RCP;
using SymEngine::factorial;
using SymEngine::Mul;
using SymEngine::Pow;
using SymEngine::Integer;

/** Compute the n:th rising factorial of (k/2).
 *
 * @param k
 * @param n
 * @return
 */
inline Expression rf_half(unsigned long k, unsigned long n) {
	long r = 1L;
	long s = (long) k;
	for (int i = 0; i < n; ++i) {
		r = r * s;
		s += 2L;
	}
	return {rational(r, 1L << n)};
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

/** Compute the binomial coefficient for (n,k) where n is a rational number.
 *
 * @param n Rational number.
 * @param k Integer. k != 1 and n != k  FIXME(JO) Check what this means.
 * @return
 */
inline Expression binomial_rational(RCP<const Basic> n, unsigned long k) {
	Expression d = Expression(n) - k;
	Expression result(1);
	for (int i = 1; i < k + 1; i++) {
		d += 1;
		result *= d;
	}
	return result / factorial(k);
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

