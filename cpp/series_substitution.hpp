#pragma once

#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/symbol.h>
#include <symengine/number.h>

#include "series.hpp"
#include "polynomials.hpp"
#include "symbols.hpp"

namespace point_to_ellipse_series {

using SymEngine::Expression;
using SymEngine::RCP;
using SymEngine::Basic;
using SymEngine::Number;
using SymEngine::Symbol;
using SymEngine::Add;
using SymEngine::Mul;
using SymEngine::Pow;
using SymEngine::rcp_static_cast;

inline static std::vector<RCP<const Basic>>
get_add_args(const RCP<const Basic> &expr) {
	if (SymEngine::is_a<Add>(*expr)) {
		return expr->get_args();
	} else {
		return {expr};
	}
}

inline static std::vector<RCP<const Basic>>
get_mul_args(const RCP<const Basic> &expr) {
	if (SymEngine::is_a<Mul>(*expr)) {
		return expr->get_args();
	} else {
		return {expr};
	}
}

static std::shared_ptr<SeriesBase> poly_bell_substitution(const Expression &p) {
	std::shared_ptr<SeriesBase> seqTot = std::make_shared<SeriesEmpty>();

	for (const auto &polyTerm: get_add_args(expand(p.get_basic()))) {
		std::shared_ptr<SeriesBase> seqTerm = std::make_shared<SeriesFactor>(
				Expression(1));

		for (const auto &termFactor: get_mul_args(polyTerm)) {
			if (SymEngine::is_a_Number(*termFactor)) {
				seqTerm = (*seqTerm) * Expression(termFactor);
			} else if (SymEngine::is_a<Pow>(*termFactor) ||
					   SymEngine::is_a<Symbol>(*termFactor)) {
				RCP<const Basic> base;
				RCP<const Basic> exp;

				if (SymEngine::is_a<Pow>(*termFactor)) {
					auto pow_ptr = rcp_static_cast<const Pow>(termFactor);
					base = pow_ptr->get_base();
					exp = pow_ptr->get_exp();
				} else {
					base = termFactor;
					exp = integer(1);
				}

				auto sym = rcp_static_cast<const Symbol>(base);
				std::string name = sym->get_name();
				std::string ix = name.substr(name.find('_')); // "_0", "_1", ...

				// Build Bell polynomial generator.
				int j = SymEngine::down_cast<const SymEngine::Integer &>(
						*exp).as_int();

				auto bellLambdaGen = [j, ix](int n) -> Expression {
					if (n >= j) {
						return partial_ordinary_bell_polynomial(n, j, "a" + ix);
					} else {
						return {0};
					}
				};

				seqTerm = (*seqTerm) * std::make_shared<Series>(bellLambdaGen);
			} else {
				throw std::runtime_error(
						"Unhandled factor in SymEngine expression. (poly_bell_substitution)");
			}
		}

		seqTot = (*seqTot) + seqTerm;
	}

	return seqTot;
}

// Coefficient of the power of a double power series where the first series
// starts from 0 and the second starts from 1.
static std::shared_ptr<SeriesBase> double_series_power_coeff(int n, int i) {
	// Polynomial for b_{n,i} in terms of {a_0, ..., a_n}
	Expression b_ni = ordinary_potential_polynomial(n, i, "a");
	// Substitute Bell polynomials for coefficients
	return poly_bell_substitution(b_ni);
}

/** Mirrors Python a_nk_ser.
 *
 * finite series sum from max(k, n+n_offset) to n+k.
 *
 * @param n n index.
 * @param k k index.
 * @param n_offset Offset for the lower bound.
 * @param d_nkl Triple-index coefficient function.
 * @return Expression for the finite series in e2.
 */
static Expression a_nk_ser(int n, int k, int n_offset,
						   const std::function<mpq_class(int, int, int)> &d_nkl) {
	Expression a_nk(0);
	for (int l = std::max(k, n + n_offset); l <= n + k; ++l)
		a_nk = a_nk + mpq_to_expr(d_nkl(n, k, l)) * pow(e2, l);
	return a_nk;
}

/** Mirrors Python a_nk_C.
 *
 * finite series sum from 1 to k, only if k >= n+1.
 *
 * @param n n index.
 * @param k k index.
 * @param c_nkl Triple-index coefficient function.
 * @return Expression for the finite series in e2.
 */
static Expression a_nk_C(int n, int k,
						 const std::function<mpq_class(int, int, int)> &c_nkl) {
	Expression a_nk(0);
	if (k < n + 1)
		return a_nk;
	for (int l = 1; l <= k; ++l)
		a_nk = a_nk + mpq_to_expr(c_nkl(n, k, l)) * pow(e2, l);
	return a_nk;
}

/** Substitute a_{n,k} symbols in polynomial p with the result of a_nk callback.
 *
 * Mirrors Python a_nk_sub. The callback takes (n, k) and returns an Expression
 * for the substitution.
 *
 * @param p Polynomial to substitute into.
 * @param a_nk Callback returning the substitution expression for given (n, k).
 * @return Polynomial with substitutions applied.
 */
static Expression a_nk_sub(const Expression &p,
						   const std::function<Expression(int, int)> &a_nk) {
	Expression pTot{0};

	// Loop over additive terms
	for (const auto &polyTerm: get_add_args(SymEngine::expand(p.get_basic()))) {
		Expression pTerm = 1;

		// Loop over multiplicative factors
		for (const auto &termFactor: get_mul_args(polyTerm)) {
			if (SymEngine::is_a_Number(*termFactor)) {
				// Numeric constant
				pTerm = pTerm * Expression(termFactor);
			} else if (SymEngine::is_a<Pow>(*termFactor) ||
					   SymEngine::is_a<Symbol>(*termFactor)) {
				// Handle x^y or x^1
				RCP<const Basic> base;
				RCP<const Basic> exp;

				if (SymEngine::is_a<Pow>(*termFactor)) {
					auto pow_ptr = rcp_static_cast<const Pow>(termFactor);
					base = pow_ptr->get_base();
					exp = pow_ptr->get_exp();
				} else {
					base = termFactor;
					exp = integer(1);
				}

				// Extract indices from variable name: e.g. "a_3_2"
				auto sym = rcp_static_cast<const Symbol>(base);
				std::string name = sym->get_name();
				// Name format is "a_n_k", so split
				auto first_us = name.find('_');
				auto second_us = name.find('_', first_us + 1);

				int n = std::stoi(
						name.substr(first_us + 1, second_us - first_us - 1));
				int k = std::stoi(name.substr(second_us + 1));

				int exp_int = SymEngine::down_cast<const SymEngine::Integer &>(
						*exp).as_int();

				pTerm = pTerm * pow(a_nk(n, k), exp_int);
			} else {
				std::cout << p << std::endl;
				std::cout << Expression(termFactor) << std::endl;
				throw std::runtime_error(
						"Unhandled factor in SymEngine expression. (a_nk_sub)");
			}
		}

		pTot = pTot + pTerm;
	}

	return pTot;
}

static std::shared_ptr<SeriesBase> poly_bell_substitution2(const APoly &poly) {
	if (poly.empty())
		return std::make_shared<SeriesFactor>(Expression{0});

	std::shared_ptr<SeriesBase> seqTot = std::make_shared<SeriesEmpty>();

	for (const auto &[monomial, coeff] : poly) {
		std::shared_ptr<SeriesBase> seqTerm = std::make_shared<SeriesFactor>(
				mpq_to_expr(coeff));

		for (const auto &[j, e] : monomial) {
			std::string base = "a_" + std::to_string(j);
			auto bellLambdaGen = [e, base](int n) -> Expression {
				return (n >= e) ? partial_ordinary_bell_polynomial2(n, e, base)
								: Expression{0};
			};
			seqTerm = (*seqTerm) * std::make_shared<Series>(bellLambdaGen);
		}

		seqTot = (*seqTot) + seqTerm;
	}

	return seqTot;
}

static std::shared_ptr<SeriesBase> double_series_power_coeff2(int n, int i) {
	APoly b_ni = ordinary_potential_polynomial2(n, i);
	return poly_bell_substitution2(b_ni);
}

static Expression a_nk_ser2(int n, int k, int n_offset,
							const std::function<mpq_class(int, int, int)> &d_nkl) {
	Expression a_nk(0);
	for (int l = std::max(k, n + n_offset); l <= n + k; ++l)
		a_nk = a_nk + mpq_to_expr(d_nkl(n, k, l)) * pow(e2, l);
	return a_nk;
}

static Expression a_nk_C2(int n, int k,
						  const std::function<mpq_class(int, int, int)> &c_nkl) {
	Expression a_nk(0);
	if (k < n + 1)
		return a_nk;
	for (int l = 1; l <= k; ++l)
		a_nk = a_nk + mpq_to_expr(c_nkl(n, k, l)) * pow(e2, l);
	return a_nk;
}

static Expression a_nk_sub2(const Expression &p,
							const std::function<Expression(int, int)> &a_nk) {
	Expression pTot{0};

	for (const auto &polyTerm: get_add_args(SymEngine::expand(p.get_basic()))) {
		Expression pTerm = 1;

		for (const auto &termFactor: get_mul_args(polyTerm)) {
			if (SymEngine::is_a_Number(*termFactor)) {
				pTerm = pTerm * Expression(termFactor);
			} else if (SymEngine::is_a<Pow>(*termFactor) ||
					   SymEngine::is_a<Symbol>(*termFactor)) {
				RCP<const Basic> base;
				RCP<const Basic> exp;

				if (SymEngine::is_a<Pow>(*termFactor)) {
					auto pow_ptr = rcp_static_cast<const Pow>(termFactor);
					base = pow_ptr->get_base();
					exp = pow_ptr->get_exp();
				} else {
					base = termFactor;
					exp = integer(1);
				}

				auto sym = rcp_static_cast<const Symbol>(base);
				std::string name = sym->get_name();
				auto first_us = name.find('_');
				auto second_us = name.find('_', first_us + 1);

				int n = std::stoi(
						name.substr(first_us + 1, second_us - first_us - 1));
				int k = std::stoi(name.substr(second_us + 1));

				int exp_int = SymEngine::down_cast<const SymEngine::Integer &>(
						*exp).as_int();

				pTerm = pTerm * pow(a_nk(n, k), exp_int);
			} else {
				std::cout << p << std::endl;
				std::cout << Expression(termFactor) << std::endl;
				throw std::runtime_error(
						"Unhandled factor in SymEngine expression. (a_nk_sub2)");
			}
		}

		pTot = pTot + pTerm;
	}

	return pTot;
}

// ---------------------------------------------------------------------------
// LExpr pipeline (GMP-only, no SymEngine).
// ---------------------------------------------------------------------------

// Evaluate an LExpr tree as a dense E2Poly by substituting each a_{j,m}
// via var_fn. DAG memoization ensures each shared sub-node is computed once.
static E2Poly lexpr_eval_e2poly(
		const LExpr& e,
		const std::function<E2Poly(int, int)>& var_fn,
		std::unordered_map<const LNode*, E2Poly>& memo) {
	if (!e) return {};
	auto it = memo.find(e.get());
	if (it != memo.end()) return it->second;

	E2Poly result = std::visit([&](const auto& node) -> E2Poly {
		using T = std::decay_t<decltype(node)>;
		if constexpr (std::is_same_v<T, LVarNode>)
			return var_fn(node.j, node.m);
		else if constexpr (std::is_same_v<T, LConstNode>)
			return E2Poly{node.c};
		else if constexpr (std::is_same_v<T, LAddNode>)
			return e2poly_add(lexpr_eval_e2poly(node.a, var_fn, memo),
							  lexpr_eval_e2poly(node.b, var_fn, memo));
		else  // LMulNode
			return e2poly_mul(lexpr_eval_e2poly(node.a, var_fn, memo),
							  lexpr_eval_e2poly(node.b, var_fn, memo));
	}, e->data);

	return memo[e.get()] = result;
}

static E2Poly lexpr_eval_e2poly(
		const LExpr& e,
		const std::function<E2Poly(int, int)>& var_fn) {
	std::unordered_map<const LNode*, E2Poly> memo;
	return lexpr_eval_e2poly(e, var_fn, memo);
}

// GMP-only counterpart of a_nk_ser2: sum_{l} d_nkl(n,k,l) * e2^l.
static E2Poly a_nk_ser_lexpr(int n, int k, int n_offset,
							 const std::function<mpq_class(int, int, int)>& d_nkl) {
	int l_lo = std::max(k, n + n_offset);
	int l_hi = n + k;
	if (l_lo > l_hi || k < 1) return {};
	E2Poly result(l_hi + 1, mpq_class(0));
	for (int l = l_lo; l <= l_hi; ++l)
		result[l] = d_nkl(n, k, l);
	return result;
}

// GMP-only counterpart of a_nk_C2: sum_{l=1}^{k} c_nkl(n,k,l) * e2^l.
static E2Poly a_nk_C_lexpr(int n, int k,
						   const std::function<mpq_class(int, int, int)>& c_nkl) {
	if (k < n + 1) return {};
	E2Poly result(k + 1, mpq_class(0));
	for (int l = 1; l <= k; ++l)
		result[l] = c_nkl(n, k, l);
	return result;
}

// Substitute a_j^e in each APoly term with a TSeries<LExpr> of Bell polynomials.
static std::shared_ptr<TSeriesBase<LExpr>>
poly_bell_substitution_lexpr(const APoly& poly) {
	if (poly.empty())
		return std::make_shared<TSeriesFactor<LExpr>>(LExpr{});  // zero series
	std::shared_ptr<TSeriesBase<LExpr>> seqTot = std::make_shared<TSeriesEmpty<LExpr>>();
	for (const auto& [monomial, coeff] : poly) {
		std::shared_ptr<TSeriesBase<LExpr>> seqTerm =
				std::make_shared<TSeriesFactor<LExpr>>(lexpr_const(coeff));
		for (const auto& [j, e] : monomial) {
			auto gen = [j=j, e=e](int k) -> LExpr {
				return (k >= e) ? partial_bell_lexpr(k, e, j) : LExpr{};
			};
			seqTerm = (*seqTerm) * std::make_shared<TSeries<LExpr>>(gen);
		}
		seqTot = (*seqTot) + seqTerm;
	}
	return seqTot;
}

static std::shared_ptr<TSeriesBase<LExpr>>
double_series_power_coeff_lexpr(int n, int i) {
	return poly_bell_substitution_lexpr(ordinary_potential_polynomial2(n, i));
}

}  // namespace point_to_ellipse_series
