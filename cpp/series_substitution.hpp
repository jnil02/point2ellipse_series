#pragma once

#include "series.hpp"
#include "polynomials.hpp"

namespace point_to_ellipse_series {

// ---------------------------------------------------------------------------
// LExpr pipeline (GMP-only).
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

// GMP-only counterpart of a_nk_ser: sum_{l} d_nkl(n,k,l) * e2^l.
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

// GMP-only counterpart of a_nk_C: sum_{l=1}^{k} c_nkl(n,k,l) * e2^l.
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
