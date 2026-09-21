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

// GMP-only counterpart of a_nk_C: sum_{n} c(k,l,n) * e2^n.
//
// The generator a_{l,k} carries the parity of the underlying c coefficients: it
// is zero unless k>=l+1 and k==l+1 (mod 2), and within a nonzero generator only
// the e2-powers n==k (mod 2) contribute. Off-parity generators/terms are
// identically zero, so pruning them here only avoids provably-zero evaluations.
static E2Poly a_nk_C_lexpr(int l, int k,
						   const std::function<mpq_class(int, int, int)>& c) {
	if (k < l + 1 || (k - l - 1) % 2 != 0) return {};
	E2Poly result(k + 1, mpq_class(0));
	for (int n = 2 - k % 2; n <= k; n += 2)  // e2-power n == k (mod 2); others vanish
		result[n] = c(k, l, n);
	return result;
}

// Substitute a_j^e in each APoly term with a TSeries<LExpr> of Bell polynomials.
//
// start(j) is the lowest power z^{start(j)} at which generator a_j begins. The
// default assumes every generator starts at z^1, giving the standard partial
// ordinary Bell nonzero condition k>=e. For generators that start at z^{j+1}
// (the inside-evolute a_l), pass start=[](int l){return l+1;}, which tightens
// the guard to k>=e*(j+1). Since a_{j,m}=0 for m<start(j), the pruned terms are
// identically zero, so this only removes provably-zero terms.
static std::shared_ptr<TSeriesBase<LExpr>>
poly_bell_substitution_lexpr(const APoly& poly,
							 const std::function<int(int)>& start = [](int){ return 1; }) {
	if (poly.empty())
		return std::make_shared<TSeriesFactor<LExpr>>(LExpr{});  // zero series
	std::shared_ptr<TSeriesBase<LExpr>> seqTot = std::make_shared<TSeriesEmpty<LExpr>>();
	for (const auto& [monomial, coeff] : poly) {
		std::shared_ptr<TSeriesBase<LExpr>> seqTerm =
				std::make_shared<TSeriesFactor<LExpr>>(lexpr_const(coeff));
		for (const auto& [j, e] : monomial) {
			// a_j^e starts at z^{e*start(j)}; below that partial_bell_lexpr vanishes.
			int thr = e * start(j);
			auto gen = [j=j, e=e, thr](int k) -> LExpr {
				return (k >= thr) ? partial_bell_lexpr(k, e, j) : LExpr{};
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

// Inside-evolute variant: the inner generator a_l starts at z^{l+1}, so the Bell
// guard tightens to k>=i*(l+1). Values are identical to the untightened version
// (pruned terms vanish after the a_{n,k} substitution); it only avoids building
// provably-zero terms.
static std::shared_ptr<TSeriesBase<LExpr>>
double_series_power_coeff_evo_lexpr(int l, int i) {
	return poly_bell_substitution_lexpr(ordinary_potential_polynomial2(l, i),
										[](int l){ return l + 1; });
}

}  // namespace point_to_ellipse_series
