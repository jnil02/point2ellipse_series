#pragma once

#include <gmpxx.h>
#include <variant>
#include <memory>
#include <unordered_map>
#include <functional>

// Forward-declare LNode so LAddNode/LMulNode can hold shared_ptr<const LNode>.
struct LNode;

struct LVarNode   { int j, m; };    // symbolic variable a_{j,m}
struct LConstNode { mpq_class c; }; // rational constant
struct LAddNode   { std::shared_ptr<const LNode> a, b; }; // lazy a + b
struct LMulNode   { std::shared_ptr<const LNode> a, b; }; // lazy a * b

// All four alternative types are complete here, so the variant is well-formed.
struct LNode {
	std::variant<LVarNode, LConstNode, LAddNode, LMulNode> data;
};

// LExpr = shared_ptr to an immutable DAG node.  nullptr represents zero.
using LExpr = std::shared_ptr<const LNode>;

inline LExpr lexpr_var(int j, int m) {
	return std::make_shared<LNode>(LVarNode{j, m});
}
inline LExpr lexpr_const(mpq_class c) {
	return std::make_shared<LNode>(LConstNode{std::move(c)});
}

// Arithmetic — all O(1).  nullptr short-circuits as the additive/absorbing zero.

inline LExpr operator+(LExpr a, const LExpr& b) {
	if (!a) return b;
	if (!b) return a;
	return std::make_shared<LNode>(LAddNode{std::move(a), b});
}

inline LExpr& operator+=(LExpr& a, const LExpr& b) {
	return a = std::move(a) + b;
}

inline LExpr operator*(const LExpr& a, const LExpr& b) {
	if (!a || !b) return nullptr;
	return std::make_shared<LNode>(LMulNode{a, b});
}

inline LExpr operator*(LExpr a, const mpq_class& s) {
	if (!a || s == 0) return nullptr;
	if (s == 1) return a;
	return std::make_shared<LNode>(LMulNode{std::move(a), lexpr_const(s)});
}

inline LExpr operator*(const mpq_class& s, LExpr a) {
	return std::move(a) * s;
}

// Evaluation: substitute each a_{j,m} via var_fn and compute the result.
// Uses a DAG memo-table so shared sub-expressions are evaluated exactly once.

inline mpq_class lexpr_eval(
		const LExpr& e,
		const std::function<mpq_class(int, int)>& var_fn,
		std::unordered_map<const LNode*, mpq_class>& memo) {
	if (!e) return mpq_class(0);

	auto it = memo.find(e.get());
	if (it != memo.end()) return it->second;

	mpq_class result = std::visit([&](const auto& node) -> mpq_class {
		using T = std::decay_t<decltype(node)>;
		if constexpr (std::is_same_v<T, LVarNode>) {
			return var_fn(node.j, node.m);
		} else if constexpr (std::is_same_v<T, LConstNode>) {
			return node.c;
		} else if constexpr (std::is_same_v<T, LAddNode>) {
			return lexpr_eval(node.a, var_fn, memo)
				   + lexpr_eval(node.b, var_fn, memo);
		} else { // LMulNode
			mpq_class lv = lexpr_eval(node.a, var_fn, memo);
			if (lv == 0) return mpq_class(0);
			return lv * lexpr_eval(node.b, var_fn, memo);
		}
	}, e->data);

	result.canonicalize();
	return memo[e.get()] = result;
}

inline mpq_class lexpr_eval(
		const LExpr& e,
		const std::function<mpq_class(int, int)>& var_fn) {
	std::unordered_map<const LNode*, mpq_class> memo;
	return lexpr_eval(e, var_fn, memo);
}
