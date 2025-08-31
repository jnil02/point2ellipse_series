#pragma once

#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/symbol.h>
#include <symengine/number.h>

#include "series.hpp"
#include "polynomials.hpp"
#include "symbols.hpp"

using SymEngine::Expression;
using SymEngine::RCP;
using SymEngine::Basic;
using SymEngine::Number;
using SymEngine::Symbol;
using SymEngine::Add;
using SymEngine::Mul;
using SymEngine::Pow;
using SymEngine::rcp_static_cast;

std::vector<RCP<const Basic>> get_add_args(const RCP<const Basic> &expr) {
	if (SymEngine::is_a<Add>(*expr)) {
		return expr->get_args();
	} else {
		return {expr};
	}
}

std::vector<RCP<const Basic>> get_mul_args(const RCP<const Basic> &expr) {
	if (SymEngine::is_a<Mul>(*expr)) {
		return expr->get_args();
	} else {
		return {expr};
	}
}

std::shared_ptr<SeriesBase> poly_bell_substitution(const Expression &p) {
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
				std::string ix = name.substr(name.find('_')); // "_0", "_1",

				// Build Bell polynomial generator
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

// Coefficient of the power of a double power series where the first series start from 0 and the second starts from 1.
std::shared_ptr<SeriesBase> double_series_power_coeff(int n, int i) {
	// Polynomial for b_{n,i} in terms of {a_0, ..., a_n}
	Expression b_ni = ordinary_potential_polynomial(n, i, "a");
	// Substitute Bell polynomials for coefficients
	return poly_bell_substitution(b_ni);
}

// C++ version of a_nk_sub
Expression a_nk_sub(const Expression &p,
					int n_offset,
					const std::function<Expression(int, int, int)> &d_nkl) {
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

				// Construct substitution sum
				Expression a_nk = 0;
				int exp_int = SymEngine::down_cast<const SymEngine::Integer &>(
						*exp).as_int();

				for (int l = std::max(k, n + n_offset); l <= n + k; ++l) {
					a_nk = a_nk + d_nkl(n, k, l) * pow(e2, l);
				}

				pTerm = pTerm * pow(a_nk, exp_int);
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
