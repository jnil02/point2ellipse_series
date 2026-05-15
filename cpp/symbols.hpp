#pragma once

#include <symengine/symbol.h>

inline const SymEngine::Expression varrho(SymEngine::symbol("varrho"));
inline const SymEngine::Expression psi(SymEngine::symbol("psi"));
inline const SymEngine::Expression sin_psi(SymEngine::symbol("sin_psi"));
inline const SymEngine::Expression cos_psi(SymEngine::symbol("cos_psi"));
inline const SymEngine::RCP<const SymEngine::Symbol> e2sym = SymEngine::symbol("e2");
inline const SymEngine::Expression e2(e2sym);  // e² variable
inline const SymEngine::Expression rho_ae2(SymEngine::symbol("rho_ae2"));  // rho/(a*e²)
inline const SymEngine::Expression b_a(SymEngine::symbol("b_a"));           // b/a = sqrt(1-e²)