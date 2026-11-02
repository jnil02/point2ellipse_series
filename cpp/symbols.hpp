#pragma once

#include <symengine/symbol.h>

inline const SymEngine::Expression varrho(SymEngine::symbol("varrho"));
inline const SymEngine::Expression psi(SymEngine::symbol("psi"));
inline const SymEngine::Expression sin_psi(SymEngine::symbol("sin_psi"));
inline const SymEngine::Expression cos_psi(SymEngine::symbol("cos_psi"));
inline const SymEngine::RCP<const SymEngine::Symbol> e2sym = SymEngine::symbol("e2");
inline const SymEngine::Expression e2(e2sym);  // e² variable
inline const SymEngine::Expression sigma(SymEngine::symbol("sigma"));  // rho/(a*e²)
inline const SymEngine::Expression epsilon(SymEngine::symbol("epsilon"));           // b/a = sqrt(1-e²)