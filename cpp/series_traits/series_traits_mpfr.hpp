#pragma once

#include "series_traits.hpp"
#include <mpreal.h>

namespace point_to_ellipse_series {

template<>
inline mpfr::mpreal to<mpfr::mpreal>(const mpq_class& c) {
	mpfr::mpreal num, den;
	mpfr_set_z(num.mpfr_ptr(), c.get_num().get_mpz_t(), mpfr::mpreal::get_default_rnd());
	mpfr_set_z(den.mpfr_ptr(), c.get_den().get_mpz_t(), mpfr::mpreal::get_default_rnd());
	return num / den;
}

template<>
inline mpfr::mpreal ipow<mpfr::mpreal>(const mpfr::mpreal& base, int exp) {
	return mpfr::pow(base, static_cast<long>(exp));
}

template<>
inline mpfr::mpreal isin_mul<mpfr::mpreal>(const mpfr::mpreal& psi_v, int n) {
	return mpfr::sin(mpfr::mpreal(2 * n) * psi_v);
}

template<>
inline mpfr::mpreal icos_mul<mpfr::mpreal>(const mpfr::mpreal& psi_v, int n) {
	return mpfr::cos(mpfr::mpreal(2 * n) * psi_v);
}

}  // namespace point_to_ellipse_series