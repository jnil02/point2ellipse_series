#pragma once

#include "series_traits.hpp"
#include <cmath>

namespace point_to_ellipse_series {

template<> inline double to<double>(const mpq_class& c) {
	return c.get_d();
}

template<> inline double ipow<double>(const double& b, int e) {
	return std::pow(b, e);
}


}  // namespace point_to_ellipse_series