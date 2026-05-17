#pragma once

/*
 * Reference ellipse parameters.
 *
 * Always define ELLIPSE_A (semi-major axis).
 * Define the shape via exactly one of:
 *   ELLIPSE_B     — semi-minor axis
 *   ELLIPSE_INV_F — inverse flattening a/(a-b)
 *
 * Macros must be consistent across all translation units in a binary.
 * The recommended way to override is via CMake target_compile_definitions,
 * which applies to the whole target at once:
 *
 *   # Unit ellipse with b/a = 0.5 (large evolute, good for convergence studies):
 *   target_compile_definitions(my_target PRIVATE
 *       ELLIPSE_A="1.0"
 *       ELLIPSE_B="0.5")
 *
 *   # WGS84 via inverse flattening:
 *   target_compile_definitions(my_target PRIVATE
 *       ELLIPSE_A="6378137.0"
 *       ELLIPSE_INV_F="298.257223563")
 *
 * Alternatively, define before including this header in a single-TU program:
 *
 *   #define ELLIPSE_A "1.0"
 *   #define ELLIPSE_B "0.5"
 *   #include "ellipse.hpp"
 *
 * Default: WGS84 (a and 1/f).
 */

#ifndef ELLIPSE_A
#define ELLIPSE_A "6378137.0"
#endif

// Default shape: WGS84 1/f (used when neither ELLIPSE_B nor ELLIPSE_INV_F is
// defined). To use semi-minor axis instead, define ELLIPSE_B.
#if !defined(ELLIPSE_B) && !defined(ELLIPSE_INV_F)
#define ELLIPSE_INV_F "298.257223563"
#endif
