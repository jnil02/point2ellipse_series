cmake_minimum_required(VERSION 3.15)

project(point_to_ellipse_Fourier_series_expansion CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

# Default to an optimized build. Pass -DCMAKE_BUILD_TYPE=Debug (or
# RelWithDebInfo) to override. Note that Release and RelWithDebInfo define
# NDEBUG, which disables the coefficient index-range asserts.
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type." FORCE)
endif()

# Repository-root directory for the Python-generated CSV fixtures/outputs.
set(PTE_TEST_DATA_DIR "${PROJECT_SOURCE_DIR}/../test_data")

# ---------------------------------------------------------------------------
# Third-party dependencies (located at configure time so that missing
# libraries/headers fail early with a clear message).
# ---------------------------------------------------------------------------

find_path(GMPXX_INCLUDE_DIR gmpxx.h)
find_library(GMP_LIBRARY gmp)
find_library(GMPXX_LIBRARY gmpxx)
find_library(MPFR_LIBRARY mpfr)
find_path(MPREAL_INCLUDE_DIR mpreal.h)
find_path(SYMENGINE_INCLUDE_DIR symengine/expression.h)
find_library(SYMENGINE_LIBRARY symengine)

foreach(dep GMPXX_INCLUDE_DIR GMP_LIBRARY GMPXX_LIBRARY MPFR_LIBRARY
            MPREAL_INCLUDE_DIR SYMENGINE_INCLUDE_DIR SYMENGINE_LIBRARY)
    if(NOT ${dep})
        message(FATAL_ERROR "Dependency not found: ${dep}. "
                "See the Dependencies section of README.md.")
    endif()
endforeach()

# GMP with its C++ interface (mpq_class).
add_library(PTE::GMP INTERFACE IMPORTED)
target_include_directories(PTE::GMP INTERFACE ${GMPXX_INCLUDE_DIR})
target_link_libraries(PTE::GMP INTERFACE ${GMPXX_LIBRARY} ${GMP_LIBRARY})

# MPFR with the header-only MPFR C++ wrapper (mpfr::mpreal).
add_library(PTE::MPFR INTERFACE IMPORTED)
target_include_directories(PTE::MPFR INTERFACE ${MPREAL_INCLUDE_DIR})
target_link_libraries(PTE::MPFR INTERFACE ${MPFR_LIBRARY} PTE::GMP)

# SymEngine (built with MPFR support).
add_library(PTE::SymEngine INTERFACE IMPORTED)
target_include_directories(PTE::SymEngine INTERFACE ${SYMENGINE_INCLUDE_DIR})
target_link_libraries(PTE::SymEngine INTERFACE ${SYMENGINE_LIBRARY} PTE::MPFR)

# ---------------------------------------------------------------------------
# Coefficient library (independent of the ellipse parameters).
# ---------------------------------------------------------------------------

set(PTECOEF_SOURCE
        coefficients.cpp
        coefficients_evo.cpp
        )

add_library(ptecoef_shared SHARED ${PTECOEF_SOURCE})
add_library(ptecoef_static STATIC ${PTECOEF_SOURCE})

foreach(lib ptecoef_shared ptecoef_static)
    # Both build as libptecoef (.so / .a).
    set_target_properties(${lib} PROPERTIES OUTPUT_NAME ptecoef)
    # PUBLIC: the coefficient headers expose mpq_class, and in-tree consumers
    # get the local headers ahead of installed (potentially old) ones.
    target_include_directories(${lib}
            PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
            )
    target_link_libraries(${lib} PUBLIC PTE::GMP)
endforeach()

include(GNUInstallDirs) # Defines CMAKE_INSTALL_LIBDIR.

# Install coefficient library. (run ldconfig subsequently.)
install(TARGETS ptecoef_shared ptecoef_static
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}  # .so
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}  # .a
        )

install(FILES
        coefficients.hpp
        coefficients_evo.hpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# ---------------------------------------------------------------------------
# Executables. All link the in-tree static coefficient library.
#
# Ellipse parameters are compile-time macros (see ellipse/ellipse_params.hpp):
# define ELLIPSE_A plus exactly one of ELLIPSE_B / ELLIPSE_INV_F per target.
# ---------------------------------------------------------------------------

# Demo.
add_executable(main demo/main.cpp)
target_link_libraries(main PRIVATE ptecoef_static PTE::SymEngine)
target_compile_definitions(main
        PRIVATE
        ELLIPSE_A="6378137.0"
        ELLIPSE_INV_F="298.257223563"
        )

# Convergence sweep binary for the inside-evolute series.
# Default ellipse: unit ellipse a=1, b=0.5 (large evolute; convenient for
# convergence/ROC studies).
add_executable(sweep_evo convergence/sweep_evo.cpp)
target_include_directories(sweep_evo PRIVATE convergence)
target_link_libraries(sweep_evo PRIVATE ptecoef_static PTE::MPFR)
target_compile_definitions(sweep_evo
        PRIVATE
        ELLIPSE_A="1.0"
        ELLIPSE_B="0.5"
        TEST_DATA_DIR="${PTE_TEST_DATA_DIR}"
        )

# Coefficient ROC diagnostic for the inside-evolute series.
add_executable(diag_coeff_evo convergence/diag_coeff_evo.cpp)
target_include_directories(diag_coeff_evo PRIVATE convergence)
target_link_libraries(diag_coeff_evo PRIVATE ptecoef_static PTE::MPFR)
target_compile_definitions(diag_coeff_evo
        PRIVATE
        ELLIPSE_A="1.0"
        ELLIPSE_B="0.25"
        TEST_DATA_DIR="${PTE_TEST_DATA_DIR}"
        )

# Coefficient micro-benchmarks.
add_executable(bench_coefficients benchmarks/bench_coefficients.cpp)
target_link_libraries(bench_coefficients PRIVATE ptecoef_static PTE::SymEngine)
target_compile_definitions(bench_coefficients
        PRIVATE
        ELLIPSE_A="1.0"
        ELLIPSE_B="0.25"
        )

# Inside-evolute series evaluation benchmarks (mpreal numeric path).
add_executable(bench_series_evo benchmarks/bench_series_evo.cpp)
target_link_libraries(bench_series_evo PRIVATE ptecoef_static PTE::SymEngine)
target_compile_definitions(bench_series_evo
        PRIVATE
        ELLIPSE_A="1.0"
        ELLIPSE_B="0.25"
        )

# ---------------------------------------------------------------------------
# Unit tests (Catch2 v2 single header <catch.hpp>). Disable with
# -DBUILD_TESTING=OFF; skipped with a warning if the header is not found.
# ---------------------------------------------------------------------------

include(CTest)

# Regenerate the tracked CSV fixtures from the Python reference (needs Python 3
# with SymPy and mpmath). Not part of the default build:
#   cmake --build <build-dir> --target generate_test_data
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    set(PTE_PYTHON_DIR "${PROJECT_SOURCE_DIR}/../python")
    add_custom_target(generate_test_data
            COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${PTE_PYTHON_DIR}"
                    ${Python3_EXECUTABLE} ${PTE_PYTHON_DIR}/tests/generate_test_data.py
            COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${PTE_PYTHON_DIR}"
                    ${Python3_EXECUTABLE} ${PTE_PYTHON_DIR}/tests/generate_series_test_data.py
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/.."
            COMMENT "Regenerating the CSV fixtures in test_data/ from the Python reference"
            VERBATIM
            )
endif()

if(BUILD_TESTING)
    find_path(CATCH2_INCLUDE_DIR catch.hpp)
    if(CATCH2_INCLUDE_DIR)
        add_subdirectory(tests)
    else()
        message(WARNING "Catch2 v2 header <catch.hpp> not found; skipping "
                "the unit tests. Set CATCH2_INCLUDE_DIR to enable them.")
    endif()
endif()
