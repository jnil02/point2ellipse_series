cmake_minimum_required(VERSION 3.15)  #

project(point_to_ellipse_Fourier_series_expansion)

set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS -O2)

add_executable(main
        main.cpp
        coefficients.cpp
        )

target_link_libraries(main
        symengine
        mpfr
        gmp
        )

target_compile_definitions(main
        PRIVATE
        ELLIPSE_A="6378137.0"
        ELLIPSE_INV_F="298.257223563"
        )

# Build coefficient library.
set(PTECOEF_SOURCE
        coefficients.cpp
        coefficients_evo.cpp
        )

add_library(ptecoef_shared SHARED ${PTECOEF_SOURCE})
add_library(ptecoef_static STATIC ${PTECOEF_SOURCE})

set_target_properties(ptecoef_shared PROPERTIES OUTPUT_NAME ptecoef)
set_target_properties(ptecoef_static PROPERTIES OUTPUT_NAME ptecoef)

target_link_libraries(ptecoef_shared PRIVATE gmpxx gmp)
target_link_libraries(ptecoef_static PRIVATE gmpxx gmp)

include(GNUInstallDirs)

# Install coefficient library. (run ldconfig subsequently.)
install(TARGETS ptecoef_shared ptecoef_static
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}  # .so
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}  # .a
        )

install(FILES
        coefficients.hpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Unit tests.
add_subdirectory(tests)

# Convergence sweep binary for the inside-evolute series.
# Default ellipse: unit ellipse a=1, b=0.25 (ae bound becomes binding at ~62°, good for ROC studies).
# Override by redefining ELLIPSE_A + ELLIPSE_B (or ELLIPSE_INV_F) below.
add_executable(sweep_evo
        sweep_evo.cpp
        coefficients.cpp
        coefficients_evo.cpp
        )

target_include_directories(sweep_evo
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        )

target_link_libraries(sweep_evo
        PRIVATE
        mpfr
        gmp
        )

target_compile_definitions(sweep_evo
        PRIVATE
        ELLIPSE_A="1.0"
        ELLIPSE_B="0.25"
        TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/../test_data"
        )

# Coefficient micro-benchmarks.
add_executable(bench_coefficients
        benchmarks/bench_coefficients.cpp
        coefficients.cpp
        coefficients_evo.cpp
        )

target_include_directories(bench_coefficients
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        )

target_link_libraries(bench_coefficients
        PRIVATE
        symengine
        mpfr
        gmp
        )

target_compile_definitions(bench_coefficients
        PRIVATE
        ELLIPSE_A="1.0"
        ELLIPSE_B="0.25"
        )

# Inside-evolute series evaluation benchmarks (mpreal numeric path).
add_executable(bench_series_evo
        benchmarks/bench_series_evo.cpp
        coefficients.cpp
        coefficients_evo.cpp
        )

target_include_directories(bench_series_evo
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        )

target_link_libraries(bench_series_evo
        PRIVATE
        symengine
        mpfr
        gmp
        )

target_compile_definitions(bench_series_evo
        PRIVATE
        ELLIPSE_A="1.0"
        ELLIPSE_B="0.25"
        )