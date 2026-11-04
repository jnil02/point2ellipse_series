# Catch2 v2 single header.
add_library(PTE::Catch2 INTERFACE IMPORTED)
target_include_directories(PTE::Catch2 INTERFACE ${CATCH2_INCLUDE_DIR})

# Add a Catch2 test executable linked against the in-tree coefficient library
# (whose PUBLIC include directory gives the local headers priority over
# installed, potentially old, ones) and register it with CTest.
#
#   pte_add_test(<name> SOURCES <src>... [LIBS <lib>...] [DEFS <def>...])
function(pte_add_test name)
    cmake_parse_arguments(ARG "" "" "SOURCES;LIBS;DEFS" ${ARGN})
    add_executable(${name} ${ARG_SOURCES})
    target_link_libraries(${name} PRIVATE ptecoef_static PTE::Catch2 ${ARG_LIBS})
    target_compile_definitions(${name} PRIVATE ${ARG_DEFS})
    add_test(NAME ${name} COMMAND ${name})
endfunction()

# Coefficients cross-checked against the Python-generated CSV fixtures.
pte_add_test(test_coefficients
        SOURCES test_coefficients.cpp test_stirling.cpp
        DEFS TEST_DATA_DIR="${PTE_TEST_DATA_DIR}"
        )

# Series values cross-checked against the Python-generated CSV fixtures.
pte_add_test(test_series_values
        SOURCES test_series_values.cpp
        DEFS TEST_DATA_DIR="${PTE_TEST_DATA_DIR}"
        )

# Far-field series against the closed-form reference (WGS84).
pte_add_test(test_series
        SOURCES test_series.cpp
        LIBS PTE::SymEngine
        DEFS ELLIPSE_A="6378137.0" ELLIPSE_INV_F="298.257223563"
        )

# Inside-evolute series against the closed-form reference (WGS84).
pte_add_test(test_series_evo
        SOURCES test_series_evo.cpp
        LIBS PTE::SymEngine
        DEFS ELLIPSE_A="6378137.0" ELLIPSE_INV_F="298.257223563"
        )

# Lazy-expression machinery.
pte_add_test(test_lexpr
        SOURCES test_lexpr.cpp
        )
