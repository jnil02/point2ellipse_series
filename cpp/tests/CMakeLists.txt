add_executable(test_coefficients
        test_coefficients.cpp
        test_stirling.cpp
        ../coefficients.cpp
        ../coefficients_evo.cpp
        )

target_include_directories(test_coefficients
        PRIVATE
        # Points to root/cpp/. This is to give local header files priority over
        # installed (potentially old) header files.
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        )

target_link_libraries(test_coefficients
        PRIVATE
        gmpxx
        gmp
        )

target_compile_definitions(test_coefficients
        PRIVATE
        TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/../test_data"
        )

include(CTest)
add_test(NAME test_coefficients COMMAND test_coefficients)

add_executable(test_series
        test_series.cpp
        ../coefficients.cpp
        ../coefficients_evo.cpp
        )

target_include_directories(test_series
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        )

target_link_libraries(test_series
        PRIVATE
        symengine
        mpfr
        gmp
        )

target_compile_definitions(test_series
        PRIVATE
        ELLIPSE_A="6378137.0"
        ELLIPSE_INV_F="298.257223563"
        )

add_test(NAME test_series COMMAND test_series)

add_executable(test_series_evo
        test_series_evo.cpp
        ../coefficients.cpp
        ../coefficients_evo.cpp
        )

target_include_directories(test_series_evo
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        )

target_link_libraries(test_series_evo
        PRIVATE
        symengine
        mpfr
        gmp
        )

target_compile_definitions(test_series_evo
        PRIVATE
        ELLIPSE_A="6378137.0"
        ELLIPSE_INV_F="298.257223563"
        )

add_test(NAME test_series_evo COMMAND test_series_evo)

add_executable(test_lexpr
        test_lexpr.cpp
        )

target_include_directories(test_lexpr
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        )

target_link_libraries(test_lexpr
        PRIVATE
        gmpxx
        gmp
        )

add_test(NAME test_lexpr COMMAND test_lexpr)