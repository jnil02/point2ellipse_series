add_executable(test_coefficients
        test_coefficients.cpp
        test_stirling.cpp
        ../coefficients.cpp
        ../coefficients_evo.cpp
        ../series.cpp
        )

target_include_directories(test_coefficients
        PRIVATE
        # Points to root/cpp/. This is to give local header files priority over
        # installed (potentially old) header files.
        ${CMAKE_CURRENT_SOURCE_DIR}/..
        )

target_link_libraries(test_coefficients
        PRIVATE
        symengine
        mpfr
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
        ../series.cpp
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

add_test(NAME test_series COMMAND test_series)

add_executable(test_series_evo
        test_series_evo.cpp
        ../coefficients.cpp
        ../coefficients_evo.cpp
        ../series.cpp
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

add_test(NAME test_series_evo COMMAND test_series_evo)