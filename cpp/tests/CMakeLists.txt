add_executable(test_coefficients
        test_coefficients.cpp
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