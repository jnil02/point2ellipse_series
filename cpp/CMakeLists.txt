cmake_minimum_required(VERSION 3.15)  #

project(point_to_ellipse_Fourier_series_expansion)

set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS -O2)

add_executable(main
        main.cpp
        series.cpp
        )

target_link_libraries(main
        symengine
        mpfr
        gmp
        )
