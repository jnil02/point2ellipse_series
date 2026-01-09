cmake_minimum_required(VERSION 3.15)  #

project(point_to_ellipse_Fourier_series_expansion)

set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS -O2)

add_executable(main
        main.cpp
        series.cpp
        coefficients.cpp
        )

target_link_libraries(main
        symengine
        mpfr
        gmp
        )

# Build coefficient library.
set(PTECOEF_SOURCE
        coefficients.cpp
        series.cpp
        )

add_library(ptecoef_shared SHARED ${PTECOEF_SOURCE})
add_library(ptecoef_static STATIC ${PTECOEF_SOURCE})

set_target_properties(ptecoef_shared PROPERTIES OUTPUT_NAME ptecoef)
set_target_properties(ptecoef_static PROPERTIES OUTPUT_NAME ptecoef)

include(GNUInstallDirs)

# Install coefficient library. (run ldconfig subsequently.)
install(TARGETS ptecoef_shared ptecoef_static
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}  # .so
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}  # .a
        )

install(FILES
        coefficients.hpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})