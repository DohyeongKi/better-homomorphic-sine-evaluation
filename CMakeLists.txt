cmake_minimum_required(VERSION 3.10)

include_directories( /usr/local/include ) # path to NTL and gmp

link_directories( /usr/local/lib ) # path to NTL and amp

set( CMAKE_CXX_STANDARD 11 )

add_compile_options( -pthread -O2 )

ADD_EXECUTABLE( find_polynomial
                find_polynomial.cpp )

TARGET_LINK_LIBRARIES( find_polynomial ntl gmp pthread)