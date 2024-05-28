cmake_minimum_required(VERSION 3.5)
project(scurvy)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-parentheses -fdiagnostics-color=always -O2")

find_package(Eigen3 REQUIRED)
message(EIGEN3_INCLUDE_DIRS="${EIGEN3_INCLUDE_DIRS}")

add_executable(stress_test
        src/stress_test.cpp
        src/solve.cpp
        src/include/basics.h
        src/include/cv.h
        src/include/maths.h
        src/include/ncv.h
        src/include/ncv_ca.h
        src/include/ncv_nca.h
        src/include/scurvy.h
        src/include/solve.h)

target_include_directories(stress_test PUBLIC src/include "${EIGEN3_INCLUDE_DIRS}")
