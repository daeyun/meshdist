cmake_minimum_required(VERSION 3.7)
project(meshdist)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)

FIND_PACKAGE(Boost 1.5 REQUIRED)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -march=native -funroll-loops -mfpmath=sse -ftree-vectorize")

set(SOURCE_FILES main.cc meshdist.h)

include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})
include_directories(third_party)
add_executable(meshdist ${SOURCE_FILES})

add_executable(tests test_meshdist.cc)
