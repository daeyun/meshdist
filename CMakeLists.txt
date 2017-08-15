cmake_minimum_required(VERSION 3.4)

project(meshdist)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(assimp REQUIRED)
find_package(CGAL REQUIRED)

SET(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.5 REQUIRED COMPONENTS log REQUIRED)

enable_testing()
find_package(GTest REQUIRED)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -march=native -funroll-loops -mfpmath=both")


set(SOURCE_FILES meshdist.h meshdist.cc mesh_io.h mesh_io.cc meshdist_cgal.cc meshdist_cgal.h meshdist_trimesh.h)

add_subdirectory(third_party)

include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${TRIMESH2_INCLUDE_DIR})
include_directories(${ASSIMP_INCLUDE_DIR})
include_directories(${CGAL_INCLUDE_DIR})
include_directories(third_party)


add_library(meshdist STATIC ${SOURCE_FILES})
target_link_libraries(meshdist ${TRIMESH2_LIBRARY})
target_link_libraries(meshdist ${ASSIMP_LIBRARIES})
target_link_libraries(meshdist ${Boost_LIBRARIES})
target_link_libraries(meshdist ${CGAL_LIBRARIES})


add_executable(main main.cc)
target_link_libraries(main meshdist)


include_directories(${GTEST_INCLUDE_DIRS})
add_executable(tests meshdist_test.cc)
target_link_libraries(tests ${GTEST_BOTH_LIBRARIES})
target_link_libraries(tests meshdist)
