# This script defines the following variables:
#
#   TRIMESH2_LIBRARY
#   TRIMESH2_INCLUDE_DIR 

if (NOT COMMAND ExternalProject_Add)
    include(ExternalProject)
endif ()

if (NOT DEFINED PROCESSOR_COUNT)
    include(ProcessorCount)
    ProcessorCount(PROCESSOR_COUNT)
endif ()
if (PROCESSOR_COUNT EQUAL 0)
    message(WARNING "Unable to determine the number of processors. Default: PROCESSOR_COUNT=1")
    set(${PROCESSOR_COUNT} 1)
endif ()

# Download and build.

ExternalProject_Add(
        trimesh2_project
        URL http://gfx.cs.princeton.edu/proj/trimesh2/src/trimesh2-2.12.tar.gz
        URL_HASH SHA1=8868efe10a27b95af9f221b7d767c5dd5f1f14ec
        SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/trimesh2
        CONFIGURE_COMMAND ""
        BUILD_COMMAND make -j${PROCESSOR_COUNT}
        INSTALL_COMMAND ""
        BUILD_IN_SOURCE 1
)

ExternalProject_Get_Property(trimesh2_project SOURCE_DIR)

add_library(trimesh2 STATIC IMPORTED GLOBAL)

if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set_property(
            TARGET trimesh2
            PROPERTY IMPORTED_LOCATION ${SOURCE_DIR}/lib.Linux64/libtrimesh.a
    )
else ()
    set_property(
            TARGET trimesh2
            PROPERTY IMPORTED_LOCATION ${SOURCE_DIR}/lib.Linux/libtrimesh.a
    )
endif ()

add_dependencies(trimesh2 trimesh2_project)

set(TRIMESH2_LIBRARY "trimesh2"
        CACHE INTERNAL "")
set(TRIMESH2_INCLUDE_DIR "${SOURCE_DIR}/include"
        CACHE INTERNAL "")

unset(SOURCE_DIR)
