cmake_minimum_required(VERSION 3.7.0)
project(pFIRE)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

if( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE Debug CACHE STRING
             "Choose the type of build, options are: None Debug Release"
        FORCE )
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Weffc++ -O0 -ggdb -march=native -fno-stack-protector")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O3 -march=native")

# Output programs to a "bin" directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${pFIRE_SOURCE_DIR}/bin)

# Ensure we can find FindPETSc and friends
list(APPEND CMAKE_MODULE_PATH ${pFIRE_SOURCE_DIR}/cmake-modules)

find_package(PETSc REQUIRED)
include_directories(${PETSC_INCLUDES})
set(EXTRA_LIBS ${EXTRA_LIBS} ${PETSC_LIBRARIES})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDES})
set(EXTRA_LIBS ${EXTRA_LIBS} ${MPI_LIBRARIES})

find_package(OpenImageIO REQUIRED)
include_directories(${OPENIMAGEIO_INCLUDE_DIR})
set(EXTRA_LIBS ${EXTRA_LIBS} ${OPENIMAGEIO_LIBRARIES})

add_subdirectory(${pFIRE_SOURCE_DIR}/src)
