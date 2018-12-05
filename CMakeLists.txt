cmake_minimum_required(VERSION 3.7.0)
project(pFIRE)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not permitted. Make a separate folder for"
          "building:\nmkdir build; cd build; cmake ..\n"
          "Before that, remove the files already created:\nrm -rf CMakeCache.txt CMakeFiles")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

if( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE Debug CACHE STRING
             "Choose the type of build, options are: None Debug Release"
        FORCE )
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Weffc++ -Wuninitialized -O0 -ggdb -march=native -fstack-protector-strong")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O3 -march=native -fno-fast-math")

option(DEBUG_VERBOSE "Verbose debug outputs")
if(DEBUG_VERBOSE)
  add_definitions(-DDEBUG_VERBOSE)
  message(STATUS "Enabled verbose debugging messages")
endif(DEBUG_VERBOSE)

option(DEBUG_OUTPUT_INTERMEDIATES "Save all intermediate frames in registration")
if(DEBUG_OUTPUT_INTERMEDIATES)
  add_definitions(-DDEBUG_DUMP_INTERMEDIATES)
endif(DEBUG_OUTPUT_INTERMEDIATES)


# Output programs to a "bin" directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${pFIRE_SOURCE_DIR}/bin)

# Ensure we can find FindPETSc and friends
list(APPEND CMAKE_MODULE_PATH ${pFIRE_SOURCE_DIR}/cmake-modules)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDES})
set(EXTRA_LIBS ${EXTRA_LIBS} ${MPI_LIBRARIES})

find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDES})
set(EXTRA_LIBS ${EXTRA_LIBS} ${HDF5_LIBRARIES})

if(NOT HDF5_IS_PARALLEL)
  message(FATAL_ERROR "Parallel HDF5 is required.")
endif(NOT HDF5_IS_PARALLEL)

find_package(PETSc REQUIRED)
include_directories(${PETSC_INCLUDES})
set(EXTRA_LIBS ${EXTRA_LIBS} ${PETSC_LIBRARIES})
add_definitions(-DPETSC_CLANGUAGE_CXX=1)

find_package(Boost REQUIRED COMPONENTS filesystem program_options)
include_directories(${BOOST_INCLUDE_DIRS})
set(EXTRA_LIBS ${EXTRA_LIBS} ${Boost_LIBRARIES})


find_package(DCMTK)
if(DCMTK_FOUND)
include_directories(${DCMTK_INCLUDES})
set(EXTRA_LIBS ${EXTRA_LIBS} ${DCMTK_LIBRARIES})
  add_definitions(-DUSE_DCMTK)
endif(DCMTK_FOUND)

find_package(OpenImageIO)
if(OPENIMAGEIO_FOUND)
include_directories(${OPENIMAGEIO_INCLUDE_DIR})
set(EXTRA_LIBS ${EXTRA_LIBS} ${OPENIMAGEIO_LIBRARIES})
  add_definitions(-DUSE_OIIO)
endif(OPENIMAGEIO_FOUND)

add_subdirectory(${pFIRE_SOURCE_DIR}/src)

enable_testing()
add_subdirectory(test)

include(FeatureSummary)
feature_summary(WHAT ALL)
