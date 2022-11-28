#
#   Copyright 2019 University of Sheffield
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http:#www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

cmake_minimum_required(VERSION 3.7.0)
project(pFIRE)
include(CMakePrintHelpers)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not permitted. Make a separate folder for "
          "building:\nmkdir build; cd build; cmake ..\n"
          "Before that, remove the files already created:\nrm -rf CMakeCache.txt CMakeFiles")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

cmake_policy(SET CMP0074 NEW)

if( NOT CMAKE_BUILD_TYPE )
    set( CMAKE_BUILD_TYPE Debug CACHE STRING
             "Choose the type of build, options are: None Debug Release"
        FORCE )
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wuninitialized -O0 -ggdb -march=native -fstack-protector-strong")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O3 -march=native -fno-fast-math")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO  "${CMAKE_CXX_FLAGS_RELEASE} -ggdb")

option(DEBUG_VERBOSE "Verbose debug outputs")
if(DEBUG_VERBOSE)
  add_definitions(-DDEBUG_VERBOSE)
  message(STATUS "Enabled verbose debugging messages")
endif(DEBUG_VERBOSE)

option(DEBUG_CHECKS "Enable internal debugging checks (slow)")
if(DEBUG_CHECKS)
  add_definitions(-DDEBUG_CHECKS)
  message(STATUS "Enabled internal debugging checks (will reduce performance)")
endif(DEBUG_CHECKS)

option(DEBUG_PETSC "Debug PETSc object creation and deletion")
if(DEBUG_PETSC)
  add_definitions(-DDEBUG_PETSC_OBJECTS)
  message(STATUS "Enabled debugging of PETSc object lifetimes")
endif(DEBUG_PETSC)

option(DEBUG_OUTPUT_INTERMEDIATES "Save all intermediate frames in registration")
if(DEBUG_OUTPUT_INTERMEDIATES)
  add_definitions(-DDEBUG_DUMP_INTERMEDIATES)
endif(DEBUG_OUTPUT_INTERMEDIATES)

set(GITCPP_PATH "src/gitstate.cpp")
include(cmake-modules/GitState.cmake)
#WriteGitStateFile("${pFIRE_SOURCE_DIR}")
#WriteGitVersionSrc("${pFIRE_SOURCE_DIR}/src/gitstate.cpp")

# Output programs to a "bin" directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${pFIRE_SOURCE_DIR}/bin)

# Ensure we can find FindPETSc and friends
list(APPEND CMAKE_MODULE_PATH ${pFIRE_SOURCE_DIR}/cmake-modules)

# Define flags to enable dependencies
option(ENABLE_OpenImageIO "Enable Linking OpenImageIO library" [FALSE])
option(ENABLE_DCMTK       "Enable Linking DCMTK library" [FALSE])

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDES})
set(EXTRA_LIBS ${EXTRA_LIBS} ${MPI_LIBRARIES})

find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
message(STATUS "HDF5_INCLUDES: ${HDF5_INCLUDE_DIRS}")
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

if(ENABLE_DCMTK)
find_package(DCMTK)
endif(ENABLE_DCMTK)

if(DCMTK_FOUND)
include_directories(${DCMTK_INCLUDES})
set(EXTRA_LIBS ${EXTRA_LIBS} ${DCMTK_LIBRARIES})
  add_definitions(-DUSE_DCMTK)
endif(DCMTK_FOUND)

#list(APPEND CMAKE_FIND_ROOT_PATH "... insigneo-pfire/docker-pfire/singularity/fix_oiio_dcmtk/oiio/dist/lib64/cmake/OpenImageIO/")


if(ENABLE_OpenImageIO)
find_package(OpenImageIO CONFIG)
endif(ENABLE_OpenImageIO)

cmake_print_variables(OpenImageIO_FOUND)



if(OpenImageIO_FOUND AND ENABLE_OpenImageIO)
	
	#cmake_print_variables(OpenImageIO_INCLUDE_DIR)
	#cmake_print_variables(OpenImageIO_LIBRARIES)
	#cmake_print_variables(OpenImageIO_LIB_DIR)
	# OpenImageIO cmake config is missing library names so declaring here until they fix it
	set (OpenImageIO_LIBRARIES "${OpenImageIO_LIB_DIR}/libOpenImageIO_Util.so;${OpenImageIO_LIB_DIR}/libOpenImageIO.so")
	include_directories(${OpenImageIO_INCLUDE_DIR})
	set(EXTRA_LIBS ${EXTRA_LIBS} ${OpenImageIO_LIBRARIES} )
	add_definitions(-DUSE_OIIO)
	
			
else(OpenImageIO_FOUND AND ENABLE_OpenImageIO)
	message(STATUS "OPENIMAGEIO_FOUND ELSE branch")
	
	include_directories(${OpenImageIO_INCLUDE_DIR})
	set(EXTRA_LIBS ${EXTRA_LIBS} ${OpenImageIO_LIBRARIES} ) #/usr/local/oiio/lib64/libOpenImageIO.a /usr/local/oiio/lib64/libOpenImageIO_Util.a)
	
	add_definitions(-DUSE_OIIO)

endif(OpenImageIO_FOUND AND ENABLE_OpenImageIO)

#find_library(JPEG_LIBRARY NAMES jpeg jpg libjpeg  REQUIRED NO_SYSTEM_ENVIRONMENT_PATH PATHS /usr/local/oiio/lib64)
#cmake_print_variables(JPEG_LIBRARY)

set(EXTRA_LIBS ${EXTRA_LIBS} ${OpenImageIO_LIBRARIES})

add_subdirectory(${pFIRE_SOURCE_DIR}/src)

option(BUILD_TESTS "Built pFIRE unit tests" ON)
if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(unit_test)
endif(BUILD_TESTS)

include(FeatureSummary)
feature_summary(WHAT ALL)
