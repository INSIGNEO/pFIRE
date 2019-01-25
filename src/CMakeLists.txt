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

file(GLOB libpFIRESOURCES
  "*.cpp")

list(REMOVE_ITEM libpFIRESOURCES "${CMAKE_CURRENT_SOURCE_DIR}/sandbox.cpp")
list(REMOVE_ITEM libpFIRESOURCES "${CMAKE_CURRENT_SOURCE_DIR}/pfire.cpp")

if(NOT OPENIMAGEIO_FOUND)
  list(REMOVE_ITEM libpFIRESOURCES "${CMAKE_CURRENT_SOURCE_DIR}/oiioloader.cpp")
  list(REMOVE_ITEM libpFIRESOURCES "${CMAKE_CURRENT_SOURCE_DIR}/oiiowriter.cpp")
endif(NOT OPENIMAGEIO_FOUND)

if(NOT DCMTK_FOUND)
  list(REMOVE_ITEM libpFIRESOURCES "${CMAKE_CURRENT_SOURCE_DIR}/dcmloader.cpp")
endif(NOT DCMTK_FOUND)

add_library(libpfire ${libpFIRESOURCES})
set_target_properties(libpfire PROPERTIES PREFIX "")
target_link_libraries(libpfire ${EXTRA_LIBS})

add_executable(pfire pfire.cpp)
set_target_properties(pfire PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(pfire libpfire)
