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

file(GLOB pfire_libSOURCES
  "*.cpp")

list(REMOVE_ITEM pfire_libSOURCES "${CMAKE_CURRENT_SOURCE_DIR}/sandbox.cpp")
list(REMOVE_ITEM pfire_libSOURCES "${CMAKE_CURRENT_SOURCE_DIR}/pfire.cpp")

if(NOT OpenImageIO_FOUND)
  list(REMOVE_ITEM pfire_libSOURCES "${CMAKE_CURRENT_SOURCE_DIR}/oiioloader.cpp")
  list(REMOVE_ITEM pfire_libSOURCES "${CMAKE_CURRENT_SOURCE_DIR}/oiiowriter.cpp")
endif(NOT OpenImageIO_FOUND)

if(NOT DCMTK_FOUND)
  list(REMOVE_ITEM pfire_libSOURCES "${CMAKE_CURRENT_SOURCE_DIR}/dcmloader.cpp")
endif(NOT DCMTK_FOUND)

cmake_print_variables(${pfire_libSOURCES})
cmake_print_variables(${EXTRA_LIBS} )
cmake_print_variables(${DCMTK_LIBRARIES})


#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_library(pfire_lib SHARED  ${pfire_libSOURCES} )

#set_target_properties(pfire_lib PROPERTIES PREFIX "")
#target_link_libraries(pfire_lib ${EXTRA_LIBS} ) # ${DCMTK_LIBRARIES})
target_link_libraries(pfire_lib PUBLIC ${EXTRA_LIBS})
target_link_libraries(pfire_lib PUBLIC OpenImageIO::OpenImageIO )#${OpenImageIO_LIBRARIES})
target_link_libraries(pfire_lib PUBLIC DCMTK::dcmtk)
# alternative approach based on modern CMAKE
#target_link_libraries(pfire_lib ${DCMTK_LIBRARIES})



add_executable(pfire pfire.cpp)
set_target_properties(pfire PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(pfire  PUBLIC pfire_lib)


install(TARGETS pfire DESTINATION bin)

install(TARGETS pfire_lib DESTINATION lib)
