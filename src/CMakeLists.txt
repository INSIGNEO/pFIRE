file(GLOB libpFIRESOURCES
  "*.cpp")

list(REMOVE_ITEM libpFIRESOURCES "${CMAKE_CURRENT_SOURCE_DIR}/sandbox.cpp")
list(REMOVE_ITEM libpFIRESOURCES "${CMAKE_CURRENT_SOURCE_DIR}/pfire.cpp")

if(NOT OIIO_FOUND)
  list(REMOVE_ITEM libpFIRESOURCES "${CMAKE_CURRENT_SOURCE_DIR}/oiioloader.cpp")
endif(NOT OIIO_FOUND)

add_library(libpfire ${libpFIRESOURCES})
target_link_libraries(libpfire ${EXTRA_LIBS})

add_executable(pfire pfire.cpp)
set_target_properties(pfire PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(pfire libpfire)
