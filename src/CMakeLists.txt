file(GLOB libpFIRESOURCES
  "*.cpp")

list(REMOVE_ITEM libpFIRESOURCES "${CMAKE_CURRENT_SOURCE_DIR}/sandbox.cpp")
list(REMOVE_ITEM libpFIRESOURCES "${CMAKE_CURRENT_SOURCE_DIR}/pfire.cpp")

if(NOT OPENIMAGEIO_FOUND)
  list(REMOVE_ITEM libpFIRESOURCES "${CMAKE_CURRENT_SOURCE_DIR}/oiioloader.cpp")
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
