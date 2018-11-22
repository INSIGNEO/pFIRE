file(GLOB libpFIRESOURCES
  "*.cpp")

list(REMOVE_ITEM libpFIRESOURCES "${CMAKE_CURRENT_SOURCE_DIR}/sandbox.cpp")
list(REMOVE_ITEM libpFIRESOURCES "${CMAKE_CURRENT_SOURCE_DIR}/pfire.cpp")

add_library(libpfire ${libpFIRESOURCES})
set_target_properties(libpfire PROPERTIES PREFIX "")
target_link_libraries(libpfire ${EXTRA_LIBS})

add_executable(pfire pfire.cpp)
set_target_properties(pfire PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(pfire libpfire)


add_executable(sandbox sandbox.cpp)
set_target_properties(sandbox PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(sandbox libpfire ${EXTRA_LIBS})
