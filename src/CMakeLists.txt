file(GLOB pFIRESOURCES
  "*.cpp")

list(REMOVE_ITEM pFIRESOURCES "${CMAKE_CURRENT_SOURCE_DIR}/sandbox.cpp")
message("Removing sandbox.cpp from pFIRE glob")

add_executable(pfire ${pFIRESOURCES})
set_target_properties(pfire PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(pfire ${EXTRA_LIBS})

#set (SANDBOXSOURCES
#  sandbox
#  )

#add_executable(sandbox ${SANDBOXSOURCES})
#target_link_libraries(sandbox ${EXTRA_LIBS})
