find_package (Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
include_directories (${pFIRE_SOURCE_DIR}/src ${Boost_INCLUDE_DIRS})
set(EXTRA_LIBS ${EXTRA_LIBS} ${Boost_LIBRARIES})

add_definitions(-DBOOST_TEST_DYN_LINK=1)
add_executable(test_indexing test_indexing.cpp)
target_link_libraries(test_indexing libpfire ${Boost_LIBRARIES})
add_test(NAME Indexing COMMAND test_indexing)

add_definitions(-DBOOST_TEST_DYN_LINK=1)
add_executable(test_workspace test_workspace.cpp)
target_link_libraries(test_workspace libpfire ${Boost_LIBRARIES})
add_test(NAME Workspace COMMAND test_workspace)

add_definitions(-DBOOST_TEST_DYN_LINK=1)
add_executable(test_gradients test_gradients.cpp)
target_link_libraries(test_gradients libpfire ${Boost_LIBRARIES})
add_test(NAME Gradients COMMAND test_gradients)

add_executable(test_map test_map.cpp)
target_link_libraries(test_map libpfire ${Boost_LIBRARIES})
add_test(Name Map COMMAND test_map)
