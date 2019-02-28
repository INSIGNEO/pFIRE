find_package (Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
include_directories (${pFIRE_SOURCE_DIR}/src ${Boost_INCLUDE_DIRS})
set(EXTRA_LIBS ${EXTRA_LIBS} ${Boost_LIBRARIES})
add_definitions(-DBOOST_TEST_DYN_LINK=1)

function(add_mpi_test name executable nprocs)
  set(test_params -np ${nprocs} $<TARGET_FILE:${executable}>)
  add_test(NAME ${name} COMMAND mpirun ${test_params}
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endfunction(add_mpi_test)

add_executable(test_indexing test_indexing.cpp)
target_link_libraries(test_indexing libpfire ${Boost_LIBRARIES})
add_mpi_test(Indexing test_indexing 4)

add_executable(test_workspace test_workspace.cpp)
target_link_libraries(test_workspace libpfire ${Boost_LIBRARIES})
add_mpi_test(Workspace test_workspace 4)

add_executable(test_gradients test_gradients.cpp)
target_link_libraries(test_gradients libpfire ${Boost_LIBRARIES})
add_mpi_test(Gradients test_gradients 4)

add_executable(test_hdfwriter test_hdfwriter.cpp)
target_link_libraries(test_hdfwriter libpfire ${Boost_LIBRARIES})
add_mpi_test(HDFWriter test_hdfwriter 4)

add_executable(test_shirtloader test_shirtloader.cpp)
target_link_libraries(test_shirtloader libpfire ${Boost_LIBRARIES})
add_mpi_test(ShIRTLoader test_shirtloader 4)
