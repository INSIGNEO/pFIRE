find_package (Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
include_directories (${pFIRE_SOURCE_DIR}/src ${Boost_INCLUDE_DIRS})
set(EXTRA_LIBS ${EXTRA_LIBS} ${Boost_LIBRARIES})
add_definitions(-DBOOST_TEST_DYN_LINK=1)

function(add_mpi_test name executable nprocs)
  add_executable(${executable} ${executable}.cpp)
  target_link_libraries(${executable} libpfire ${Boost_LIBRARIES})
  set(mpi_test_params -np ${nprocs} $<TARGET_FILE:${executable}>)
  add_test(NAME ${name}-mpi COMMAND mpirun ${mpi_test_params}
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(${name}-mpi PROPERTIES LABELS mpi)
  add_test(NAME ${name}-serial COMMAND $<TARGET_FILE:${executable}>
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(${name}-serial PROPERTIES LABELS serial)
endfunction(add_mpi_test)

#add_mpi_test(Indexing test_indexing 4)
#add_mpi_test(Gradients test_gradients 4)
#add_mpi_test(HDFWriter test_hdfwriter 4)
#add_mpi_test(ShIRTLoader test_shirtloader 4)

add_mpi_test(warp test_warp 4)
add_mpi_test(hdfwriter test_hdfwriter 4)
