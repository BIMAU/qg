find_package(Matlab COMPONENTS MAIN_PROGRAM)

include_directories(${PROJECT_SOURCE_DIR}/src)

set(QG_MEX_FILES
  QG_init
  QG_free
  QG_rhs
  QG_bilin
  QG_jacob
  QG_jacobian
  QG_mass
  QG_apply
  QG_compute_precon
  QG_solve
  QG_set_par
  QG_get_par)

foreach(i ${QG_MEX_FILES})
  matlab_add_mex(NAME "${i}" SRC "${i}.cpp")

  # Matlab 2018b hack
  target_compile_options(${i} PRIVATE "-fvisibility=default")

  target_link_libraries(${i} qg)

  install(TARGETS ${i} DESTINATION matlab)
endforeach()

set(QG_MATLAB_FILES
  QG.m
  plotQG.m
  crs2sp.m
  stochforcing.m)

install(FILES ${QG_MATLAB_FILES} DESTINATION matlab)

enable_testing()

set(MATLAB_TEST_FILES
  test_basic_things
  test_jacobian
  test_periodic)

foreach(i ${MATLAB_TEST_FILES})
  matlab_add_unit_test(
    NAME "${i}"
    TIMEOUT 30
    UNITTEST_FILE "${CMAKE_CURRENT_SOURCE_DIR}/MLQG_tests/${i}.m"
    ADDITIONAL_PATH "${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_CURRENT_SOURCE_DIR}"
    )
endforeach()

set(CTEST_OPTIONS --force-new-ctest-process)
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_OPTIONS})
