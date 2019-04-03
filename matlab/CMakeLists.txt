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
  QG.m)

install(FILES ${QG_MATLAB_FILES} DESTINATION matlab)

# enable_testing()
# matlab_add_unit_test(
#     NAME stokes0
#     TIMEOUT 30
#     UNITTEST_FILE ${CMAKE_CURRENT_SOURCE_DIR}/test/test_stokes0.m
#     ADDITIONAL_PATH "${CMAKE_CURRENT_BINARY_DIR}:${CMAKE_CURRENT_SOURCE_DIR}"
# )
