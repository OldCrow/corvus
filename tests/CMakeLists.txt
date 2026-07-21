add_executable(test_erf test_erf.cpp)
target_link_libraries(test_erf PRIVATE corvus::corvus)
add_test(NAME erf COMMAND test_erf)

add_executable(bench_erf bench_erf.cpp)
target_link_libraries(bench_erf PRIVATE corvus::corvus)
# Deliberately not registered with ctest: run manually on a quiet machine.

add_executable(test_erf_ulp test_erf_ulp.cpp)
target_link_libraries(test_erf_ulp PRIVATE corvus::corvus)
add_test(NAME erf_ulp
  COMMAND test_erf_ulp ${CMAKE_CURRENT_SOURCE_DIR}/data/erf_reference.txt)
