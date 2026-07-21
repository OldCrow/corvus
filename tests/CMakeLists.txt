# Apply the same dev warnings as the library (CORVUS_WARN_FLAGS from the
# top-level lists; empty unless CORVUS_DEV_WARNINGS).
function(corvus_test_target name)
  target_link_libraries(${name} PRIVATE corvus::corvus)
  target_compile_options(${name} PRIVATE ${CORVUS_WARN_FLAGS})
endfunction()

add_executable(test_erf test_erf.cpp)
corvus_test_target(test_erf)
add_test(NAME erf COMMAND test_erf)

add_executable(bench_erf bench_erf.cpp)
corvus_test_target(bench_erf)
# Deliberately not registered with ctest: run manually on a quiet machine.

add_executable(test_erf_ulp test_erf_ulp.cpp)
corvus_test_target(test_erf_ulp)
add_test(NAME erf_ulp
  COMMAND test_erf_ulp ${CMAKE_CURRENT_SOURCE_DIR}/data/erf_reference.txt)

add_executable(test_erfc test_erfc.cpp)
corvus_test_target(test_erfc)
add_test(NAME erfc COMMAND test_erfc)

add_executable(test_erfc_ulp test_erfc_ulp.cpp)
corvus_test_target(test_erfc_ulp)
add_test(NAME erfc_ulp
  COMMAND test_erfc_ulp ${CMAKE_CURRENT_SOURCE_DIR}/data/erfc_reference.txt)

add_executable(bench_erfc bench_erfc.cpp)
corvus_test_target(bench_erfc)
# Deliberately not registered with ctest: run manually on a quiet machine.
