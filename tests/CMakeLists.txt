# Apply the same dev warnings as the library (CORVUS_WARN_FLAGS from the
# top-level lists; empty unless CORVUS_DEV_WARNINGS).
function(corvus_test_target name)
  target_link_libraries(${name} PRIVATE corvus::corvus)
  target_compile_options(${name} PRIVATE ${CORVUS_WARN_FLAGS})
endfunction()

# A test that compiles an INTERNAL kernel header instead of linking the public
# API: it needs Highway and the source root directly, plus the same
# target-selection defines the library gets (see CORVUS_HWY_TARGET_DEFS).
function(corvus_kernel_test_target name)
  # corvus::corvus for expect_target.h (active_target); hwy::hwy because the
  # test compiles kernel headers itself.
  target_link_libraries(${name} PRIVATE corvus::corvus hwy::hwy)
  target_include_directories(${name} PRIVATE
    ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_definitions(${name} PRIVATE ${CORVUS_HWY_TARGET_DEFS})
  target_compile_features(${name} PRIVATE cxx_std_20)
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

add_executable(test_exp_dd test_exp_dd.cpp)
corvus_kernel_test_target(test_exp_dd)
add_test(NAME exp_dd
  COMMAND test_exp_dd ${CMAKE_CURRENT_SOURCE_DIR}/data/exp_dd_reference.txt)

add_executable(test_log_dd test_log_dd.cpp)
corvus_kernel_test_target(test_log_dd)
add_test(NAME log_dd
  COMMAND test_log_dd ${CMAKE_CURRENT_SOURCE_DIR}/data/log_dd_reference.txt)
