add_executable(test_erf test_erf.cpp)
target_link_libraries(test_erf PRIVATE corvus::corvus)
add_test(NAME erf COMMAND test_erf)
