# Examples are user-facing programs, not tests — but they verify their own
# output and return non-zero on failure, so CI can run them as smoke tests and
# they cannot rot into printing wrong numbers convincingly.
#
# Deliberately NOT registered with ctest: the gates in tests/ are the accuracy
# contract, and mixing demonstrations into that list would blur what a green
# ctest run means.
function(corvus_example name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE corvus::corvus)
  target_compile_options(${name} PRIVATE ${CORVUS_WARN_FLAGS})
  set_target_properties(${name} PROPERTIES FOLDER examples)
endfunction()

corvus_example(normal_distribution)
