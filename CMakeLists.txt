cmake_minimum_required(VERSION 3.24)
project(corvus VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(CORVUS_BUILD_TESTS "Build corvus tests" ${PROJECT_IS_TOP_LEVEL})

# Compile-time target restriction, e.g. -DCORVUS_DISABLED_TARGETS="HWY_AVX3|HWY_AVX2"
# forces dispatch to stop at the next-lower tier. Used for native per-tier
# validation on a machine whose best tier is higher (Ryzen AVX-512 -> AVX2/SSE2).
set(CORVUS_DISABLED_TARGETS "" CACHE STRING
    "Pipe-separated HWY_* target macros to exclude from compilation/dispatch")

include(cmake/FindOrFetchHighway.cmake)

add_library(corvus STATIC
  src/erf.cpp
)
add_library(corvus::corvus ALIAS corvus)

target_include_directories(corvus
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(corvus PRIVATE hwy::hwy)

if(CORVUS_DISABLED_TARGETS)
  target_compile_definitions(corvus PRIVATE
    "HWY_DISABLED_TARGETS=(${CORVUS_DISABLED_TARGETS})")
endif()

if(CORVUS_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# Install/export is only supported against a system (find_package) Highway.
# A FetchContent-built Highway is not in our export set; installing corvus
# would then produce a broken package. See PLAN.md open items.
if(CORVUS_HWY_PROVIDER STREQUAL "system")
  include(GNUInstallDirs)
  include(CMakePackageConfigHelpers)

  install(TARGETS corvus EXPORT corvus-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY include/corvus DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(EXPORT corvus-targets
    NAMESPACE corvus::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/corvus)

  configure_package_config_file(
    cmake/corvus-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/corvus-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/corvus)
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/corvus-config-version.cmake
    COMPATIBILITY SameMajorVersion)
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/corvus-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/corvus-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/corvus)
else()
  message(STATUS "corvus: Highway was fetched, not found — 'install' target disabled")
endif()
