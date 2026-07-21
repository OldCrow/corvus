cmake_minimum_required(VERSION 3.25)  # 3.25: FetchContent SYSTEM keyword
project(corvus VERSION 0.1.0 LANGUAGES CXX
  DESCRIPTION "SIMD-vectorized statistical special functions"
  HOMEPAGE_URL "https://github.com/OldCrow/corvus")

# --- Build type ---------------------------------------------------------
# Single-config generators default to Release: this is a numerical library
# and Debug-by-accident produces misleading benchmarks (house rule: perf
# numbers from optimized builds only). Intended roles:
#   Release        -O3, no debug info: benchmarks, accuracy validation, dist
#   RelWithDebInfo -O2 -g: profiling (Instruments/perf need symbols)
#   Debug          -O0 -g: debugger sessions; pair with CORVUS_SANITIZE
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
  STRINGS Release RelWithDebInfo Debug)

# compile_commands.json for clangd/IDEs; harmless for consumers.
if(PROJECT_IS_TOP_LEVEL)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# --- Options (all CORVUS_-prefixed) -------------------------------------
option(CORVUS_BUILD_TESTS "Build corvus tests" ${PROJECT_IS_TOP_LEVEL})
option(CORVUS_DEV_WARNINGS "Enable -Wall -Wextra -Wpedantic (dev builds)"
       ${PROJECT_IS_TOP_LEVEL})
option(CORVUS_WERROR "Treat warnings as errors (CI)" OFF)

# Compile-time target restriction, e.g. -DCORVUS_DISABLED_TARGETS="HWY_AVX3|HWY_AVX2"
# forces dispatch to stop at the next-lower tier. Used for native per-tier
# validation on a machine whose best tier is higher (Ryzen AVX-512 -> AVX2/SSE2).
set(CORVUS_DISABLED_TARGETS "" CACHE STRING
    "Pipe-separated HWY_* target macros to exclude from compilation/dispatch")

# Sanitizers for Debug-type dev builds, e.g. -DCORVUS_SANITIZE=address;undefined
# (UBSan is cheap insurance for numerical kernel code; not for benchmarks).
set(CORVUS_SANITIZE "" CACHE STRING
    "Semicolon-separated sanitizers (address;undefined;thread)")

include(cmake/FindOrFetchHighway.cmake)

add_library(corvus STATIC
  src/erf.cpp
  src/erf_data.cpp
  src/erfc.cpp
)
add_library(corvus::corvus ALIAS corvus)

# PUBLIC: the header uses std::span, so consumers need C++20 too — the
# feature requirement travels through the installed package config.
target_compile_features(corvus PUBLIC cxx_std_20)
set_target_properties(corvus PROPERTIES
  CXX_EXTENSIONS OFF
  # PIC so the static lib can later be linked into shared objects
  # (pybind11 bindings per the pylibhmm/pylibstats pattern).
  POSITION_INDEPENDENT_CODE ON)

target_include_directories(corvus
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(corvus PRIVATE hwy::hwy)

# Dev warnings: PRIVATE only — never exported to consumers. CORVUS_WARN_FLAGS
# is also applied to test/bench targets in tests/CMakeLists.txt.
set(CORVUS_WARN_FLAGS "")
if(CORVUS_DEV_WARNINGS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    list(APPEND CORVUS_WARN_FLAGS -Wall -Wextra -Wpedantic)
    if(CORVUS_WERROR)
      list(APPEND CORVUS_WARN_FLAGS -Werror)
    endif()
  elseif(MSVC)
    list(APPEND CORVUS_WARN_FLAGS /W4)
    if(CORVUS_WERROR)
      list(APPEND CORVUS_WARN_FLAGS /WX)
    endif()
  endif()
endif()
target_compile_options(corvus PRIVATE ${CORVUS_WARN_FLAGS})

if(CORVUS_SANITIZE)
  list(JOIN CORVUS_SANITIZE "," _corvus_san)
  target_compile_options(corvus PRIVATE -fsanitize=${_corvus_san})
  target_link_options(corvus PUBLIC -fsanitize=${_corvus_san})
endif()

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
