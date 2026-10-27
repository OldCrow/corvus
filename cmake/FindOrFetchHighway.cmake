# Prefer an installed Highway; fall back to fetching a pinned version.
# Sets CORVUS_HWY_PROVIDER to "system" or "fetched".
#
# Version floor matches the accuracy-audit pin below: every accuracy
# claim was validated against Highway 1.4.x, and older system Highway
# (e.g. Ubuntu's libhwy-dev 1.0.x) predates parts of the surface the
# facade has needed over time (the floor was originally established on
# hn::ReduceMax, removed from the facade at #10; the audited pin, not
# any single op, is what the floor protects). Better to fall back to
# the fetched pin than to accept any system version. Bump the floor
# together with the pin and a revalidation pass, here and in
# cmake/corvus-config.cmake.in's find_dependency.

find_package(hwy 1.4 CONFIG QUIET)

if(hwy_FOUND)
  set(CORVUS_HWY_PROVIDER "system")
  message(STATUS "corvus: using system Highway ${hwy_VERSION}")

  # macOS, Apple Clang, CMAKE_OSX_SYSROOT empty (CMake's default when the
  # selected Xcode's SDK is the compiler's own, so no -isysroot is passed):
  # clang then adds /usr/local/include to its builtin search list as a
  # USER directory, and CMake, seeing the same directory among the
  # compiler's implicit include directories, emits no -isystem for the
  # imported hwy::hwy target at all. Highway's headers are therefore
  # compiled as our own code and hwy/base.h trips -Wpedantic
  # (-Wgnu-zero-variadic-macro-arguments) under -Werror (#62). With any
  # explicit -isysroot clang drops its builtin /usr/local/include, the
  # directory is no longer implicit, and CMake's -isystem does the job.
  # Restoring the flag by hand makes clang discard its own duplicate
  # ("ignoring duplicate directory") and keep the system one. Scoped to
  # Apple Clang on macOS: GCC already treats its builtin dirs as system,
  # and -isystem on a libstdc++ host's /usr/include breaks #include_next.
  if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    get_target_property(_corvus_hwy_incs hwy::hwy INTERFACE_INCLUDE_DIRECTORIES)
    foreach(_dir IN LISTS _corvus_hwy_incs)
      if(_dir IN_LIST CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES)
        set_property(TARGET hwy::hwy APPEND PROPERTY
          INTERFACE_COMPILE_OPTIONS "$<$<COMPILE_LANGUAGE:CXX>:-isystem${_dir}>")
        message(STATUS "corvus: Highway include ${_dir} is implicit for this "
          "compiler; passing -isystem explicitly so its headers stay system headers (#62)")
      endif()
    endforeach()
    unset(_corvus_hwy_incs)
  endif()
else()
  include(FetchContent)
  set(HWY_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
  set(HWY_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(HWY_ENABLE_CONTRIB ON CACHE BOOL "" FORCE)
  set(HWY_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
  set(HWY_FORCE_STATIC_LIBS ON CACHE BOOL "" FORCE)
  set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

  # Pin matches the version the accuracy audit was validated against
  # (docs/ACCURACY.md); bump only together with a revalidation pass.
  # SYSTEM: fetched Highway headers must not trip our -Wpedantic -Werror
  # (installed Highway is an imported target and gets this implicitly).
  FetchContent_Declare(highway
    GIT_REPOSITORY https://github.com/google/highway.git
    GIT_TAG        1.4.0
    GIT_SHALLOW    TRUE
    SYSTEM)
  FetchContent_MakeAvailable(highway)
  set(CORVUS_HWY_PROVIDER "fetched")
  message(STATUS "corvus: fetched Highway 1.4.0 via FetchContent")
endif()

if(TARGET hwy AND NOT TARGET hwy::hwy)
  add_library(hwy::hwy ALIAS hwy)
endif()
