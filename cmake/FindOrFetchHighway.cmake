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
# together with the pin and a revalidation pass.

find_package(hwy 1.4 CONFIG QUIET)

if(hwy_FOUND)
  set(CORVUS_HWY_PROVIDER "system")
  message(STATUS "corvus: using system Highway ${hwy_VERSION}")
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
