# Prefer an installed Highway; fall back to fetching a pinned version.
# Sets CORVUS_HWY_PROVIDER to "system" or "fetched".

find_package(hwy CONFIG QUIET)

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
  FetchContent_Declare(highway
    GIT_REPOSITORY https://github.com/google/highway.git
    GIT_TAG        1.4.0
    GIT_SHALLOW    TRUE)
  FetchContent_MakeAvailable(highway)
  set(CORVUS_HWY_PROVIDER "fetched")
  message(STATUS "corvus: fetched Highway 1.4.0 via FetchContent")
endif()

if(TARGET hwy AND NOT TARGET hwy::hwy)
  add_library(hwy::hwy ALIAS hwy)
endif()
