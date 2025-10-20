# Purpose: pick a vcpkg location early, then chain-load the real vcpkg toolchain

# prioritize BITLOOP_ROOT for a single souce of truth (if installed)
if (DEFINED ENV{VCPKG_ROOT} AND EXISTS "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake")
  set(_vcpkg_dir "$ENV{VCPKG_ROOT}")
  message(STATUS "vcpkg toolchain detected at: ${_vcpkg_dir}")
elseif (DEFINED ENV{BITLOOP_ROOT} AND EXISTS "$ENV{BITLOOP_ROOT}/vcpkg/scripts/buildsystems/vcpkg.cmake")
  set(_vcpkg_dir "$ENV{BITLOOP_ROOT}/vcpkg")
  message(STATUS "vcpkg toolchain detected at: ${_vcpkg_dir}")
else()
  message(FATAL_ERROR
    "vcpkg not found.\n"
    "Tried:\n"
    "  - $ENV{VCPKG_ROOT} (if set)\n"
    "  - $ENV{BITLOOP_ROOT}/vcpkg (if set)\n")
endif()

# 5) Export and chain-load
set(ENV{VCPKG_ROOT} "${_vcpkg_dir}")
set(VCPKG_ROOT "${_vcpkg_dir}" CACHE PATH "" FORCE)
set(_vcpkg_toolchain "${_vcpkg_dir}/scripts/buildsystems/vcpkg.cmake")

#if (DEFINED ENV{BITLOOP_ROOT})
#  set(OVERLAY_PORTS_PATH "$ENV{BITLOOP_ROOT}/vcpkg-ports/ports")
#  list(APPEND VCPKG_OVERLAY_PORTS ${OVERLAY_PORTS_PATH})
#endif()


message(STATUS "Using vcpkg at: ${_vcpkg_dir}")
include("${_vcpkg_toolchain}")
