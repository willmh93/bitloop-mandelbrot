# discovery_toolchain.cmake
#
# Toolchain used ONLY for Bitloop's project-discovery configure pass.
# This pass must be host-native and robust even when the outer build is a
# cross-compile (e.g. Emscripten). It should NOT inherit emcc/em++ variables.

# Mark discovery mode
set(BITLOOP_DISCOVERY ON CACHE BOOL "" FORCE)
set(BITLOOP_INTERNAL_DISCOVERY_RUN 1 CACHE BOOL "" FORCE)

# Make CMake feature/ABI checks less invasive.
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY CACHE STRING "" FORCE)

# Force host compilers so the discovery configure never tries to use em++.
# (The discovery pass should not compile/link target binaries.)
if (WIN32)
  # In VS CMake/Ninja environments, cl.exe is typically on PATH.
  set(CMAKE_C_COMPILER   cl CACHE FILEPATH "" FORCE)
  set(CMAKE_CXX_COMPILER cl CACHE FILEPATH "" FORCE)
endif()

# For non-Windows hosts, prefer the default toolchain on PATH.
if (NOT WIN32)
  if (NOT DEFINED CMAKE_C_COMPILER)
    set(CMAKE_C_COMPILER   cc CACHE FILEPATH "" FORCE)
  endif()
  if (NOT DEFINED CMAKE_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER c++ CACHE FILEPATH "" FORCE)
  endif()
endif()

# Keep discovery fast and avoid spurious failures.
set(CMAKE_C_COMPILER_WORKS   TRUE CACHE BOOL "" FORCE)
set(CMAKE_CXX_COMPILER_WORKS TRUE CACHE BOOL "" FORCE)

message(STATUS "discovery_toolchain.cmake <<<<<<<<<<<<<<<<<<<<<<<<")