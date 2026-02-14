# tooling/findBitloop.cmake
include_guard(GLOBAL)

get_property(bl_in_try_compile GLOBAL PROPERTY IN_TRY_COMPILE)
if(bl_in_try_compile)
  return()
endif()

# tooling -> project root
get_filename_component(bl_project_root "${CMAKE_CURRENT_LIST_DIR}/.." REALPATH)

if(BITLOOP_DISCOVERY)
  include("${CMAKE_CURRENT_LIST_DIR}/toolchain/discovery_shims.cmake")
  return()
endif()

if(NOT TARGET bitloop::bitloop)
  if(EXISTS "${bl_project_root}/bitloop/CMakeLists.txt")
    add_subdirectory("${bl_project_root}/bitloop" "${CMAKE_BINARY_DIR}/_bitloop")
  else()
    find_package(bitloop CONFIG REQUIRED)
  endif()
endif()
