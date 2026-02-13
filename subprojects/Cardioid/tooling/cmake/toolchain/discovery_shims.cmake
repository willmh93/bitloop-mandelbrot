# bitloop_discovery_shims.cmake

# _ChatGPT_ record root once (GLOBAL PROPERTY, not a normal variable)
get_property(bl_root_set GLOBAL PROPERTY BL_ROOT_PROJECT SET)
if(NOT bl_root_set)
  set_property(GLOBAL PROPERTY BL_ROOT_PROJECT "${CMAKE_SOURCE_DIR}")
endif()

macro(bitloop_new_project)
endmacro()

macro(bitloop_show)
endmacro()

macro(bitloop_hide)
endmacro()

macro(bitloop_add_dependency proj relpath)
  get_property(m GLOBAL PROPERTY BL_CHILD_MANIFESTS)
  if(NOT m)
    set(m "")
  endif()

  get_filename_component(child_dir "${relpath}" REALPATH BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  set(child_manifest "${child_dir}/vcpkg.json")

  if(EXISTS "${child_manifest}")
    list(APPEND m "${child_manifest}")
    list(REMOVE_DUPLICATES m)
    set_property(GLOBAL PROPERTY BL_CHILD_MANIFESTS "${m}")
  endif()

  # _ChatGPT_ descend so children can register their own children
  if(EXISTS "${child_dir}/CMakeLists.txt")
    string(MD5 child_hash "${child_dir}")
    add_subdirectory("${child_dir}" "${CMAKE_BINARY_DIR}/_bl_discovery/${child_hash}" EXCLUDE_FROM_ALL)
  endif()
endmacro()

macro(bitloop_finalize)
  get_property(root GLOBAL PROPERTY BL_ROOT_PROJECT)
  if(CMAKE_CURRENT_SOURCE_DIR STREQUAL root)
    if(DEFINED BITLOOP_DISCOVERY_OUT AND NOT BITLOOP_DISCOVERY_OUT STREQUAL "")
      get_property(m GLOBAL PROPERTY BL_CHILD_MANIFESTS)
      if(NOT m)
        set(m "")
      endif()

      file(WRITE  "${BITLOOP_DISCOVERY_OUT}" "set(BL_CHILD_MANIFESTS\n")
      foreach(p IN LISTS m)
        file(APPEND "${BITLOOP_DISCOVERY_OUT}" "  \"${p}\"\n")
      endforeach()
      file(APPEND "${BITLOOP_DISCOVERY_OUT}" ")\n")
    endif()
  endif()

  return()
endmacro()
