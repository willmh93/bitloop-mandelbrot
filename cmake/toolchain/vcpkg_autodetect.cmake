# vcpkg_autodetect.cmake
#    - Searches for an existing vcpkg in the current directory/workspace
#    - If one isn't found, it clones vcpkg at a pinned <sha>
#    - Also picks a location for shared binary sources/cache to avoid unnecessary rebuilds
#    - then chain-load that vcpkg toolchain

# Helper: Walks up directory level looking for a ".bitloop-workspace" file
function(find_bitloop_workspace start_dir out_dir is_root_ws)
  set(limit 20)

  get_filename_component(dir "${start_dir}" REALPATH)
  set(found "")
  set(found_is_root OFF)

  while(limit GREATER 0)
    set(marker "${dir}/.bitloop-workspace")
    if(EXISTS "${marker}")
      # Read + normalize + lowercase for case-insensitive parse
      file(READ "${marker}" _mk)
      string(REPLACE "\r\n" "\n" _mk "${_mk}")
      string(REPLACE "\r"   "\n" _mk "${_mk}")
      string(TOLOWER "${_mk}" _mk_lower)

      # Match "root = <token>" anywhere (no ^/$ anchors -> search whole text)
      # CMake regex: no \s; use [ \t]
      string(REGEX MATCH "[ \t]*root[ \t]*=[ \t]*([^ \t\r\n]+)" _m "${_mk_lower}")

      set(is_root OFF)
      if(_m)
        # capture group (the value) is in CMAKE_MATCH_1
        set(_val "${CMAKE_MATCH_1}")
        if(_val STREQUAL "true" OR _val STREQUAL "1" OR _val STREQUAL "on" OR _val STREQUAL "yes")
          set(is_root ON)
        endif()
      endif()

      if(is_root)
        # Hard stop: this is the authoritative root
        set(${out_dir}    "${dir}" PARENT_SCOPE)
        set(${is_root_ws} ON       PARENT_SCOPE)
        return()
      endif()

      # Nearest non-hard marker: remember it, but keep walking in case a hard root appears above
      if(found STREQUAL "")
        set(found "${dir}")
        set(found_is_root OFF)
      endif()
    endif()

    # climb
    get_filename_component(parent "${dir}" DIRECTORY)
    if(parent STREQUAL dir)
      break()
    endif()
    set(dir "${parent}")
    math(EXPR limit "${limit}-1")
  endwhile()

  # Return nearest (if any), and whether it was a hard root (it wasn't if we got here)
  set(${out_dir}    "${found}"       PARENT_SCOPE)
  set(${is_root_ws} "${found_is_root}" PARENT_SCOPE)
endfunction()

# Helper: Clones vcpkg to 'dest_dir' at the pinned <sha>
function(clone_vcpkg dest_dir pinned_sha)
  get_filename_component(_dst "${dest_dir}" REALPATH)
  get_filename_component(_dst_parent "${_dst}" DIRECTORY)
  file(MAKE_DIRECTORY "${_dst_parent}")

  if(EXISTS "${_dst}/scripts/buildsystems/vcpkg.cmake")
    set(VCPKG_ROOT "${_dst}" PARENT_SCOPE)
    return()
  endif()

  find_package(Git QUIET)
  set(_ok OFF)
  if(GIT_FOUND)
    execute_process(COMMAND "${GIT_EXECUTABLE}" clone https://github.com/microsoft/vcpkg "${_dst}"
                    RESULT_VARIABLE rv)
    if(rv EQUAL 0)
      execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${_dst}" checkout ${pinned_sha}
                      RESULT_VARIABLE rv2)
      if(rv2 EQUAL 0)
        set(_ok ON)
      else()
        file(REMOVE_RECURSE "${_dst}")
      endif()
    endif()
  endif()

  if(NOT _ok)
    set(_zip "${_dst_parent}/vcpkg-${pinned_sha}.zip")
    file(DOWNLOAD "https://github.com/microsoft/vcpkg/archive/${pinned_sha}.zip" "${_zip}" TLS_VERIFY ON)
    file(ARCHIVE_EXTRACT INPUT "${_zip}" DESTINATION "${_dst_parent}")
    file(RENAME "${_dst_parent}/vcpkg-${pinned_sha}" "${_dst}")
  endif()

  if(NOT EXISTS "${_dst}/scripts/buildsystems/vcpkg.cmake")
    message(FATAL_ERROR "[vcpkg] Expected toolchain not found under ${_dst}")
  endif()

  set(VCPKG_ROOT "${_dst}" PARENT_SCOPE)
endfunction()



# Look for workspace root
find_bitloop_workspace(${CMAKE_SOURCE_DIR} WORKSPACE_DIR IS_ROOT_WORKSPACE)

set(VCPKG_PINNED_SHA "74e6536215718009aae747d86d84b78376bf9e09" CACHE STRING "Pinned vcpkg commit SHA")


# Potential vcpkg directory paths
set(LOCAL_VCPKG_DIR       "${CMAKE_SOURCE_DIR}/.vcpkg")
set(WORKSPACE_VCPKG_DIR   "${WORKSPACE_DIR}/.vcpkg")

# Potential vcpkg toolchain paths
set(LOCAL_VCPKG_PATH      "${LOCAL_VCPKG_DIR}/scripts/buildsystems/vcpkg.cmake")
set(WORKSPACE_VCPKG_PATH  "${WORKSPACE_VCPKG_DIR}/scripts/buildsystems/vcpkg.cmake")

# -------------- Set flags --------------
set(FOUND_WORKSPACE          FALSE)
set(FOUND_LOCAL_VCPKG        FALSE)
set(FOUND_WORKSPACE_VCPKG    FALSE)
set(FOUND_UNUSED_LOCAL_VCPKG FALSE)

# Found workspace?
if (EXISTS ${WORKSPACE_DIR})
  set(FOUND_WORKSPACE TRUE)
endif()

# Found local vcpkg?
if (EXISTS ${LOCAL_VCPKG_PATH})
  set(FOUND_LOCAL_VCPKG TRUE)
endif()

# Found workspace vcpkg?
if (FOUND_WORKSPACE AND EXISTS ${WORKSPACE_VCPKG_PATH})
  set(FOUND_WORKSPACE_VCPKG TRUE)
endif()

# Both workspace and local vcpkg found?
if (FOUND_WORKSPACE AND FOUND_LOCAL_VCPKG)
    # Are they the same? Treat as workspace, not local
    if (FOUND_WORKSPACE STREQUAL FOUND_LOCAL_VCPKG)
        # They're the same, disable local
        set(FOUND_LOCAL_VCPKG FALSE)
    else()
        # Different - workspace vcpkg takes priority
        set(FOUND_UNUSED_LOCAL_VCPKG TRUE)
    endif()
endif()

# ---------------------------------------


message(STATUS "")
message(STATUS "-------------- Searching for vcpkg --------------")
message(STATUS "WORKSPACE_DIR:          ${WORKSPACE_DIR}")
message(STATUS "IS_ROOT_WORKSPACE:      ${IS_ROOT_WORKSPACE}")
message(STATUS "FOUND_WORKSPACE:        ${FOUND_WORKSPACE}")
message(STATUS "FOUND_WORKSPACE_VCPKG:  ${FOUND_WORKSPACE_VCPKG}")
message(STATUS "FOUND_LOCAL_VCPKG:      ${FOUND_LOCAL_VCPKG}")
message(STATUS "-------------------------------------------------")
message(STATUS "")

if (FOUND_WORKSPACE)
    if (FOUND_WORKSPACE_VCPKG)
        message(STATUS "Found existing workspace vcpkg: ${WORKSPACE_VCPKG_DIR}")
    else()
        message(STATUS "Cloning vcpkg in to workspace")
        clone_vcpkg(${WORKSPACE_VCPKG_DIR} ${VCPKG_PINNED_SHA})
    endif()
    set(_vcpkg_dir ${WORKSPACE_VCPKG_DIR})

    if (FOUND_UNUSED_LOCAL_VCPKG)
        message(STATUS "Found local vcpkg which will be ignored as shared workspace vcpkg takes priority.")
        message(STATUS "Consider removing the local vcpkg/cache")
    endif()

elseif (FOUND_LOCAL_VCPKG)
    message(STATUS "Found existing local vcpkg: ${LOCAL_VCPKG_DIR}")

else()
    message(STATUS "No vcpkg found - cloning local vcpkg")
    clone_vcpkg(${LOCAL_VCPKG_DIR} ${VCPKG_PINNED_SHA})
    set(_vcpkg_dir ${LOCAL_VCPKG_DIR})

endif()


# Override VCPKG_ROOT with established vcpkg dir
set(ENV{VCPKG_ROOT} "${_vcpkg_dir}")
set(VCPKG_ROOT "${_vcpkg_dir}" CACHE PATH "" FORCE)

# Pick a root location relative to the determined VCPKG_ROOT
get_filename_component(_root_dir "${_vcpkg_dir}/.." REALPATH)

set(_cache_dir      "${_root_dir}/.vcpkg-cache")
#set(_installed_dir  "${_root_dir}/.vcpkg-installed")

file(MAKE_DIRECTORY "${_cache_dir}")
#file(MAKE_DIRECTORY "${_installed_dir}")

# Shared binary cache / vcpkg-installed dir
#set(ENV{VCPKG_INSTALLED_DIR}         "${_installed_dir}")
set(ENV{VCPKG_DEFAULT_BINARY_CACHE}  "${_cache_dir}")
set(ENV{VCPKG_BINARY_SOURCES}        "clear;files,${_cache_dir},readwrite")
set(ENV{VCPKG_DEFAULT_BINARY_CACHE}  "${_vcpkg_dir}/../.vcpkg-cache")

# Finally, load vcpkg toolchain
message(STATUS "Using vcpkg at: ${_vcpkg_dir}")
include("${_vcpkg_dir}/scripts/buildsystems/vcpkg.cmake")
