# vcpkg_autodetect.cmake
#    - Searches for an existing vcpkg in the current directory/workspace
#    - If one isn't found, it clones vcpkg at a pinned <sha>
#    - Also picks a location for shared binary sources/cache to avoid unnecessary rebuilds
#    - then chain-load that vcpkg toolchain

# Helper: Walks up directory level looking for a ".bitloop-workspace" file
function(find_workspace_root start_dir out_var)
  set(dir "${start_dir}")
  set(levels 0)
  set(limit 15)
  if (ARGC GREATER 3)
    set(limit "${ARGV3}")
  endif()

  # Normalize the start path (in case of relative paths / symlinks)
  get_filename_component(dir "${dir}" REALPATH)

  while (TRUE)
    if (EXISTS "${dir}/.bitloop-workspace")
      set(${out_var} "${dir}" PARENT_SCOPE)
      return()
    endif()

    # Stop if we've hit the filesystem root (parent == self)
    get_filename_component(parent "${dir}" DIRECTORY)
    if (parent STREQUAL dir)
      break()
    endif()

    # Increment level, exit if we exceed limit
    math(EXPR levels "${levels}+1")
    if (levels GREATER limit)
      break()
    endif()

    set(dir "${parent}")
  endwhile()

  set(${out_var} "" PARENT_SCOPE)
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
find_workspace_root(${CMAKE_SOURCE_DIR} WORKSPACE_DIR 20)

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
set(_installed_dir  "${_root_dir}/.vcpkg-installed")

file(MAKE_DIRECTORY "${_cache_dir}")
file(MAKE_DIRECTORY "${_installed_dir}")

# Shared binary cache / vcpkg-installed dir
set(ENV{VCPKG_INSTALLED_DIR}         "${_installed_dir}")
set(ENV{VCPKG_DEFAULT_BINARY_CACHE}  "${_cache_dir}")
set(ENV{VCPKG_BINARY_SOURCES}        "clear;files,${_cache_dir},readwrite")
set(ENV{VCPKG_DEFAULT_BINARY_CACHE}  "${_vcpkg_dir}/../.vcpkg-cache")

# Finally, load vcpkg toolchain
message(STATUS "Using vcpkg at: ${_vcpkg_dir}")
include("${_vcpkg_dir}/scripts/buildsystems/vcpkg.cmake")
