# tooling/merge_vcpkg_manifests.cmake
#
# Usage:
#   include("${_root_dir}/tooling/merge_vcpkg_manifests.cmake")
#   bl_merge_vcpkg_manifests(
#     ROOT_MANIFEST      "${_root_dir}/vcpkg.json"
#     CHILD_MANIFESTS    "${BL_CHILD_MANIFESTS}"   # list of .../vcpkg.json files
#     OUT_DIR            "${BL_MERGED_MANIFEST_DIR}"
#     [VERBOSE]
#   )
#

include_guard(GLOBAL)

# -------------------------
# JSON helpers
# -------------------------

function(_bl_read_text_no_bom out _path)
  file(READ "${_path}" _txt)
  string(SUBSTRING "${_txt}" 0 1 _c0)
  if(_c0 STREQUAL "?")
    string(SUBSTRING "${_txt}" 1 -1 _txt)
  endif()
  string(STRIP "${_txt}" _txt)
  set(${out} "${_txt}" PARENT_SCOPE)
endfunction()

function(_bl_json_get out json key)
  string(JSON _v ERROR_VARIABLE _err GET "${json}" "${key}")
  if(_err)
    set(${out} "" PARENT_SCOPE)
  else()
    set(${out} "${_v}" PARENT_SCOPE)
  endif()
endfunction()

function(_bl_json_has out_bool json key)
  string(JSON _t ERROR_VARIABLE _err TYPE "${json}" "${key}")
  if(_err)
    set(${out_bool} FALSE PARENT_SCOPE)
  else()
    set(${out_bool} TRUE PARENT_SCOPE)
  endif()
endfunction()

function(_bl_json_escape out s)
  # Minimal JSON string escaping
  set(_x "${s}")
  string(REPLACE "\\" "\\\\" _x "${_x}")
  string(REPLACE "\"" "\\\"" _x "${_x}")
  string(REPLACE "\n" "\\n" _x "${_x}")
  string(REPLACE "\r" "\\r" _x "${_x}")
  string(REPLACE "\t" "\\t" _x "${_x}")
  set(${out} "${_x}" PARENT_SCOPE)
endfunction()

# -------------------------
# vcpkg-configuration.json merging
# -------------------------

function(_bl_load_manifest_configuration out_config_json manifest_vcpkg_json_path)
  get_filename_component(_dir "${manifest_vcpkg_json_path}" DIRECTORY)
  set(_cfg_file "${_dir}/vcpkg-configuration.json")

  if(EXISTS "${_cfg_file}")
    _bl_read_text_no_bom(_cfg "${_cfg_file}")
    set(${out_config_json} "${_cfg}" PARENT_SCOPE)
    return()
  endif()

  _bl_read_text_no_bom(_manifest "${manifest_vcpkg_json_path}")

  _bl_json_has(_has_conf "${_manifest}" "configuration")
  if(_has_conf)
    _bl_json_get(_cfg "${_manifest}" "configuration")
    if(NOT "${_cfg}" STREQUAL "")
      set(${out_config_json} "${_cfg}" PARENT_SCOPE)
      return()
    endif()
  endif()

  _bl_json_has(_has_conf2 "${_manifest}" "vcpkg-configuration")
  if(_has_conf2)
    _bl_json_get(_cfg "${_manifest}" "vcpkg-configuration")
    if(NOT "${_cfg}" STREQUAL "")
      set(${out_config_json} "${_cfg}" PARENT_SCOPE)
      return()
    endif()
  endif()

  set(${out_config_json} "" PARENT_SCOPE)
endfunction()

function(_bl_parse_default_registry json out_kind out_repo out_baseline)
  _bl_json_get(_k "${json}" "kind")
  _bl_json_get(_r "${json}" "repository")
  _bl_json_get(_b "${json}" "baseline")
  set(${out_kind} "${_k}" PARENT_SCOPE)
  set(${out_repo} "${_r}" PARENT_SCOPE)
  set(${out_baseline} "${_b}" PARENT_SCOPE)
endfunction()

function(_bl_merge_registry_entry kind repo baseline packages_list reg_keys reg_kind reg_repo reg_base reg_packages)
  set(_keys  "${${reg_keys}}")
  set(_kinds "${${reg_kind}}")
  set(_repos "${${reg_repo}}")
  set(_bases "${${reg_base}}")
  set(_pkgs  "${${reg_packages}}")

  set(_key "${kind}|${repo}")
  list(FIND _keys "${_key}" _idx)

  # normalize incoming packages list
  set(_new_pkgs ${packages_list})
  list(REMOVE_DUPLICATES _new_pkgs)
  list(SORT _new_pkgs)
  list(JOIN _new_pkgs "|" _new_joined_pipe) # <-- IMPORTANT (not ';')

  if(_idx EQUAL -1)
    list(APPEND _keys  "${_key}")
    list(APPEND _kinds "${kind}")
    list(APPEND _repos "${repo}")
    list(APPEND _bases "${baseline}")
    list(APPEND _pkgs  "${_new_joined_pipe}") # <-- single element per registry

    set(${reg_keys}      "${_keys}"  PARENT_SCOPE)
    set(${reg_kind}      "${_kinds}" PARENT_SCOPE)
    set(${reg_repo}      "${_repos}" PARENT_SCOPE)
    set(${reg_base}      "${_bases}" PARENT_SCOPE)
    set(${reg_packages}  "${_pkgs}"  PARENT_SCOPE)
    return()
  endif()

  list(GET _bases ${_idx} _old_base)
  if(NOT "${_old_base}" STREQUAL "${baseline}")
    message(FATAL_ERROR
      "vcpkg registry baseline conflict for repository:\n"
      "  repo     : ${repo}\n"
      "  baseline : ${_old_base}\n"
      "  baseline : ${baseline}"
    )
  endif()

  # Merge old + new packages (split old pipe-string back into a list)
  list(GET _pkgs ${_idx} _old_joined_pipe)

  set(_old_list "")
  if(NOT "${_old_joined_pipe}" STREQUAL "")
    string(REPLACE "|" ";" _old_list "${_old_joined_pipe}")
  endif()

  set(_merged ${_old_list})
  foreach(p IN LISTS _new_pkgs)
    list(APPEND _merged "${p}")
  endforeach()

  list(REMOVE_DUPLICATES _merged)
  list(SORT _merged)
  list(JOIN _merged "|" _merged_pipe)

  list(REMOVE_AT _pkgs ${_idx})
  list(INSERT _pkgs ${_idx} "${_merged_pipe}")

  set(${reg_packages} "${_pkgs}" PARENT_SCOPE)
endfunction()

function(_bl_emit_vcpkg_configuration_json out_json)
  # Inputs expected in caller scope:
  #   _schema (optional)
  #   _default_kind/_default_repo/_default_base (optional)
  #   _reg_kind/_reg_repo/_reg_base/_reg_pkgs (parallel lists; _reg_pkgs entries are pipe-joined)

  set(_out "{")
  set(_need_comma FALSE)

  # $schema
  if(NOT "${_schema}" STREQUAL "")
    _bl_json_escape(_sch "${_schema}")
    string(APPEND _out "\"$schema\":\"${_sch}\"")
    set(_need_comma TRUE)
  endif()

  # default-registry
  if(NOT "${_default_kind}" STREQUAL "")
    if(_need_comma)
      string(APPEND _out ",")
    endif()

    _bl_json_escape(_k "${_default_kind}")
    _bl_json_escape(_b "${_default_base}")

    string(APPEND _out "\"default-registry\":{")
    string(APPEND _out "\"kind\":\"${_k}\",")

    if(NOT "${_default_repo}" STREQUAL "")
      _bl_json_escape(_r "${_default_repo}")
      string(APPEND _out "\"repository\":\"${_r}\",")
    endif()

    string(APPEND _out "\"baseline\":\"${_b}\"")
    string(APPEND _out "}")

    set(_need_comma TRUE)
  endif()

  # registries
  list(LENGTH _reg_kind _n_reg_kind)
  list(LENGTH _reg_repo _n_reg_repo)
  list(LENGTH _reg_base _n_reg_base)
  list(LENGTH _reg_pkgs _n_reg_pkgs)

  if(NOT _n_reg_kind EQUAL _n_reg_repo OR
     NOT _n_reg_kind EQUAL _n_reg_base OR
     NOT _n_reg_kind EQUAL _n_reg_pkgs)
    message(FATAL_ERROR
      "Internal error: registry arrays misaligned:\n"
      "  kind=${_n_reg_kind} repo=${_n_reg_repo} base=${_n_reg_base} pkgs=${_n_reg_pkgs}\n"
      "  _reg_pkgs='${_reg_pkgs}'"
    )
  endif()

  if(_n_reg_kind GREATER 0)
    if(_need_comma)
      string(APPEND _out ",")
    endif()

    string(APPEND _out "\"registries\":[")
    math(EXPR _last "${_n_reg_kind} - 1")

    foreach(i RANGE 0 ${_last})
      if(NOT i EQUAL 0)
        string(APPEND _out ",")
      endif()

      list(GET _reg_kind ${i} _k)
      list(GET _reg_repo ${i} _r)
      list(GET _reg_base ${i} _b)
      list(GET _reg_pkgs ${i} _pkgs_pipe)

      _bl_json_escape(_k_e "${_k}")
      _bl_json_escape(_r_e "${_r}")
      _bl_json_escape(_b_e "${_b}")

      string(APPEND _out "{")
      string(APPEND _out "\"kind\":\"${_k_e}\",")

      if(NOT "${_r}" STREQUAL "")
        string(APPEND _out "\"repository\":\"${_r_e}\",")
      endif()

      string(APPEND _out "\"baseline\":\"${_b_e}\",")

      # packages: split pipe-joined string back into a CMake list
      set(_pkgs_list "")
      if(NOT "${_pkgs_pipe}" STREQUAL "")
        string(REPLACE "|" ";" _pkgs_list "${_pkgs_pipe}")
      endif()

      string(APPEND _out "\"packages\":[")
      set(_pfirst TRUE)
      foreach(p IN LISTS _pkgs_list)
        if(NOT _pfirst)
          string(APPEND _out ",")
        endif()
        set(_pfirst FALSE)
        _bl_json_escape(_p_e "${p}")
        string(APPEND _out "\"${_p_e}\"")
      endforeach()
      string(APPEND _out "]")

      string(APPEND _out "}")
    endforeach()

    string(APPEND _out "]")
  endif()

  string(APPEND _out "}")

  # Validate generated JSON parses at all (parse error detection)
  string(JSON _dummy ERROR_VARIABLE _json_err GET "${_out}" "default-registry")
  if(_json_err AND NOT _json_err STREQUAL "NOTFOUND" AND NOT _json_err STREQUAL "NOT_FOUND")
    message(FATAL_ERROR
      "Generated vcpkg configuration JSON is invalid:\n"
      "${_json_err}\n"
      "JSON:\n${_out}"
    )
  endif()

  set(${out_json} "${_out}" PARENT_SCOPE)
endfunction()


function(bl_merge_vcpkg_configuration out_config_json root_vcpkg_json child_vcpkg_jsons)
  set(_schema "")
  set(_default_kind "")
  set(_default_repo "")
  set(_default_base "")

  set(_reg_keys "")
  set(_reg_kind "")
  set(_reg_repo "")
  set(_reg_base "")
  set(_reg_pkgs "")

  function(_incorporate_config cfg_json source_label)
    if("${cfg_json}" STREQUAL "")
      return()
    endif()

    _bl_json_get(_sch "${cfg_json}" "\$schema")
    if(NOT "${_sch}" STREQUAL "")
      if("${_schema}" STREQUAL "")
        set(_schema "${_sch}" PARENT_SCOPE)
      elseif(NOT "${_schema}" STREQUAL "${_sch}")
      endif()
    endif()

    _bl_json_has(_has_def "${cfg_json}" "default-registry")
    if(_has_def)
      _bl_json_get(_def "${cfg_json}" "default-registry")
      _bl_parse_default_registry("${_def}" _k _r _b)

      if(NOT "${_k}" STREQUAL "" OR NOT "${_b}" STREQUAL "")
        if("${_default_kind}" STREQUAL "")
          set(_default_kind "${_k}" PARENT_SCOPE)
          set(_default_repo "${_r}" PARENT_SCOPE)
          set(_default_base "${_b}" PARENT_SCOPE)
        else()
          if(NOT "${_default_kind}" STREQUAL "${_k}" OR
             NOT "${_default_repo}" STREQUAL "${_r}" OR
             NOT "${_default_base}" STREQUAL "${_b}")
            message(FATAL_ERROR
              "default-registry conflict across the project tree.\n"
              "Existing: kind=${_default_kind}, repo=${_default_repo}, baseline=${_default_base}\n"
              "New     : kind=${_k}, repo=${_r}, baseline=${_b}\n"
              "Source  : ${source_label}\n"
              "A single merged environment requires a single default-registry definition."
            )
          endif()
        endif()
      endif()
    endif()

    # Merge registries
    string(JSON _len ERROR_VARIABLE _err LENGTH "${cfg_json}" "registries")
    if(_err AND NOT _err STREQUAL "NOTFOUND" AND NOT _err STREQUAL "NOT_FOUND")
      message(FATAL_ERROR "Error reading registries[] from ${source_label}: ${_err}")
    endif()

    if(NOT _len STREQUAL "" AND NOT _len STREQUAL "NOTFOUND" AND NOT _len STREQUAL "NOT_FOUND" AND NOT _len EQUAL 0)
      math(EXPR _last "${_len} - 1")
      foreach(i RANGE 0 ${_last})
        string(JSON _entry ERROR_VARIABLE _e0 GET "${cfg_json}" "registries" ${i})
        if(_e0 AND NOT _e0 STREQUAL "NOTFOUND" AND NOT _e0 STREQUAL "NOT_FOUND")
          message(FATAL_ERROR "Error reading registries[${i}] from ${source_label}: ${_e0}")
        endif()

        _bl_json_get(_k "${_entry}" "kind")
        _bl_json_get(_repo "${_entry}" "repository")
        _bl_json_get(_base "${_entry}" "baseline")

        if("${_k}" STREQUAL "" OR "${_base}" STREQUAL "")
          message(FATAL_ERROR "Registry entry missing 'kind' or 'baseline' in ${source_label}")
        endif()

        # packages array -> list
        string(JSON _plen ERROR_VARIABLE _pe LENGTH "${_entry}" "packages")
        if(_pe AND NOT _pe STREQUAL "NOTFOUND" AND NOT _pe STREQUAL "NOT_FOUND")
          message(FATAL_ERROR "Error reading registry packages in ${source_label}: ${_pe}")
        endif()
        if(_plen STREQUAL "" OR _plen STREQUAL "NOTFOUND" OR _plen STREQUAL "NOT_FOUND")
          message(FATAL_ERROR "Registry entry missing 'packages' list in ${source_label} (required for non-default registries).")
        endif()

        set(_plist "")
        if(NOT _plen EQUAL 0)
          math(EXPR _plast "${_plen} - 1")
          foreach(j RANGE 0 ${_plast})
            string(JSON _p ERROR_VARIABLE _pe2 GET "${_entry}" "packages" ${j})
            if(_pe2 AND NOT _pe2 STREQUAL "NOTFOUND" AND NOT _pe2 STREQUAL "NOT_FOUND")
              message(FATAL_ERROR "Error reading registry packages[${j}] in ${source_label}: ${_pe2}")
            endif()
            list(APPEND _plist "${_p}")
          endforeach()
        endif()

        _bl_merge_registry_entry("${_k}" "${_repo}" "${_base}" "${_plist}"
          _reg_keys _reg_kind _reg_repo _reg_base _reg_pkgs
        )

        set(_reg_keys "${_reg_keys}" PARENT_SCOPE)
        set(_reg_kind "${_reg_kind}" PARENT_SCOPE)
        set(_reg_repo "${_reg_repo}" PARENT_SCOPE)
        set(_reg_base "${_reg_base}" PARENT_SCOPE)
        set(_reg_pkgs "${_reg_pkgs}" PARENT_SCOPE)
      endforeach()
    endif()
  endfunction()

  # Root config
  _bl_load_manifest_configuration(_root_cfg "${root_vcpkg_json}")
  _incorporate_config("${_root_cfg}" "root:${root_vcpkg_json}")

  # Child configs
  foreach(m IN LISTS child_vcpkg_jsons)
    _bl_load_manifest_configuration(_cfg "${m}")
    _incorporate_config("${_cfg}" "child:${m}")
  endforeach()

  # Emit merged configuration JSON object
  _bl_emit_vcpkg_configuration_json(_out)
  set(${out_config_json} "${_out}" PARENT_SCOPE)
endfunction()

# -------------------------
# Dependency map (names/frags/kinds) with conflict policy
# -------------------------

function(_bl_dep_put names_var frags_var kinds_var name frag kind)
  set(_names "${${names_var}}")
  set(_frags "${${frags_var}}")
  set(_kinds "${${kinds_var}}")

  list(FIND _names "${name}" _idx)
  if(_idx EQUAL -1)
    list(APPEND _names "${name}")
    list(APPEND _frags "${frag}")
    list(APPEND _kinds "${kind}")

    set(${names_var} "${_names}" PARENT_SCOPE)
    set(${frags_var} "${_frags}" PARENT_SCOPE)
    set(${kinds_var} "${_kinds}" PARENT_SCOPE)
    return()
  endif()

  list(GET _kinds ${_idx} _old_kind)
  list(GET _frags ${_idx} _old_frag)

  # OBJECT beats STRING
  if(_old_kind STREQUAL "STRING" AND kind STREQUAL "OBJECT")
    list(REMOVE_AT _frags ${_idx})
    list(INSERT _frags ${_idx} "${frag}")
    list(REMOVE_AT _kinds ${_idx})
    list(INSERT _kinds ${_idx} "OBJECT")

    set(${frags_var} "${_frags}" PARENT_SCOPE)
    set(${kinds_var} "${_kinds}" PARENT_SCOPE)
    return()
  endif()

  # If both are OBJECT but different, fail (safe default).
  if(_old_kind STREQUAL "OBJECT" AND kind STREQUAL "OBJECT" AND NOT _old_frag STREQUAL frag)
    message(FATAL_ERROR
      "vcpkg dependency object conflict for port '${name}'.\n"
      "Existing: ${_old_frag}\n"
      "New     : ${frag}\n"
      "Unify these across manifests or move the more-specific entry into the root manifest."
    )
  endif()

  # Otherwise keep existing (STRING vs STRING duplicates are ignored)
  set(${frags_var} "${_frags}" PARENT_SCOPE)
  set(${kinds_var} "${_kinds}" PARENT_SCOPE)
endfunction()

function(_bl_collect_deps manifest_path manifest_json names_var frags_var kinds_var)
  # Pull current lists from caller
  set(_names "${${names_var}}")
  set(_frags "${${frags_var}}")
  set(_kinds "${${kinds_var}}")

  # Query length of dependencies
  set(_err "")
  string(JSON _len ERROR_VARIABLE _err LENGTH "${manifest_json}" "dependencies")

  # Real error only if err is not the success sentinel
  if(NOT _err STREQUAL "" AND NOT _err STREQUAL "NOTFOUND" AND NOT _err STREQUAL "NOT_FOUND")
    message(FATAL_ERROR
      "Failed to read 'dependencies' in: ${manifest_path}\n"
      "CMake JSON error: ${_err}"
    )
  endif()

  # If dependencies missing, LENGTH generally can't be evaluated; _len may be empty/NOTFOUND.
  if(_len STREQUAL "" OR _len STREQUAL "NOTFOUND" OR _len STREQUAL "NOT_FOUND")
    set(${names_var} "${_names}" PARENT_SCOPE)
    set(${frags_var} "${_frags}" PARENT_SCOPE)
    set(${kinds_var} "${_kinds}" PARENT_SCOPE)
    return()
  endif()

  if(_len EQUAL 0)
    set(${names_var} "${_names}" PARENT_SCOPE)
    set(${frags_var} "${_frags}" PARENT_SCOPE)
    set(${kinds_var} "${_kinds}" PARENT_SCOPE)
    return()
  endif()

  math(EXPR _last "${_len} - 1")
  foreach(i RANGE 0 ${_last})
    set(_err "")
    string(JSON _type ERROR_VARIABLE _err TYPE "${manifest_json}" "dependencies" ${i})
    if(NOT _err STREQUAL "" AND NOT _err STREQUAL "NOTFOUND" AND NOT _err STREQUAL "NOT_FOUND")
      message(FATAL_ERROR "Failed to read dependencies[${i}] type in ${manifest_path}\n${_err}")
    endif()

    if(_type STREQUAL "STRING")
      set(_err "")
      string(JSON _name ERROR_VARIABLE _err GET "${manifest_json}" "dependencies" ${i})
      if(NOT _err STREQUAL "" AND NOT _err STREQUAL "NOTFOUND" AND NOT _err STREQUAL "NOT_FOUND")
        message(FATAL_ERROR "Failed to read dependencies[${i}] in ${manifest_path}\n${_err}")
      endif()

      # Store as JSON string fragment
      _bl_json_escape(_esc "${_name}")
      set(_frag "\"${_esc}\"")
      set(_kind "STRING")
      set(_dep_name "${_name}")

    elseif(_type STREQUAL "OBJECT")
      set(_err "")
      string(JSON _dep_name ERROR_VARIABLE _err GET "${manifest_json}" "dependencies" ${i} "name")
      if(NOT _err STREQUAL "" AND NOT _err STREQUAL "NOTFOUND" AND NOT _err STREQUAL "NOT_FOUND")
        message(FATAL_ERROR "Dependency object missing name at index ${i} in ${manifest_path}\n${_err}")
      endif()

      set(_err "")
      string(JSON _frag ERROR_VARIABLE _err GET "${manifest_json}" "dependencies" ${i})
      if(NOT _err STREQUAL "" AND NOT _err STREQUAL "NOTFOUND" AND NOT _err STREQUAL "NOT_FOUND")
        message(FATAL_ERROR "Failed to read dependency object at index ${i} in ${manifest_path}\n${_err}")
      endif()
      set(_kind "OBJECT")

    else()
      message(FATAL_ERROR "Unsupported dependency entry type '${_type}' in ${manifest_path}")
    endif()

    # Insert / merge policy
    list(FIND _names "${_dep_name}" _idx)
    if(_idx EQUAL -1)
      list(APPEND _names "${_dep_name}")
      list(APPEND _frags "${_frag}")
      list(APPEND _kinds "${_kind}")
    else()
      list(GET _kinds ${_idx} _old_kind)
      list(GET _frags ${_idx} _old_frag)

      # OBJECT beats STRING
      if(_old_kind STREQUAL "STRING" AND _kind STREQUAL "OBJECT")
        list(REMOVE_AT _frags ${_idx})
        list(INSERT _frags ${_idx} "${_frag}")
        list(REMOVE_AT _kinds ${_idx})
        list(INSERT _kinds ${_idx} "OBJECT")
      elseif(_old_kind STREQUAL "OBJECT" AND _kind STREQUAL "OBJECT" AND NOT _old_frag STREQUAL _frag)
        message(FATAL_ERROR
          "vcpkg dependency object conflict for '${_dep_name}'\n"
          "  Existing: ${_old_frag}\n"
          "  New     : ${_frag}\n"
          "  File    : ${manifest_path}"
        )
      endif()
    endif()
  endforeach()

  # Export to caller scope (this is the crucial part)
  set(${names_var} "${_names}" PARENT_SCOPE)
  set(${frags_var} "${_frags}" PARENT_SCOPE)
  set(${kinds_var} "${_kinds}" PARENT_SCOPE)
endfunction()


function(_bl_build_deps_array_json names frags out_json)
  # Produce a stable ordering by name without breaking alignment with frags:
  # build "name|idx" pairs, sort pairs, then output frags by idx.
  list(LENGTH names _n)
  if(_n EQUAL 0)
    set(${out_json} "[]" PARENT_SCOPE)
    return()
  endif()

  math(EXPR _last "${_n} - 1")
  set(_pairs "")
  foreach(i RANGE 0 ${_last})
    list(GET names ${i} _nm)
    list(APPEND _pairs "${_nm}|${i}")
  endforeach()
  list(SORT _pairs)

  set(_json "[")
  set(_first TRUE)
  foreach(p IN LISTS _pairs)
    string(REPLACE "|" ";" _tmp "${p}")
    list(GET _tmp 1 _idx)
    list(GET frags ${_idx} _frag)

    if(_first)
      set(_first FALSE)
    else()
      string(APPEND _json ",")
    endif()
    string(APPEND _json "${_frag}")
  endforeach()
  string(APPEND _json "]")

  set(${out_json} "${_json}" PARENT_SCOPE)
endfunction()


# -------------------------
# Main entry
# -------------------------

function(bl_merge_vcpkg_manifests)
  cmake_parse_arguments(BL
    "ALLOW_CHILD_OVERRIDES;VERBOSE"
    "ROOT_MANIFEST;OUT_DIR;MODE"
    "CHILD_MANIFESTS"
    ${ARGV}
  )

  if(NOT BL_ROOT_MANIFEST OR NOT EXISTS "${BL_ROOT_MANIFEST}")
    message(FATAL_ERROR "bl_merge_vcpkg_manifests: ROOT_MANIFEST missing/not found: ${BL_ROOT_MANIFEST}")
  endif()
  if(NOT BL_OUT_DIR)
    message(FATAL_ERROR "bl_merge_vcpkg_manifests: OUT_DIR is required")
  endif()

  file(MAKE_DIRECTORY "${BL_OUT_DIR}")
  set(_out_manifest "${BL_OUT_DIR}/vcpkg.json")

  # Normalize child list: unique + sorted
  set(_children ${BL_CHILD_MANIFESTS})
  list(REMOVE_DUPLICATES _children)
  list(SORT _children)

  #message(STATUS "ROOT_MANIFEST: ${BL_ROOT_MANIFEST}")

  #foreach(m IN LISTS _children)
  #  message(STATUS "Detected child manifest: ${m}")
  #endforeach()

  # Build merged configuration JSON (object string)
  bl_merge_vcpkg_configuration(_merged_cfg_json
      "${BL_ROOT_MANIFEST}"
      "${_children}"
  )

  # Load root JSON
  file(READ "${BL_ROOT_MANIFEST}" _root_json)

  _bl_json_get(_root_name      "${_root_json}" "name")
  _bl_json_get(_root_version   "${_root_json}" "version-string")
  _bl_json_get(_root_baseline  "${_root_json}" "builtin-baseline")
  _bl_json_get(_root_overrides "${_root_json}" "overrides")
  _bl_json_get(_root_features  "${_root_json}" "features")
  _bl_json_get(_root_deps_raw  "${_root_json}" "dependencies")
  if("${_root_deps_raw}" STREQUAL "")
    set(_root_deps_raw "[]")
  endif()

  # Validate/collect children
  set(_names "")
  set(_frags "")
  set(_kinds "")
  _bl_collect_deps("${BL_ROOT_MANIFEST}" "${_root_json}" _names _frags _kinds)

  foreach(m IN LISTS _children)
    file(READ "${m}" _child_json)

    _bl_json_get(_child_baseline "${_child_json}" "builtin-baseline")
    if(_root_baseline AND _child_baseline AND NOT _root_baseline STREQUAL _child_baseline)
      message(FATAL_ERROR
        "builtin-baseline mismatch:\n"
        "  root : ${_root_baseline}\n"
        "  child: ${_child_baseline}\n"
        "  file : ${m}\n"
        "Unify baselines to use a single merged install plan."
      )
    endif()

    _bl_json_has(_child_has_overrides "${_child_json}" "overrides")
    if(_child_has_overrides AND NOT BL_ALLOW_CHILD_OVERRIDES)
      message(FATAL_ERROR
        "Child manifest uses 'overrides' (not merged by default).\n"
        "  file: ${m}\n"
        "Move overrides to the root manifest, or pass ALLOW_CHILD_OVERRIDES to ignore this check."
      )
    endif()

    _bl_collect_deps("${m}" "${_child_json}" _names _frags _kinds)
  endforeach()

  _bl_build_deps_array_json("${_names}" "${_frags}" _deps_union_json)

  # Emit merged manifest JSON (minimal fields preserved)
  set(_out "{")
  set(_first TRUE)

  if(_root_name)
    _bl_json_escape(_n "${_root_name}")
    if(NOT _first)
      string(APPEND _out ",") 
    else()
      set(_first FALSE) 
    endif()
    string(APPEND _out "\"name\":\"${_n}\"")
  endif()

  if(_root_version)
    _bl_json_escape(_v "${_root_version}")
    if(NOT _first) 
      string(APPEND _out ",") 
    else()
      set(_first FALSE)
    endif()
    string(APPEND _out "\"version-string\":\"${_v}\"")
  endif()

  if(_root_baseline)
    _bl_json_escape(_b "${_root_baseline}")
    if(NOT _first) 
      string(APPEND _out ",")
    else()
      set(_first FALSE) 
    endif()
    string(APPEND _out "\"builtin-baseline\":\"${_b}\"")
  endif()

  if(NOT "${_root_overrides}" STREQUAL "")
    if(NOT _first) 
      string(APPEND _out ",")
    else() 
      set(_first FALSE) 
    endif()
    string(APPEND _out "\"overrides\":${_root_overrides}")
  endif()

  if(NOT _first) 
      string(APPEND _out ",")
  else() 
      set(_first FALSE) 
  endif()
  string(APPEND _out "\"dependencies\":${_deps_union_json}")

  string(APPEND _out ",\"configuration\":${_merged_cfg_json}")

  string(APPEND _out "}")

  file(WRITE "${_out_manifest}" "${_out}\n")

  if(BL_VERBOSE)
    message(STATUS "[vcpkg-merge] wrote UNION manifest: ${_out_manifest}")
  endif()

endfunction()
