# tooling/cmake/toolchain/bitloop_discovery_toolchain.cmake
set(BITLOOP_INTERNAL_DISCOVERY_RUN 1 CACHE INTERNAL "")

# If you include vcpkg.cmake here for compiler integration, ensure it does not install.
set(VCPKG_MANIFEST_INSTALL OFF CACHE BOOL "" FORCE)
