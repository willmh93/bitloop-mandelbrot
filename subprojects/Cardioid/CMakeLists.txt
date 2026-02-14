cmake_minimum_required(VERSION 3.21)
project(Cardioid LANGUAGES C CXX)
include(CheckCXXCompilerFlag)
include("${CMAKE_CURRENT_LIST_DIR}/tooling/findBitloop.cmake")

# ---------------------------------  flags  ----------------------------------
# if this Cmake project open as root, BITLOOP_DEV_MODE overrides child projects

# set(BITLOOP_DEV_MODE TRUE)

# ------------------  configure project hierarchy/sources  -------------------

file(GLOB SIM_SOURCES CONFIGURE_DEPENDS
	"Cardioid/*.cpp"
	"Cardioid/*.h")

# 2. Create new project
bitloop_new_project(Cardioid ${SIM_SOURCES})
bitloop_finalize()
