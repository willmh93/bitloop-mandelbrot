cmake_minimum_required(VERSION 3.21)
project(Mandelbrot LANGUAGES C CXX)
include(CheckCXXCompilerFlag)

if (NOT TARGET bitloop::bitloop)
  find_package(bitloop CONFIG REQUIRED)
endif()

# ---------------------------------  flags  ----------------------------------

# if this Cmake project open as root, BITLOOP_DEV_MODE overrides child projects,
# otherwise it inherits parent project value.

# set(BITLOOP_DEV_MODE TRUE)


# ------------------  configure project hierarchy/sources  -------------------

set(SRC_DIR "${CMAKE_CURRENT_LIST_DIR}/Mandelbrot")
file(GLOB_RECURSE SIM_SOURCES CONFIGURE_DEPENDS 
	"${SRC_DIR}/*.cpp" 
	"${SRC_DIR}/*.h" 
	"${SRC_DIR}/*.hpp"
)

bitloop_new_project(Mandelbrot ${SIM_SOURCES})
bitloop_add_dependency(Mandelbrot  "subprojects/Cardioid")
bitloop_hide(Cardioid)
bitloop_finalize()

# ---------  configure your target (compile options, packages, etc)  ---------

#find_package(zfp CONFIG REQUIRED)
#target_link_libraries(Mandelbrot PRIVATE zfp::zfp)

if (MSVC)
  target_compile_options(Mandelbrot PRIVATE /bigobj)
endif()

# fast math
if (NOT CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
	if (MSVC)
		target_compile_options(Mandelbrot PRIVATE /arch:AVX2 /fp:fast)
	else() # GCC/Clang
		target_compile_options(Mandelbrot PRIVATE -mavx2 -mfma -ffp-contract=fast -ffast-math -fno-math-errno)
	endif()
else()
	target_compile_options(Mandelbrot PRIVATE -O3 -flto -msimd128 -mrelaxed-simd -ffast-math -ffp-contract=fast -fno-math-errno)
	target_link_options(Mandelbrot PRIVATE -O3 -flto)
endif()

