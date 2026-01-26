cmake_minimum_required(VERSION 3.21)
project(Mandelbrot LANGUAGES C CXX)
include(CheckCXXCompilerFlag)

if (NOT TARGET bitloop::bitloop)
  find_package(bitloop CONFIG REQUIRED)
endif()

# ----- configure project hierarchy (discovery) -----

# Gather source files for this project
file(GLOB SIM_SOURCES CONFIGURE_DEPENDS
	"Mandelbrot/external/*.cpp"
	"Mandelbrot/*.cpp"
	"Mandelbrot/*.hpp"
	"Mandelbrot/*.h"
)

# flags (applied only if project open as root)
#set(BITLOOP_DEV_MODE TRUE)

bitloop_new_project(Mandelbrot  ${SIM_SOURCES})
bitloop_add_dependency(Mandelbrot  "subprojects/Cardioid")
bitloop_hide(Cardioid)
bitloop_finalize()

# ----- configure project packages/compile/link options -----

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