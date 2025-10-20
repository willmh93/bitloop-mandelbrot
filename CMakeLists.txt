cmake_minimum_required(VERSION 3.21)
project(Mandelbrot  LANGUAGES CXX)
include(CheckCXXCompilerFlag)

if (NOT TARGET bitloop::bitloop)
  find_package(bitloop CONFIG REQUIRED)
endif()

# Gather source files for this project
file(GLOB SIM_SOURCES CONFIGURE_DEPENDS
	"Mandelbrot/*.cpp"
	"Mandelbrot/*.h"
)

bitloop_new_project(Mandelbrot  ${SIM_SOURCES})

# MSVC
if (MSVC)
	target_compile_options(Mandelbrot PRIVATE /arch:AVX2 /fp:fast)
else()
	# GCC/Clang
	target_compile_options(Mandelbrot PRIVATE -mavx2 -mfma -ffp-contract=fast -ffast-math -fno-math-errno)
endif()


bitloop_add_dependency(Mandelbrot  "subprojects/Cardioid")

bitloop_hide(Cardioid)

message(STATUS "Actually where we should finalize: Mandelbrot")
bitloop_finalize()
