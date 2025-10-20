cmake_minimum_required(VERSION 3.21)
project(Cardioid  LANGUAGES CXX)

# Link to bitloop
if (NOT TARGET bitloop::bitloop)
  find_package(bitloop CONFIG REQUIRED)
endif()

file(GLOB SIM_SOURCES CONFIGURE_DEPENDS
	"Cardioid/*.cpp"
	"Cardioid/*.h")

# 2. Create new project
bitloop_new_project(Cardioid  ${SIM_SOURCES})

bitloop_finalize()
