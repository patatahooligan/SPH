cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)

project(sph)

find_package(VTK REQUIRED)

include_directories(${VTK_INCLUDE_DIRS})

add_subdirectory(
	${CMAKE_CURRENT_SOURCE_DIR}/thesis
	${CMAKE_CURRENT_SOURCE_DIR}/thesis/build
)
add_subdirectory(
	${CMAKE_CURRENT_SOURCE_DIR}/surface-from-vtp
	${CMAKE_CURRENT_SOURCE_DIR}/surface-from-vtp/build
)
