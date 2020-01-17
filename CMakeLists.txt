cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

project(sph)

add_subdirectory(
	${CMAKE_CURRENT_SOURCE_DIR}/thesis
	${CMAKE_CURRENT_SOURCE_DIR}/thesis/build
)
add_subdirectory(
	${CMAKE_CURRENT_SOURCE_DIR}/surface-from-vtp
	${CMAKE_CURRENT_SOURCE_DIR}/surface-from-vtp/build
)
