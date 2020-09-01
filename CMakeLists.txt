cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)

project(sph)

find_package(VTK REQUIRED)
find_package(TBB REQUIRED)
find_package(tinyxml2 REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(VTK_INCLUDE_DIRS})
link_libraries(
	${VTK_LIBRARIES}
	${TBB_LIBRARIES}
	${Boost_LIBRARIES}
	tinyxml2
)

add_subdirectory(
	${CMAKE_CURRENT_SOURCE_DIR}/thesis
	${CMAKE_CURRENT_SOURCE_DIR}/thesis/build
)
add_subdirectory(
	${CMAKE_CURRENT_SOURCE_DIR}/surface-from-vtp
	${CMAKE_CURRENT_SOURCE_DIR}/surface-from-vtp/build
)
