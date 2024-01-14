cmake_minimum_required(VERSION 3.10)

# set the project name
project(mesh_generator)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set (SOURCES
    src/VolumeRemesher/makePolyhedralMesh.cpp
    src/VolumeRemesher/delaunay.cpp
    src/VolumeRemesher/embed.cpp
    src/VolumeRemesher/conforming_mesh.cpp
    src/VolumeRemesher/extended_predicates.cpp
    src/VolumeRemesher/BSP.cpp
    src/VolumeRemesher/inOutPartition.cpp
)

try_compile(HAS_SIMD_INSTRUCTIONS ${CMAKE_BINARY_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/test_simd.cpp")

# make a static library for use within external programs
add_library(${PROJECT_NAME}_lib
    ${SOURCES}
)


if(NOT ${HAS_SIMD_INSTRUCTIONS})
	MESSAGE(WARNING "No simd supported")
	target_compile_definitions(${PROJECT_NAME}_lib PUBLIC VR_DISABLE_SIMD)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME}_lib PRIVATE "/fp:strict")
	# use intrinsic functions
	target_compile_options(${PROJECT_NAME}_lib PRIVATE "/Oi")
	# reserve enough stack size
	# target_compile_options(${PROJECT_NAME}_lib PRIVATE "/STACK:8421376")
	set_target_properties(${PROJECT_NAME}_lib PROPERTIES LINK_FLAGS --STACK,8421376)

	# turn off annoying warnings
	target_compile_definitions(${PROJECT_NAME}_lib PRIVATE "_CRT_SECURE_NO_WARNINGS")
	target_compile_definitions(${PROJECT_NAME}_lib PRIVATE "_CRT_NONSTDC_NO_WARNINGS")
	# workaround for known bug in C++17 on MSVC
	target_compile_definitions(${PROJECT_NAME}_lib PRIVATE "_HAS_STD_BYTE=0")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	# set standard optimization level
	target_compile_options(${PROJECT_NAME}_lib PRIVATE -O2)
	# grant IEEE 754 compliance
	target_compile_options(${PROJECT_NAME}_lib PRIVATE -frounding-math)
	# reserve enough stack size
	target_compile_options(${PROJECT_NAME}_lib PRIVATE -Wl,-z,stacksize=8421376)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	#Disable all optimizations
	target_compile_options(${PROJECT_NAME}_lib PRIVATE -O0)
	# reserve enough stack size
	target_compile_options(${PROJECT_NAME}_lib PRIVATE -Wl,-z,stacksize=8421376)
endif()


include(gmp)


# Public include directory
target_include_directories(${PROJECT_NAME}_lib PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(${PROJECT_NAME}_lib PUBLIC gmp::gmp)
target_compile_features(${PROJECT_NAME}_lib PUBLIC cxx_std_17)


add_library(VolumeRemesher::VolumeRemesher ALIAS ${PROJECT_NAME}_lib)


# add the executable
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(NOT HAS_PARENT)
	add_executable(${PROJECT_NAME}
    	src/VolumeRemesher/main.cpp
	)
	target_link_libraries(${PROJECT_NAME} PUBLIC VolumeRemesher::VolumeRemesher)
	target_compile_definitions(${PROJECT_NAME} PRIVATE VOL_REM_MODELS_PATH=\"${CMAKE_SOURCE_DIR}/models\")
endif()

target_compile_definitions(${PROJECT_NAME}_lib PUBLIC USE_GNU_GMP_CLASSES)
