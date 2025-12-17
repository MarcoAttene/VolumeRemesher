cmake_minimum_required(VERSION 3.10)

# Set the target architecture.
# All modern x86/x64 processors support AVX2.
# Older x86/x64 processors may support SSE2 but not AVX2.
# Very old x86/x64 processors, or non x86/x64 
# processors, do not support any of the two.
set(ENABLE_SSE2 True)
set(ENABLE_AVX2 True)

#####################


# set the project name
project(mesh_generator)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set (SOURCES
    src/makePolyhedralMesh.cpp
    src/delaunay.cpp
    src/embed.cpp
    src/conforming_mesh.cpp
    src/extended_predicates.cpp
    src/BSP.cpp
    src/inOutPartition.cpp
)

# add the executable
add_executable(${PROJECT_NAME} 
    src/main.cpp
    ${SOURCES}
)

# make a static library for use within external programs
add_library(${PROJECT_NAME}_lib 
    ${SOURCES}
)

set(ALL_TARGETS ${PROJECT_NAME} ${PROJECT_NAME}_lib)

foreach (TARGET ${ALL_TARGETS})
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	set(GMP_HOME_DIR "D:/SYNC_DATA/Sviluppo_Software/Third_Party/mpir-master/msvc/vs19")
	set(GMPC_HOME_DIR ${GMP_HOME_DIR}/lib_mpir_gc/x64/Release)
	set(GMPXX_HOME_DIR ${GMP_HOME_DIR}/lib_mpir_cxx/x64/Release)
	set(GMP_INCLUDE_DIR ${GMPXX_HOME_DIR})
	set(GMP_LIB_DIR ${GMPC_HOME_DIR} ${GMPXX_HOME_DIR})
	set(GMP_LIB_NAME "mpir")

# grant IEEE 754 compliance
	target_compile_options(${TARGET} PUBLIC "/fp:strict")
# use intrinsic functions
	target_compile_options(${TARGET} PUBLIC "/Oi")
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:AVX2")
	elseif(ENABLE_SSE2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:SSE2")
	endif()
# turn off annoying warnings
	target_compile_options(${TARGET} PUBLIC "/D _CRT_SECURE_NO_WARNINGS") 
	target_compile_options(${TARGET} PUBLIC "/D _CRT_NONSTDC_NO_WARNINGS") 
# workaround for known bug in C++17 on MSVC
	target_compile_options(${TARGET} PUBLIC "/D _HAS_STD_BYTE=0")
else()
	set(GMP_HOME_DIR "gmp/6.2.1_1")
	set(GMP_LIB_NAME "gmp")
	set(GMP_INCLUDE_DIR ${GMP_HOME_DIR}/include)
	set(GMP_LIB_DIR ${GMP_HOME_DIR}/lib)

# set standard optimization level
	target_compile_options(${TARGET} PUBLIC -O2)
# grant IEEE 754 compliance
	target_compile_options(${TARGET} PUBLIC -frounding-math)
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME} PUBLIC "-mavx2")
		target_compile_options(${PROJECT_NAME} PUBLIC "-mfma")
	elseif(ENABLE_SSE2)
		target_compile_options(${PROJECT_NAME} PUBLIC "-msse2")
	endif()
endif()

# Public include directory
target_include_directories(${TARGET} PUBLIC
	src
	include
	${GMP_INCLUDE_DIR}
)

target_link_directories(${TARGET} PUBLIC ${GMP_LIB_DIR})
target_link_libraries(${TARGET} PUBLIC ${GMP_LIB_NAME})

endforeach(TARGET)

