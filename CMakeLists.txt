cmake_minimum_required(VERSION 3.10)

# use gcc and g++ instead of clang
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
set(CMAKE_CXX_COMPILER "g++-10")
endif()

# set the project name
project(mesh_generator)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
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
# reserve enough stack size
	target_link_options(${TARGET} PUBLIC "/STACK:8421376")
# turn off annoying warnings
	target_compile_options(${TARGET} PUBLIC "/D _CRT_SECURE_NO_WARNINGS") 
	target_compile_options(${TARGET} PUBLIC "/D _CRT_NONSTDC_NO_WARNINGS") 
# workaround for known bug in C++17 on MSVC
	target_compile_options(${TARGET} PUBLIC "/D _HAS_STD_BYTE=0")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set(GMP_HOME_DIR "gmp/6.2.1_1")
	set(GMP_LIB_NAME "gmp")
	set(GMP_INCLUDE_DIR ${GMP_HOME_DIR}/include)
	set(GMP_LIB_DIR ${GMP_HOME_DIR}/lib)

# set standard optimization level
	target_compile_options(${TARGET} PUBLIC -O2)
# grant IEEE 754 compliance
	target_compile_options(${TARGET} PUBLIC -frounding-math)
# reserve enough stack size	
	target_compile_options(${TARGET} PUBLIC -Wl,-z,stacksize=8421376)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	set(GMP_HOME_DIR "gmp/6.2.1_1")
	set(GMP_LIB_NAME "gmp")
	set(GMP_INCLUDE_DIR ${GMP_HOME_DIR}/include)
	set(GMP_LIB_DIR ${GMP_HOME_DIR}/lib)

#Disable all optimizations
	target_compile_options(${TARGET} PUBLIC -O0)
# reserve enough stack size	
	target_compile_options(${TARGET} PUBLIC -Wl,-z,stacksize=8421376)
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

