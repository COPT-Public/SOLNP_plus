project(SOLNP C)
cmake_minimum_required(VERSION 3.5.1)

if (CMAKE_HOST_WIN32)
    set(WINDOWS 1)
    message("Compile SOLNP+ on Windows")
elseif(CMAKE_HOST_APPLE)
    set(MACOS 1)
    message("Compile SOLNP+ on MacOS")
elseif(CMAKE_HOST_UNIX)
    message("Compile SOLNP+ on Linux")
    set(LINUX 1)
endif()
#set(SOLNP_COMPILER GCC)

#set(CMAKE_OSX_ARCHITECTURES arm64)
set(CMAKE_OSX_ARCHITECTURES x86_64)
set(CMAKE_C_STANDARD 11)

SET(CMAKE_BUILD_TYPE "Release")

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

find_package(Python3)

# Configure MKL and OSQP roots
set(ENV{MKL_ROOT} YOUR_MKL_ROOT_PATH) 

set(ENV{OSQP_ROOT} YOUR_OSQP_ROOT_PATH)

# Configure MKL paths
set(ENV{MKL_LIB_PATH} $ENV{MKL_ROOT}/lib/intel64)
set(ENV{MKL_include_PATH} $ENV{MKL_ROOT}/include)

# Configure osqp paths
set(ENV{OSQP_LIB_PATH} $ENV{OSQP_ROOT}/lib)
set(ENV{OSQP_include_PATH} $ENV{OSQP_ROOT}/include)

# Include project directory
include_directories(./include)
include_directories($ENV{OSQP_include_PATH})
include_directories($ENV{MKL_include_PATH})

# file(GLOB_RECURSE mysrc ./source/*.c)
file(GLOB_RECURSE mysrc 
                  ./source/der_info.c
                  ./source/linalg.c
                  ./source/linsys.c
                  ./source/qp_solver.c
                  ./source/Rescue.c
                  ./source/solnp_py.c
                  ./source/solnp.c
                  ./source/subnp_qp.c
                  ./source/util.c
                  )

file(GLOB_RECURSE myinclude ./include/*.h)

add_library(solnp SHARED ${mysrc} ${myinclude})

if(WINDOWS)
    target_link_libraries(solnp $ENV{MKL_LIB_PATH}/mkl_core.lib)
    target_link_libraries(solnp $ENV{MKL_LIB_PATH}/mkl_intel_lp64.lib)
    target_link_libraries(solnp $ENV{MKL_LIB_PATH}/mkl_sequential.lib)
else()
    target_link_libraries(solnp -Wl,--start-group
                                $ENV{MKL_LIB_PATH}/libmkl_core.a
                                $ENV{MKL_LIB_PATH}/libmkl_intel_lp64.a
                                $ENV{MKL_LIB_PATH}/libmkl_sequential.a
                                -Wl,--end-group)
endif()
target_link_libraries(solnp $ENV{OSQP_LIB_PATH}/libosqp.a)

