project(SOLNP C)
cmake_minimum_required(VERSION 3.5.1)
set(LINK_MKL 0)
if (LINK_MKL)
    message("Linking MKL")
    set(BLA_VENDOR Intel10_64lp)
else()
    message("Linking LAPACK")
endif()

if (CMAKE_HOST_WIN32)
    set(WINDOWS 1)
    set(CMAKE_OSX_ARCHITECTURES x86_64)
    message("Compile SOLNP+ on Windows")
    # set(BLA_STATIC ON)
    find_package(LAPACK REQUIRED)
elseif(CMAKE_HOST_APPLE)
    set(MACOS 1)
    set(CMAKE_OSX_ARCHITECTURES arm64)
    message("Compile SOLNP+ on MacOS")
    find_library(ACCELERATE_FRAMEWORK Accelerate)
elseif(CMAKE_HOST_UNIX)
    set(LINUX 1)
    set(CMAKE_OSX_ARCHITECTURES x86_64)
    message("Compile SOLNP+ on Linux")
    # set(BLA_STATIC ON)
    find_package(LAPACK REQUIRED)
endif()
#set(SOLNP_COMPILER GCC)

set(CMAKE_C_STANDARD 11)

SET(CMAKE_BUILD_TYPE "Release")

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

find_package(Python3)

# # Configure MKL and OSQP roots
# set(ENV{MKL_ROOT} YOUR_MKL_PATH)
# set(ENV{LAPACK_ROOT} YOUR_LAPACK_PATH) 

# # Configure MKL paths
# set(ENV{MKL_LIB_PATH} $ENV{MKL_ROOT}/lib/intel64)
# set(ENV{MKL_include_PATH} $ENV{MKL_ROOT}/include)

# Configure osqp paths
include(FindOSQPConf.cmake)
if (OSQP_LIBRARY-NOTFOUND)
    message(FATAL_ERROR "OSQP not found")
endif()

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
                  ./source/solnp_util.c
                  ./source/solnp_c_interface.c
                  )

file(GLOB_RECURSE myinclude ./include/*.h)

add_library(solnp SHARED ${mysrc} ${myinclude})


target_link_libraries(solnp PRIVATE $ENV{OSQP_LIB_PATH}/libosqp.a)
if(MACOS)
    target_link_libraries(solnp PRIVATE "${ACCELERATE_FRAMEWORK}")
else()
    target_link_libraries(solnp PRIVATE $ENV{LAPACK_LIBRARIES})
endif()

# if(WINDOWS)
#     if(LINK_MKL)
#         target_link_libraries(solnp $ENV{MKL_LIB_PATH}/mkl_core.lib)
#         target_link_libraries(solnp $ENV{MKL_LIB_PATH}/mkl_intel_lp64.lib)
#         target_link_libraries(solnp $ENV{MKL_LIB_PATH}/mkl_sequential.lib)
#     else()
#         target_link_libraries(solnp $ENV{LAPACK_ROOT}/lib/liblapack.a)
#         target_link_libraries(solnp $ENV{LAPACK_ROOT}/lib/libblas.a)
#         target_link_libraries(solnp gfortran)
#         target_link_libraries(solnp quadmath)
#         # Add definition for name manging
#         add_definitions("-DADD_")
#     endif()
#     target_link_libraries(solnp $ENV{OSQP_LIB_PATH}/libosqp.a)

# elseif(MACOS)
#     target_link_libraries(solnp PRIVATE $ENV{OSQP_LIB_PATH}/libosqp.a)
#     if(LINK_MKL)
#         message( FATAL_ERROR "MKL not supported on MACOS" )
#     else()
#         target_link_libraries(solnp PRIVATE "${ACCELERATE_FRAMEWORK}")
#     endif()
# else()
#     if(LINK_MKL)
#         target_link_libraries(solnp -Wl,--start-group
#                                     $ENV{MKL_LIB_PATH}/libmkl_core.a
#                                     $ENV{MKL_LIB_PATH}/libmkl_intel_lp64.a
#                                     $ENV{MKL_LIB_PATH}/libmkl_sequential.a
#                                     -Wl,--end-group)
#     else()
#         # target_link_libraries(solnp $ENV{LAPACK_ROOT}/lib/liblapack.a)
#         target_link_libraries(solnp $ENV{LAPACK_LIBRARIES})
#     endif()
#     target_link_libraries(solnp $ENV{OSQP_LIB_PATH}/libosqp.a)
# endif()

