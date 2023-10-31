project(SOLNP_demo C)
cmake_minimum_required(VERSION 3.5.1)

set(LINK_MKL 0)

if (CMAKE_HOST_WIN32)
    set(WINDOWS 1)
    set(CMAKE_OSX_ARCHITECTURES x86_64)
    message("Compile SOLNP+ on Windows")
elseif(CMAKE_HOST_APPLE)
    set(MACOS 1)
    set(CMAKE_OSX_ARCHITECTURES arm64)
    message("Compile SOLNP+ on MacOS")
elseif(CMAKE_HOST_UNIX)
    set(LINUX 1)
    set(CMAKE_OSX_ARCHITECTURES x86_64)
    message("Compile SOLNP+ on Linux")
endif()
#set(SOLNP_COMPILER GCC)

set(CMAKE_C_STANDARD 11)

SET(CMAKE_BUILD_TYPE "Release")

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)


# Include project directory
include_directories(../../include)

# file(GLOB_RECURSE mysrc ./source/*.c ./P78_example.c)
# file(GLOB_RECURSE mysrc ./P78_example.c)

file(GLOB_RECURSE myinclude ../../include/*.h)

add_executable(solnp_demo ./demo.c ${myinclude})
add_executable(solnp_p78 ./P78_example.c ${myinclude})
add_executable(solnp_p93 ./P93_example.c ${myinclude})
add_executable(solnp_rosenbrock ./rosenbrock_example.c ${myinclude})



if(WINDOWS)
    target_link_libraries(solnp_demo ${CMAKE_CURRENT_SOURCE_DIR}/../../build/libsolnp.dll)

elseif(MACOS)
    target_link_libraries(solnp_demo ${CMAKE_CURRENT_SOURCE_DIR}/../../build/libsolnp.dylib)
    target_link_libraries(solnp_p78 ${CMAKE_CURRENT_SOURCE_DIR}/../../build/libsolnp.dylib)
    target_link_libraries(solnp_p93 ${CMAKE_CURRENT_SOURCE_DIR}/../../build/libsolnp.dylib)
    target_link_libraries(solnp_rosenbrock ${CMAKE_CURRENT_SOURCE_DIR}/../../build/libsolnp.dylib)
else()
    target_link_libraries(solnp_demo ${CMAKE_CURRENT_SOURCE_DIR}/../../build/libsolnp.so)
endif()

