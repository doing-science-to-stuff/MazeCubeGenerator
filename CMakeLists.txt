#
# CMakeLists.txt
# MazeCubeGen: maze cube generator
#
# Copyright (c) 2020 Bryan Franklin. All rights reserved.
#
cmake_minimum_required(VERSION 2.8)
project (mcg)

# The version number.
set (MCG_VERSION_MAJOR 0)
set (MCG_VERSION_MINOR 0)
set (MCG_VERSION_PATCH 1)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/mcg_config.h.in"
  "${PROJECT_BINARY_DIR}/mcg_config.h"
  )

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# add the binary tree to the search path for include files
# so that we will find mcg_config.h
include_directories("${PROJECT_BINARY_DIR}")

AUX_SOURCE_DIRECTORY(. SOURCE_FILES)
add_executable(mcg ${SOURCE_FILES})
set_property(TARGET mcg PROPERTY C_STANDARD 99)

# add a debug build target
add_custom_target(debug-build
                  COMMAND bash -c 'cmake -DCMAKE_BUILD_TYPE=Debug . && make'
                  )

# add a valgrind testing target
add_custom_target(valgrind
                  COMMAND valgrind --leak-check=full --show-reachable=yes --vgdb-error=0 -v ./mcg
                  DEPENDS debug-build
                  )

# add common system libraries
target_link_libraries(ndt m)
