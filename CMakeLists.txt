#
# CMakeLists.txt
# MazeCubeGen: maze cube generator
#
# Copyright (c) 2020 Bryan Franklin. All rights reserved.
#
cmake_minimum_required(VERSION 2.9)
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

# add a target to remove files cmake creates
add_custom_target(full-clean
                  COMMAND bash -c 'rm -frv CMakeCache.txt CMakeFiles Makefile'
                  )

# add a debug build target
add_custom_target(debug-build
                  COMMAND bash -c 'cmake -DCMAKE_BUILD_TYPE=Debug . && make'
                  )

# add a valgrind testing targets
add_custom_target(valgrind1
                  COMMAND valgrind --leak-check=full --show-reachable=yes --vgdb-error=0 -v ./mcg --help || true
                  DEPENDS debug-build
                  )
add_custom_target(valgrind2
                  COMMAND valgrind --leak-check=full --show-reachable=yes --vgdb-error=0 -v ./mcg -d 9,11,13 -r 12345 -k 1 -o vg.txt -f vg_flat.stl -m vg.stl -s -p vg_sol.stl -g vg.gv
                  DEPENDS debug-build
                  )
add_custom_target(valgrind3
                  COMMAND valgrind --leak-check=full --show-reachable=yes --vgdb-error=0 -v ./mcg -i no_such_file.txt -o no_such_output.txt || true
                  DEPENDS debug-build
                  )
add_custom_target(valgrind4
                  COMMAND valgrind --leak-check=full --show-reachable=yes --vgdb-error=0 -v ./mcg -i vg.txt -s -o vg_out.txt
                  DEPENDS debug-build
                  )
add_custom_target(valgrind5
                  COMMAND valgrind --leak-check=full --show-reachable=yes --vgdb-error=0 -v ./mcg -o vg2.txt || true
                  DEPENDS debug-build
                  )
add_custom_target(valgrind
                  DEPENDS valgrind1 valgrind2 valgrind3 valgrind4 valgrind5
                  )

# add common system libraries
target_link_libraries(mcg m)
