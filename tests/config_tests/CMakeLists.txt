
#
# This file tests the cmake config files being exported/installed by rawtoaces.
# It simulates a project using the rawtoaces library as a dependency.
#
# DO NOT include it into a parent CMakeLists.txt.
#
# Usage: 
#    cmake -S . -B build -D rawtoaces_DIR=PATH_TO_INSTALL_DIR 

cmake_minimum_required(VERSION 3.10)
project( rawtoaces_config_test )

find_package(rawtoaces CONFIG REQUIRED )
find_package(OpenImageIO CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)

enable_testing()
add_subdirectory(core)
add_subdirectory(util)
