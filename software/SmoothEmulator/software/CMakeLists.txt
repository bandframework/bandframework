cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(msu_smooth LANGUAGES CXX)

message("Will Compile with " ${CMAKE_CXX_COMPILER})
message("For different compiler enter: cmake -D CMAKE_CXX_COMPILER=YOUR_PREFERRED_COMPILER .")

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++20 -O3")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Eigen3 REQUIRED)

message("eigen3_include_dir is $ENV{EIGEN3_INCLUDE_DIR}")

include_directories(SYSTEM include $ENV{EIGEN3_INCLUDE_DIR} ${PYBIND11_INCLUDE_DIR} ${Python_INCLUDE_DIR})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")

file(GLOB smooth_src src/smooth/*.cc)
add_library(smooth STATIC ${smooth_src})

file(GLOB smoothutils_src src/smoothutils/*/*.cc src/smoothutils/*/*/*.cc)
add_library(msu_smoothutils STATIC ${smoothutils_src})
