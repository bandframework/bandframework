cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
set(GITHOME_BAND_SMOOTH $ENV{GITHOME_BAND_SMOOTH})
set(GITHOME_MSU_PRATT $ENV{GITHOME_MSU_PRATT})
project(msu_smooth CXX)

message("Will Compile with " ${CMAKE_CXX_COMPILER})
message("For different compiler enter: cmake -D YOUR_PREFERRED_COMPILER .")

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Eigen3 REQUIRED)

message("eigen3_include_dir is ${EIGEN3_INCLUDE_DIR}")

include_directories(SYSTEM include ${GITHOME_BAND_SMOOTH}/software/include ${EIGEN3_INCLUDE_DIR} $ENV{EIGEN3_INCLUDE_DIR})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")

file(GLOB smooth_src src/smooth/*.cc)
add_library(smooth STATIC ${smooth_src})

file(GLOB smoothutils_src src/smoothutils/*/*.cc src/smoothutils/*/*/*.cc)
add_library(msu_smoothutils STATIC ${smoothutils_src})