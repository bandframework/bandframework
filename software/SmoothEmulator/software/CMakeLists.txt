cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(GITHOME_BAND $ENV{GITHOME_BAND})
project(msu_smooth CXX C)
message("CXX is ${CXX}")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -std=c++17")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)

message("eigen3_include_dir is ${EIGEN3_INCLUDE_DIR}")

include_directories(SYSTEM include ${GITHOME_BAND}/SmoothEmulator/software/include ${GITHOME_BAND}/SmoothEmulator/SmoothUtilities/software/include ${GSL_INCLUDE_DIR} $ENV{GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} $ENV{EIGEN3_INCLUDE_DIR} include)

message("${GITHOME_MSU}/SmoothEmulator/SmoothUtilities/software/include")

file(GLOB smooth_src src/*.cc)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")
add_library(smooth STATIC ${smooth_src})
