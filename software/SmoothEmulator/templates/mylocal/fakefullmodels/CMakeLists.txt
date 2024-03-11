cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
set(CMAKE_VERBOSE_MAKEFILE on)
project(smooth_fakemodels)
set(GITHOME_BAND_SMOOTH $ENV{GITHOME_BAND_SMOOTH})
set (MAIN_PROGRAM_DIR ".")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")

add_custom_target(extern_smooth COMMAND make WORKING_DIRECTORY ${GITHOME_BAND_SMOOTH}/software)

find_package(Eigen3 REQUIRED)
message("eigen3_include_dir is ${EIGEN3_INCLUDE_DIR}")

set (CXX $ENV{CXX})
message("CXX is ${CXX}")
set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2")
set(CMAKE_MODULE_PATH "${MAIN_PROGRAM_DIR}/cmake")
include_directories(SYSTEM  ${EIGEN3_INCLUDE_DIR} ${GITHOME_BAND_SMOOTH}/software/include $ENV{EIGEN3_INCLUDE_DIR})      

find_library(${GITHOME_BAND_SMOOTH}/software/lib msu_smooth)
find_library(${GITHOME_BAND_SMOOTH}/software/lib msu_smoothutils)

# fakerhic
#################################################
set(SOURCE_FILES fakerhic.cc)
add_executable(fakerhic fakerhic.cc)
add_dependencies(fakerhic extern_smooth)
target_link_libraries(fakerhic ${GITHOME_BAND_SMOOTH}/software/lib/libmsu_smoothutils.a
 ${EIGEN3_LIBRARY})
#################################################

