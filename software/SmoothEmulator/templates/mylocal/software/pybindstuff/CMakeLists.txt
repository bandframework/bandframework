cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
set(CMAKE_VERBOSE_MAKEFILE on)
project(emulator_smooth CXX)
set(GITHOME_BAND_SMOOTH $ENV{GITHOME_BAND_SMOOTH})
set (MAIN_PROGRAM_DIR "../main_programs")
set (FAKEFULLMODEL_DIR "../fakefullmodels")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "$ENV{PWD}/../bin")
set (PYBIND11_FINDPYTHON ON)

add_custom_target(extern_smooth COMMAND make WORKING_DIRECTORY ${GITHOME_BAND_SMOOTH}/software)

message("Will Compile with " ${CMAKE_CXX_COMPILER})
message("For different compiler enter: cmake -D CMAKE_CXX_COMPILER=YOUR_PREFERRED_COMPILER .")
set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2 -fPIC")
set(CMAKE_MODULE_PATH "${MAIN_PROGRAM_DIR}/cmake")
find_package(Eigen3 REQUIRED)
message("eigen3_include_dir is ${EIGEN3_INCLUDE_DIR}")
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR} ${GITHOME_BAND_SMOOTH}/software/include $ENV{EIGEN3_INCLUDE_DIR})      
link_directories(${GITHOME_BAND_SMOOTH}/software/lib)

find_library(${GITHOME_BAND_SMOOTH}/software/lib msu_smooth)
find_library(${GITHOME_BAND_SMOOTH}/software/lib msu_smoothutils)
find_package(pybind11 REQUIRED)


add_library(emulator_smooth MODULE pybind_main.cc)
add_dependencies(emulator_smooth extern_smooth)


target_link_libraries(emulator_smooth PRIVATE pybind11::module pybind11::lto pybind11::windows_extras ${GITHOME_BAND_SMOOTH}/software/lib/libsmooth.a  ${GITHOME_BAND_SMOOTH}/software/lib/libmsu_smoothutils.a   ${EIGEN3_LIBRARY} )
pybind11_extension(emulator_smooth)
if(NOT MSVC AND NOT ${CMAKE_BUILD_TYPE} MATCHES Debug|RelWithDebInfo)
    # Strip unnecessary sections of the binary on Linux/macOS
    pybind11_strip(emulator_smooth)
endif()

set_target_properties(emulator_smooth PROPERTIES CXX_VISIBILITY_PRESET "hidden"
                                         CUDA_VISIBILITY_PRESET "hidden")