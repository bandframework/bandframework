cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
set(CMAKE_VERBOSE_MAKEFILE on)
project(emulator_smooth CXX)
set(SMOOTH_HOME $ENV{SMOOTH_HOME})
set (MAIN_PROGRAM_DIR "../main_programs")
set (FAKEFULLMODEL_DIR "../fakefullmodels")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "$ENV{PWD}/../bin")
set (PYBIND11_FINDPYTHON ON)
SET(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()

add_custom_target(extern_smooth COMMAND make WORKING_DIRECTORY ${SMOOTH_HOME}/software)

message("Will Compile with " ${CMAKE_CXX_COMPILER})
message("For different compiler enter: cmake -D CMAKE_CXX_COMPILER=YOUR_PREFERRED_COMPILER .")
set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++17 -O2 -fPIC -fvisibility=hidden")
set(CMAKE_MODULE_PATH "${MAIN_PROGRAM_DIR}/cmake")
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${SMOOTH_HOME}/software/include $ENV{EIGEN3_INCLUDE_DIR})
link_directories(${SMOOTH_HOME}/software/lib)

find_library(${SMOOTH_HOME}/software/lib msu_smooth)
find_library(${SMOOTH_HOME}/software/lib msu_smoothutils)
find_package(pybind11 REQUIRED)


add_library(emulator_smooth MODULE pybind_main.cc)
add_dependencies(emulator_smooth extern_smooth)

if(APPLE)
  set_target_properties(${MODULE} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(APPLE)

target_link_libraries(emulator_smooth PRIVATE pybind11::module pybind11::lto pybind11::windows_extras ${SMOOTH_HOME}/software/lib/libsmooth.a  ${SMOOTH_HOME}/software/lib/libmsu_smoothutils.a   ${EIGEN3_LIBRARY} )
pybind11_extension(emulator_smooth)
if(NOT MSVC AND NOT ${CMAKE_BUILD_TYPE} MATCHES Debug|RelWithDebInfo)
    # Strip unnecessary sections of the binary on Linux/macOS
    pybind11_strip(emulator_smooth)
endif()

set_target_properties(emulator_smooth PROPERTIES CXX_VISIBILITY_PRESET "hidden"
                                         CUDA_VISIBILITY_PRESET "hidden")
