cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(msu_smooth LANGUAGES CXX)

#find_package(Python3 3.9 EXACT)
#set(PYBIND11 LANGUAGES CXX)
#set(pybind11_DIR /Users/scottpratt/Library/Python/3.9/lib/python/site-packages/pybind11/share/cmake/pybind11)
#set(PYBIND11_INCLUDE_DIR /Users/scottpratt/Library/Python/3.9/lib/python/site-packages/pybind11/include)
#set(Python_INCLUDE_DIR /usr/local/Cellar/python@3.9/3.9.19_1/Frameworks/Python.framework/Versions/3.9/include/python3.9)
#find_package(pybind11 CONFIG REQUIRED)
#set(PYBIND11 FINDPYTHON ON)
#pybind11_add_module(msu_smooth include/msu_smooth/master.h)
#install(TARGETS msu_smooth DESTINATION .)

message("Will Compile with " ${CMAKE_CXX_COMPILER})
message("For different compiler enter: cmake -D YOUR_PREFERRED_COMPILER .")

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++20 -O2")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Eigen3 REQUIRED)

message("eigen3_include_dir is ${EIGEN3_INCLUDE_DIR}")

include_directories(SYSTEM include include ${EIGEN3_INCLUDE_DIR} $ENV{EIGEN3_INCLUDE_DIR} ${PYBIND11_INCLUDE_DIR} ${Python_INCLUDE_DIR})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")

file(GLOB smooth_src src/smooth/*.cc)
add_library(smooth STATIC ${smooth_src})

file(GLOB smoothutils_src src/smoothutils/*/*.cc src/smoothutils/*/*/*.cc)
add_library(msu_smoothutils STATIC ${smoothutils_src})

#set_target_properties(msu_smooth PROPERTIES LINKER_LANGUAGE CXX)