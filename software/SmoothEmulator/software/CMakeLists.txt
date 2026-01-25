cmake_minimum_required(VERSION 4.0 FATAL_ERROR)
project(msu_smooth LANGUAGES CXX)

message("Will Compile with " ${CMAKE_CXX_COMPILER})
message("For different compiler enter: cmake -D CMAKE_CXX_COMPILER=YOUR_PREFERRED_COMPILER .")

set(CMAKE_CXX_FLAGS "-W -Wall -march=native -std=c++20 -O3 -fPIC -Wno-unused-result")

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Eigen3 REQUIRED)

message("eigen3_include_dir is $ENV{EIGEN3_INCLUDE_DIR}")

include_directories(SYSTEM include $ENV{EIGEN3_INCLUDE_DIR} ${PYBIND11_INCLUDE_DIR} ${Python_INCLUDE_DIR})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")

file(GLOB smooth_src src/smooth/*.cc)
add_library(smooth STATIC ${smooth_src})

file(GLOB smoothutils_src src/smoothutils/*/*.cc src/smoothutils/*/*/*.cc)
add_library(smoothutils STATIC ${smoothutils_src})

##########################################

project(smooth_mainprograms CXX)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "../bin")
#################################################
link_directories(lib)
find_library(lib libsmoothutils.a)
find_library(lib libsmooth.a)

# trainingpoint_optimizer
#################################################
set(SOURCE_FILES MainPrograms/trainingpoint_optimizer_main.cc)
add_executable(trainingpoint_optimizer MainPrograms/trainingpoint_optimizer_main.cc)

get_target_property(TARGET_INCLUDES trainingpoint_optimizer INCLUDE_DIRECTORIES)

add_dependencies(trainingpoint_optimizer smooth)
add_dependencies(trainingpoint_optimizer smoothutils)
target_link_libraries(trainingpoint_optimizer libsmooth.a  libsmoothutils.a   ${EIGEN3_LIBRARY})

#################################################
# smoothy_testattrainingpts
#################################################
set(SOURCE_FILES MainPrograms/smoothy_testattrainingpts_main.cc)
add_executable(smoothy_testattrainingpts MainPrograms/smoothy_testattrainingpts_main.cc)
add_dependencies(smoothy_testattrainingpts smooth)
add_dependencies(smoothy_testattrainingpts smoothutils)
target_link_libraries(smoothy_testattrainingpts libsmooth.a  libsmoothutils.a   ${EIGEN3_LIBRARY})

#################################################
# smoothy_testvsfullmodel
#################################################
set(SOURCE_FILES MainPrograms/smoothy_testvsfullmodel_main.cc)
add_executable(smoothy_testvsfullmodel MainPrograms/smoothy_testvsfullmodel_main.cc)
add_dependencies(smoothy_testvsfullmodel smooth)
add_dependencies(smoothy_testvsfullmodel smoothutils)
target_link_libraries(smoothy_testvsfullmodel libsmooth.a  libsmoothutils.a   ${EIGEN3_LIBRARY})

#################################################
# smoothy_calcobs
#################################################
set(SOURCE_FILES MainPrograms/smoothy_calcobs_main.cc)
add_executable(smoothy_calcobs MainPrograms/smoothy_calcobs_main.cc)
add_dependencies(smoothy_calcobs smooth)
add_dependencies(smoothy_calcobs smoothutils)
target_link_libraries(smoothy_calcobs libsmooth.a libsmoothutils.a   ${EIGEN3_LIBRARY})

#################################################
# mcmc
#################################################
set(SOURCE_FILES MainPrograms/smoothy_mcmc_main.cc)
add_executable(smoothy_mcmc MainPrograms/smoothy_mcmc_main.cc)
add_dependencies(smoothy_mcmc smooth)
add_dependencies(smoothy_mcmc smoothutils)
target_link_libraries(smoothy_mcmc libsmooth.a libsmoothutils.a   ${EIGEN3_LIBRARY})

#################################################
# fake full model
#################################################
set(SOURCE_FILES TutorialPrograms/fakefullmodel.cc)
add_executable(fakefullmodel TutorialPrograms/fakefullmodel.cc)
add_dependencies(fakefullmodel smooth)
add_dependencies(fakefullmodel smoothutils)
target_link_libraries(fakefullmodel
libsmooth.a libsmoothutils.a ${EIGEN3_LIBRARY})

#################################################
# create fake info files
#################################################
set(SOURCE_FILES TutorialPrograms/fakeinfo.cc)
add_executable(fakeinfo TutorialPrograms/fakeinfo.cc)
add_dependencies(fakeinfo smooth)
add_dependencies(fakeinfo smooth smoothutils)
target_link_libraries(fakeinfo libsmooth.a libsmoothutils.a ${EIGEN3_LIBRARY})
#################################################


