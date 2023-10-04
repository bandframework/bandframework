cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(msu_commonutils)
set(CXX $ENV{CXX})
message("CXX is ${CXX}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fpie")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wmissing-declarations -Wpointer-arith -Wshadow -Wuninitialized -Winit-self -Wundef -Wcast-align -Wformat=2 -Wold-style-cast -Werror=switch -std=c++17 -mfpmath=sse")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall  -std=c++17")
message(STATUS "CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} include)
message("Includes: ${GSL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}")

file(GLOB msu_commonutils_src ./src/*/*.cc ./src/*/*/*.cc)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "lib")
add_library(msu_commonutils STATIC ${msu_commonutils_src})
target_link_libraries(msu_commonutils ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${Boost_LIBRARIES})


