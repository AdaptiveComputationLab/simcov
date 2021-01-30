cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(LungModel)
#----------------------------------------------------------------------------------------#
#   project specification
#----------------------------------------------------------------------------------------#

# use PROJECT_SOURCE_DIR, not CMAKE_SOURCE_DIR
set(LUNGMODEL_SOURCE_DIR ${PROJECT_SOURCE_DIR})
set(LUNGMODEL_BINARY_DIR ${PROJECT_BINARY_DIR})

message("Building ${CMAKE_BUILD_TYPE} version")

# check compiler and set appropriate linker flags
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
include(CheckCXXSourceCompiles)
include(CheckCXXSourceRuns)

# all c++ must be at least c++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

if (NOT APPLE)
  add_definitions(-D_XOPEN_SOURCE=700)
endif()

include(CMakeDefinitions.txt)

include_directories(../include include)

add_subdirectory(src)
