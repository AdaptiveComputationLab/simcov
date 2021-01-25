cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-sign-compare -Wno-unused-function -Wno-catch-value")

add_definitions(-DALL_DEFNS=\"${all_defs}\")

add_executable(lungmodel main.cpp)

install(TARGETS lungmodel DESTINATION ${CMAKE_INSTALL_PREFIX}/)
