cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

# Cmake versions after 3.12 complain the ROOT path policy is not set.
# The following sets the policy to respect ROOT package paths.
if (POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

message(STATUS "finding UPCXX for SimCov build")
#find_package(UPCXX 2019.9.0.1 REQUIRED)
find_package(UPCXX REQUIRED)
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.13 AND DEFINED UPCXX_LIBRARIES)
  message(STATUS "UPCXX_UTILS is using the UPCXX::upcxx library interface: ${UPCXX_LIBRARIES}")
  find_package(UPCXX)
else()
  find_program(UPCXX_EXEC upcxx)
  set(CMAKE_CXX_COMPILER ${UPCXX_EXEC})
  message(STATUS "UPCXX_UTILS is using upcxx directly as the UPCXX::upcxx library interface is not available in this low version of cmake: ${CMAKE_VERSION}"    )
  unset(UPCXX_LIBRARIES)
endif()

message("Building ${CMAKE_BUILD_TYPE} version")

find_package(Threads REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder -Wno-sign-compare -Wno-unused-function -Wno-catch-value -pthread")

# a bit of jiggery-pokery to get a list of all compile time definitions in the code
get_directory_property(dir_defs COMPILE_DEFINITIONS)
set(all_defs "")
message("Definitions:")
foreach(d ${dir_defs})
    message("  " ${d})
    set(all_defs "${all_defs} ${d}")
endforeach()

add_definitions(-DALL_DEFNS=\"${all_defs}\")

add_executable(simcov main.cpp tissue.cpp utils.cpp)

target_link_libraries(simcov Threads::Threads ${UPCXX_LIBRARIES} ${UPCXX_UTILS_LIBRARIES} SIMCOV_VERSION_LIB)

#install(PROGRAMS simcov.py DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
install(TARGETS simcov DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
