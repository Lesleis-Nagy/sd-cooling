cmake_minimum_required(VERSION 3.17)
project(sd_cooling)

set(CMAKE_CXX_STANDARD 20)

# Project directories
set(sd_cooling_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
set(sd_cooling_SRC_CMDLINE_DIR "${CMAKE_SOURCE_DIR}/src_cmdline")
set(sd_cooling_THIRD_PARTY_DIR "${CMAKE_SOURCE_DIR}/third_party")

set(sd_cooling_INCLUDE_BUILD_DIR "${CMAKE_BINARY_DIR}/include")
file(MAKE_DIRECTORY ${sd_cooling_INCLUDE_BUILD_DIR})

include_directories(${sd_cooling_INCLUDE_DIR})
include_directories(${sd_cooling_INCLUDE_BUILD_DIR})

#
# Set up boost
#
find_package(Boost COMPONENTS program_options REQUIRED)

#
# Set up eigen library
#
set(EIGEN_INCLUDE_DIRECTORIES "${sd_cooling_THIRD_PARTY_DIR}/eigen")
include_directories(${EIGEN_INCLUDE_DIRECTORIES})
message(STATUS "Eigen directories: ${EIGEN_INCLUDE_DIRECTORIES}")

#
# Set up csv library
#
set(CSV_INCLUDE_DIRECTORIES "${sd_cooling_THIRD_PARTY_DIR}/csv")
include_directories(${CSV_INCLUDE_DIRECTORIES})
message(STATUS "CSV directories: ${CSV_INCLUDE_DIRECTORIES}")

#
# Set up JSON library
#
set(JSON_INCLUDE_DIRECTORIES "${sd_cooling_THIRD_PARTY_DIR}/json/include")
include_directories(${JSON_INCLUDE_DIRECTORIES})
message(STATUS "JSON directories: ${JSON_INCLUDE_DIRECTORIES}")

#
# Set up RYAML library.
#
set(RYAML_INCLUDE_DIRECTORIES "${sd_cooling_THIRD_PARTY_DIR}/ryaml/include")
include_directories(${RYAML_INCLUDE_DIRECTORIES})
message(STATUS "RYAML directories: ${RYAML_INCLUDE_DIRECTORIES}")

#
# Build command line tools.
#
add_subdirectory(${sd_cooling_SRC_CMDLINE_DIR})

