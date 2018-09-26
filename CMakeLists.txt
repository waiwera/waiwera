cmake_minimum_required(VERSION 3.1...3.12)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project (waiwera VERSION 0.5.0
  LANGUAGES Fortran C)

# PETSc
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")
find_package(PETSc REQUIRED)
include_directories(${PETSC_INCLUDES})
add_definitions(${PETSC_DEFINITIONS})

# FSON
find_package(FSON REQUIRED)
include_directories(${FSON_MODULES})

file(GLOB SOURCES "src/*.F90")

add_executable(waiwera ${SOURCES})
target_link_libraries(waiwera ${PETSC_LIBRARIES} ${FSON_LIBRARY})
