# Waiwera CMake script

cmake_minimum_required(VERSION 3.1...3.12)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

include(ExternalProject)
find_package(Git REQUIRED)

project (waiwera
  VERSION 0.5.0
  LANGUAGES Fortran C)

set(EXENAME waiwera)
set(TEST_EXENAME test_all)

# external libraries:
set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_HOME_DIRECTORY}/cmake-modules")
find_package(PETSc 3.10 REQUIRED)
find_package(FSON)
if (NOT ${FSON_LIBRARY})
  include(install_FSON)
endif()
find_package(FRUIT)

# compiler options:
set(COMPILER_OPTIONS_COMMON -fPIC -ffree-line-length-none)
set(COMPILER_OPTIONS_DEBUG -g -O0 -Wall -Wno-unused-dummy-argument -Wno-unused-function
  -Wno-return-type -Wno-maybe-uninitialized -fcheck=all)
set(COMPILER_OPTIONS_RELEASE -O3 -march=native -mtune=native)

set(WAIWERA_MAIN_SOURCE src/waiwera.F90)
set(WAIWERA_SOURCES
  src/capillary_pressure.F90
  src/cell.F90
  src/dag.F90
  src/dictionary.F90
  src/dm_utils.F90
  src/eos.F90
  src/eos_setup.F90
  src/eos_wae.F90
  src/eos_wce.F90
  src/eos_we.F90
  src/eos_w.F90
  src/eos_wge.F90
  src/face.F90
  src/flow_simulation.F90
  src/fluid.F90
  src/fson_mpi.F90
  src/hdf5io.F90
  src/IAPWS.F90
  src/IFC67.F90
  src/initial.F90
  src/interpolation.F90
  src/kinds.F90
  src/list.F90
  src/logfile.F90
  src/mesh.F90
  src/minc.F90
  src/mpi_utils.F90
  src/ncg_air_thermodynamics.F90
  src/ncg_co2_thermodynamics.F90
  src/ncg_thermodynamics.F90
  src/ode.F90
  src/powertable.F90
  src/profiling.F90
  src/relative_permeability.F90
  src/rock.F90
  src/root_finder.F90
  src/source_control.F90
  src/source.F90
  src/source_setup.F90
  src/thermodynamics.F90
  src/thermodynamics_setup.F90
  src/timestepper.F90
  src/utils.F90
  src/version.F90
  src/zone.F90
  src/zone_label.F90
  )
add_executable(${EXENAME} ${WAIWERA_MAIN_SOURCE} ${WAIWERA_SOURCES})

target_compile_definitions(${EXENAME} PRIVATE ${PETSC_DEFINITIONS})
target_include_directories(${EXENAME} PRIVATE
  ${PETSC_INCLUDES} ${FSON_MODULES})
target_compile_options(${EXENAME} PRIVATE
  ${COMPILER_OPTIONS_COMMON}
  $<$<CONFIG:Debug>:${COMPILER_OPTIONS_DEBUG}>
  $<$<CONFIG:Release>:${COMPILER_OPTIONS_RELEASE}>)
target_link_libraries(${EXENAME} ${PETSC_LIBRARIES} ${FSON_LIBRARY})
set_target_properties(${EXENAME} PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/dist/
)

# unit tests:

set(TEST_SOURCES
  test/unit/src/test_all.F90
  test/unit/src/capillary_pressure_test.F90
  test/unit/src/cell_test.F90
  test/unit/src/dag_test.F90
  test/unit/src/dictionary_test.F90
  test/unit/src/dm_utils_test.F90
  test/unit/src/eos_test.F90
  test/unit/src/eos_we_test.F90
  test/unit/src/eos_wge_test.F90
  test/unit/src/eos_w_test.F90
  test/unit/src/face_test.F90
  test/unit/src/flow_simulation_test.F90
  test/unit/src/fluid_test.F90
  test/unit/src/fson_mpi_test.F90
  test/unit/src/IAPWS_test.F90
  test/unit/src/IFC67_test.F90
  test/unit/src/initial_test.F90
  test/unit/src/interpolation_test.F90
  test/unit/src/list_test.F90
  test/unit/src/mesh_test.F90
  test/unit/src/minc_test.F90
  test/unit/src/ncg_air_thermodynamics_test.F90
  test/unit/src/ncg_co2_thermodynamics_test.F90
  test/unit/src/powertable_test.F90
  test/unit/src/relative_permeability_test.F90
  test/unit/src/rock_test.F90
  test/unit/src/root_finder_test.F90
  test/unit/src/setup_test.F90
  test/unit/src/source_control_test.F90
  test/unit/src/source_setup_test.F90
  test/unit/src/source_test.F90
  test/unit/src/timestepper_test.F90
  test/unit/src/utils_test.F90
  test/unit/src/zone_test.F90
  )
add_executable(${TEST_EXENAME} ${TEST_SOURCES} ${WAIWERA_SOURCES})

target_compile_definitions(${TEST_EXENAME} PRIVATE ${PETSC_DEFINITIONS})
target_include_directories(${TEST_EXENAME} PRIVATE
  ${PETSC_INCLUDES} ${FSON_MODULES} ${FRUIT_MODULES})
target_compile_options(${TEST_EXENAME} PRIVATE
  ${COMPILER_OPTIONS_COMMON}
  $<$<CONFIG:Debug>:${COMPILER_OPTIONS_DEBUG}>
  $<$<CONFIG:Release>:${COMPILER_OPTIONS_RELEASE}>)
target_link_libraries(${TEST_EXENAME}
  ${PETSC_LIBRARIES} ${FSON_LIBRARY} ${FRUIT_LIBRARY})
set_target_properties(${TEST_EXENAME} PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test/unit/dist/
)
