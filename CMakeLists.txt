# Waiwera CMake script

cmake_minimum_required(VERSION 3.1...3.12)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project (waiwera
  VERSION 0.5.0
  LANGUAGES Fortran C)

# external libraries:
set(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/cmake-modules)
find_package(PETSc REQUIRED)
find_package(FSON REQUIRED)

# compiler options:
set(COMPILER_OPTIONS_COMMON -fPIC -ffree-line-length-none)
set(COMPILER_OPTIONS_DEBUG -g -O0 -Wall -Wno-unused-dummy-argument -Wno-unused-function
  -Wno-return-type -Wno-maybe-uninitialized -fcheck=all)
set(COMPILER_OPTIONS_RELEASE -O3 -march=native -mtune=native)

set(SOURCES
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
  src/waiwera.F90
  src/zone.F90
  src/zone_label.F90
  )
add_executable(waiwera ${SOURCES})

target_compile_definitions(waiwera PRIVATE ${PETSC_DEFINITIONS})
target_include_directories(waiwera PRIVATE ${PETSC_INCLUDES} ${FSON_MODULES})
target_compile_options(waiwera PRIVATE
  ${COMPILER_OPTIONS_COMMON}
  $<$<CONFIG:Debug>:${COMPILER_OPTIONS_DEBUG}>
  $<$<CONFIG:Release>:${COMPILER_OPTIONS_RELEASE}>)
target_link_libraries(waiwera ${PETSC_LIBRARIES} ${FSON_LIBRARY})
set_target_properties(waiwera PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/dist/
)
