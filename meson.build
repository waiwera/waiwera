# Waiwera Meson build script

project('waiwera', ['c', 'fortran'],
        version: '0.5.0',
        license: 'GPL3',
        meson_version: '>=0.43.0')

PETSc = dependency('PETSc', version: '>=3.10.1')
FSON  = dependency('FSON',  version: '>=1.0.0')

waiwera_sources = [
  'src/capillary_pressure.F90',
  'src/cell.F90',
  'src/dag.F90',
  'src/dictionary.F90',
  'src/dm_utils.F90',
  'src/eos.F90',
  'src/eos_setup.F90',
  'src/eos_wae.F90',
  'src/eos_wce.F90',
  'src/eos_we.F90',
  'src/eos_w.F90',
  'src/eos_wge.F90',
  'src/face.F90',
  'src/flow_simulation.F90',
  'src/fluid.F90',
  'src/fson_mpi.F90',
  'src/hdf5io.F90',
  'src/IAPWS.F90',
  'src/IFC67.F90',
  'src/initial.F90',
  'src/interpolation.F90',
  'src/kinds.F90',
  'src/list.F90',
  'src/logfile.F90',
  'src/mesh.F90',
  'src/minc.F90',
  'src/mpi_utils.F90',
  'src/ncg_air_thermodynamics.F90',
  'src/ncg_co2_thermodynamics.F90',
  'src/ncg_thermodynamics.F90',
  'src/ode.F90',
  'src/powertable.F90',
  'src/profiling.F90',
  'src/relative_permeability.F90',
  'src/rock.F90',
  'src/root_finder.F90',
  'src/source_control.F90',
  'src/source.F90',
  'src/source_setup.F90',
  'src/thermodynamics.F90',
  'src/thermodynamics_setup.F90',
  'src/timestepper.F90',
  'src/utils.F90',
  'src/version.F90',
  'src/zone.F90',
  'src/zone_label.F90']

waiwera_executable_sources = waiwera_sources + ['src/waiwera.F90']

waiwera_deps = [PETSc, FSON]
waiwera_rpath = ''
foreach dep: waiwera_deps
  waiwera_rpath += dep.get_pkgconfig_variable('libdir') + ':'
endforeach

executable('waiwera',
           waiwera_executable_sources,
           dependencies: waiwera_deps,
           install_rpath: waiwera_rpath,
           install: true)
