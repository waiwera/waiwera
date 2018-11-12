# Waiwera Meson build script

project('waiwera', ['c', 'fortran'],
        version: '0.5.0',
        license: 'GPL3',
        meson_version: '>=0.43.0')

PETSc = dependency('PETSc', version: '>=3.10.1')
fson  = dependency('fson',  version: '>=1.0.0', fallback: ['fson', 'fson_dep'])
zofu  = dependency('zofu',  version: '>=0.1.0')

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
waiwera_deps = [PETSc, fson]

if meson.get_compiler('fortran').get_id() == 'gcc'
  f_args = ['-ffree-line-length-none',
            '-Wno-unused-dummy-argument',
            '-Wno-unused-function',
            '-Wno-return-type',
            '-Wno-maybe-uninitialized']
else
  f_args = []
endif

if get_option('set_rpath')
  waiwera_rpath = join_paths(get_option('prefix'), get_option('libdir'))
  waiwera_rpath += ':' + PETSc.get_pkgconfig_variable('libdir')
else
  waiwera_rpath = ''
endif

executable('waiwera',
           waiwera_executable_sources,
           dependencies: waiwera_deps,
           fortran_args: f_args,
           install_rpath: waiwera_rpath,
           install: true)

waiwera = shared_library('waiwera',
           waiwera_executable_sources,
           dependencies: waiwera_deps,
           fortran_args: f_args,
           install: false)

test_sources = ['capillary_pressure_test',
                'IFC67_test',
                'IAPWS_test']
foreach test_name: test_sources
  test_src = join_paths(meson.current_source_dir(),
                        'test', 'unit', 'src', test_name + '.F90')
  driver_src_name = test_name + '_driver.F90'

  test_driver_src = configure_file(
    output: driver_src_name,
    command: ['zofu-driver', test_src, driver_src_name, '--mpi'])

  test_exe = executable(test_name,
                      [test_driver_src, test_src],
                      link_with: waiwera,
                      dependencies: [zofu, PETSc, fson])
  test(test_name, test_exe, is_parallel: false)
endforeach
