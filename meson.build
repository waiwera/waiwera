# Waiwera Meson build script

project('waiwera', ['c', 'fortran'],
        version: '1.1.0',
        license: 'GPL3',
        meson_version: '>=0.46.0')

if not get_option('mpi_wrapper_compiler')
  mpi = dependency('mpi', language: 'fortran')
endif

PETSc = dependency('PETSc', version: '>=3.12.0')
fson  = dependency('fson',  version: '>=1.0.0', fallback: ['fson', 'fson_dep'])
zofu  = dependency('zofu',  version: '>=0.1.0', fallback: ['zofu', 'zofu_dep'])

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
  'src/fson_utils.F90',
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
if not get_option('mpi_wrapper_compiler')
  waiwera_deps += mpi
endif

if meson.get_compiler('fortran').get_id() == 'gcc'
  f_args = ['-ffree-line-length-none',
            '-Wno-unused-dummy-argument',
            '-Wno-unused-function',
            '-Wno-return-type',
            '-Wno-maybe-uninitialized']
  if not get_option('buildtype').startswith('release')
    f_args += ['-fcheck=all']
  endif
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

unit_test_deps = [zofu, PETSc]
if not get_option('mpi_wrapper_compiler')
  unit_test_deps += mpi
endif

unit_test_utils = shared_library('unit_test_utils',
                            [join_paths(meson.current_source_dir(),
                                        'test', 'unit', 'src', 'unit_test_utils.F90')],
                            link_with: waiwera,
                            dependencies: unit_test_deps)

test_sources = ['capillary_pressure_test',
                'cell_test',
                'dag_test',
                'dictionary_test',
                'dm_utils_test',
                'eos_test',
                'eos_we_test',
                'eos_wge_test',
                'eos_w_test',
                'face_test',
                'flow_simulation_test',
                'fluid_test',
                'fson_mpi_test',
                'IAPWS_test',
                'IFC67_test',
                'initial_test',
                'interpolation_test',
                'list_test',
                'mesh_test',
                'minc_test',
                'ncg_air_thermodynamics_test',
                'ncg_co2_thermodynamics_test',
                'powertable_test',
                'relative_permeability_test',
                'rock_test',
                'root_finder_test',
                'source_control_test',
                'source_setup_test',
                'source_test',
                'timestepper_test',
                'utils_test',
                'zone_test']

zofu_driver = find_program('zofu-driver', required: false)
if not zofu_driver.found()
  zofu_subproject = subproject('zofu')
  zofu_driver = zofu_subproject.get_variable('driver_exe')
endif

test_env = environment()
test_env.set('WAIWERA_TEST_DATA_PATH',
             join_paths(meson.current_source_dir(),
                        'test', 'unit', 'data', ''))

foreach test_name: test_sources

  test_src = join_paths(meson.current_source_dir(),
                        'test', 'unit', 'src', test_name + '.F90')

  test_lib = static_library(test_name + '_lib', test_src,
                            link_with: [waiwera, unit_test_utils],
                            dependencies: [waiwera_deps, zofu],
                            fortran_args: f_args)

  test_driver_src = custom_target(
    test_name + '_driver',
    input: test_src,
    output: test_name + '_driver.F90',
    command: [zofu_driver, '@INPUT@', '@OUTPUT@', '--mpi'])

  test_exe = executable(test_name,
                        [test_src, test_driver_src],
                        link_with: [waiwera, unit_test_utils, test_lib],
                        dependencies: [waiwera_deps, zofu],
                        fortran_args: f_args)

  test(test_name, test_exe, env: test_env, is_parallel: false)

endforeach
