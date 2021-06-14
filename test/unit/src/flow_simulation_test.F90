module flow_simulation_test

  ! Tests for flow_simulation module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use fson
  use flow_simulation_module
  use unit_test_utils_module

  implicit none
  private

  character(len = 512) :: data_path

  public :: setup, teardown
  public :: test_flow_simulation_init, &
       test_flow_simulation_fluid_properties, test_flow_simulation_lhs
  public :: test_setup_gravity, test_setup_flux

contains

!------------------------------------------------------------------------

  subroutine setup()

    use profiling_module, only: init_profiling

    ! Locals:
    PetscErrorCode :: ierr
    PetscInt :: ios

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call init_profiling()

    call get_environment_variable('WAIWERA_TEST_DATA_PATH', &
         data_path, status = ios)
    if (ios /= 0) data_path = ''

  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown()

    PetscErrorCode :: ierr

    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

  subroutine flow_simulation_basic_test(test, sim, title, thermo, eos, &
       dim, dof)

    ! Tests basic flow simulation parameters.

    class(unit_test_type), intent(in out) :: test
    type(flow_simulation_type), intent(in) :: sim
    character(*), intent(in) :: title, thermo, eos
    PetscInt, intent(in) :: dim, dof
    ! Locals:
    PetscInt :: sim_dim, sim_dof
    Vec :: x
    PetscErrorCode :: ierr
    PetscMPIInt :: rank
    
    call DMGetDimension(sim%mesh%dm, sim_dim, ierr); CHKERRQ(ierr)
    call DMGetGlobalVector(sim%mesh%dm, x, ierr); CHKERRQ(ierr)
    call VecGetSize(x, sim_dof, ierr); CHKERRQ(ierr)
    call DMRestoreGlobalVector(sim%mesh%dm, x, ierr); CHKERRQ(ierr)

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       call test%assert(title, sim%title, "Flow simulation title")
       call test%assert(thermo, sim%thermo%name, "Flow simulation thermodynamics")
       call test%assert(eos, sim%eos%name, "Flow simulation EOS")
       call test%assert(dim, sim_dim, "Flow simulation mesh dimension")
       call test%assert(dof, sim_dof, "Flow simulation mesh dof")
    end if

  end subroutine flow_simulation_basic_test

!------------------------------------------------------------------------
! Unit test routines:
!------------------------------------------------------------------------

  subroutine test_flow_simulation_init(test)

    ! Test flow_simulation init() method
    ! This uses a simple problem with a 12-cell rectangular mesh and two rock types.

    use fson_mpi_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    character(:), allocatable :: path
    type(flow_simulation_type) :: sim
    type(fson_value), pointer :: json
    PetscErrorCode :: err

    path = trim(adjustl(data_path)) // "flow_simulation/init/"
    json => fson_parse_mpi(path // "test_init.json")

    call sim%init(json, err = err)
    call fson_destroy_mpi(json)

    call test%assert(0, err, "flow_simulation init() error")

    call flow_simulation_basic_test(test, sim, title = "Test flow simulation init", &
         thermo = "IAPWS-97", eos = "w", dim = 3, dof = 12)

    call vec_diff_test(test, sim%solution, "primary", path, sim%mesh%cell_index)

    call vec_diff_test(test, sim%rock, "rock", path, sim%mesh%cell_index)

    call sim%destroy()

  end subroutine test_flow_simulation_init

!------------------------------------------------------------------------

  subroutine test_flow_simulation_fluid_properties(test)
    ! Test flow_simulation fluid_properties() method

    use fson_mpi_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(flow_simulation_type) :: sim
    type(fson_value), pointer :: json
    character(:), allocatable :: path
    PetscReal :: time = 0._dp
    PetscErrorCode :: err

    path = trim(adjustl(data_path)) // "flow_simulation/fluid_properties/"
    json => fson_parse_mpi(trim(path) // "test_fluid_properties.json")

    call sim%init(json, err = err)
    call fson_destroy_mpi(json)

    call test%assert(0, err, "fluid_properties error")
    call sim%pre_solve(time, sim%solution, err = err)
    call vec_diff_test(test, sim%fluid, "fluid", path, sim%mesh%cell_index)
    
    call sim%destroy()

  end subroutine test_flow_simulation_fluid_properties

!------------------------------------------------------------------------

  subroutine test_flow_simulation_lhs(test)
    ! Test LHS function

    use fson_mpi_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(flow_simulation_type) :: sim
    type(fson_value), pointer :: json
    character(:), allocatable :: path
    PetscReal, parameter :: time = 0._dp
    PetscReal :: interval(2) = [time, time]
    Vec :: lhs
    PetscErrorCode :: ierr, err

    path = trim(adjustl(data_path)) // "flow_simulation/lhs/"
    json => fson_parse_mpi(trim(path) // "test_lhs.json")

    call sim%init(json, err = err)
    call fson_destroy_mpi(json)

    call test%assert(0, err, "lhs error")
    call DMGetGlobalVector(sim%mesh%dm, lhs, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(lhs, "lhs", ierr); CHKERRQ(ierr)

    call sim%pre_solve(time, sim%solution, err = err)
    call sim%lhs(time, interval, sim%solution, lhs, err)
    call vec_diff_test(test, lhs, "lhs", path, sim%mesh%cell_index)

    call DMRestoreGlobalVector(sim%mesh%dm, lhs, ierr); CHKERRQ(ierr)
    call sim%destroy()

  end subroutine test_flow_simulation_lhs

!------------------------------------------------------------------------

  subroutine test_setup_gravity(test)
    ! Test gravity setup

    use fson_mpi_module
    use IAPWS_module
    use eos_we_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/2D.msh"}, ' // &
         '"logfile": {"filename": "", "echo": false}}')
    call gravity_test(json, "2D default gravity", [0._dp, 0._dp, 0._dp])
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/2D.msh"}, ' // &
         '"gravity": null, ' // &
         '"logfile": {"filename": "", "echo": false}}')
    call gravity_test(json, "2D null gravity", [0._dp, 0._dp, 0._dp])
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/2D.msh"}, ' // &
         '"gravity": 9.81, "logfile": {"filename": "", "echo": false}}')
    call gravity_test(json, "2D scalar gravity", [0._dp, -9.81_dp, 0._dp])
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/2D.msh"}, ' // &
         '"gravity": [-9.8, 0.0], "logfile": {"filename": "", "echo": false}}')
    call gravity_test(json, "2D vector gravity", [-9.8_dp, 0._dp, 0._dp])
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/block3.exo"}, ' // &
         '"logfile": {"filename": "", "echo": false}}')
    call gravity_test(json, "3D default gravity", [0._dp, 0._dp, -9.8_dp])
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/block3.exo"}, ' // &
         '"gravity": 9.80665, "logfile": {"filename": "", "echo": false}}')
    call gravity_test(json, "3D scalar gravity", [0._dp, 0._dp, -9.80665_dp])
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/block3.exo"}, ' // &
         '"gravity": [0., 0., -9.81], "logfile": {"filename": "", "echo": false}}')
    call gravity_test(json, "3D vector gravity", [0._dp, 0._dp, -9.81_dp])
    call fson_destroy_mpi(json)

  contains

    subroutine gravity_test(json, title, expected_gravity)

      type(fson_value), pointer, intent(in out) :: json
      character(*), intent(in) :: title
      PetscReal, intent(in) :: expected_gravity(:)
      ! Locals:
      PetscMPIInt :: rank
      type(flow_simulation_type) :: sim
      PetscErrorCode :: err, ierr

      call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
      call sim%init(json = json, err = err)
      if (rank == 0) then
         call test%assert(expected_gravity, sim%gravity, title)
      end if
      call sim%destroy()

    end subroutine gravity_test

  end subroutine test_setup_gravity

!------------------------------------------------------------------------

  subroutine test_setup_flux(test)
    ! Test flux setup

    use fson_mpi_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/2D.msh"}, ' // &
         '"output": {"fields": {"flux": ["water"]}}, ' // &
         '"logfile": {"filename": "", "echo": false}}')
    call flux_test(json, "2D", 4, 172)
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/2D.msh"}, ' // &
         '"boundaries": [{"faces": {"cells": [0, 12, 24, 36, 48, 60, 72, 84], ' // &
         '  "normal": [-1, 0]}}], ' // &
         '"output": {"fields": {"flux": ["water"]}}, ' // &
         '"logfile": {"filename": "", "echo": false}}')
    call flux_test(json, "2D -x bdy", 4, 180)
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/2D.msh"}, ' // &
         '"boundaries": [{"faces": ' // &
         '  {"cells": [84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95], ' // &
         '  "normal": [0, 1]}}], ' // &
         '"output": {"fields": {"flux": ["water"]}}, ' // &
         '"logfile": {"filename": "", "echo": false}}')
    call flux_test(json, "2D +y bdy", 4, 184)
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/col10.exo"}, ' // &
         '"output": {"fields": {"flux": ["water"]}}, ' // &
         '"logfile": {"filename": "", "echo": false}}')
    call flux_test(json, "column", 4, 9)
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/col10.exo"}, ' // &
         '"boundaries": [{"faces": ' // &
         '  {"cells": [0], "normal": [0, 0, 1]}}], ' // &
         '"output": {"fields": {"flux": ["water"]}}, ' // &
         '"logfile": {"filename": "", "echo": false}}')
    call flux_test(json, "column bdy", 4, 10)
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/3D.exo"}, ' // &
         '"output": {"fields": {"flux": ["water"]}}, ' // &
         '"logfile": {"filename": "", "echo": false}}')
    call flux_test(json, "3D", 4, 75)
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/3D.exo"}, ' // &
         '"boundaries": [{"faces": ' // &
         '  {"cells": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], ' // &
         '  "normal": [0, 0, 1]}}], ' // &
         '"output": {"fields": {"flux": ["water"]}}, ' // &
         '"logfile": {"filename": "", "echo": false}}')
    call flux_test(json, "3D bdy", 4, 87)
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '  "filename": "' // trim(adjustl(data_path)) // 'mesh/3D.exo", ' // &
         '  "zones": {"all": {"-": null}}, ' // &
         '  "minc": {"rock": {"zones": "all", ' // &
         '                    "fracture": {"type": "rock"}, ' // &
         '                    "matrix": {"type": "rock"}}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1}}}}, ' // &
         '"rock": {"types": [{"name": "rock", "zones": "all"}]}, ' // &
         '"output": {"fields": {"flux": ["water"]}}, ' // &
         '"logfile": {"filename": "", "echo": false}}')
    call flux_test(json, "3D MINC", 4, 75 + 36)
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '  "filename": "' // trim(adjustl(data_path)) // 'mesh/3D.exo", ' // &
         '  "zones": {"all": {"-": null}, "top": {"z": [-125, 0]}}, ' // &
         '  "minc": {"rock": {"zones": "top", ' // &
         '                    "fracture": {"type": "rock"}, ' // &
         '                    "matrix": {"type": "rock"}}, ' // &
         '           "geometry": {"fracture": {"volume": 0.1}}}}, ' // &
         '"rock": {"types": [{"name": "rock", "zones": "all"}]}, ' // &
         '"output": {"fields": {"flux": ["water"]}}, ' // &
         '"logfile": {"filename": "", "echo": false}}')
    call flux_test(json, "3D partial MINC", 4, 75 + 24)
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/hybrid10.msh"}, ' // &
         '"output": {"fields": {"flux": ["water"]}}, ' // &
         '"logfile": {"filename": "", "echo": false}}')
    call flux_test(json, "hybrid", 4, 12)
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/hybrid10.msh"}, ' // &
         '"boundaries": [{"faces": {"cells": [6, 9], ' // &
         '  "normal": [0, 0, 1]}}], ' // &
         '"output": {"fields": {"flux": ["water"]}}, ' // &
         '"logfile": {"filename": "", "echo": false}}')
    call flux_test(json, "hybrid bdy", 4, 14)
    call fson_destroy_mpi(json)

  contains

    subroutine flux_test(json, title, expected_vector_blocksize, &
         expected_vector_size)

      type(fson_value), pointer, intent(in out) :: json
      character(*), intent(in) :: title
      PetscInt, intent(in) :: expected_vector_blocksize
      PetscInt, intent(in) :: expected_vector_size ! number of blocks
      ! Locals:
      PetscMPIInt :: rank
      type(flow_simulation_type) :: sim
      DM :: flux_dm
      Vec :: global_flux
      PetscInt :: bs, vec_size
      PetscErrorCode :: err, ierr

      call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

      call sim%init(json = json, err = err)
      if (rank == 0) then
         call test%assert(0, err, title // ': init error')
      end if

      call VecGetDM(sim%flux, flux_dm, ierr); CHKERRQ(ierr)
      call DMGetGlobalVector(flux_dm, global_flux, ierr); CHKERRQ(ierr)
      call DMLocalToGlobal(flux_dm, sim%flux, &
           INSERT_VALUES, global_flux, ierr); CHKERRQ(ierr)

      call VecGetBlockSize(global_flux, bs, ierr); CHKERRQ(ierr)
      call VecGetSize(global_flux, vec_size, ierr); CHKERRQ(ierr)
      call DMRestoreGlobalVector(flux_dm, global_flux, ierr); CHKERRQ(ierr)

      if (rank == 0) then
         call test%assert(expected_vector_blocksize, bs, &
              title // ': vector blocksize')
         call test%assert(expected_vector_size, vec_size / bs, &
              title // ': vector size')
      end if

      call sim%destroy()

    end subroutine flux_test

  end subroutine test_setup_flux

!------------------------------------------------------------------------

end module flow_simulation_test
