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
  public :: test_setup_gravity

  type(flow_simulation_type) :: sim

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

  subroutine flow_simulation_basic_test(test, title, thermo, eos, dim, dof)

    ! Tests basic flow simulation parameters.

    class(unit_test_type), intent(in out) :: test
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
    type(fson_value), pointer :: json
    PetscErrorCode :: err

    path = trim(adjustl(data_path)) // "flow_simulation/init/"
    json => fson_parse_mpi(path // "test_init.json")

    call sim%init(json, err = err)
    call fson_destroy_mpi(json)

    call test%assert(0, err, "flow_simulation init() error")

    call flow_simulation_basic_test(test, title = "Test flow simulation init", &
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
    PetscMPIInt :: rank
    PetscViewer :: viewer
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    PetscErrorCode :: ierr, err

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    viewer = PETSC_NULL_VIEWER
    call thermo%init()

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/2D.msh"}}')
    call eos%init(json, thermo)
    call sim%mesh%init(json)
    call sim%setup_gravity(json)
    call sim%mesh%configure(eos, sim%gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)
    if (rank == 0) then
       call test%assert([0._dp, 0._dp, 0._dp], sim%gravity, "2D default gravity")
    end if
    call sim%mesh%destroy()

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/2D.msh"}, "gravity": null}')
    call sim%mesh%init(json)
    call sim%setup_gravity(json)
    call sim%mesh%configure(eos, sim%gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)
    if (rank == 0) then
       call test%assert([0._dp, 0._dp, 0._dp], sim%gravity, "2D null gravity")
    end if
    call sim%mesh%destroy()

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/2D.msh"}, ' // &
         '"gravity": 9.81}')
    call sim%mesh%init(json)
    call sim%setup_gravity(json)
    call sim%mesh%configure(eos, sim%gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)
    if (rank == 0) then
       call test%assert([0._dp, -9.81_dp, 0._dp], sim%gravity, "2D scalar gravity")
    end if
    call sim%mesh%destroy()

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/2D.msh"}, ' // &
         '"gravity": [-9.8, 0.0]}')
    call sim%mesh%init(json)
    call sim%setup_gravity(json)
    call sim%mesh%configure(eos, sim%gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)
    if (rank == 0) then
       call test%assert([-9.8_dp, 0._dp, 0._dp], sim%gravity, "2D vector gravity")
    end if
    call sim%mesh%destroy()

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/block3.exo"}}')
    call sim%mesh%init(json)
    call sim%setup_gravity(json)
    call sim%mesh%configure(eos, sim%gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)
    if (rank == 0) then
       call test%assert([0._dp, 0._dp, -9.8_dp], sim%gravity, "3D default gravity")
    end if
    call sim%mesh%destroy()

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/block3.exo"}, ' // &
         '"gravity": 9.80665}')
    call sim%mesh%init(json)
    call sim%setup_gravity(json)
    call sim%mesh%configure(eos, sim%gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)
    if (rank == 0) then
       call test%assert([0._dp, 0._dp, -9.80665_dp], sim%gravity, "3D scalar gravity")
    end if
    call sim%mesh%destroy()

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "' // trim(adjustl(data_path)) // 'mesh/block3.exo"}, ' // &
         '"gravity": [0., 0., -9.81]}')
    call sim%mesh%init(json)
    call sim%setup_gravity(json)
    call sim%mesh%configure(eos, sim%gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)
    if (rank == 0) then
       call test%assert([0._dp, 0._dp, -9.81_dp], sim%gravity, "3D vector gravity")
    end if
    call sim%mesh%destroy()

    call eos%destroy()
    call thermo%destroy()

  end subroutine test_setup_gravity

!------------------------------------------------------------------------

end module flow_simulation_test
