module flow_simulation_test

  ! Tests for flow_simulation module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use fruit
  use fson
  use flow_simulation_module

  implicit none
  private

public :: test_flow_simulation_init, &
     test_flow_simulation_fluid_properties, test_flow_simulation_lhs
public :: test_setup_gravity
public :: vec_write, vec_diff_test

PetscReal, parameter :: tol = 1.e-6_dp
type(flow_simulation_type) :: sim

contains

!------------------------------------------------------------------------
! Utility routines:
!------------------------------------------------------------------------

  subroutine vec_write(v, name, path, cell_index, field_indices, field_group)

    ! Writes vec v to HDF file with specified name and path (for
    ! generating reference values to test against).

    use hdf5io_module, only: vec_sequence_view_hdf5

    Vec, intent(in) :: v
    character(*), intent(in) :: name, path
    IS, intent(in) :: cell_index
    PetscInt, intent(in) :: field_indices(:)
    character(*), intent(in) :: field_group
    ! Locals:
    PetscErrorCode :: ierr
    PetscViewer :: viewer
    PetscInt :: time_index
    PetscReal :: time

    time_index = 0
    time = 0._dp

    call PetscViewerHDF5Open(PETSC_COMM_WORLD, &
         trim(path) // trim(name) // ".h5", &
         FILE_MODE_WRITE, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5PushGroup(viewer, "/", ierr); CHKERRQ(ierr)

    call ISView(cell_index, viewer, ierr); CHKERRQ(ierr)
    call vec_sequence_view_hdf5(v, field_indices, field_group, time_index, &
         time, viewer)

    call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

  end subroutine vec_write

!------------------------------------------------------------------------

  subroutine vec_diff_test(v, name, path, cell_index)

    ! Tests vec v against values from HDF5 file with specified base name,
    ! at the given path.

    use dm_utils_module, only: vec_reorder

    Vec, intent(in) :: v
    character(*), intent(in) :: name, path
    IS, intent(in) :: cell_index
    ! Locals:
    DM :: dm
    Vec :: vread, diff
    IS :: output_cell_index
    PetscViewer :: viewer
    PetscReal :: diffnorm
    PetscErrorCode :: ierr
    PetscMPIInt :: rank

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMGetGlobalVector(dm, vread, ierr); CHKERRQ(ierr)
    call DMGetGlobalVector(dm, diff, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(vread, name, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5Open(PETSC_COMM_WORLD, trim(path) // trim(name) // ".h5", &
         FILE_MODE_READ, viewer, ierr)
    CHKERRQ(ierr)
    call PetscViewerHDF5PushGroup(viewer, "/", ierr); CHKERRQ(ierr)
    call ISDuplicate(cell_index, output_cell_index, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(output_cell_index, "cell_index", ierr)
    CHKERRQ(ierr)
    call ISLoad(output_cell_index, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5PushGroup(viewer, "fields", ierr); CHKERRQ(ierr)
    call VecLoad(vread, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
    call vec_reorder(vread, output_cell_index, cell_index)
    call VecDuplicate(vread, diff, ierr); CHKERRQ(ierr)
    call VecCopy(vread, diff, ierr); CHKERRQ(ierr)
    call VecAXPY(diff, -1._dp, v, ierr); CHKERRQ(ierr)
    call VecNorm(diff, NORM_2, diffnorm, ierr); CHKERRQ(ierr)
    if (rank == 0) then
       call assert_equals(0._dp, diffnorm, tol, &
            "Flow simulation " // trim(name) // " vector")
    end if
    call DMRestoreGlobalVector(dm, diff, ierr); CHKERRQ(ierr)
    call DMRestoreGlobalVector(dm, vread, ierr); CHKERRQ(ierr)

  end subroutine vec_diff_test

!------------------------------------------------------------------------

  subroutine flow_simulation_basic_test(title, thermo, eos, dim, dof)

    ! Tests basic flow simulation parameters.

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
       call assert_equals(title, sim%title, "Flow simulation title")
       call assert_equals(thermo, sim%thermo%name, "Flow simulation thermodynamics")
       call assert_equals(eos, sim%eos%name, "Flow simulation EOS")
       call assert_equals(dim, sim_dim, "Flow simulation mesh dimension")
       call assert_equals(dof, sim_dof, "Flow simulation mesh dof")
    end if

  end subroutine flow_simulation_basic_test

!------------------------------------------------------------------------
! Unit test routines:
!------------------------------------------------------------------------

  subroutine test_flow_simulation_init

    ! Test flow_simulation init() method
    ! This uses a simple problem with a 12-cell rectangular mesh and two rock types.

    use fson_mpi_module

    ! Locals:
    character(26), parameter :: path = "data/flow_simulation/init/"
    type(fson_value), pointer :: json
    PetscErrorCode :: err

    json => fson_parse_mpi(trim(path) // "test_init.json")

    call sim%init(json, err = err)
    call fson_destroy_mpi(json)

    call assert_equals(0, err, "flow_simulation init() error")

    call flow_simulation_basic_test(title = "Test flow simulation init", &
         thermo = "IAPWS-97", eos = "w", dim = 3, dof = 12)

    call vec_diff_test(sim%solution, "primary", path, sim%mesh%cell_index)

    call vec_diff_test(sim%rock, "rock", path, sim%mesh%cell_index)

    call sim%destroy()

  end subroutine test_flow_simulation_init

!------------------------------------------------------------------------

  subroutine test_flow_simulation_fluid_properties
    ! Test flow_simulation fluid_properties() method

    use fson_mpi_module

    ! Locals:
    type(fson_value), pointer :: json
    character(64), parameter :: path = "data/flow_simulation/fluid_properties/"
    PetscReal :: time = 0._dp
    PetscErrorCode :: err

    json => fson_parse_mpi(trim(path) // "test_fluid_properties.json")

    call sim%init(json, err = err)
    call fson_destroy_mpi(json)

    call assert_equals(0, err, "fluid_properties error")
    call sim%pre_solve(time, sim%solution, err = err)
    call vec_diff_test(sim%fluid, "fluid", path, sim%mesh%cell_index)
    
    call sim%destroy()

  end subroutine test_flow_simulation_fluid_properties

!------------------------------------------------------------------------

  subroutine test_flow_simulation_lhs
    ! Test LHS function

    use fson_mpi_module

    ! Locals:
    type(fson_value), pointer :: json
    character(64), parameter :: path = "data/flow_simulation/lhs/"
    PetscReal, parameter :: time = 0._dp
    PetscReal :: interval(2) = [time, time]
    Vec :: lhs
    PetscErrorCode :: ierr, err

    json => fson_parse_mpi(trim(path) // "test_lhs.json")

    call sim%init(json, err = err)
    call fson_destroy_mpi(json)

    call assert_equals(0, err, "lhs error")
    call DMGetGlobalVector(sim%mesh%dm, lhs, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(lhs, "lhs", ierr); CHKERRQ(ierr)

    call sim%pre_solve(time, sim%solution, err = err)
    call sim%lhs(time, interval, sim%solution, lhs, err)
    call vec_diff_test(lhs, "lhs", path, sim%mesh%cell_index)

    call DMRestoreGlobalVector(sim%mesh%dm, lhs, ierr); CHKERRQ(ierr)
    call sim%destroy()

  end subroutine test_flow_simulation_lhs

!------------------------------------------------------------------------

  subroutine test_setup_gravity
    ! Test gravity setup

    use fson_mpi_module
    use IAPWS_module
    use eos_we_module

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
         '"filename": "data/mesh/2D.msh"}}')
    call eos%init(json, thermo)
    call sim%mesh%init(eos, json)
    call sim%setup_gravity(json)
    call sim%mesh%configure(sim%gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)
    call sim%mesh%destroy_distribution_data()
    if (rank == 0) then
       call assert_equals([0._dp, 0._dp, 0._dp], sim%gravity, 3, tol, &
            "2D default gravity")
    end if
    call sim%mesh%destroy()

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "data/mesh/2D.msh"}, "gravity": null}')
    call sim%mesh%init(eos, json)
    call sim%setup_gravity(json)
    call sim%mesh%configure(sim%gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)
    call sim%mesh%destroy_distribution_data()
    if (rank == 0) then
       call assert_equals([0._dp, 0._dp, 0._dp], sim%gravity, 3, tol, &
            "2D null gravity")
    end if
    call sim%mesh%destroy()

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "data/mesh/2D.msh"}, ' // &
         '"gravity": 9.81}')
    call sim%mesh%init(eos, json)
    call sim%setup_gravity(json)
    call sim%mesh%configure(sim%gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)
    call sim%mesh%destroy_distribution_data()
    if (rank == 0) then
       call assert_equals([0._dp, -9.81_dp, 0._dp], sim%gravity, 3, tol, &
         "2D scalar gravity")
    end if
    call sim%mesh%destroy()

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "data/mesh/2D.msh"}, ' // &
         '"gravity": [-9.8, 0.0]}')
    call sim%mesh%init(eos, json)
    call sim%setup_gravity(json)
    call sim%mesh%configure(sim%gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)
    call sim%mesh%destroy_distribution_data()
    if (rank == 0) then
       call assert_equals([-9.8_dp, 0._dp, 0._dp], sim%gravity, 3, tol, &
         "2D vector gravity")
    end if
    call sim%mesh%destroy()

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "data/mesh/block3.exo"}}')
    call sim%mesh%init(eos, json)
    call sim%setup_gravity(json)
    call sim%mesh%configure(sim%gravity, json, viewer = viewer, err = err)
    call sim%mesh%destroy_distribution_data()
    call fson_destroy_mpi(json)
    if (rank == 0) then
       call assert_equals([0._dp, 0._dp, -9.8_dp], sim%gravity, 3, tol, &
            "3D default gravity")
    end if
    call sim%mesh%destroy()

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "data/mesh/block3.exo"}, ' // &
         '"gravity": 9.80665}')
    call sim%mesh%init(eos, json)
    call sim%setup_gravity(json)
    call sim%mesh%configure(sim%gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)
    call sim%mesh%destroy_distribution_data()
    if (rank == 0) then
       call assert_equals([0._dp, 0._dp, -9.80665_dp], sim%gravity, 3, tol, &
            "3D scalar gravity")
    end if
    call sim%mesh%destroy()

    json => fson_parse_mpi(str = '{"mesh": {' // &
         '"filename": "data/mesh/block3.exo"}, ' // &
         '"gravity": [0., 0., -9.81]}')
    call sim%mesh%init(eos, json)
    call sim%setup_gravity(json)
    call sim%mesh%configure(sim%gravity, json, viewer = viewer, err = err)
    call fson_destroy_mpi(json)
    call sim%mesh%destroy_distribution_data()
    if (rank == 0) then
       call assert_equals([0._dp, 0._dp, -9.81_dp], sim%gravity, 3, tol, &
            "3D vector gravity")
    end if
    call sim%mesh%destroy()

    call eos%destroy()
    call thermo%destroy()

  end subroutine test_setup_gravity

!------------------------------------------------------------------------

end module flow_simulation_test
