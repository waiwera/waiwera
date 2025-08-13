module unit_test_utils_module
  !! Unit testing utilities

#include <petsc/finclude/petsc.h>

  use petsc
  use zofu
  use kinds_module

  implicit none
  private

  public :: transition_compare, fluid_compare
  public :: vec_write, vec_diff_test

contains

!------------------------------------------------------------------------

  subroutine transition_compare(test, expected_primary, expected_region, &
       expected_transition, expected_err, primary, fluid, transition, &
       err, message)

    use fluid_module

    ! Runs assertions to test EOS transition

    class(unit_test_type), intent(in out) :: test
    PetscReal, intent(in) :: expected_primary(:), primary(:)
    PetscInt, intent(in) :: expected_region, expected_err
    type(fluid_type), intent(in) :: fluid
    PetscBool, intent(in) :: expected_transition, transition
    PetscErrorCode, intent(in) :: err
    character(*), intent(in) :: message
    ! Locals:
    PetscReal, parameter :: tol = 1.e-6_dp

    call test%assert(expected_err, err, trim(message) // " error")
    if (err == 0) then
       call test%assert(expected_primary, primary, &
            trim(message) // " primary", tol = tol)
       call test%assert(expected_region, nint(fluid%region), &
            trim(message) // " region")
       call test%assert(expected_transition, transition, &
            trim(message) // " transition")
    end if

  end subroutine transition_compare

!------------------------------------------------------------------------

  subroutine fluid_compare(test, fluid1, fluid2, message, tol)

    use fluid_module

    ! Runs assertions to see if two fluids have the same properties

    class(unit_test_type), intent(in out) :: test
    type(fluid_type), intent(in) :: fluid1, fluid2
    character(*), intent(in) :: message
    PetscReal, intent(in), optional :: tol
    ! Locals:
    PetscInt :: p
    character(1) :: phase_str
    PetscReal :: effective_tol
    PetscReal, parameter :: default_tol = 1.e-10_dp

    if (present(tol)) then
       effective_tol = tol
    else
       effective_tol = default_tol
    end if

    call test%assert(fluid1%pressure, fluid2%pressure, &
         trim(message) // " pressure", tol = effective_tol)
    call test%assert(fluid1%temperature, fluid2%temperature, &
         trim(message) // " temperature", tol = effective_tol)
    call test%assert(fluid1%liquidlike_fraction, fluid2%liquidlike_fraction, &
         trim(message) // " liquid-like fraction", tol = effective_tol)
    call test%assert(fluid1%partial_pressure, fluid2%partial_pressure, &
         trim(message) // " partial pressure", tol = effective_tol)

    do p = 1, fluid1%num_phases
       associate (phase1 => fluid1%phase(p), phase2 => fluid2%phase(p))
         if ((phase1%saturation > 0._dp) .or. (phase2%saturation > 0._dp)) then
            write(phase_str, '(i1)') p
            call test%assert(phase1%data, phase2%data, &
                 trim(message) // " phase " // phase_str, tol = effective_tol)
         end if
       end associate
    end do

  end subroutine fluid_compare

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

  subroutine vec_diff_test(test, v, name, path, cell_index)

    ! Tests vec v against values from HDF5 file with specified base name,
    ! at the given path.

    use dm_utils_module, only: vec_reorder

    class(unit_test_type), intent(in out) :: test
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
       call test%assert(0._dp, diffnorm, &
            "Flow simulation " // trim(name) // " vector")
    end if
    call DMRestoreGlobalVector(dm, diff, ierr); CHKERRQ(ierr)
    call DMRestoreGlobalVector(dm, vread, ierr); CHKERRQ(ierr)

  end subroutine vec_diff_test

!------------------------------------------------------------------------

end module unit_test_utils_module
