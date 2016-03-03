module flow_simulation_test

  ! Tests for flow_simulation module

  use kinds_module
  use mpi_module
  use fruit
  use fson
  use flow_simulation_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

public :: test_flow_simulation_init, &
     test_flow_simulation_fluid_properties, test_flow_simulation_lhs
public :: vec_write, vec_diff_test

PetscReal, parameter :: tol = 1.e-6_dp
type(flow_simulation_type) :: sim

contains

!------------------------------------------------------------------------
! Utility routines:
!------------------------------------------------------------------------

  subroutine get_rocktype_counts(rock_count)

    ! Returns an array on the output rank containing the total numbers
    ! (summed over all processors) of non-ghost cells with each rock
    ! type index.

    use rock_module, only : rocktype_label_name

    PetscInt, allocatable, intent(out) :: rock_count(:)
    ! Locals:
    PetscInt :: ir, ic, c, IS_size, ghost
    PetscInt :: num_rocktypes, total_num_rocktypes
    IS :: rock_IS
    DMLabel :: ghost_label
    PetscInt, allocatable :: nrc(:)
    PetscInt, pointer :: rock_cells(:)
    PetscErrorCode :: ierr
    
    call DMGetLabelSize(sim%mesh%dm, rocktype_label_name, &
         num_rocktypes, ierr); CHKERRQ(ierr)
    call MPI_allreduce(num_rocktypes, total_num_rocktypes, 1, MPI_INTEGER, &
         MPI_MAX, mpi%comm, ierr)
    call DMGetLabel(sim%mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    allocate(nrc(total_num_rocktypes), rock_count(total_num_rocktypes))
    do ir = 1, total_num_rocktypes
       nrc(ir) = 0 ! number of rocktype cells on processor
       call DMGetStratumIS(sim%mesh%dm, rocktype_label_name, ir, &
            rock_IS, ierr); CHKERRQ(ierr)
       if (rock_IS /= 0) then
          call ISGetIndicesF90(rock_IS, rock_cells, ierr); CHKERRQ(ierr)
          IS_size = size(rock_cells)
          do ic = 1, IS_size
             c = rock_cells(ic)
             call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
             if (ghost < 0) then
                nrc(ir) = nrc(ir) + 1
             end if
          end do
          call ISRestoreIndicesF90(rock_IS, rock_cells, ierr); CHKERRQ(ierr)
       end if
       call ISDestroy(rock_IS, ierr); CHKERRQ(ierr)
       call MPI_reduce(nrc(ir), rock_count(ir), 1, MPI_INTEGER, MPI_SUM, &
            mpi%output_rank, mpi%comm, ierr)
    end do
    deallocate(nrc)

  end subroutine get_rocktype_counts

!------------------------------------------------------------------------

  subroutine vec_write(v, name, path, cell_index)

    ! Writes vec v to HDF file with specified name and path (for
    ! generating reference values to test against).

    Vec, intent(in) :: v
    character(*), intent(in) :: name, path
    IS, intent(in) :: cell_index
    ! Locals:
    PetscErrorCode :: ierr
    PetscViewer :: viewer

    call PetscViewerHDF5Open(mpi%comm, &
         trim(path) // trim(name) // ".h5", &
         FILE_MODE_WRITE, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5PushGroup(viewer, "/", ierr); CHKERRQ(ierr)
    call VecView(v, viewer, ierr); CHKERRQ(ierr)
    call ISView(cell_index, viewer, ierr); CHKERRQ(ierr)
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

    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMGetGlobalVector(dm, vread, ierr); CHKERRQ(ierr)
    call DMGetGlobalVector(dm, diff, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(vread, name, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5Open(mpi%comm, trim(path) // trim(name) // ".h5", &
         FILE_MODE_READ, viewer, ierr)
    CHKERRQ(ierr)
    call PetscViewerHDF5PushGroup(viewer, "/", ierr); CHKERRQ(ierr)
    call ISDuplicate(cell_index, output_cell_index, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(output_cell_index, "cell_order", ierr)
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
    if (mpi%rank == mpi%output_rank) then
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
    
    call DMGetDimension(sim%mesh%dm, sim_dim, ierr); CHKERRQ(ierr)
    call DMGetGlobalVector(sim%mesh%dm, x, ierr); CHKERRQ(ierr)
    call VecGetSize(x, sim_dof, ierr); CHKERRQ(ierr)
    call DMRestoreGlobalVector(sim%mesh%dm, x, ierr); CHKERRQ(ierr)

    if (mpi%rank == mpi%output_rank) then
       call assert_equals(title, sim%title, "Flow simulation title")
       call assert_equals(thermo, sim%thermo%name, "Flow simulation thermodynamics")
       call assert_equals(eos, sim%eos%name, "Flow simulation EOS")
       call assert_equals(dim, sim_dim, "Flow simulation mesh dimension")
       call assert_equals(dof, sim_dof, "Flow simulation mesh dof")
    end if

  end subroutine flow_simulation_basic_test

!------------------------------------------------------------------------

  subroutine flow_simulation_label_test(rock_cells)

    ! Tests flow simulation labels.

    use boundary_module, only : open_boundary_label_name
    use rock_module, only : rocktype_label_name

    PetscInt, intent(in) :: rock_cells(:)
    ! Locals:
    PetscBool :: open_bdy, has_rock_label
    PetscInt, allocatable :: sim_rock_cells(:)
    PetscErrorCode :: ierr

    ! Open boundary label:
    call DMHasLabel(sim%mesh%dm, open_boundary_label_name, open_bdy, &
         ierr); CHKERRQ(ierr)
    if (mpi%rank == mpi%output_rank) then
       call assert_equals(.true., open_bdy, "Flow simulation open boundary label")
    end if

    ! Rock type label:
    call DMHasLabel(sim%mesh%dm, rocktype_label_name, has_rock_label, &
         ierr); CHKERRQ(ierr)
    if (mpi%rank == mpi%output_rank) then
       call assert_equals(.true., has_rock_label, "Flow simulation rocktype label")
    end if
    if (has_rock_label) then
       call get_rocktype_counts(sim_rock_cells)
       if (mpi%rank == mpi%output_rank) then
          call assert_equals(size(rock_cells), size(sim_rock_cells), &
               "Flow simulation num rocktypes")
          call assert_equals(rock_cells, sim_rock_cells, &
               size(rock_cells), "Flow simulation num rocktype cells")
       end if
       deallocate(sim_rock_cells)
    end if

  end subroutine flow_simulation_label_test

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

    json => fson_parse_mpi(trim(path) // "test_init.json")

    call sim%init(json)
    call fson_destroy_mpi(json)

    call flow_simulation_basic_test(title = "Test flow simulation init", &
         thermo = "IAPWS-97", eos = "w", dim = 3, dof = 12)

    call flow_simulation_label_test(rock_cells = [8, 4])

    call vec_diff_test(sim%solution, "primary", path, sim%mesh%cell_order_inv)

    call vec_diff_test(sim%rock, "rock", path, sim%mesh%cell_order_inv)

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

    call sim%init(json)
    call fson_destroy_mpi(json)

    call sim%pre_solve(time, sim%solution, err)
    call vec_diff_test(sim%fluid, "fluid", path, sim%mesh%cell_order_inv)
    
    call sim%destroy()

  end subroutine test_flow_simulation_fluid_properties

!------------------------------------------------------------------------

  subroutine test_flow_simulation_lhs
    ! Test LHS function

    use fson_mpi_module

    ! Locals:
    type(fson_value), pointer :: json
    character(64), parameter :: path = "data/flow_simulation/lhs/"
    PetscReal :: time = 0._dp
    Vec :: lhs
    PetscErrorCode :: ierr, err

    json => fson_parse_mpi(trim(path) // "test_lhs.json")

    call sim%init(json)
    call fson_destroy_mpi(json)

    call DMGetGlobalVector(sim%mesh%dm, lhs, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(lhs, "lhs", ierr); CHKERRQ(ierr)

    call sim%pre_solve(time, sim%solution, err)
    call sim%lhs(time, sim%solution, lhs, err)
    call vec_diff_test(lhs, "lhs", path, sim%mesh%cell_order_inv)

    call DMRestoreGlobalVector(sim%mesh%dm, lhs, ierr); CHKERRQ(ierr)
    call sim%destroy()

  end subroutine test_flow_simulation_lhs

!------------------------------------------------------------------------

end module flow_simulation_test
