module simulation_test

  ! Tests for simulation module

  use kinds_module
  use mpi_module
  use fruit
  use fson
  use simulation_module

  implicit none
  private

#include <petsc-finclude/petsc.h90>

public :: test_simulation_init

PetscReal, parameter :: tol = 1.e-6_dp

contains

!------------------------------------------------------------------------
! Utility routines:
!------------------------------------------------------------------------

  subroutine get_rocktype_counts(sim, rock_count)

    ! Returns an array on the output rank containing the total numbers
    ! (summed over all processors) of non-ghost cells with each rock
    ! type index.

    use rock_module, only : rocktype_label_name

    type(simulation_type), intent(in) :: sim
    PetscInt, allocatable, intent(out) :: rock_count(:)
    ! Locals:
    PetscInt :: ir, ic, c, IS_size, ghost
    PetscInt :: num_rocktypes, total_num_rocktypes
    IS :: rock_IS
    DMLabel :: ghost_label
    PetscInt, allocatable :: nrc(:)
    PetscInt, pointer :: rock_cells(:)
    PetscErrorCode :: ierr
    
    call DMPlexGetLabelSize(sim%mesh%dm, rocktype_label_name, &
         num_rocktypes, ierr); CHKERRQ(ierr)
    call MPI_allreduce(num_rocktypes, total_num_rocktypes, 1, MPI_INTEGER, &
         MPI_MAX, mpi%comm, ierr)
    call DMPlexGetLabel(sim%mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    allocate(nrc(total_num_rocktypes), rock_count(total_num_rocktypes))
    do ir = 1, total_num_rocktypes
       nrc(ir) = 0 ! number of rocktype cells on processor
       call DMPlexGetStratumIS(sim%mesh%dm, rocktype_label_name, ir, &
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
       call MPI_reduce(nrc(ir), rock_count(ir), 1, MPI_INTEGER, MPI_SUM, &
            mpi%output_rank, mpi%comm, ierr)
    end do
    deallocate(nrc)

  end subroutine get_rocktype_counts

!------------------------------------------------------------------------

  subroutine vec_diff_test(v, name, path)

    ! Tests vec v against values from HDF5 file with specified base name,
    ! at the given path.

    Vec, intent(in) :: v
    character(*), intent(in) :: name, path
    ! Locals:
    Vec :: vread, diff
    PetscViewer :: viewer
    PetscReal :: diffnorm
    PetscErrorCode :: ierr

    call VecCreate(mpi%comm, vread, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(vread, name, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5Open(mpi%comm, path // trim(name) // ".h5", &
         FILE_MODE_READ, viewer, ierr)
    CHKERRQ(ierr)
    call PetscViewerHDF5PushGroup(viewer, "/", ierr); CHKERRQ(ierr)
    call VecLoad(vread, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
    call VecDuplicate(vread, diff, ierr); CHKERRQ(ierr)
    call VecCopy(vread, diff, ierr); CHKERRQ(ierr)
    call VecAXPY(diff, -1._dp, v, ierr); CHKERRQ(ierr)
    call VecNorm(diff, NORM_2, diffnorm, ierr); CHKERRQ(ierr)
    call assert_equals(0._dp, diffnorm, tol, &
         "Simulation " // trim(name) // " vector")
    call VecDestroy(vread, ierr); CHKERRQ(ierr)
    call VecDestroy(diff, ierr); CHKERRQ(ierr)

  end subroutine vec_diff_test

!------------------------------------------------------------------------

  subroutine simulation_basic_test(sim, title, thermo, eos, dim, dof)

    ! Tests basic simulation parameters.

    type(simulation_type), intent(in) :: sim
    character(*), intent(in) :: title, thermo, eos
    PetscInt, intent(in) :: dim, dof
    ! Locals:
    PetscInt :: sim_dim, sim_dof
    Vec :: x
    PetscErrorCode :: ierr
    
    call DMGetDimension(sim%mesh%dm, sim_dim, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(sim%mesh%dm, x, ierr); CHKERRQ(ierr)
    call VecGetSize(x, sim_dof, ierr); CHKERRQ(ierr)
    call VecDestroy(x, ierr); CHKERRQ(ierr)

    if (mpi%rank == mpi%output_rank) then
       call assert_equals(title, sim%title, "Simulation title")
       call assert_equals(thermo, sim%thermo%name, "Simulation thermodynamics")
       call assert_equals(eos, sim%eos%name, "Simulation EOS")
       call assert_equals(dim, sim_dim, "Simulation mesh dimension")
       call assert_equals(dof, sim_dof, "Simulation mesh dof")
    end if

  end subroutine simulation_basic_test

!------------------------------------------------------------------------

  subroutine simulation_label_test(sim, rock_cells)

    ! Tests simulation labels.

    use mesh_module, only : open_boundary_label_name
    use rock_module, only : rocktype_label_name

    type(simulation_type), intent(in) :: sim
    PetscInt, intent(in) :: rock_cells(:)
    ! Locals:
    PetscBool :: open_bdy, has_rock_label
    PetscInt, allocatable :: sim_rock_cells(:)
    PetscErrorCode :: ierr

    ! Open boundary label:
    call DMPlexHasLabel(sim%mesh%dm, open_boundary_label_name, open_bdy, &
         ierr); CHKERRQ(ierr)
    if (mpi%rank == mpi%output_rank) then
       call assert_equals(.true., open_bdy, "Simulation open boundary label")
    end if

    ! Rock type label:
    call DMPlexHasLabel(sim%mesh%dm, rocktype_label_name, has_rock_label, &
         ierr); CHKERRQ(ierr)
    if (mpi%rank == mpi%output_rank) then
       call assert_equals(.true., has_rock_label, "Simulation rocktype label")
    end if
    if (has_rock_label) then
       call get_rocktype_counts(sim, sim_rock_cells)
       if (mpi%rank == mpi%output_rank) then
          call assert_equals(size(rock_cells), size(sim_rock_cells), &
               "Simulation num rocktypes")
          call assert_equals(rock_cells, sim_rock_cells, &
               size(rock_cells), "Simulation num rocktype cells")
       end if
       deallocate(sim_rock_cells)
    end if

  end subroutine simulation_label_test

!------------------------------------------------------------------------
! Unit test routines:
!------------------------------------------------------------------------

  subroutine test_simulation_init

    ! Test simulation init() method
    ! This uses a simple problem with a 12-cell rectangular mesh and two rock types.

    type(simulation_type) :: sim
    ! Locals:
    character(:), allocatable :: path
    PetscReal, parameter :: expected_initial = 2.e5_dp
    PetscInt :: initial_size
    PetscReal, pointer :: initial(:)
    PetscReal, allocatable :: expected_initial_array(:)
    PetscErrorCode :: ierr

    path = "data/simulation/init/"

    call sim%init(path // "test_init.json")

    call simulation_basic_test(sim, title = "Test simulation init", &
         thermo = "IAPWS-97", eos = "W", dim = 3, dof = 12)

    call simulation_label_test(sim, rock_cells = [9, 3])

    call VecGetArrayF90(sim%initial, initial, ierr); CHKERRQ(ierr)
    initial_size = size(initial)
    allocate(expected_initial_array(initial_size))
    expected_initial_array = expected_initial
    call assert_equals(expected_initial_array, initial, initial_size, &
         "Simulation initial conditions")
    deallocate(expected_initial_array)
    call VecRestoreArrayF90(sim%initial, initial, ierr); CHKERRQ(ierr)

    ! test fluid vector initialized

    call vec_diff_test(sim%rock, "rock", path)

    ! test timestepper initialized

    call sim%destroy()

  end subroutine test_simulation_init

!------------------------------------------------------------------------

end module simulation_test
