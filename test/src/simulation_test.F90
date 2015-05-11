module simulation_test

  ! Tests for simulation module

  use kinds_module
  use mpi_module
  use fruit
  use fson
  use simulation_module
  use timestepping_module

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
    ! NB: at present this will only work if the number of processors is the
    ! same as that used to write the HDF5 file.
    ! This should be fixed when DMPlex supports a 'natural ordering' similar
    ! to that used by DMDA. At that point the code here will need to be altered.

    Vec, intent(in) :: v
    character(*), intent(in) :: name, path
    ! Locals:
    DM :: dm
    Vec :: vread, diff
    PetscViewer :: viewer
    PetscReal :: diffnorm
    PetscErrorCode :: ierr

    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMGetGlobalVector(dm, vread, ierr); CHKERRQ(ierr)
    call DMGetGlobalVector(dm, diff, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(vread, name, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5Open(mpi%comm, path // trim(name) // ".h5", &
         FILE_MODE_READ, viewer, ierr)
    CHKERRQ(ierr)
    call PetscViewerHDF5PushGroup(viewer, "/", ierr); CHKERRQ(ierr)
    call PetscViewerHDF5PushGroup(viewer, "fields", ierr); CHKERRQ(ierr)
    call VecLoad(vread, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
    call VecDuplicate(vread, diff, ierr); CHKERRQ(ierr)
    call VecCopy(vread, diff, ierr); CHKERRQ(ierr)
    call VecAXPY(diff, -1._dp, v, ierr); CHKERRQ(ierr)
    call VecNorm(diff, NORM_2, diffnorm, ierr); CHKERRQ(ierr)
    call assert_equals(0._dp, diffnorm, tol, &
         "Simulation " // trim(name) // " vector")
    call DMRestoreGlobalVector(dm, diff, ierr); CHKERRQ(ierr)
    call DMRestoreGlobalVector(dm, vread, ierr); CHKERRQ(ierr)

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

subroutine timestepper_test(ts, method, initial_time, initial_stepsize, &
     final_time, max_num_steps, max_stepsize)

  ! Tests timestepper initialization.

  type(timestepper_type), intent(in) :: ts
  character(*), intent(in) :: method
  PetscReal, intent(in) :: initial_time, initial_stepsize
  PetscReal, intent(in) :: final_time, max_stepsize
  PetscInt, intent(in)  :: max_num_steps

  if (mpi%rank == mpi%output_rank) then
     call assert_equals(method, ts%method%name, "Simulation timestepper method")
     call assert_equals(initial_time, ts%steps%current%time, "Simulation initial time")
     call assert_equals(initial_stepsize, ts%steps%next_stepsize, "Simulation initial stepsize")
     call assert_equals(final_time, ts%steps%final_time, "Simulation final time")
     call assert_equals(max_num_steps, ts%steps%max_num, "Simulation max num steps")
     call assert_equals(max_stepsize, ts%steps%adaptor%max_stepsize, "Simulation max stepsize")
  end if

end subroutine timestepper_test

!------------------------------------------------------------------------
! Unit test routines:
!------------------------------------------------------------------------

  subroutine test_simulation_init

    ! Test simulation init() method
    ! This uses a simple problem with a 12-cell rectangular mesh and two rock types.

    type(simulation_type) :: sim
    ! Locals:
    character(:), allocatable :: path

    path = "data/simulation/init/"

    call sim%init(path // "test_init.json")

    call simulation_basic_test(sim, title = "Test simulation init", &
         thermo = "IAPWS-97", eos = "W", dim = 3, dof = 12)

    call simulation_label_test(sim, rock_cells = [9, 3])

    call vec_diff_test(sim%initial, "initial", path)

    ! test fluid vector initialized

    call vec_diff_test(sim%rock, "rock", path)

    call timestepper_test(sim%timestepper, "Backward Euler", 0.0_dp, &
         100._dp, 1000._dp, 15, 0.0_dp)

    call sim%destroy()

  end subroutine test_simulation_init

!------------------------------------------------------------------------

end module simulation_test
