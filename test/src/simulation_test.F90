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
       end if
       call MPI_reduce(nrc(ir), rock_count(ir), 1, MPI_INTEGER, MPI_SUM, &
            mpi%output_rank, mpi%comm, ierr)
    end do
    deallocate(nrc)

  end subroutine get_rocktype_counts

!------------------------------------------------------------------------

  subroutine test_simulation_init

    ! Test simulation init() method.
    ! This uses a simple problem with a 12-cell rectangular mesh and two rock types.

    use mesh_module, only : open_boundary_label_name
    use rock_module, only : rocktype_label_name

    type(simulation_type) :: sim
    ! Locals:
    character(max_filename_length), parameter :: filename = "data/simulation/init/test_init.json"
    character(20), parameter :: expected_title = "Test simulation init"
    character(16), parameter :: expected_thermo = "IAPWS-97"
    character(1), parameter  :: expected_eos = "W"
    PetscInt, parameter :: expected_dim = 3, num_cells = 12, num_primary = 1
    PetscInt, parameter :: expected_num_rocktypes = 2
    PetscInt, parameter :: expected_num_rocktype_cells(2) = [9, 3]
    PetscInt :: dim, global_solution_dof, num_rocktypes
    PetscInt, allocatable :: rock_count(:)
    Vec :: x
    PetscBool :: open_bdy, has_rock_label
    PetscErrorCode :: ierr

    call sim%init(filename)

    call DMGetDimension(sim%mesh%dm, dim, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(sim%mesh%dm, x, ierr); CHKERRQ(ierr)
    call VecGetSize(x, global_solution_dof, ierr); CHKERRQ(ierr)
    call VecDestroy(x, ierr); CHKERRQ(ierr)

    if (mpi%rank == mpi%output_rank) then
       call assert_equals(expected_title, sim%title, "Simulation title")
       call assert_equals(expected_thermo, sim%thermo%name, &
            "Simulation thermodynamics")
       call assert_equals(expected_eos, sim%eos%name, "Simulation EOS")
       call assert_equals(expected_dim, dim, "Simulation mesh dimension")
       call assert_equals(num_cells * num_primary, global_solution_dof, &
            "Simulation mesh dof")
    end if

    call DMPlexHasLabel(sim%mesh%dm, open_boundary_label_name, open_bdy, &
         ierr); CHKERRQ(ierr)
    if (mpi%rank == mpi%output_rank) then
       call assert_equals(.true., open_bdy, "Simulation open boundary label")
    end if

    call DMPlexHasLabel(sim%mesh%dm, rocktype_label_name, has_rock_label, &
         ierr); CHKERRQ(ierr)
    if (mpi%rank == mpi%output_rank) then
       call assert_equals(.true., has_rock_label, "Simulation rocktype label")
    end if
    if (has_rock_label) then
       call get_rocktype_counts(sim, rock_count)
       if (mpi%rank == mpi%output_rank) then
          call assert_equals(expected_num_rocktypes, size(rock_count), &
               "Simulation num rocktypes")
          call assert_equals(expected_num_rocktype_cells, rock_count, &
               expected_num_rocktypes, "Simulation num rocktype cells")
       end if
       deallocate(rock_count)
    end if

    call sim%destroy()

  end subroutine test_simulation_init

!------------------------------------------------------------------------

end module simulation_test
