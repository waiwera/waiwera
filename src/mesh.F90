module mesh_module
  !! Module for mesh type.

  use kinds_module
  use mpi_module
  use fson

  implicit none

  private

#include <petsc-finclude/petscsys.h>
#include <petsc-finclude/petscdef.h>
#include <petsc-finclude/petscdm.h>

  integer, parameter, public :: max_mesh_filename_length = 200

  type, public :: mesh_type
     !! Mesh type.
     private
     character(max_mesh_filename_length), public :: filename
     DM, public :: dm
     PetscInt, public :: start_cell, end_cell, end_interior_cells
     PetscInt, public :: num_cells, num_interior_cells
   contains
     procedure :: distribute => mesh_distribute
     procedure :: construct_ghost_cells => mesh_construct_ghost_cells
     procedure :: setup_discretization => mesh_setup_discretization
     procedure, public :: init => mesh_init
     procedure, public :: destroy => mesh_destroy
  end type mesh_type

contains

!------------------------------------------------------------------------

  subroutine mesh_distribute(self)
    !! Distributes mesh over processors.
    
    class(mesh_type), intent(in out) :: self
    ! Locals:
    DM :: dist_dm
    PetscErrorCode :: ierr

    call DMPlexDistribute(self%dm, 0, PETSC_NULL_OBJECT, dist_dm, ierr)
    CHKERRQ(ierr)
    if (dist_dm /= 0) then
       call DMDestroy(self%dm, ierr); CHKERRQ(ierr)
       self%dm = dist_dm
    end if
    
  end subroutine mesh_distribute

!------------------------------------------------------------------------

  subroutine mesh_construct_ghost_cells(self)
    !! Constructs ghost cells on the mesh.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    DM :: ghost_dm
    PetscErrorCode :: ierr

    call DMPlexConstructGhostCells(self%dm, PETSC_NULL_CHARACTER, &
         PETSC_NULL_INTEGER, ghost_dm, ierr); CHKERRQ(ierr)
    if (ghost_dm /= 0) then
       call DMDestroy(self%dm, ierr); CHKERRQ(ierr);
       self%dm = ghost_dm
    end if

  end subroutine mesh_construct_ghost_cells

!------------------------------------------------------------------------

  subroutine mesh_setup_discretization(self, dof, dim)
    !! Sets up finite volume discretization for the mesh.

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: dof !! Degrees of freedom
    PetscInt, intent(in) :: dim !! Spatial dimension
    ! Locals:
    PetscFV :: fvm
    PetscDS :: ds
    PetscErrorCode :: ierr

    call PetscFVCreate(mpi%comm, fvm, ierr); CHKERRQ(ierr)
    call PetscFVSetNumComponents(fvm, dof, ierr); CHKERRQ(ierr)
    call PetscFVSetSpatialDimension(fvm, dim, ierr); CHKERRQ(ierr)

    call DMGetDS(self%dm, ds, ierr); CHKERRQ(ierr)
    call PetscDSAddDiscretization(ds, fvm, ierr); CHKERRQ(ierr)

  end subroutine mesh_setup_discretization

!------------------------------------------------------------------------

  subroutine mesh_get_bounds(self)
    !! Gets cell bounds on current processor.

    class(mesh_type), intent(in out) :: self

    call DMPlexGetHybridBounds(self%dm, self%end_interior_cell, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         ierr)
    CHKERRQ(ierr)

    call DMPlexGetHeightStratum(self%dm, 0, self%start_cell, &
         self%end_cell, ierr)
    CHKERRQ(ierr)

    self%num_cells = self%end_cell - self%start_cell
    self%num_interior_cells = self%end_interior_cell - self%start_cell

  end subroutine mesh_get_bounds

!------------------------------------------------------------------------

  subroutine mesh_init(self, json, dof)
    !! Initializes mesh, reading filename from JSON input file.
    !! If the filename is not present, an error is raised.
    !! Otherwise, the PETSc DM is read in.

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    PetscInt, intent(in) :: dof !! Degrees of freedom on cells
    ! Locals:
    PetscErrorCode :: ierr
    type(fson_value), pointer :: mesh
    PetscInt, parameter :: dim = 3

    if (mpi%rank == mpi%input_rank) then
       call fson_get(json, "mesh", mesh)
       if (associated(mesh)) then
          call fson_get(mesh, ".", self%filename)
       else
          write (*,'(a,a)') "Stopping: mesh not found in input."
          stop
       end if
    end if

    call MPI_bcast(self%filename, max_mesh_filename_length, MPI_CHARACTER, &
         mpi%input_rank, mpi%comm, ierr)

    ! Read in DM:
    call DMPlexCreateFromFile(mpi%comm, self%filename, PETSC_TRUE, self%dm, ierr)
    CHKERRQ(ierr)

    ! Set up adjacency for finite volume mesh:
    call DMPlexSetAdjacencyUseCone(self%dm, PETSC_TRUE, ierr); CHKERRQ(ierr)
    call DMPlexSetAdjacencyUseClosure(self%dm, PETSC_FALSE, ierr); CHKERRQ(ierr)

    call self%distribute()

    call self%construct_ghost_cells()

    call self%setup_discretization(dof, dim)

    call self%get_bounds()

  end subroutine mesh_init

!------------------------------------------------------------------------

  subroutine mesh_destroy(self)
    !! Destroys mesh object.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    
    call DMDestroy(self%dm, ierr); CHKERRQ(ierr)

  end subroutine mesh_destroy

!------------------------------------------------------------------------

end module mesh_module
