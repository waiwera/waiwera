module mesh_module
  !! Module for mesh type.

  use kinds_module
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
   contains
     procedure, public :: init => mesh_init
     procedure, public :: destroy => mesh_destroy
  end type mesh_type

contains

!------------------------------------------------------------------------

  subroutine mesh_init(self, comm, rank, io_rank, json)
    !! Initializes mesh, reading filename from JSON input file.
    !! If the filename is not present, an error is raised.
    !! Otherwise, the PETSc DM is read in.

    class(mesh_type), intent(in out) :: self
    MPI_Comm, intent(in) :: comm !! MPI communicator
    PetscMPIInt, intent(in) :: rank !! MPI rank
    PetscMPIInt, intent(in) :: io_rank !! MPI rank of I/O process
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    ! Locals:
    PetscErrorCode :: ierr
    type(fson_value), pointer :: mesh

    if (rank == io_rank) then
       call fson_get(json, "mesh", mesh)
       if (associated(mesh)) then
          call fson_get(mesh, ".", self%filename)
       else
          write (*,'(a,a)') "Stopping: mesh not found in input."
          stop
       end if
    end if

    call MPI_bcast(self%filename, max_mesh_filename_length, MPI_CHARACTER, &
         io_rank, comm, ierr)

    ! Read in DM:
    call DMPlexCreateFromFile(comm, self%filename, PETSC_TRUE, self%dm, ierr)
    CHKERRQ(ierr)
    
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
