module rock_mechanics_module
  implicit none
  private
#include <petsc/finclude/petsc.h90>

  PetscInt, parameter :: max_primary_variable_name_length = 12
  PetscInt, parameter :: max_mesh_filename_length = 200
  integer, parameter :: dp=kind(0.d0)
  
  type, public :: rock_mechanics_type
     private
     DM :: dm_rock
     PetscInt :: ndof
     character(max_primary_variable_name_length), allocatable :: primary_variable_names(:)
     character(max_mesh_filename_length) :: filename
     
   contains
     private
     procedure :: mesh_init => rock_mechanics_mesh_init
     procedure, public :: init => rock_mechanics_init

  end type rock_mechanics_type

contains
  !------------------------------------------------------------------------
  ! public routines
  !------------------------------------------------------------------------
  subroutine rock_mechanics_init(self,json,logfile)
    use fson
    use logfile_module
    !use rock_parameters_module, only : setup_rocktype_labels
    
    class(rock_mechanics_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    type(logfile_type), intent(in out), optional :: logfile
    
    print *,"in rock mechanics init"
    call self%mesh_init(json,logfile)

    !call self%setup_rocktype_labels()
    
  end subroutine rock_mechanics_init
  !------------------------------------------------------------------------

  
  !------------------------------------------------------------------------
  ! private routines
  !------------------------------------------------------------------------
   subroutine rock_mechanics_mesh_init(self,json,logfile)
    use fson
    use mpi_module
    use logfile_module
    
    class(rock_mechanics_type), intent (in out) :: self
    type(fson_value), pointer, intent(in) :: json
    type(logfile_type), intent(in out), optional :: logfile
    
    ! locals
    PetscErrorCode :: ierr
    type(fson_value), pointer :: mesh

    if (mpi%rank == mpi%input_rank) then
       call fson_get(json, "mesh", mesh)
       if (associated(mesh)) then
          call fson_get(mesh, ".", self%filename)
       else
          if (present(logfile)) then
             call logfile%write(LOG_LEVEL_ERR, 'mesh', 'init', &
                  str_key = 'stop           ', &
                  str_value = 'mesh not found in input.', &
                  rank = mpi%input_rank)
          end if
          stop
       end if
    end if

    print *,"mesh filename is:",self%filename
    
    call MPI_bcast(self%filename, max_mesh_filename_length, MPI_CHARACTER, &
         mpi%input_rank, mpi%comm, ierr)

    ! Read in DM:
    call DMPlexCreateFromFile(mpi%comm, self%filename, PETSC_TRUE, self%dm_rock, ierr)
    CHKERRQ(ierr)    
    call DMGetDimension(self%dm_rock, self%ndof, ierr); CHKERRQ(ierr)
    self%primary_variable_names = ["displacement"]
    ! this is from mesh.F90 -- not sure what it does
    call DMPlexSetAdjacencyUseCone(self%dm_rock, PETSC_TRUE, ierr); CHKERRQ(ierr)
    call DMPlexSetAdjacencyUseClosure(self%dm_rock, PETSC_FALSE, ierr); CHKERRQ(ierr)
    
  end subroutine rock_mechanics_mesh_init
  !------------------------------------------------------------------------
end module rock_mechanics_module
