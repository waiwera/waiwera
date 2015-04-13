module mesh_test

  ! Tests for mesh module

  use kinds_module
  use mpi_module
  use fruit
  use fson
  use mesh_module

  implicit none
  private

#include <petsc-finclude/petscsys.h>
#include <petsc-finclude/petscdef.h>

public :: test_init

contains

!------------------------------------------------------------------------

  subroutine test_init

    ! Mesh init test

    use eos_module, only: max_primary_variable_name_length
    use face_module

    character(max_mesh_filename_length) :: filename = "data/mesh_test_init.json"
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    character(:), allocatable :: primary_variable_names(:)
    Vec :: v
    type(face_type) :: face
    PetscInt :: celldof, facedof, num_primary
    PetscErrorCode :: ierr
    PetscInt, parameter :: num_cells = 3, num_faces = 16
    
    primary_variable_names = [character(max_primary_variable_name_length) :: &
         "Pressure", "Temperature"]

    if (mpi%rank == mpi%input_rank) then
       json => fson_parse(filename)
    end if

    call mesh%init(json, primary_variable_names)

    call DMCreateGlobalVector(mesh%dm, v, ierr)
    call VecGetSize(v, celldof, ierr)
    if (mpi%rank == mpi%input_rank) then
       num_primary = size(primary_variable_names)
       call assert_equals(num_cells * num_primary, celldof, "global cell dof")
    end if
    call VecDestroy(v, ierr)

    call VecGetSize(mesh%face_geom, facedof, ierr)
    if (mpi%rank == mpi%input_rank) then
       call assert_equals(num_faces * face%dof(), facedof, "global face dof")
    end if

    call mesh%destroy()
    deallocate(primary_variable_names)

  end subroutine test_init

!------------------------------------------------------------------------

end module mesh_test
