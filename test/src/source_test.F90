module source_test

  ! Test for source module

  use kinds_module
  use mpi_module
  use fruit
  use source_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

public :: test_setup_source_vector

contains
  
!------------------------------------------------------------------------

  subroutine test_setup_source_vector

    ! setup_source_vector() test

    use fson
    use fson_mpi_module
    use mesh_module
    use flow_simulation_test, only: vec_write, vec_diff_test
    use eos_module, only: max_primary_variable_name_length
    use boundary_module, only: open_boundary_label_name

    character(16), parameter :: path = "data/source/"
    PetscInt, parameter :: nc = 2
    PetscInt :: np, range_start, range_end
    character(max_primary_variable_name_length), allocatable :: &
         primary_variable_names(:)
    PetscBool :: isothermal
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    Vec :: source
    PetscSection :: section
    PetscLayout :: layout
    PetscErrorCode :: ierr

    primary_variable_names = &
         [character(max_primary_variable_name_length) :: &
         "P1", "P2", "T"]
    np = size(primary_variable_names)

    isothermal = (nc == np)

    json => fson_parse_mpi(trim(path) // "test_source.json")

    call mesh%init(json)
    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(primary_variable_names)

    call DMGetDefaultGlobalSection(mesh%dm, section, ierr); CHKERRQ(ierr)
    call PetscSectionGetValueLayout(mpi%comm, section, layout, ierr)
    CHKERRQ(ierr)
    call PetscLayoutGetRange(layout, range_start, range_end, ierr)
    CHKERRQ(ierr)
    call PetscLayoutDestroy(layout, ierr)
    CHKERRQ(ierr)

    call setup_source_vector(json, mesh%dm, np, isothermal, source, &
         range_start)
    call vec_diff_test(source, "source", path, mesh%cell_order)

    call VecDestroy(source, ierr); CHKERRQ(ierr)
    call mesh%destroy()
    call fson_destroy_mpi(json)

  end subroutine test_setup_source_vector

!------------------------------------------------------------------------

end module source_test

