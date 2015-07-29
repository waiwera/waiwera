module source_module
 !! Module for handling sinks and sources.

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  public :: setup_source_vector

contains

!------------------------------------------------------------------------

  subroutine setup_source_vector(json, dm, num_primary, source)
    !! Sets up sinks and sources. Source strengths are stored (for
    !! now) in the source vector, with values for all components in
    !! all cells.

    use kinds_module
    use fson
    use fson_mpi_module

    type(fson_value), pointer, intent(in) :: json
    DM, intent(in) :: dm
    PetscInt, intent(in) :: num_primary
    Vec, intent(in out) :: source
    ! Locals:
    PetscErrorCode :: ierr
    PetscInt :: c, isrc, i
    type(fson_value), pointer :: sources, src
    PetscInt :: num_sources, cell, count
    PetscReal :: val
    PetscReal, allocatable :: values(:)
    PetscInt, allocatable :: indices(:)
    PetscInt, parameter :: default_component = 1

    call DMCreateGlobalVector(dm, source, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(source, "source", ierr); CHKERRQ(ierr)

    call VecSet(source, 0._dp, ierr); CHKERRQ(ierr)

    if (fson_has_mpi(json, "source")) then

       call VecGetSize(source, count, ierr); CHKERRQ(ierr)
       allocate(values(count), indices(count))
       call fson_get_mpi(json, "source", sources)
       num_sources = fson_value_count_mpi(sources, ".")
       do isrc = 1, num_sources
          src => fson_value_get_mpi(sources, isrc)
          call fson_get_mpi(src, "cell", val = cell)
          call fson_get_mpi(src, "component", default_component, c)
          call fson_get_mpi(src, "value", val = val)
          values(cell * num_primary + c) = val
       end do

       do i = 1, count
          indices(i) = i-1
       end do
       call VecSetValues(source, count, indices, &
            values, INSERT_VALUES, ierr); CHKERRQ(ierr)
       call VecAssemblyBegin(source, ierr); CHKERRQ(ierr)
       call VecAssemblyEnd(source, ierr); CHKERRQ(ierr)
       deallocate(values, indices)

    end if

  end subroutine setup_source_vector

!------------------------------------------------------------------------

end module source_module
