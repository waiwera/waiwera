module dm_utils_module
  !! Module for PETSc DM utilities.

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  public :: set_dm_data_layout, section_offset, vec_section

contains

!------------------------------------------------------------------------

  subroutine set_dm_data_layout(dm, num_components, field_dim, &
       field_name)
    !! Sets data layout on default section of the given DM.

    DM, intent(in out) :: dm
    PetscInt, target, intent(in) :: num_components(:) !! Number of components in each field
    PetscInt, intent(in) :: field_dim(:)  !! Dimension each field is defined on (0 = nodes, etc.)
    character(*), intent(in), optional :: field_name(:) !! Name of each field
    ! Locals:
    PetscInt :: dim
    PetscSection :: section
    PetscInt :: num_fields, i, num_bc
    PetscInt, allocatable, target :: num_dof(:)
    PetscInt, target :: bc_field(1)
    IS, target :: bc_comps(1), bc_points(1)
    PetscInt, pointer :: pnum_components(:), pnum_dof(:), pbc_field(:)
    IS, pointer :: pbc_comps(:), pbc_points(:)
    PetscErrorCode :: ierr

    call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)
    num_fields = size(num_components)
    allocate(num_dof(num_fields*(dim+1)))
    num_dof = 0
    do i = 1, num_fields
       num_dof((i-1) * (dim+1) + field_dim(i) + 1) = num_components(i)
    end do

    ! Boundary conditions (none):
    num_bc = 0
    bc_field(1) = 0

    pnum_components => num_components
    pnum_dof => num_dof
    pbc_field => bc_field
    pbc_comps => bc_comps
    pbc_points => bc_points

    call DMPlexCreateSection(dm, dim, num_fields, pnum_components, &
         pnum_dof, num_bc, pbc_field, pbc_comps, pbc_points, &
         PETSC_NULL_OBJECT, section, ierr); CHKERRQ(ierr)

    if (present(field_name)) then
       do i = 1, num_fields
          call PetscSectionSetFieldName(section, i-1, field_name(i), ierr)
          CHKERRQ(ierr)
       end do
    end if

    call DMSetDefaultSection(dm, section, ierr); CHKERRQ(ierr)
    call PetscSectionDestroy(section, ierr); CHKERRQ(ierr)
    ! Create the global section:
    call DMGetDefaultGlobalSection(dm, section, ierr); CHKERRQ(ierr)
    deallocate(num_dof)

  end subroutine set_dm_data_layout

!------------------------------------------------------------------------

  subroutine section_offset(section, p, offset, ierr)
    !! Wrapper for PetscSectionGetOffset(), adding one to the result for
    !! Fortran 1-based indexing.

    PetscSection, intent(in) :: section !! PETSc section
    PetscInt, intent(in) :: p !! Mesh point
    PetscInt, intent(out) :: offset
    PetscErrorCode, intent(out) :: ierr

    call PetscSectionGetOffset(section, p, offset, ierr)
    offset = offset + 1

  end subroutine section_offset

!------------------------------------------------------------------------

  subroutine vec_section(v, section)
    !! Gets default PETSc section from DM of a vector v.

    Vec, intent(in) :: v
    PetscSection, intent(out) :: section
    ! Locals:
    DM :: dm
    PetscErrorCode :: ierr

    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMGetDefaultSection(dm, section, ierr); CHKERRQ(ierr)

  end subroutine vec_section

!------------------------------------------------------------------------

end module dm_utils_module
