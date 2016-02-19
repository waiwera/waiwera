module dm_utils_module
  !! Module for PETSc DM utilities.

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  public :: set_dm_data_layout, section_offset, global_section_offset
  public :: global_vec_section, local_vec_section
  public :: global_to_local_vec_section, restore_dm_local_vec
  public :: global_vec_range_start, vec_reorder

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

  subroutine global_vec_range_start(v, range_start)
    !! Gets global start of global range from PetscLayout of global section
    !! on DM that vector v is defined on.

    use mpi_module

    Vec, intent(in) :: v
    PetscInt, intent(out) :: range_start
    ! Locals:
    PetscSection :: section
    PetscLayout :: layout
    PetscInt :: range_end
    PetscErrorCode :: ierr

    call global_vec_section(v, section)

    call PetscSectionGetValueLayout(mpi%comm, section, layout, ierr)
    CHKERRQ(ierr)
    call PetscLayoutGetRange(layout, range_start, range_end, ierr)
    CHKERRQ(ierr)
    call PetscLayoutDestroy(layout, ierr)
    CHKERRQ(ierr)

  end subroutine global_vec_range_start

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

  subroutine global_section_offset(section, p, range_start, offset, ierr)
    !! Wrapper for PetscSectionGetOffset(), adding one to the result for
    !! Fortran 1-based indexing. For global sections, we also need to
    !! subtract the layout range start to get indices suitable for
    !! indexing into arrays.

    PetscSection, intent(in) :: section !! PETSc section
    PetscInt, intent(in) :: p !! Mesh point
    PetscInt, intent(in) :: range_start !! Start of PetscLayout range
    PetscInt, intent(out) :: offset
    PetscErrorCode, intent(out) :: ierr
    ! Locals:

    call PetscSectionGetOffset(section, p, offset, ierr)
    offset = offset + 1 - range_start

  end subroutine global_section_offset

!------------------------------------------------------------------------

  subroutine global_vec_section(v, section)
    !! Gets default global PETSc section from DM of a vector v, and
    !! range_start from the PetscLayout of the section.

    Vec, intent(in) :: v
    PetscSection, intent(out) :: section
    ! Locals:
    DM :: dm
    PetscErrorCode :: ierr

    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMGetDefaultGlobalSection(dm, section, ierr); CHKERRQ(ierr)

  end subroutine global_vec_section

!------------------------------------------------------------------------

  subroutine local_vec_section(v, section)
    !! Gets default local PETSc section from DM of a vector v.

    Vec, intent(in) :: v
    PetscSection, intent(out) :: section
    ! Locals:
    DM :: dm
    PetscErrorCode :: ierr

    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMGetDefaultSection(dm, section, ierr); CHKERRQ(ierr)

  end subroutine local_vec_section

!------------------------------------------------------------------------

  subroutine global_to_local_vec_section(v, local_v, section)
    !! Takes a global vector v and returns a local vector, with values
    !! scattered from the global vector, and the default local PETSc
    !! section from the DM of the global vector.

    Vec, intent(in) :: v
    Vec, intent(out) :: local_v
    PetscSection, intent(out) :: section
    ! Locals:
    DM :: dm
    PetscErrorCode :: ierr

    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMGetDefaultSection(dm, section, ierr); CHKERRQ(ierr)

    call DMGetLocalVector(dm, local_v, ierr); CHKERRQ(ierr)
    call DMGlobalToLocalBegin(dm, v, INSERT_VALUES, local_v, ierr)
    CHKERRQ(ierr)
    call DMGlobalToLocalEnd(dm, v, INSERT_VALUES, local_v, ierr)
    CHKERRQ(ierr)

  end subroutine global_to_local_vec_section

!------------------------------------------------------------------------

  subroutine restore_dm_local_vec(local_v)
    !! Restores local vector to its DM.

    Vec, intent(out) :: local_v
    ! Locals:
    DM :: dm
    PetscErrorCode :: ierr

    call VecGetDM(local_v, dm, ierr); CHKERRQ(ierr)
    call DMRestoreLocalVector(dm, local_v, ierr); CHKERRQ(ierr)

  end subroutine restore_dm_local_vec

!------------------------------------------------------------------------

  subroutine write_vec_vtk(v, filename)
    !! Writes vector v to VTK file.

    use mpi_module

    Vec, intent(in) :: v
    character(len = *), intent(in) :: filename
    ! Locals:
    PetscViewer :: viewer
    PetscErrorCode :: ierr

    call PetscViewerCreate(mpi%comm, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerSetType(viewer, PETSCVIEWERVTK, ierr); CHKERRQ(ierr)
    call PetscViewerFileSetName(viewer, filename, ierr); CHKERRQ(ierr)
    call VecView(v, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

  end subroutine write_vec_vtk

!------------------------------------------------------------------------

  subroutine vec_reorder(v, old_order, new_order)
    !! Reorders a vector from the specified old ordering to new ordering.

    use mpi_module

    Vec, intent(in out) :: v
    IS, intent(in), optional :: old_order
    IS, intent(in) :: new_order
    ! Locals:
    Vec :: vinitial
    VecScatter :: scatter
    PetscInt :: blocksize
    IS :: old_order_block, new_order_block
    PetscInt, pointer :: indices(:)
    PetscErrorCode :: ierr

    call VecDuplicate(v, vinitial, ierr); CHKERRQ(ierr)
    call VecCopy(v, vinitial, ierr); CHKERRQ(ierr)
    call VecGetBlockSize(v, blocksize, ierr); CHKERRQ(ierr)

    ! Create index sets with the appropriate block size:
    call ISGetIndicesF90(old_order, indices, ierr); CHKERRQ(ierr)
    call ISCreateBlock(mpi%comm, blocksize, size(indices), indices, &
         PETSC_COPY_VALUES, old_order_block, ierr); CHKERRQ(ierr)
    call ISRestoreIndicesF90(old_order, indices, ierr); CHKERRQ(ierr)
    call ISGetIndicesF90(new_order, indices, ierr); CHKERRQ(ierr)
    call ISCreateBlock(mpi%comm, blocksize, size(indices), indices, &
         PETSC_COPY_VALUES, new_order_block, ierr); CHKERRQ(ierr)
    call ISRestoreIndicesF90(new_order, indices, ierr); CHKERRQ(ierr)

    call VecScatterCreate(vinitial, old_order_block, v, new_order_block, &
         scatter, ierr); CHKERRQ(ierr)
    call ISDestroy(old_order_block, ierr); CHKERRQ(ierr)
    call ISDestroy(new_order_block, ierr); CHKERRQ(ierr)

    call VecScatterBegin(scatter, vinitial, v, INSERT_VALUES, &
         SCATTER_FORWARD, ierr); CHKERRQ(ierr)
    call VecScatterEnd(scatter, vinitial, v, INSERT_VALUES, &
         SCATTER_FORWARD, ierr); CHKERRQ(ierr)

    call VecScatterDestroy(scatter, ierr); CHKERRQ(ierr)
    call VecDestroy(vinitial, ierr); CHKERRQ(ierr)

  end subroutine vec_reorder

!------------------------------------------------------------------------

end module dm_utils_module
