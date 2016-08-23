module dm_utils_module
  !! Module for utilities related to PETSc DMs (which handle data management on parallel meshes), and Vecs.

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  public :: set_dm_data_layout, section_offset, global_section_offset
  public :: global_vec_section, local_vec_section
  public :: global_to_local_vec_section, restore_dm_local_vec
  public :: global_vec_range_start, vec_reorder
  public :: dm_cell_normal_face
  public :: write_vec_vtk
  public :: vec_max_pointwise_abs_scale

contains

!------------------------------------------------------------------------

  subroutine set_dm_data_layout(dm, num_components, field_dim, &
       field_name)
    !! Sets data layout on default section of the given DM.

    DM, intent(in out) :: dm !! DM object
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
    ! Create the global section:
    call DMGetDefaultGlobalSection(dm, section, ierr); CHKERRQ(ierr)
    deallocate(num_dof)

  end subroutine set_dm_data_layout

!------------------------------------------------------------------------

  subroutine global_vec_range_start(v, range_start)
    !! Gets global start of global range from PetscLayout of global section
    !! on DM that vector v is defined on.

    use mpi_module

    Vec, intent(in) :: v !! Global vector
    PetscInt, intent(out) :: range_start !! Range start for global section
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

    PetscSection, intent(in) :: section !! Local section
    PetscInt, intent(in) :: p !! Mesh point in DM
    PetscInt, intent(out) :: offset !! Offset value
    PetscErrorCode, intent(out) :: ierr !! Error flag

    call PetscSectionGetOffset(section, p, offset, ierr)
    offset = offset + 1

  end subroutine section_offset

!------------------------------------------------------------------------

  subroutine global_section_offset(section, p, range_start, offset, ierr)
    !! Wrapper for PetscSectionGetOffset(), adding one to the result for
    !! Fortran 1-based indexing. For global sections, we also need to
    !! subtract the layout range start to get indices suitable for
    !! indexing into arrays.

    PetscSection, intent(in) :: section !! Global section
    PetscInt, intent(in) :: p !! Mesh point
    PetscInt, intent(in) :: range_start !! Start of PetscLayout range
    PetscInt, intent(out) :: offset !! Offset value
    PetscErrorCode, intent(out) :: ierr !! Error flag
    ! Locals:

    call PetscSectionGetOffset(section, p, offset, ierr)
    offset = offset + 1 - range_start

  end subroutine global_section_offset

!------------------------------------------------------------------------

  subroutine global_vec_section(v, section)
    !! Gets default global PETSc section from DM of a vector v, and
    !! range_start from the PetscLayout of the section.

    Vec, intent(in) :: v !! Global vector
    PetscSection, intent(out) :: section !! Default global section
    ! Locals:
    DM :: dm
    PetscErrorCode :: ierr

    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMGetDefaultGlobalSection(dm, section, ierr); CHKERRQ(ierr)

  end subroutine global_vec_section

!------------------------------------------------------------------------

  subroutine local_vec_section(v, section)
    !! Gets default local PETSc section from DM of a vector v.

    Vec, intent(in) :: v !! Local vector
    PetscSection, intent(out) :: section !! Default local section
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

    Vec, intent(in) :: v !! Global vector
    Vec, intent(out) :: local_v !! Local vector
    PetscSection, intent(out) :: section !! Default local section
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

    Vec, intent(out) :: local_v !! Local vector
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

    Vec, intent(in) :: v !! Vector
    character(len = *), intent(in) :: filename !! VTK output filename
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

  subroutine vec_reorder(v, old_index, new_index)
    !! Reorders a vector from the specified old ordering to new
    !! ordering.

    use mpi_module

    Vec, intent(in out) :: v !! Global vector to be reordered
    IS, intent(in) :: old_index !! Index set for old ordering
    IS, intent(in) :: new_index !! Index set for new ordering
    ! Locals:
    Vec :: vinitial
    VecScatter :: scatter
    PetscInt :: blocksize
    IS :: old_index_block, new_index_block
    PetscInt, pointer :: indices(:)
    PetscErrorCode :: ierr

    call VecDuplicate(v, vinitial, ierr); CHKERRQ(ierr)
    call VecCopy(v, vinitial, ierr); CHKERRQ(ierr)
    call VecGetBlockSize(v, blocksize, ierr); CHKERRQ(ierr)

    call ISGetIndicesF90(old_index, indices, ierr); CHKERRQ(ierr)
    call ISCreateBlock(mpi%comm, blocksize, size(indices), indices, &
         PETSC_COPY_VALUES, old_index_block, ierr); CHKERRQ(ierr)
    call ISRestoreIndicesF90(old_index, indices, ierr); CHKERRQ(ierr)

    call ISGetIndicesF90(new_index, indices, ierr); CHKERRQ(ierr)
    call ISCreateBlock(mpi%comm, blocksize, size(indices), indices, &
         PETSC_COPY_VALUES, new_index_block, ierr); CHKERRQ(ierr)
    call ISRestoreIndicesF90(new_index, indices, ierr); CHKERRQ(ierr)

    call VecScatterCreate(vinitial, old_index_block, v, &
         new_index_block, scatter, ierr); CHKERRQ(ierr)
    call ISDestroy(old_index_block, ierr); CHKERRQ(ierr)
    call ISDestroy(new_index_block, ierr); CHKERRQ(ierr)

    call VecScatterBegin(scatter, vinitial, v, INSERT_VALUES, &
         SCATTER_FORWARD, ierr); CHKERRQ(ierr)
    call VecScatterEnd(scatter, vinitial, v, INSERT_VALUES, &
         SCATTER_FORWARD, ierr); CHKERRQ(ierr)

    call VecScatterDestroy(scatter, ierr); CHKERRQ(ierr)
    call VecDestroy(vinitial, ierr); CHKERRQ(ierr)

  end subroutine vec_reorder

!------------------------------------------------------------------------

  subroutine dm_cell_normal_face(dm, c, normal, f)
    !! Returns index of DM mesh face on the specified cell c, with
    !! outward normal vector closest to the specified one.

    DM, intent(in) :: dm !! DM
    PetscInt, intent(in) :: c !! Cell mesh point in DM
    PetscReal, intent(in) :: normal(:) !! Normal vector
    PetscInt, intent(out) :: f !! Face mesh point in DM
    ! Locals:
    PetscInt :: i, num_faces, imax
    PetscInt :: dim, start_cell, end_cell
    PetscErrorCode :: ierr
    PetscInt, pointer :: faces(:)
    PetscReal, allocatable, target :: centroid(:), face_normal(:)
    PetscReal, pointer :: pcentroid(:), pface_normal(:)
    PetscReal, allocatable :: cos_theta(:)
    PetscReal :: area, normal_norm, face_normal_norm

    call DMPlexGetHeightStratum(dm, 0, start_cell, &
         end_cell, ierr); CHKERRQ(ierr)

    if ((start_cell <= c) .and. (c < end_cell)) then

       call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)
       allocate(centroid(dim), face_normal(dim))
       pcentroid => centroid
       pface_normal => face_normal
       call DMPlexGetConeSize(dm, c, num_faces, ierr); CHKERRQ(ierr)
       call DMPlexGetCone(dm, c, faces, ierr); CHKERRQ(ierr)
       allocate(cos_theta(num_faces))
       normal_norm = norm2(normal)

       do i = 1, num_faces
          call DMPlexComputeCellGeometryFVM(dm, faces(i), area, &
               pcentroid, pface_normal, ierr); CHKERRQ(ierr)
          face_normal_norm = norm2(face_normal)
          cos_theta(i) = dot_product(normal, face_normal) / &
               (normal_norm * face_normal_norm)
       end do

       imax = maxloc(cos_theta, 1)
       f = faces(imax)

       nullify(pcentroid, pface_normal)
       deallocate(cos_theta, centroid, face_normal)

    else
       f = -1
    end if

  end subroutine dm_cell_normal_face

!------------------------------------------------------------------------

  subroutine vec_max_pointwise_abs_scale(v, scale, tol, maxval, maxloc)
    !! Returns pointwise max absolute value of v / abs(scale), and the
    !! index at which the maximum occurs. Where the value of scale is
    !! less than tol, the value of tol is used instead.

    Vec, intent(in) :: v !! Value vector
    Vec, intent(in) :: scale !! Scale vector
    PetscReal, intent(in) :: tol !! Tolerance for small values
    PetscReal, intent(out) :: maxval !! Maximum value
    PetscInt, intent(out) :: maxloc !! Index of maximum value
    ! Locals:
    DM :: dm
    Vec :: scaled_v
    PetscReal, pointer, contiguous :: v_array(:), scale_array(:)
    PetscReal, pointer, contiguous :: scaled_v_array(:)
    PetscInt :: i, low, hi
    PetscReal :: scale_i
    PetscErrorCode :: ierr

    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMGetGlobalVector(dm, scaled_v, ierr); CHKERRQ(ierr)

    call VecGetArrayReadF90(v, v_array, ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(scale, scale_array, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(scaled_v, scaled_v_array, ierr); CHKERRQ(ierr)

    call VecGetOwnershipRange(v, low, hi, ierr); CHKERRQ(ierr)
    do i = 1, hi - low
       scale_i = max(abs(scale_array(i)), tol)
       scaled_v_array(i) = abs(v_array(i)) / scale_i
    end do

    call VecRestoreArrayReadF90(v, v_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(scale, scale_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(scaled_v, scaled_v_array, ierr)
    CHKERRQ(ierr)

    call VecMax(scaled_v, maxloc, maxval, ierr); CHKERRQ(ierr)

    call DMRestoreGlobalVector(dm, scaled_v, ierr); CHKERRQ(ierr)

  end subroutine vec_max_pointwise_abs_scale

!------------------------------------------------------------------------

end module dm_utils_module
