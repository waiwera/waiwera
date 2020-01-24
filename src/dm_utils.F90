!   Copyright 2016 University of Auckland.

!   This file is part of Waiwera.

!   Waiwera is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.

!   Waiwera is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.

!   You should have received a copy of the GNU Lesser General Public License
!   along with Waiwera.  If not, see <http://www.gnu.org/licenses/>.

module dm_utils_module
  !! Module for utilities related to PETSc DMs (which handle data management on parallel meshes), and Vecs.

#include <petsc/finclude/petsc.h>

  use petsc

  implicit none
  private

  type, public :: dm_stratum_type
     !! Type for point stratum (cells, faces, edges or vertices) in a
     !! DM. This includes functionality for MINC DM point calculations.
     private
     PetscInt, public :: start, end, end_interior, end_non_ghost
     PetscInt, public, allocatable :: minc_shift(:)
     PetscInt, public :: num_minc_points, num_partition_ghosts
   contains
     private
     procedure, public :: size => dm_stratum_size
     procedure, public :: contains_point => dm_stratum_contains_point
     procedure, public :: destroy => dm_stratum_destroy
     procedure :: minc_point_single => dm_stratum_minc_point_single
     procedure :: minc_point_array => dm_stratum_minc_point_array
     generic, public :: minc_point => minc_point_single, minc_point_array
  end type dm_stratum_type

  interface natural_to_local_cell_index
     module procedure natural_to_local_cell_index_single
     module procedure natural_to_local_cell_index_array
  end interface natural_to_local_cell_index

  interface local_to_natural_cell_index
     module procedure local_to_natural_cell_index_single
     module procedure local_to_natural_cell_index_array
  end interface local_to_natural_cell_index

  public :: dm_get_strata, dm_point_stratum_height
  public :: dm_create_section, dm_set_data_layout, dm_set_default_data_layout
  public :: dm_setup_global_section
  public :: section_offset, global_section_offset
  public :: global_vec_section, local_vec_section
  public :: global_to_local_vec_section, restore_dm_local_vec
  public :: global_vec_range_start, vec_reorder
  public :: dm_cell_normal_face
  public :: write_vec_vtk
  public :: vec_max_pointwise_abs_scale
  public :: dm_set_fv_adjacency
  public :: dm_get_num_partition_ghost_points, dm_get_bdy_cell_shift
  public :: dm_get_end_interior_cell
  public :: dm_get_natural_to_global_ao, dm_get_cell_index
  public :: natural_to_local_cell_index, local_to_natural_cell_index
  public :: dm_natural_order_IS
  public :: create_path_dm
  public :: get_field_subvector, section_get_field_names
  public :: dm_global_cell_field_dof, dm_check_create_label
  public :: dm_label_partition_ghosts, dm_label_boundary_ghosts
  public :: dm_distribute_local_vec, dm_distribute_global_vec
  public :: dm_distribute_index_set
  public :: vec_copy_common_local
  public :: mat_type_is_block, mat_coloring_perturbed_columns
  public :: dm_copy_cone_orientation, dm_cell_counts

contains

!------------------------------------------------------------------------
! DM stratum type:
!------------------------------------------------------------------------

  PetscInt function dm_stratum_size(self) result(n)
    !! Returns size of DM stratum.

    class(dm_stratum_type), intent(in) :: self

    n = self%end - self%start

  end function dm_stratum_size

!------------------------------------------------------------------------

  PetscBool function dm_stratum_contains_point(self, p) result(has)
    !! Returns whether specified DM point p is inside the stratum.

    class(dm_stratum_type), intent(in) :: self
    PetscInt, intent(in) :: p

    has = ((self%start <= p) .and. (p < self%end))

  end function dm_stratum_contains_point

!------------------------------------------------------------------------

  subroutine dm_stratum_destroy(self)
    !! Destroys DM stratum.

    class(dm_stratum_type), intent(in out) :: self

    if (allocated(self%minc_shift)) then
       deallocate(self%minc_shift)
    end if

  end subroutine dm_stratum_destroy

!------------------------------------------------------------------------

  PetscInt function dm_stratum_minc_point_single(self, p, m) &
       result(minc_p)
    !! Returns point for MINC level m in MINC DM corresponding to
    !! point p in original DM. Partition ghost cells are shifted up by
    !! the number of MINC points in the stratum. For MINC points (m >
    !! 0), p should be the index of the fracture cell in the list of
    !! fracture cells for the given MINC level.

    class(dm_stratum_type), intent(in) :: self
    PetscInt, intent(in) :: p, m

    minc_p = p + self%minc_shift(m)
    if ((m == 0) .and. (p >= self%end_non_ghost)) then
       minc_p = minc_p + self%num_minc_points
    end if

  end function dm_stratum_minc_point_single

!........................................................................

  function dm_stratum_minc_point_array(self, p_array, m) &
       result(minc_p_array)
    !! Returns array of MINC points for array of original DM points.

    class(dm_stratum_type), intent(in) :: self
    PetscInt, intent(in) :: p_array(:), m
    PetscInt, allocatable :: minc_p_array(:)
    ! Locals:
    PetscInt :: p, i

    allocate(minc_p_array(size(p_array)))
    do i = 1, size(p_array)
       p = p_array(i)
       minc_p_array(i) = self%minc_point_single(p, m)
    end do

  end function dm_stratum_minc_point_array

!------------------------------------------------------------------------

  subroutine dm_get_strata(dm, depth, strata)
    !! Gets mesh depth and array of strata for the given DM.

    DM, intent(in) :: dm
    PetscInt, intent(out) :: depth
    type(dm_stratum_type), allocatable, intent(out) :: strata(:)
    ! Locals:
    PetscMPIInt :: np
    PetscInt :: h, dummy, end_interior_cell
    PetscErrorCode :: ierr

    call MPI_comm_size(PETSC_COMM_WORLD, np, ierr)

    call DMPlexGetDepth(dm, depth, ierr); CHKERRQ(ierr)
    allocate(strata(0: depth))

    do h = 0, depth
       call DMPlexGetHeightStratum(dm, h, strata(h)%start, &
            strata(h)%end, ierr); CHKERRQ(ierr)
       strata(h)%end_interior = strata(h)%end
    end do

    call DMPlexGetGhostCellStratum(dm, end_interior_cell, &
         dummy, ierr); CHKERRQ(ierr)
    if (end_interior_cell >= 0) then
       strata(0)%end_interior = end_interior_cell
    end if

    if (np == 1) then
       strata%end_non_ghost = strata%end
    else
       do h = 0, depth
          strata(h)%num_partition_ghosts = dm_get_num_partition_ghost_points(dm, h)
          strata(h)%end_non_ghost = strata(h)%end - strata(h)%num_partition_ghosts
       end do
    end if

  end subroutine dm_get_strata

!------------------------------------------------------------------------

  PetscInt function dm_point_stratum_height(strata, p) result(height)
    !! Returns height of stratum containing point p, given array of
    !! strata.

    class(dm_stratum_type), intent(in) :: strata(0:)
    PetscInt, intent(in) :: p
    ! Locals:
    PetscInt :: h

    height = -1
    associate(depth => size(strata) - 1)
      do h = 0, depth
         if (strata(h)%contains_point(p)) then
            height = h
            exit
         end if
      end do
    end associate

  end function dm_point_stratum_height

!------------------------------------------------------------------------
! DM utilities:
!------------------------------------------------------------------------

  subroutine dm_set_fields(dm, num_components)
    !! Sets fields in the DM.

    DM, intent(in out) :: dm
    PetscInt, intent(in) :: num_components(:) !! Number of components in each field
    ! Locals:
    PetscInt :: dim, f
    PetscFV :: fvm
    PetscErrorCode :: ierr

    call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)
    call DMClearFields(dm, ierr); CHKERRQ(ierr)

    associate(num_fields => size(num_components))
      do f = 1, num_fields
         call PetscFVCreate(PETSC_COMM_WORLD, fvm, ierr); CHKERRQ(ierr)
         call PetscFVSetFromOptions(fvm, ierr); CHKERRQ(ierr)
         call PetscFVSetNumComponents(fvm, num_components(f), ierr); CHKERRQ(ierr)
         call PetscFVSetSpatialDimension(fvm, dim, ierr); CHKERRQ(ierr)
         call DMAddField(dm, PETSC_NULL_DMLABEL, fvm, ierr); CHKERRQ(ierr)
         call PetscFVDestroy(fvm, ierr); CHKERRQ(ierr)
      end do
    end associate

    call DMCreateDS(dm, ierr); CHKERRQ(ierr)

  end subroutine dm_set_fields

!------------------------------------------------------------------------

  PetscSection function dm_create_section(dm, num_components, field_dim, &
       field_name) result(section)
    !! Creates section from the given DM and data layout parameters.

    DM, intent(in) :: dm !! DM object
    PetscInt, target, intent(in) :: num_components(:) !! Number of components in each field
    PetscInt, intent(in) :: field_dim(:)  !! Dimension each field is defined on (0 = nodes, etc.)
    character(*), intent(in), optional :: field_name(:) !! Name of each field
    ! Locals:
    PetscInt :: dim
    PetscInt :: i, num_bc
    PetscInt, allocatable, target :: num_dof(:)
    PetscInt, target :: bc_field(1)
    IS, target :: bc_comps(1), bc_points(1)
    PetscInt, pointer :: pnum_components(:), pnum_dof(:), pbc_field(:)
    IS, pointer :: pbc_comps(:), pbc_points(:)
    DMLabel, pointer :: label(:)
    PetscErrorCode :: ierr

    call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)
    associate (num_fields => size(num_components))

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
      label => NULL()

      call DMPlexCreateSection(dm, label, pnum_components, &
           pnum_dof, num_bc, pbc_field, pbc_comps, pbc_points, &
           PETSC_NULL_IS, section, ierr); CHKERRQ(ierr)

      if (present(field_name)) then
         do i = 1, num_fields
            call PetscSectionSetFieldName(section, i-1, field_name(i), ierr)
            CHKERRQ(ierr)
         end do
      end if

      deallocate(num_dof)

    end associate

  end function dm_create_section

!------------------------------------------------------------------------

  subroutine dm_set_data_layout(dm, num_components, field_dim, &
       field_name)
    !! Sets data layout on default section of the given DM.

    DM, intent(in out) :: dm !! DM object
    PetscInt, target, intent(in) :: num_components(:) !! Number of components in each field
    PetscInt, intent(in) :: field_dim(:)  !! Dimension each field is defined on (0 = nodes, etc.)
    character(*), intent(in), optional :: field_name(:) !! Name of each field
    ! Locals:
    PetscSection :: section
    PetscErrorCode :: ierr

    call dm_set_fields(dm, num_components)
    section = dm_create_section(dm, num_components, field_dim, field_name)
    call DMSetSection(dm, section, ierr); CHKERRQ(ierr)

  end subroutine dm_set_data_layout

!------------------------------------------------------------------------

  subroutine dm_set_default_data_layout(dm, dof)
    !! Sets default data layout on DM, for primary variable vector
    !! with specified number of degrees of freedom.

    DM, intent(in out) :: dm
    PetscInt, intent(in) :: dof
    ! Locals:
    PetscInt :: num_components(1), dim, field_dim(1)
    character(7) :: field_names(1)
    PetscErrorCode :: ierr

    call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)
    num_components = dof
    field_dim = dim
    field_names(1) = "Primary"

    call dm_set_data_layout(dm, num_components, field_dim, &
         field_names)

  end subroutine dm_set_default_data_layout

!------------------------------------------------------------------------

  subroutine dm_setup_global_section(dm)
    !! Sets up global section on DM, based on local section and point
    !! SF.

    DM, intent(in out) :: dm
    ! Locals:
    PetscSF :: sf
    PetscSection :: local_section, global_section
    PetscErrorCode :: ierr

    call DMGetPointSF(dm, sf, ierr); CHKERRQ(ierr)
    call DMGetSection(dm, local_section, ierr); CHKERRQ(ierr)
    call PetscSectionCreateGlobalSection(local_section, sf, PETSC_FALSE, &
         PETSC_FALSE, global_section, ierr); CHKERRQ(ierr)
    call DMSetGlobalSection(dm, global_section, ierr)
    CHKERRQ(ierr)

  end subroutine dm_setup_global_section

!------------------------------------------------------------------------

  subroutine global_vec_range_start(v, range_start)
    !! Gets global start of global range from PetscLayout of global section
    !! on DM that vector v is defined on.

    Vec, intent(in) :: v !! Global vector
    PetscInt, intent(out) :: range_start !! Range start for global section
    ! Locals:
    PetscSection :: section
    PetscLayout :: layout
    PetscInt :: range_end
    PetscErrorCode :: ierr

    call global_vec_section(v, section)

    call PetscSectionGetValueLayout(PETSC_COMM_WORLD, section, layout, &
         ierr); CHKERRQ(ierr)
    call PetscLayoutGetRange(layout, range_start, range_end, ierr)
    CHKERRQ(ierr)
    call PetscLayoutDestroy(layout, ierr)
    CHKERRQ(ierr)

  end subroutine global_vec_range_start

!------------------------------------------------------------------------

  PetscInt function section_offset(section, p) result(offset)
    !! Wrapper for PetscSectionGetOffset(), adding one to the result for
    !! Fortran 1-based indexing.

    PetscSection, intent(in) :: section !! Local section
    PetscInt, intent(in) :: p !! Mesh point in DM
    ! Locals:
    PetscErrorCode :: ierr

    call PetscSectionGetOffset(section, p, offset, ierr); CHKERRQ(ierr)
    offset = offset + 1

  end function section_offset

!------------------------------------------------------------------------

  PetscInt function global_section_offset(section, p, &
       range_start) result(offset)
    !! Wrapper for PetscSectionGetOffset(), adding one to the result for
    !! Fortran 1-based indexing. For global sections, we also need to
    !! subtract the layout range start to get indices suitable for
    !! indexing into arrays.

    PetscSection, intent(in) :: section !! Global section
    PetscInt, intent(in) :: p !! Mesh point
    PetscInt, intent(in) :: range_start !! Start of PetscLayout range
    ! Locals:
    PetscErrorCode :: ierr

    call PetscSectionGetOffset(section, p, offset, ierr); CHKERRQ(ierr)
    offset = offset + 1 - range_start

  end function global_section_offset

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
    call DMGetGlobalSection(dm, section, ierr); CHKERRQ(ierr)

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
    call DMGetSection(dm, section, ierr); CHKERRQ(ierr)

  end subroutine local_vec_section

!------------------------------------------------------------------------

  subroutine global_to_local_vec_section(v, local_v, section)
    !! Takes a global vector v and returns a local vector, with values
    !! scattered from the global vector, and the default local PETSc
    !! section from the DM of the global vector.

    Vec, intent(in) :: v !! Global vector
    Vec, intent(in out) :: local_v !! Local vector
    PetscSection, intent(out) :: section !! Default local section
    ! Locals:
    DM :: dm
    PetscErrorCode :: ierr

    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMGetSection(dm, section, ierr); CHKERRQ(ierr)

    call DMGetLocalVector(dm, local_v, ierr); CHKERRQ(ierr)
    call DMGlobalToLocal(dm, v, INSERT_VALUES, local_v, ierr); CHKERRQ(ierr)

  end subroutine global_to_local_vec_section

!------------------------------------------------------------------------

  subroutine restore_dm_local_vec(local_v)
    !! Restores local vector to its DM.

    Vec, intent(in out) :: local_v !! Local vector
    ! Locals:
    DM :: dm
    PetscErrorCode :: ierr

    call VecGetDM(local_v, dm, ierr); CHKERRQ(ierr)
    call DMRestoreLocalVector(dm, local_v, ierr); CHKERRQ(ierr)

  end subroutine restore_dm_local_vec

!------------------------------------------------------------------------

  subroutine write_vec_vtk(v, filename)
    !! Writes vector v to VTK file.

    Vec, intent(in) :: v !! Vector
    character(len = *), intent(in) :: filename !! VTK output filename
    ! Locals:
    PetscViewer :: viewer
    PetscErrorCode :: ierr

    call PetscViewerCreate(PETSC_COMM_WORLD, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerSetType(viewer, PETSCVIEWERVTK, ierr); CHKERRQ(ierr)
    call PetscViewerFileSetName(viewer, filename, ierr); CHKERRQ(ierr)
    call VecView(v, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

  end subroutine write_vec_vtk

!------------------------------------------------------------------------

  subroutine vec_reorder(v, old_index, new_index)
    !! Reorders a vector from the specified old ordering to new
    !! ordering.

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
    call ISCreateBlock(PETSC_COMM_WORLD, blocksize, size(indices), &
         indices, PETSC_COPY_VALUES, old_index_block, ierr); CHKERRQ(ierr)
    call ISRestoreIndicesF90(old_index, indices, ierr); CHKERRQ(ierr)

    call ISGetIndicesF90(new_index, indices, ierr); CHKERRQ(ierr)
    call ISCreateBlock(PETSC_COMM_WORLD, blocksize, size(indices), &
         indices, PETSC_COPY_VALUES, new_index_block, ierr); CHKERRQ(ierr)
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

    use kinds_module

    DM, intent(in) :: dm !! DM
    PetscInt, intent(in) :: c !! Cell mesh point in DM
    PetscReal, intent(in) :: normal(:) !! Normal vector
    PetscInt, intent(out) :: f !! Face mesh point in DM
    ! Locals:
    PetscInt :: i, num_faces, imax, dim, num_cells
    PetscInt :: start_cell, end_cell
    PetscErrorCode :: ierr
    PetscInt, pointer :: faces(:)
    PetscReal, allocatable, target :: centroid(:), face_normal(:)
    PetscReal, pointer :: pcentroid(:), pface_normal(:)
    PetscReal, allocatable :: cos_theta(:)
    PetscReal :: area, normal_norm, face_normal_norm

    call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(dm, 0, start_cell, &
         end_cell, ierr); CHKERRQ(ierr)

    if ((start_cell <= c) .and. (c < end_cell)) then

       allocate(centroid(3), face_normal(3))
       pcentroid => centroid
       pface_normal => face_normal
       call DMPlexGetConeSize(dm, c, num_faces, ierr); CHKERRQ(ierr)
       call DMPlexGetCone(dm, c, faces, ierr); CHKERRQ(ierr)
       allocate(cos_theta(num_faces))
       cos_theta = -1._dp
       normal_norm = norm2(normal(1: dim))

       do i = 1, num_faces
          call DMPlexGetSupportSize(dm, faces(i), num_cells, ierr)
          CHKERRQ(ierr)
          if (num_cells < 2) then
             call DMPlexComputeCellGeometryFVM(dm, faces(i), area, &
                  pcentroid, pface_normal, ierr); CHKERRQ(ierr)
             face_normal_norm = norm2(face_normal(1: dim))
             cos_theta(i) = dot_product(normal(1: dim), &
                  face_normal(1: dim)) / &
                  (normal_norm * face_normal_norm)
          end if
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

  PetscInt function dm_get_num_partition_ghost_points(dm, h) result(n)
    !! Returns number of DM partition ghost points in stratum h on
    !! current process.

    DM, intent(in) :: dm !! DM
    PetscInt, intent(in) :: h !! stratum height
    ! Locals:
    PetscMPIInt :: np
    PetscSF :: point_sf
    PetscInt :: start_point, end_point
    PetscInt :: num_roots, num_leaves
    PetscInt, pointer :: local(:)
    type(PetscSFNode), pointer :: remote(:)
    PetscErrorCode :: ierr

    call MPI_comm_size(PETSC_COMM_WORLD, np, ierr)
    if (np > 1) then
       call DMPlexGetHeightStratum(dm, h, start_point, end_point, ierr)
       call DMGetPointSF(dm, point_sf, ierr); CHKERRQ(ierr)
       call PetscSFGetGraph(point_sf, num_roots, num_leaves, &
            local, remote, ierr); CHKERRQ(ierr)
       n = count((start_point <= local) .and. (local < end_point))
    else
       n = 0
    end if

  end function dm_get_num_partition_ghost_points

!------------------------------------------------------------------------

  PetscInt function dm_get_bdy_cell_shift(dm) result(shift)
    !! Returns number of boundary cells on processes of rank lower
    !! than the current process.

    use utils_module, only: get_mpi_int_gather_array, array_cumulative_sum

    DM, intent(in) :: dm
    ! Locals:
    PetscInt :: start_cell, end_cell, end_interior_cell
    PetscInt :: num_bdy_cells
    PetscMPIInt :: rank, np
    PetscInt, allocatable :: proc_num_bdy_cells(:), &
         proc_sum_bdy_cells(:)
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call MPI_comm_size(PETSC_COMM_WORLD, np, ierr)

    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    end_interior_cell = dm_get_end_interior_cell(dm, end_cell)
    num_bdy_cells = end_cell - end_interior_cell

    proc_num_bdy_cells = get_mpi_int_gather_array()
    proc_sum_bdy_cells = get_mpi_int_gather_array()

    call MPI_gather(num_bdy_cells, 1, MPI_INTEGER, proc_num_bdy_cells, &
         1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
       proc_sum_bdy_cells = [[0], &
            array_cumulative_sum(proc_num_bdy_cells(1: np - 1))]
    end if
    call MPI_scatter(proc_sum_bdy_cells, 1, MPI_INTEGER, shift, &
         1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)

    deallocate(proc_num_bdy_cells, proc_sum_bdy_cells)

  end function dm_get_bdy_cell_shift

!------------------------------------------------------------------------

  PetscInt function dm_get_end_interior_cell(dm, end_cell) &
       result(end_interior_cell)
    !! Returns index (+1) of last interior (i.e. non-boundary) cell of
    !! the DM. In the serial case the result from
    !! DMPlexGetGhostCellStratum() is -1, so here this is corrected to
    !! end_cell.

    DM, intent(in) :: dm
    PetscInt, intent(in) :: end_cell
    ! Locals:
    PetscInt :: dummy
    PetscErrorCode :: ierr

    call DMPlexGetGhostCellStratum(dm, end_interior_cell, dummy, ierr)
    CHKERRQ(ierr)
    if (end_interior_cell < 0) end_interior_cell = end_cell

  end function dm_get_end_interior_cell

!------------------------------------------------------------------------

  AO function dm_get_natural_to_global_ao(dm, cell_natural) result(ao)
    !! Returns application ordering for natural to global mapping on
    !! the DM, given the IS local-to-natural cell mapping.

    DM, intent(in) :: dm
    IS, intent(in) :: cell_natural

    PetscMPIInt :: np
    PetscInt :: c, start_cell, end_cell, end_interior_cell, end_non_ghost_cell
    PetscInt :: num_ghost_cells, num_non_ghost_cells
    PetscInt, allocatable :: global(:), natural(:)
    PetscInt, pointer :: cell_natural_array(:)
    ISLocalToGlobalMapping :: l2g
    PetscErrorCode :: ierr

    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, &
         ierr); CHKERRQ(ierr)
    end_interior_cell = dm_get_end_interior_cell(dm, end_cell)
    call MPI_comm_size(PETSC_COMM_WORLD, np, ierr)
    if (np > 1) then
       num_ghost_cells = dm_get_num_partition_ghost_points(dm, 0)
       num_non_ghost_cells = end_interior_cell - start_cell - num_ghost_cells
       end_non_ghost_cell = start_cell + num_non_ghost_cells
       call DMGetLocalToGlobalMapping(dm, l2g, ierr); CHKERRQ(ierr)
       allocate(natural(start_cell: end_non_ghost_cell - 1), &
            global(start_cell: end_non_ghost_cell - 1))
       call ISGetIndicesF90(cell_natural, cell_natural_array, ierr)
       CHKERRQ(ierr)
       natural = cell_natural_array(1: num_non_ghost_cells)
       call ISRestoreIndicesF90(cell_natural, cell_natural_array, ierr)
       CHKERRQ(ierr)
       call ISLocalToGlobalMappingApplyBlock(l2g, num_non_ghost_cells, &
            [(c, c = start_cell, end_non_ghost_cell - 1)], global, ierr)
       CHKERRQ(ierr)
    else ! serial:
       num_non_ghost_cells = end_interior_cell - start_cell
       end_non_ghost_cell = start_cell + num_non_ghost_cells
       natural = [(c, c = start_cell, end_non_ghost_cell - 1)]
       global = natural
    end if
    call AOCreateMapping(PETSC_COMM_WORLD, num_non_ghost_cells, natural, &
         global, ao, ierr); CHKERRQ(ierr)
    deallocate(natural, global)

  end function dm_get_natural_to_global_ao

!------------------------------------------------------------------------

  PetscInt function natural_to_local_cell_index_single(ao, l2g, &
       natural) result(local)
    !! Returns local cell index corresponding to natural cell
    !! index. Off-process cells are given index values of -1.

    AO, intent(in) :: ao !! Application ordering mapping natural to global cell indices
    ISLocalToGlobalMapping, intent(in) :: l2g !! DM local to global mapping
    PetscInt, intent(in) :: natural !! Natural cell index
    ! Locals:
    PetscInt :: n, idx(1), local_array(1)
    PetscErrorCode :: ierr

    idx(1) = natural
    call AOApplicationToPetsc(ao, 1, idx, ierr); CHKERRQ(ierr)
    call ISGlobalToLocalMappingApplyBlock(l2g, IS_GTOLM_MASK, 1, &
         idx, n, local_array, ierr); CHKERRQ(ierr)
    local = local_array(1)

  end function natural_to_local_cell_index_single

!------------------------------------------------------------------------

  function natural_to_local_cell_index_array(ao, l2g, natural) &
       result(local)
    !! Returns array of local cell indices corresponding to array
    !! of natural cell indices. Any off-process cells are given index
    !! values of -1.

    AO, intent(in) :: ao !! Application ordering mapping natural to global cell indices
    ISLocalToGlobalMapping, intent(in) :: l2g !! DM local to global mapping
    PetscInt, intent(in) :: natural(:) !! Natural cell indices
    PetscInt :: local(size(natural))
    ! Locals:
    PetscInt :: n, idx(size(natural))
    PetscErrorCode :: ierr

    associate(num_cells => size(natural))
      idx = natural
      call AOApplicationToPetsc(ao, num_cells, idx, ierr); CHKERRQ(ierr)
      call ISGlobalToLocalMappingApplyBlock(l2g, IS_GTOLM_MASK, num_cells, &
           idx, n, local, ierr); CHKERRQ(ierr)
    end associate

  end function natural_to_local_cell_index_array

!------------------------------------------------------------------------

  function local_to_natural_cell_index_array(ao, l2g, local) &
       result(natural)
    !! Returns array of natural cell indices corresponding to array of
    !! local cell indices.

    AO, intent(in) :: ao !! Application ordering mapping natural to global cell indices
    ISLocalToGlobalMapping, intent(in) :: l2g !! DM local to global mapping
    PetscInt, intent(in) :: local(:) !! Local cell indices
    PetscInt :: natural(size(local))
    ! Locals:
    PetscInt :: idx(size(local))
    PetscErrorCode :: ierr

    associate(num_cells => size(local))
      call ISLocalToGlobalMappingApplyBlock(l2g, num_cells, local, idx, &
           ierr); CHKERRQ(ierr)
      call AOPetscToApplication(ao, num_cells, idx, ierr); CHKERRQ(ierr)
      natural = idx
    end associate

  end function local_to_natural_cell_index_array

!------------------------------------------------------------------------

  PetscInt function local_to_natural_cell_index_single(ao, &
       l2g, local) result(natural)
    !! Returns natural cell index corresponding to local cell index.

    AO, intent(in) :: ao !! Application ordering mapping natural to global cell indices
    ISLocalToGlobalMapping, intent(in) :: l2g !! DM local to global mapping
    PetscInt, intent(in) :: local !! Local cell index
    ! Locals:
    PetscInt :: idx(1), local_array(1)
    PetscErrorCode :: ierr

    local_array(1) = local
    call ISLocalToGlobalMappingApplyBlock(l2g, 1, local_array, idx, ierr)
    CHKERRQ(ierr)
    call AOPetscToApplication(ao, 1, idx, ierr); CHKERRQ(ierr)
    natural = idx(1)

  end function local_to_natural_cell_index_single

!------------------------------------------------------------------------

  IS function dm_natural_order_IS(dm, ao) result(natural_IS)
    !! Returns index set containing natural order for each cell,
    !! according to the specified natural-to-global AO.  Ghost cells
    !! are included but given natural order -1.

    DM, intent(in) :: dm
    AO, intent(in) :: ao
    ! Locals:
    ISLocalToGlobalMapping :: l2g
    DMLabel :: ghost_label
    PetscInt :: start_cell, end_cell, num_cells
    PetscInt :: c, ghost
    PetscInt, allocatable :: natural(:)
    PetscErrorCode :: ierr

    call DMGetLocalToGlobalMapping(dm, l2g, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, &
         ierr); CHKERRQ(ierr)
    num_cells = end_cell - start_cell
    call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    allocate(natural(0: num_cells - 1))

    do c = start_cell, end_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr)
       if (ghost < 0) then
          natural(c) = local_to_natural_cell_index(ao, l2g, c)
       else
          natural(c) = -1
       end if
    end do

    call ISCreateGeneral(PETSC_COMM_WORLD, num_cells, &
         natural, PETSC_COPY_VALUES, natural_IS, ierr)
    CHKERRQ(ierr)

  end function dm_natural_order_IS

!------------------------------------------------------------------------

  subroutine dm_get_cell_index(dm, ao, cell_index)
    !! Returns cell_index IS mapping natural cell indices to global
    !! indices of a global vector (without boundary data included).

    DM, intent(in) :: dm !! Input DM
    AO, intent(in) :: ao !! Natural-to-global AO for the DM
    IS, intent(out) :: cell_index !! Natural-to-global IS (without boundary data)
    ! Locals:
    ISLocalToGlobalMapping :: l2g
    DMLabel :: ghost_label
    PetscInt :: start_cell, end_cell, end_interior_cell
    PetscInt :: c, ic, ghost, carray(1)
    PetscInt :: num_ghost_cells, num_non_ghost_cells, bdy_cell_shift
    PetscInt, allocatable :: global(:), natural(:), global_interior(:)
    PetscInt, allocatable :: index_natural(:), index_global(:)
    AO :: ao_interior
    PetscErrorCode :: ierr

    call DMGetLocalToGlobalMapping(dm, l2g, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, &
         ierr); CHKERRQ(ierr)
    end_interior_cell = dm_get_end_interior_cell(dm, end_cell)
    num_ghost_cells = dm_get_num_partition_ghost_points(dm, 0)
    num_non_ghost_cells = end_interior_cell - start_cell - num_ghost_cells
    bdy_cell_shift = dm_get_bdy_cell_shift(dm)

    allocate(global(0: num_non_ghost_cells - 1))
    ic = 0
    call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    do c = start_cell, end_interior_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr)
       if (ghost < 0) then
          carray = c
          call ISLocalToGlobalMappingApplyBlock(l2g, 1, carray, &
               global(ic:ic), ierr); CHKERRQ(ierr)
          ic = ic + 1
       end if
    end do

    ! The index_natural array does not represent natural indices
    ! corresponding to local indices- it is just a set of natural
    ! indices owned by the current process, for the purpose of writing
    ! out the corresponding global indices:
    index_natural = global - bdy_cell_shift
    natural = global
    call AOPetscToApplication(ao, num_non_ghost_cells, natural, &
         ierr); CHKERRQ(ierr)
    global_interior = global - bdy_cell_shift
    call AOCreateMapping(PETSC_COMM_WORLD, num_non_ghost_cells, &
         natural, global_interior, ao_interior, ierr); CHKERRQ(ierr)
    deallocate(natural, global_interior)
    index_global = index_natural
    call AOApplicationToPetsc(ao_interior, num_non_ghost_cells, &
         index_global, ierr); CHKERRQ(ierr)
    call AODestroy(ao_interior, ierr); CHKERRQ(ierr)
    call ISCreateGeneral(PETSC_COMM_WORLD, num_non_ghost_cells, &
         index_global, PETSC_COPY_VALUES, cell_index, ierr)
    CHKERRQ(ierr)
    call PetscObjectSetName(cell_index, "cell_index", &
         ierr); CHKERRQ(ierr)

    deallocate(global, index_natural, index_global)

  end subroutine dm_get_cell_index

!------------------------------------------------------------------------

  subroutine dm_set_fv_adjacency(dm)
    !! Sets finite-volume adjacency for DM.

    DM, intent(in out) :: dm
    ! Locals:
    PetscErrorCode :: ierr

    call DMSetBasicAdjacency(dm, PETSC_TRUE, PETSC_FALSE, ierr)
    CHKERRQ(ierr)

  end subroutine dm_set_fv_adjacency

!------------------------------------------------------------------------

  subroutine create_path_dm(num_nodes, dm)
    !! Creates a 1-D DMPlex representing a path. Points 0 .. num_nodes
    !! -1 are the nodes; points num_nodes ... 2 * num_nodes - 2 are
    !! the edges.

    PetscInt, intent(in) :: num_nodes !! Number of nodes on the path
    DM, intent(out) :: dm !! Output DMPlex
    ! Locals:
    PetscInt :: num_edges, num_points, p
    PetscErrorCode :: ierr

    call DMPlexCreate(PETSC_COMM_WORLD, dm, ierr); CHKERRQ(ierr)
    call DMSetDimension(dm, 1, ierr); CHKERRQ(ierr)

    if (num_nodes > 0) then
       num_edges = num_nodes - 1
    else
       num_edges = 0
    end if
    num_points = num_nodes + num_edges

    call DMPlexSetChart(dm, 0, num_points, ierr)
    CHKERRQ(ierr)

    do p = num_nodes, num_points - 1
       call DMPlexSetConeSize(dm, p, 2, ierr); CHKERRQ(ierr)
    end do

    call DMSetUp(dm, ierr); CHKERRQ(ierr)

    do p = num_nodes, num_points - 1
       associate (p_node => p - num_nodes)
         call DMPlexSetCone(dm, p, [p_node, p_node + 1], ierr)
         CHKERRQ(ierr)
       end associate
    end do

    call DMPlexSymmetrize(dm, ierr); CHKERRQ(ierr)
    call DMPlexStratify(dm, ierr); CHKERRQ(ierr)

  end subroutine create_path_dm

!------------------------------------------------------------------------

  subroutine get_field_subvector(v, field, index_set, subv)
    !! Gets subvector of v for the specified field. Based on
    !! PetscSectionGetField_Internal().

    Vec, intent(in) :: v
    PetscInt, intent(in) :: field
    IS, intent(in out) :: index_set
    Vec, intent(in out) :: subv
    ! Locals:
    DM :: dm
    PetscSection :: section, global_section
    PetscInt :: pstart, pend, p, f, fc, fdof
    PetscInt :: start_cell, end_cell, end_interior_cell
    PetscInt :: num_components, gdof, poff, goff, suboff
    PetscInt :: subsize
    PetscInt, allocatable :: subindices(:)
    PetscErrorCode :: ierr

    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMGetSection(dm, section, ierr); CHKERRQ(ierr)
    call DMGetGlobalSection(dm, global_section, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr); CHKERRQ(ierr)
    end_interior_cell = dm_get_end_interior_cell(dm, end_cell)

    call PetscSectionGetChart(section, pstart, pend, ierr); CHKERRQ(ierr)
    ! Modify point range for cells- to exclude boundary ghosts:
    if ((end_cell > start_cell) .and. &
         (start_cell >= pstart) .and. (start_cell < pend)) then
       pend = end_interior_cell
    end if

    call PetscSectionGetFieldComponents(section, field, &
         num_components, ierr); CHKERRQ(ierr)
    subsize = 0
    do p = pstart, pend - 1
       call PetscSectionGetDof(global_section, p, gdof, ierr); CHKERRQ(ierr)
       if (gdof > 0) then
          call PetscSectionGetFieldDof(section, p, field, fdof, ierr)
          CHKERRQ(ierr)
          subsize = subsize + fdof
       end if
    end do
    allocate(subindices(0: subsize - 1))

    suboff = 0
    do p = pstart, pend - 1
       call PetscSectionGetDof(global_section, p, gdof, ierr); CHKERRQ(ierr)
       if (gdof > 0) then
          call PetscSectionGetOffset(global_section, p, goff, ierr)
          CHKERRQ(ierr)
          poff = 0
          do f = 0, field - 1
             call PetscSectionGetFieldDof(section, p, f, fdof, ierr)
             CHKERRQ(ierr)
             poff = poff + fdof
          end do
          call PetscSectionGetFieldDof(section, p, field, fdof, ierr)
          CHKERRQ(ierr)
          subindices(suboff: suboff + fdof - 1) = &
               goff + poff + [(fc, fc = 0, fdof - 1)]
          suboff = suboff + fdof
       end if
    end do

    call ISCreateGeneral(PETSC_COMM_WORLD, subsize, subindices, &
         PETSC_COPY_VALUES, index_set, ierr); CHKERRQ(ierr)
    deallocate(subindices)
    call VecGetSubVector(v, index_set, subv, ierr); CHKERRQ(ierr)
    call VecSetBlockSize(subv, num_components, ierr); CHKERRQ(ierr)

  end subroutine get_field_subvector

!------------------------------------------------------------------------

  subroutine section_get_field_names(section, lowercase, field_names)
    !! Returns array of section field names.

    use utils_module, only: str_to_lower

    PetscSection, intent(in) :: section
    PetscBool, intent(in) :: lowercase
    character(*), allocatable :: field_names(:)
    ! Locals:
    PetscInt :: num_fields, f
    PetscErrorCode :: ierr

    call PetscSectionGetNumFields(section, num_fields, ierr); CHKERRQ(ierr)
    allocate(field_names(num_fields))
    do f = 1, num_fields
       call PetscSectionGetFieldName(section, f - 1, field_names(f), ierr)
       if (lowercase) then
          field_names(f) = str_to_lower(field_names(f))
       end if
    end do

  end subroutine section_get_field_names

!------------------------------------------------------------------------

  PetscInt function dm_global_cell_field_dof(dm, field) result(dof)
    !! Returns degrees of freedom associated with specified DM section
    !! cell field. Checks dof on all processes in case current process
    !! has no cells.

    DM, intent(in) :: dm
    PetscInt, intent(in) :: field
    ! Locals:
    PetscSection :: section
    PetscInt :: start_cell, end_cell, field_dof
    PetscErrorCode :: ierr

    call DMGetSection(dm, section, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
    if (end_cell > start_cell) then
       call PetscSectionGetFieldDof(section, start_cell, field, &
            field_dof, ierr); CHKERRQ(ierr)
    else
       field_dof = -1
    end if
    call MPI_Allreduce(field_dof, dof, 1, MPI_INT, MPI_MAX, &
         PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)

  end function dm_global_cell_field_dof

!------------------------------------------------------------------------

  subroutine dm_check_create_label(dm, label_name)
    !! Creates label on DM with specified name. If a label of that
    !! name already exists it is removed and re-created.

    DM, intent(in out) :: dm
    character(*), intent(in) :: label_name
    ! Locals:
    PetscBool :: has_label
    DMLabel :: label
    PetscErrorCode :: ierr

    call DMHasLabel(dm, label_name, has_label, ierr); CHKERRQ(ierr)
    if (has_label) then
       call DMRemoveLabel(dm, label_name, label, ierr); CHKERRQ(ierr)
    end if
    call DMCreateLabel(dm, label_name, ierr); CHKERRQ(ierr)

  end subroutine dm_check_create_label

!------------------------------------------------------------------------

  subroutine dm_label_partition_ghosts(dm)
    !! Sets ghost (and VTK) label on partition ghost cells and
    !! faces. Based on code from PETSc DMPlexShiftLabels_Internal(),
    !! which is called from DMPlexConstructGhostCells().

    DM, intent(in out) :: dm
    ! Locals:
    PetscMPIInt :: rank
    PetscSF :: point_sf
    PetscInt :: start_cell, end_cell, start_face, end_face
    PetscInt :: c, l, f, va, vb
    PetscInt :: num_roots, num_leaves, num_cells
    PetscInt, pointer :: local(:), cells(:)
    type(PetscSFNode), pointer :: remote(:)
    DMLabel :: ghost_label, vtk_label
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call DMGetPointSF(dm, point_sf, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
    call PetscSFGetGraph(point_sf, num_roots, num_leaves, &
         local, remote, ierr); CHKERRQ(ierr)
    call dm_check_create_label(dm, "ghost")
    call dm_check_create_label(dm, "vtk")
    call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    call DMGetLabel(dm, "vtk", vtk_label, ierr); CHKERRQ(ierr)

    l = 0
    c = start_cell
    do while ((l < num_leaves) .and. (c < end_cell))
       do while ((c < local(l + 1)) .and. (c < end_cell))
          call DMLabelSetValue(vtk_label, c, 1, ierr); CHKERRQ(ierr)
          c = c + 1
       end do
       if (local(l + 1) >= end_cell) exit
       if (remote(l + 1)%rank == rank) then
          call DMLabelSetValue(vtk_label, c, 1, ierr); CHKERRQ(ierr)
       else ! partition ghost cells:
          call DMLabelSetValue(ghost_label, c, 2, ierr); CHKERRQ(ierr)
       end if
       l = l + 1
       c = c + 1
    end do

    do while (c < end_cell)
       call DMLabelSetValue(vtk_label, c, 1, ierr); CHKERRQ(ierr)
       c = c + 1
    end do

    ! Label ghost faces:
    call DMPlexGetHeightStratum(dm, 1, start_face, end_face, ierr)
    CHKERRQ(ierr)
    do f = start_face, end_face - 1
       call DMPlexGetSupportSize(dm, f, num_cells, ierr); CHKERRQ(ierr)
       if (num_cells < 2) then
          call DMLabelSetValue(ghost_label, f, 1, ierr); CHKERRQ(ierr)
       else
          call DMPlexGetSupport(dm, f, cells, ierr); CHKERRQ(ierr)
          call DMLabelGetValue(vtk_label, cells(1), va, ierr); CHKERRQ(ierr)
          call DMLabelGetValue(vtk_label, cells(2), vb, ierr); CHKERRQ(ierr)
          if ((va /= 1) .and. (vb /= 1)) then
             call DMLabelSetValue(ghost_label, f, 1, ierr); CHKERRQ(ierr)
          end if
       end if
    end do

  end subroutine dm_label_partition_ghosts

!------------------------------------------------------------------------

  subroutine dm_label_boundary_ghosts(dm, label_name)
    !! Labels boundary ghost cells with the specified label and value 1.

    DM, intent(in out) :: dm
    character(*), intent(in) :: label_name
    ! Locals:
    PetscInt :: start_cell, end_cell, end_interior_cell, dummy, c
    PetscErrorCode :: ierr

    call dm_check_create_label(dm, label_name)
    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    call DMPlexGetGhostCellStratum(dm, end_interior_cell, dummy, ierr)
    CHKERRQ(ierr)
    if (end_interior_cell >= 0) then
       do c = end_interior_cell, end_cell - 1
          call DMSetLabelValue(dm, label_name, c, 1, ierr)
          CHKERRQ(ierr)
       end do
    end if

  end subroutine dm_label_boundary_ghosts

!------------------------------------------------------------------------

  subroutine dm_distribute_local_vec(dm, sf, v)
    !! Distributes local vector v from its original DM to the
    !! specified one, via the supplied distribution star forest. The
    !! original vector v is overwritten.

    DM, intent(in) :: dm
    PetscSF, intent(in) :: sf
    Vec, intent(in out) :: v
    ! Locals:
    character(80) :: name
    DM :: v_dm, dist_v_dm
    PetscSection :: section, dist_section
    Vec :: dist_v
    PetscErrorCode :: ierr

    call VecGetDM(v, v_dm, ierr); CHKERRQ(ierr)
    call DMGetSection(v_dm, section, ierr); CHKERRQ(ierr)

    call VecCreate(PETSC_COMM_WORLD, dist_v, ierr); CHKERRQ(ierr)
    call PetscObjectGetName(v, name, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(dist_v, name, ierr); CHKERRQ(ierr)

    call DMClone(dm, dist_v_dm, ierr); CHKERRQ(ierr)
    call PetscSectionCreate(PETSC_COMM_WORLD, dist_section, ierr)
    CHKERRQ(ierr)
    call DMSetSection(dist_v_dm, dist_section, ierr)
    call VecSetDM(dist_v, dist_v_dm, ierr); CHKERRQ(ierr)

    call DMPlexDistributeField(v_dm, sf, section, v, &
         dist_section, dist_v, ierr); CHKERRQ(ierr)

    call PetscSectionDestroy(dist_section, ierr); CHKERRQ(ierr)
    call VecDestroy(v, ierr); CHKERRQ(ierr)
    v = dist_v

  end subroutine dm_distribute_local_vec

!------------------------------------------------------------------------

  subroutine dm_distribute_global_vec(dm, sf, v)
    !! Distributes global vector v from its original DM to the
    !! specified one, via the supplied distribution star forest. The
    !! original vector v is overwritten.

    DM, intent(in) :: dm
    PetscSF, intent(in) :: sf
    Vec, intent(in out) :: v
    ! Locals:
    DM :: v_dm, dist_v_dm
    Vec :: local_v, global_v
    character(80) :: name
    PetscErrorCode :: ierr

    call VecGetDM(v, v_dm, ierr); CHKERRQ(ierr)
    call DMCreateLocalVector(v_dm, local_v, ierr); CHKERRQ(ierr)
    call DMGlobalToLocal(v_dm, v, INSERT_VALUES, local_v, ierr); CHKERRQ(ierr)

    call dm_distribute_local_vec(dm, sf, local_v)

    call VecGetDM(local_v, dist_v_dm, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(dist_v_dm, global_v, ierr); CHKERRQ(ierr)
    call PetscObjectGetName(v, name, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(global_v, name, ierr); CHKERRQ(ierr)

    call DMLocalToGlobal(dist_v_dm, local_v, INSERT_VALUES, global_v, ierr)
    CHKERRQ(ierr)

    call VecDestroy(local_v, ierr); CHKERRQ(ierr)
    call VecDestroy(v, ierr); CHKERRQ(ierr)
    v = global_v

  end subroutine dm_distribute_global_vec

!------------------------------------------------------------------------

  subroutine dm_distribute_index_set(dm, sf, section, index_set)
    !! Distributes IS according to the specified distribution star
    !! forest.

    DM, intent(in) :: dm
    PetscSF, intent(in) :: sf !! Distribution star forest
    PetscSection, intent(in) :: section !! Section for existing IS
    IS, intent(in out) :: index_set
    ! Locals:
    PetscSection :: dist_section
    IS :: dist_index_set
    PetscErrorCode :: ierr

    call PetscSectionCreate(PETSC_COMM_WORLD, dist_section, ierr)
    CHKERRQ(ierr)
    call ISCreate(PETSC_COMM_WORLD, dist_index_set, ierr)
    CHKERRQ(ierr)
    call DMPlexDistributeFieldIS(dm, sf, section, &
         index_set, dist_section, &
         dist_index_set, ierr); CHKERRQ(ierr)
    call PetscSectionDestroy(dist_section, ierr); CHKERRQ(ierr)
    call ISDestroy(index_set, ierr); CHKERRQ(ierr)
    index_set = dist_index_set

  end subroutine dm_distribute_index_set

!------------------------------------------------------------------------

  subroutine vec_copy_common_local(v, w)
    !! Copies data from Vec v to w, up to the minimum of the local
    !! sizes of the two vectors.

    Vec, intent(in) :: v
    Vec, intent(in out) :: w
    ! Locals:
    PetscInt :: vsize, wsize, n
    PetscReal, pointer :: v_array(:), w_array(:)
    PetscErrorCode :: ierr

    call VecGetLocalSize(v, vsize, ierr); CHKERRQ(ierr)
    call VecGetLocalSize(w, wsize, ierr); CHKERRQ(ierr)
    n = min(vsize, wsize)

    call VecGetArrayReadF90(v, v_array, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(w, w_array, ierr); CHKERRQ(ierr)
    w_array(1:n) = v_array(1:n)
    call VecRestoreArrayF90(w, w_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(v, v_array, ierr); CHKERRQ(ierr)

  end subroutine vec_copy_common_local

!------------------------------------------------------------------------

  PetscBool function mat_type_is_block(M) result(isblock)
    !! Returns true if matrix type is a block type.

    Mat, intent(in) :: M
    ! Locals:
    MatType :: mat_type
    PetscErrorCode :: ierr

    call MatGetType(M, mat_type, ierr); CHKERRQ(ierr)
    select case (mat_type)
    case (MATBAIJ, MATSEQBAIJ, MATMPIBAIJ, MATSBAIJ, &
         MATSEQSBAIJ, MATMPISBAIJ)
       isblock = PETSC_TRUE
    case default
       isblock = PETSC_FALSE
    end select

  end function mat_type_is_block

!------------------------------------------------------------------------

  function mat_coloring_perturbed_columns(M, fd_coloring) result(columns)
    !! Returns perturbed columns from a finite difference matrix
    !! colouring, corrected if necessary to block column indices. For
    !! block-type matrices, the perturbed columns from
    !! MatFDColoringGetPerturbedColumnsF90() are already block column
    !! indices. But for other matrix types they are not.

    Mat, intent(in) :: M
    MatFDColoring, intent(in) :: fd_coloring
    PetscInt, allocatable :: columns(:)
    ! Locals:
    PetscInt :: bs
    PetscInt, pointer :: perturbed_columns(:)
    PetscErrorCode :: ierr

    call MatFDColoringGetPerturbedColumnsF90(fd_coloring, &
         perturbed_columns, ierr); CHKERRQ(ierr)
    if (mat_type_is_block(M)) then
       columns = perturbed_columns
    else
       call MatGetBlockSize(M, bs, ierr); CHKERRQ(ierr)
       columns = perturbed_columns / bs
    end if
    call MatFDColoringRestorePerturbedColumnsF90(fd_coloring, &
         perturbed_columns, ierr); CHKERRQ(ierr)

  end function mat_coloring_perturbed_columns

!------------------------------------------------------------------------

  subroutine dm_copy_cone_orientation(dm_source, p_source, dm_dest, p_dest)
    !! Copies DMPlex cone orientation from point p_source in dm_source to point
    !! p_dest in dm_dest.

    DM, intent(in) :: dm_source
    PetscInt, intent(in) :: p_source, p_dest
    DM, intent(in out) :: dm_dest
    ! Locals:
    PetscInt, pointer :: orientation(:)
    PetscErrorCode :: ierr

    call DMPlexGetConeOrientation(dm_source, p_source, orientation, &
         ierr); CHKERRQ(ierr)
    call DMPlexSetConeOrientation(dm_dest, p_dest, orientation, &
         ierr); CHKERRQ(ierr)
    call DMPlexRestoreConeOrientation(dm_source, p_source, orientation, &
         ierr); CHKERRQ(ierr)

  end subroutine dm_copy_cone_orientation

!------------------------------------------------------------------------

  subroutine dm_cell_counts(dm, cells_total, cells_min, cells_max)
    !! Returns total DM cell count over all processes (excluding ghost
    !! cells), and minimum and maximum count per process, on rank 0.

    DM, intent(in) :: dm
    PetscInt, intent(out) :: cells_total, cells_min, cells_max
    ! Locals:
    PetscInt :: start_cell, end_cell, end_interior_cell
    PetscInt :: num_ghost_cells, cells_local
    PetscErrorCode :: ierr

    call DMPlexGetHeightStratum(dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    end_interior_cell = dm_get_end_interior_cell(dm, end_cell)
    num_ghost_cells = dm_get_num_partition_ghost_points(dm, 0)
    cells_local = end_interior_cell - start_cell - num_ghost_cells

    call MPI_reduce(cells_local, cells_total, 1, MPI_INTEGER, &
           MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
    call MPI_reduce(cells_local, cells_min, 1, MPI_INTEGER, &
           MPI_MIN, 0, PETSC_COMM_WORLD, ierr)
    call MPI_reduce(cells_local, cells_max, 1, MPI_INTEGER, &
           MPI_MAX, 0, PETSC_COMM_WORLD, ierr)

  end subroutine dm_cell_counts

!------------------------------------------------------------------------

end module dm_utils_module
