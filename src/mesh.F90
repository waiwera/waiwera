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

module mesh_module
  !! Module for mesh type.

#include <petsc/finclude/petsc.h>

  use petsc
  use mpi_utils_module
  use fson
  use minc_module

  implicit none

  private

  PetscInt, parameter, public :: max_mesh_filename_length = 200
  character(len = 16), public :: open_boundary_label_name = "open_boundary" !! Name of DMLabel for identifying open boundaries
  character(len = 13), parameter, public :: &
       cell_order_label_name = "cell_order" !! Name of DMLabel for generating cell index set

  type, public :: mesh_type
     !! Mesh type.
     private
     character(max_mesh_filename_length), public :: filename !! Mesh file name
     DM, public :: original_dm, minc_dm
     DM, pointer, public :: dm !! DM representing the mesh topology
     Vec, public :: cell_geom !! Vector containing cell geometry data
     Vec, public :: face_geom !! Vector containing face geometry data
     PetscInt, public :: start_cell !! DM point containing first cell on this process
     PetscInt, public :: end_cell !! DM point one greater than last cell on this process
     PetscInt, public :: end_interior_cell !! DM point one greater than last interior (non-ghost) cell on this process
     PetscInt, public :: start_face !! DM point containing first face on this process
     PetscInt, public :: end_face !! DM point one greater than last face on this process
     PetscReal, allocatable, public :: bcs(:,:) !! Array containing boundary conditions
     IS, public :: cell_index !! Index set defining natural cell ordering
     PetscInt, public, allocatable :: ghost_cell(:), ghost_face(:) !! Ghost label values for cells and faces
     type(minc_type), allocatable, public :: minc(:)
     PetscReal, public :: permeability_rotation(3, 3) !! Rotation matrix of permeability axes
     PetscReal, public :: thickness !! Mesh thickness (for dimension < 3)
     PetscBool, public :: radial !! If mesh coordinate system is radial or Cartesian
     PetscBool, public :: has_minc !! If mesh has any MINC cells
   contains
     procedure :: assign_dm => mesh_assign_dm
     procedure :: setup_cell_order_label => mesh_setup_cell_order_label
     procedure :: setup_cell_index => mesh_setup_cell_index
     procedure :: distribute => mesh_distribute
     procedure :: construct_ghost_cells => mesh_construct_ghost_cells
     procedure :: setup_geometry => mesh_setup_geometry
     procedure :: setup_ghost_arrays => mesh_setup_ghost_arrays
     procedure :: get_bounds => mesh_get_bounds
     procedure :: setup_coordinate_parameters => mesh_setup_coordinate_parameters
     procedure :: set_boundary_face_distances => mesh_set_boundary_face_distances
     procedure :: set_permeability_rotation => mesh_set_permeability_rotation
     procedure :: modify_geometry => mesh_modify_geometry
     procedure :: override_face_properties => mesh_override_face_properties
     procedure :: setup_minc => mesh_setup_minc
     procedure :: setup_minc_dm => mesh_setup_minc_dm
     procedure, public :: init => mesh_init
     procedure, public :: configure => mesh_configure
     procedure, public :: setup_boundaries => mesh_setup_boundaries
     procedure, public :: set_boundary_values => mesh_set_boundary_values
     procedure, public :: order_vector => mesh_order_vector
     procedure, public :: destroy => mesh_destroy
  end type mesh_type

contains

!------------------------------------------------------------------------

  subroutine mesh_assign_dm(self, dm)
    !! Assigns self%dm pointer to specified DM.

    class(mesh_type), intent(in out) :: self
    DM, target, intent(in) :: dm

    self%dm => dm

  end subroutine mesh_assign_dm

!------------------------------------------------------------------------

  subroutine mesh_setup_cell_order_label(self)
    !! Sets up cell ordering label on mesh cells. This assumes the mesh
    !! has not yet been distributed.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscBool :: has_label
    PetscInt :: start_cell, end_cell, c
    PetscErrorCode :: ierr

    call DMHasLabel(self%original_dm, cell_order_label_name, has_label, &
         ierr); CHKERRQ(ierr)

    if (.not. has_label) then
       call DMCreateLabel(self%original_dm, cell_order_label_name, ierr)
       CHKERRQ(ierr)
       call DMPlexGetHeightStratum(self%original_dm, 0, start_cell, &
            end_cell, ierr); CHKERRQ(ierr)
       do c = start_cell, end_cell - 1
          call DMSetLabelValue(self%original_dm, cell_order_label_name, c, &
               c, ierr); CHKERRQ(ierr)
       end do
    end if

  end subroutine mesh_setup_cell_order_label

!------------------------------------------------------------------------

  subroutine mesh_setup_cell_index(self, viewer)
    !! Sets up cell index set from cell order label on DM.  This index
    !! set corresponds to a block size of 1.
    !! Also writes the cell interior index set to HDF5 output. This is
    !! used for post-processing and is similar to the cell index set
    !! except that the global indices apply to vectors containing only
    !! interior cells (not boundary ghost cells).

    class(mesh_type), intent(in out) :: self
    PetscViewer, intent(in out), optional :: viewer
    ! Locals:
    PetscInt :: total_count, local_count
    PetscInt :: start_cell, end_cell, end_interior_cell
    PetscInt :: cmax, fmax, emax, vmax
    PetscInt :: total_allocate_count, allocate_size
    PetscInt :: c, i, ghost, order
    DMLabel :: ghost_label, order_label
    PetscInt, allocatable :: global_index(:), natural_index(:)
    PetscInt, allocatable :: global_index_all(:), natural_index_all(:)
    PetscInt, allocatable :: index_array_all(:), index_array(:)
    PetscInt, allocatable :: local_counts(:), displacements(:)
    Vec :: v
    IS :: cell_interior_index
    PetscInt :: blocksize, start_global_index, i_global
    PetscMPIInt :: rank, num_procs
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD, num_procs, ierr)
    call DMGetLabel(self%original_dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(self%original_dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    call DMPlexGetHybridBounds(self%original_dm, cmax, fmax, emax, vmax, ierr)
    CHKERRQ(ierr)
    end_interior_cell = cmax

    ! Count interior cells:
    local_count = 0
    do c = start_cell, end_interior_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) local_count = local_count + 1
    end do
    allocate(global_index(local_count), natural_index(local_count))

    ! Get starting global index for each process:
    call DMGetGlobalVector(self%original_dm, v, ierr); CHKERRQ(ierr)
    call VecGetOwnershipRange(v, start_global_index, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    call VecGetBlockSize(v, blocksize, ierr); CHKERRQ(ierr)
    call DMRestoreGlobalVector(self%original_dm, v, ierr); CHKERRQ(ierr)
    start_global_index = start_global_index / blocksize

    ! Set up global and natural index arrays on each process:
    call DMGetLabel(self%original_dm, cell_order_label_name, order_label, &
         ierr); CHKERRQ(ierr)
    i = 1
    do c = start_cell, end_interior_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          call DMLabelGetValue(order_label, c, order, ierr); CHKERRQ(ierr)
          i_global = start_global_index + c - self%start_cell
          global_index(i) = i_global
          natural_index(i) = order
          i = i + 1
       end if
    end do

    ! Gather arrays to root process:
    if (rank == 0) then
       allocate_size = num_procs
    else ! have to allocate non-zero size, even if not actually used:
       allocate_size = 1
    end if
    allocate(local_counts(allocate_size), displacements(allocate_size))
    call MPI_gather(local_count, 1, MPI_INTEGER, local_counts, 1, &
         MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
       total_count = sum(local_counts)
       displacements(1) = 0
       do i = 2, num_procs
          displacements(i) = displacements(i-1) + local_counts(i-1)
       end do
       total_allocate_count = total_count
    else
       total_allocate_count = 1
    end if
    allocate(global_index_all(total_allocate_count), &
         natural_index_all(total_allocate_count))
    call MPI_gatherv(global_index, local_count, MPI_INTEGER, &
         global_index_all, local_counts, displacements, MPI_INTEGER, &
         0, PETSC_COMM_WORLD, ierr)
    call MPI_gatherv(natural_index, local_count, MPI_INTEGER, &
         natural_index_all, local_counts, displacements, MPI_INTEGER, &
         0, PETSC_COMM_WORLD, ierr)
    deallocate(global_index, natural_index)

    ! Set up index array on root process, and scatter:
    allocate(index_array_all(total_allocate_count))
    if (rank == 0) then
       do i = 1, total_count
          index_array_all(natural_index_all(i) + 1) = global_index_all(i)
       end do
    end if
    deallocate(global_index_all)
    allocate(index_array(local_count))
    call MPI_scatterv(index_array_all, local_counts, displacements, &
         MPI_INTEGER, index_array, local_count, MPI_INTEGER, &
         0, PETSC_COMM_WORLD, ierr)
    call ISCreateGeneral(PETSC_COMM_WORLD, local_count, index_array, &
         PETSC_COPY_VALUES, self%cell_index, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(self%cell_index, "cell_index", ierr)

    ! Set up cell interior index set:
    if (rank == 0) then
       do i = 1, total_count
          index_array_all(natural_index_all(i) + 1) = i - 1
       end do
    end if
    deallocate(natural_index_all)
    call MPI_scatterv(index_array_all, local_counts, displacements, &
         MPI_INTEGER, index_array, local_count, MPI_INTEGER, &
         0, PETSC_COMM_WORLD, ierr)
    deallocate(index_array_all, local_counts, displacements)

    call ISCreateGeneral(PETSC_COMM_WORLD, local_count, index_array, &
         PETSC_COPY_VALUES, cell_interior_index, ierr); CHKERRQ(ierr)
    deallocate(index_array)
    call PetscObjectSetName(cell_interior_index, &
         "cell_interior_index", ierr)

    if (present(viewer)) then
       call ISView(self%cell_index, viewer, ierr); CHKERRQ(ierr)
       call ISView(cell_interior_index, viewer, ierr); CHKERRQ(ierr)
    end if

    call ISDestroy(cell_interior_index, ierr); CHKERRQ(ierr)

  end subroutine mesh_setup_cell_index

!------------------------------------------------------------------------

  subroutine mesh_distribute(self)
    !! Distributes mesh over processors.
    
    class(mesh_type), intent(in out) :: self
    ! Locals:
    DM :: dist_dm
    PetscErrorCode :: ierr
    PetscInt, parameter :: overlap = 1

    call DMPlexDistribute(self%dm, overlap, PETSC_NULL_SF, &
         dist_dm, ierr); CHKERRQ(ierr)
    if (dist_dm .ne. PETSC_NULL_DM) then
       call DMDestroy(self%dm, ierr); CHKERRQ(ierr)
       self%dm = dist_dm
    end if
    
  end subroutine mesh_distribute

!------------------------------------------------------------------------

  subroutine mesh_construct_ghost_cells(self)
    !! Constructs ghost cells on open boundary faces.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    DM :: ghost_dm
    PetscErrorCode :: ierr

    call DMPlexConstructGhostCells(self%dm, open_boundary_label_name, &
         PETSC_NULL_INTEGER, ghost_dm, ierr); CHKERRQ(ierr)
    if (ghost_dm .ne. PETSC_NULL_DM) then
       call DMDestroy(self%dm, ierr); CHKERRQ(ierr);
       self%dm = ghost_dm
    end if

  end subroutine mesh_construct_ghost_cells

!------------------------------------------------------------------------

  subroutine mesh_setup_coordinate_parameters(self, json, logfile)
    !! Sets up mesh coordinate system parameters.

    use kinds_module
    use fson_mpi_module
    use logfile_module

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    PetscInt :: dim
    PetscErrorCode :: ierr
    PetscBool, parameter :: default_radial = PETSC_FALSE
    PetscReal, parameter :: default_thickness = 1._dp

    self%thickness = default_thickness
    call DMGetDimension(self%original_dm, dim, ierr); CHKERRQ(ierr)
    select case (dim)
    case(3)
       self%radial = PETSC_FALSE
    case(2)
       call fson_get_mpi(json, "mesh.radial", default_radial, self%radial, &
               logfile)
       if (.not. self%radial) then
          call fson_get_mpi(json, "mesh.thickness", default_thickness, &
               self%thickness, logfile)
       end if
    case default
       if (present(logfile)) then
          call logfile%write(LOG_LEVEL_ERR, 'mesh', 'init', &
               str_key = 'stop            ', &
               str_value = '1-D mesh not supported.', &
               rank = 0)
       end if
       stop
    end select

  end subroutine mesh_setup_coordinate_parameters

!------------------------------------------------------------------------

  subroutine mesh_set_permeability_rotation(self, json, logfile)
    !! Sets up rotation matrix for permeability axes.  Here the
    !! rotation is assumed to be only in the horizontal and defined by
    !! a single angle (anti-clockwise from x-axis). Input angle should
    !! be in degrees (converted here to radians).

    use kinds_module
    use fson_mpi_module
    use logfile_module
    use utils_module, only: degrees_to_radians, rotation_matrix_2d

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    type(fson_value), pointer :: mesh_json
    PetscReal :: angle
    PetscReal, parameter :: default_permeability_angle = 0._dp

    call fson_get_mpi(json, "mesh", mesh_json)
    call fson_get_mpi(mesh_json, "permeability_angle", &
         default_permeability_angle, angle, logfile)
    angle = degrees_to_radians(angle)

    self%permeability_rotation = 0._dp
    self%permeability_rotation(1:2, 1:2) = rotation_matrix_2d(angle)
    self%permeability_rotation(3, 3) = 1._dp

  end subroutine mesh_set_permeability_rotation

!------------------------------------------------------------------------

  subroutine mesh_setup_geometry(self, gravity)
    !! Sets up global vectors containing geometry data (e.g. cell volumes,
    !! cell centroids, face areas, face-to-centroid distances) for the mesh.

    class(mesh_type), intent(in out) :: self
    PetscReal, intent(in) :: gravity(:)
    ! Locals:
    PetscErrorCode :: ierr
    Vec :: petsc_face_geom

    call DMPlexComputeGeometryFVM(self%original_dm, self%cell_geom, &
         petsc_face_geom, ierr); CHKERRQ(ierr)

    call self%modify_geometry(petsc_face_geom, gravity)
    call self%set_boundary_face_distances()

    call PetscObjectSetName(self%cell_geom, "cell_geometry", ierr)
    CHKERRQ(ierr)
    call PetscObjectSetName(self%face_geom, "face_geometry", ierr)
    CHKERRQ(ierr)

  end subroutine mesh_setup_geometry

!------------------------------------------------------------------------

  subroutine mesh_modify_geometry(self, petsc_face_geom, gravity)
    !! Modifies cell and face geometry vectors, given the face
    !! geometry vector produced by DMPlexComputeGeometryFVM(). For 2D
    !! Cartesian meshes, cell volumes and face areas are modified to
    !! account for the mesh thickness. For radial meshes, cell volumes
    !! and face areas are computed for the solids of revolution of
    !! cells and faces in the 2D mesh, using Pappus' centroid
    !! theorem. In both cases, additional face geometry parameters are
    !! computed (e.g. distances, gravity normal).

    use kinds_module
    use cell_module
    use face_module
    use utils_module, only: rotation_matrix_2d
    use dm_utils_module, only: section_offset, local_vec_section, set_dm_data_layout

    class(mesh_type), intent(in out) :: self
    Vec, intent(in) :: petsc_face_geom
    PetscReal, intent(in) :: gravity(:)
    ! Locals:
    DM :: dm_face
    PetscSection :: face_section, petsc_face_section, cell_section
    PetscInt :: c, f, ghost_cell, ghost_face, i
    PetscInt :: start_cell, end_cell, start_face, end_face
    PetscInt :: face_offset, petsc_face_offset
    PetscInt :: cell_offset(2), offset
    type(cell_type) :: cell
    type(face_type) :: face
    type(petsc_face_type) :: petsc_face
    PetscReal, pointer, contiguous :: face_geom_array(:), petsc_face_geom_array(:)
    PetscReal, pointer, contiguous :: cell_geom_array(:)
    DMLabel :: ghost_label
    PetscInt, pointer :: cells(:)
    PetscInt :: dim, face_variable_dim(num_face_variables)
    PetscErrorCode :: ierr

    interface

       subroutine modify_cell_volume_routine(cell)
         import :: cell_type
         type(cell_type), intent(in out) :: cell
       end subroutine modify_cell_volume_routine

       subroutine modify_face_area_routine(face)
         import :: face_type
         type(face_type), intent(in out) :: face
       end subroutine modify_face_area_routine

    end interface
    procedure(modify_cell_volume_routine), pointer :: modify_cell_volume
    procedure(modify_face_area_routine), pointer :: modify_face_area

    call DMGetDimension(self%original_dm, dim, ierr); CHKERRQ(ierr)
    call DMGetLabel(self%original_dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(self%original_dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    call DMPlexGetHeightStratum(self%original_dm, 1, start_face, end_face, ierr)
    CHKERRQ(ierr)

    ! Set up cell geometry vector:
    call local_vec_section(self%cell_geom, cell_section)
    call VecGetArrayF90(self%cell_geom, cell_geom_array, ierr); CHKERRQ(ierr)

    if (dim == 2) then
       ! Adjust cell volumes:
       if (self%radial) then
          modify_cell_volume => modify_cell_volume_2d_radial
       else
          modify_cell_volume => modify_cell_volume_2d_cartesian
       end if
       do c = start_cell, end_cell - 1
          call DMLabelGetValue(ghost_label, c, ghost_cell, ierr); CHKERRQ(ierr)
          if (ghost_cell < 0) then
             call section_offset(cell_section, c, offset, ierr); CHKERRQ(ierr)
             call cell%assign_geometry(cell_geom_array, offset)
             call modify_cell_volume(cell)
          end if
       end do
    end if

    ! Set up face geometry vector:
    call DMClone(self%original_dm, dm_face, ierr); CHKERRQ(ierr)
    face_variable_dim = dim - 1
    call set_dm_data_layout(dm_face, face_variable_num_components, &
         face_variable_dim, face_variable_names)
    call DMCreateLocalVector(dm_face, self%face_geom, ierr); CHKERRQ(ierr)
    call local_vec_section(self%face_geom, face_section)
    call VecGetArrayF90(self%face_geom, face_geom_array, ierr); CHKERRQ(ierr)
    call local_vec_section(petsc_face_geom, petsc_face_section)
    call VecGetArrayF90(petsc_face_geom, petsc_face_geom_array, ierr)
    CHKERRQ(ierr)
    call face%init()

    select case (dim)
    case (3)
       modify_face_area => modify_face_area_null
    case (2)
       if (self%radial) then
          modify_face_area => modify_face_area_2d_radial
       else
          modify_face_area => modify_face_area_2d_cartesian
       end if
    end select

    do f = start_face, end_face - 1

       call DMLabelGetValue(ghost_label, f, ghost_face, ierr); CHKERRQ(ierr)

       if (ghost_face < 0) then

          call section_offset(face_section, f, face_offset, ierr); CHKERRQ(ierr)
          call section_offset(petsc_face_section, f, petsc_face_offset, ierr)
          CHKERRQ(ierr)

          call DMPlexGetSupport(self%original_dm, f, cells, ierr); CHKERRQ(ierr)
          do i = 1, 2
             call section_offset(cell_section, cells(i), cell_offset(i), ierr)
             CHKERRQ(ierr)
          end do

          call petsc_face%assign_geometry(petsc_face_geom_array, petsc_face_offset)
          call face%assign_geometry(face_geom_array, face_offset)
          call face%assign_cell_geometry(cell_geom_array, cell_offset)

          face%centroid = petsc_face%centroid
          face%area = norm2(petsc_face%area_normal)
          face%normal = petsc_face%area_normal / face%area
          face%gravity_normal = dot_product(gravity, face%normal)
          call modify_face_area(face)
          do i = 1, 2
             face%distance(i) = norm2(face%centroid - face%cell(i)%centroid)
          end do
          call face%calculate_permeability_direction(self%permeability_rotation)

       end if

    end do

    call face%destroy()
    call petsc_face%destroy()
    call VecRestoreArrayF90(self%cell_geom, cell_geom_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(self%face_geom, face_geom_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(petsc_face_geom, petsc_face_geom_array, ierr)
    CHKERRQ(ierr)

    call PetscSectionDestroy(face_section, ierr); CHKERRQ(ierr)
    call DMDestroy(dm_face, ierr); CHKERRQ(ierr)

  contains

    subroutine modify_cell_volume_2d_cartesian(cell)
      ! Volume modification for 2D Cartesian cells.
      type(cell_type), intent(in out) :: cell

      cell%volume = cell%volume * self%thickness

    end subroutine modify_cell_volume_2d_cartesian

!........................................................................

    subroutine modify_cell_volume_2d_radial(cell)
      ! Volume modification for 2D radial cells- via Pappus' centroid
      ! theorem.
      use utils_module, only: pi
      type(cell_type), intent(in out) :: cell
      ! Locals:
      PetscReal :: r

      r = cell%centroid(1)
      cell%volume = cell%volume * 2._dp * pi * r

    end subroutine modify_cell_volume_2d_radial

!........................................................................

    subroutine modify_face_area_null(face)
      ! Do-nothing area modification- for 3D cells.
      type(face_type), intent(in out) :: face

      continue

    end subroutine modify_face_area_null

!........................................................................

    subroutine modify_face_area_2d_cartesian(face)
      ! Area modification for 2D Cartesian faces.
      type(face_type), intent(in out) :: face

      face%area = face%area * self%thickness

    end subroutine modify_face_area_2d_cartesian

!........................................................................

    subroutine modify_face_area_2d_radial(face)
      ! Area modification for 2D radial faces- via Pappus' centroid
      ! theorem.
      use utils_module, only: pi
      type(face_type), intent(in out) :: face
      ! Locals:
      PetscReal :: r

      r = face%centroid(1)
      face%area = face%area * 2._dp * pi * r

    end subroutine modify_face_area_2d_radial

  end subroutine mesh_modify_geometry

!------------------------------------------------------------------------

  subroutine mesh_setup_ghost_arrays(self)
    !! Sets up arrays of ghost label values on cells and faces. This
    !! is just for faster access to these values in the rest of the
    !! code.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscInt :: c, f
    DMLabel :: ghost_label
    PetscErrorCode :: ierr

    allocate(self%ghost_cell(self%start_cell: self%end_cell - 1))
    allocate(self%ghost_face(self%start_face: self%end_face - 1))

    call DMGetLabel(self%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    do c = self%start_cell, self%end_cell - 1
       call DMLabelGetValue(ghost_label, c, self%ghost_cell(c), ierr)
       CHKERRQ(ierr)
    end do

    do f = self%start_face, self%end_face - 1
       call DMLabelGetValue(ghost_label, f, self%ghost_face(f), ierr)
       CHKERRQ(ierr)
    end do

  end subroutine mesh_setup_ghost_arrays

!------------------------------------------------------------------------

  subroutine mesh_get_bounds(self)
    !! Gets cell and face bounds  on current processor.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    PetscInt :: cmax, fmax, emax, vmax

    call DMPlexGetHybridBounds(self%dm, cmax, fmax, emax, vmax, ierr)
    CHKERRQ(ierr)
    self%end_interior_cell = cmax

    call DMPlexGetHeightStratum(self%dm, 0, self%start_cell, &
         self%end_cell, ierr)
    CHKERRQ(ierr)

    call DMPlexGetHeightStratum(self%dm, 1, self%start_face, &
         self%end_face, ierr); CHKERRQ(ierr)

  end subroutine mesh_get_bounds

!------------------------------------------------------------------------

  subroutine mesh_init(self, json, logfile)
    !! Initializes mesh, reading filename from JSON input file.
    !! If the filename is not present, an error is raised.
    !! Otherwise, the PETSc DM is read in.

    use logfile_module
    use fson_mpi_module
    use fson_value_m, only: TYPE_STRING, TYPE_OBJECT
    use dm_utils_module, only: dm_set_fv_adjacency

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscInt :: mesh_type
    PetscErrorCode :: ierr

    if (fson_has_mpi(json, "mesh")) then
       mesh_type = fson_type_mpi(json, "mesh")
       select case (mesh_type)
       case(TYPE_STRING)
          call fson_get_mpi(json, "mesh", val = self%filename)
       case (TYPE_OBJECT)
          call fson_get_mpi(json, "mesh.filename", "", self%filename)
       case default
          self%filename = ""
       end select
    end if

    if (self%filename == "") then
       if (present(logfile)) then
          call logfile%write(LOG_LEVEL_ERR, 'mesh', 'init', &
               str_key = 'stop            ', &
               str_value = 'filename not found in input.', &
               rank = 0)
       end if
       stop
    else
       ! Read in original DM:
       call DMPlexCreateFromFile(PETSC_COMM_WORLD, self%filename, PETSC_TRUE, &
            self%original_dm, ierr); CHKERRQ(ierr)
       call dm_set_fv_adjacency(self%original_dm)
       call self%assign_dm(self%original_dm)
       call self%setup_coordinate_parameters(json, logfile)
       call self%set_permeability_rotation(json, logfile)
       self%has_minc = PETSC_FALSE
    end if

  end subroutine mesh_init

!------------------------------------------------------------------------

  subroutine mesh_configure(self, primary_variable_names, gravity, viewer)
    !! Configures mesh, including distribution over processes and
    !! construction of ghost cells, setup of data layout, geometry and
    !! cell index set.

    use dm_utils_module, only: dm_setup_fv_discretization, &
         set_dm_default_data_layout

    class(mesh_type), intent(in out) :: self
    character(*), intent(in) :: primary_variable_names(:) !! Names of primary thermodynamic variables
    PetscReal, intent(in) :: gravity(:)
    PetscViewer, intent(in out), optional :: viewer !! PetscViewer for output of cell index set to HDF5 file
    ! Locals:
    PetscInt :: dof

    dof = size(primary_variable_names)

    call self%setup_cell_order_label()
    call self%distribute()
    call self%construct_ghost_cells()
    call dm_setup_fv_discretization(self%original_dm, dof)
    call set_dm_default_data_layout(self%original_dm, dof)

    if (self%has_minc) then
       call self%setup_minc_dm(dof)
    end if

    call self%get_bounds()

    call self%setup_geometry(gravity)

    call self%setup_ghost_arrays()

    call self%setup_cell_index(viewer)

  end subroutine mesh_configure

!------------------------------------------------------------------------

  subroutine mesh_destroy(self)
    !! Destroys mesh object.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscInt :: i
    PetscErrorCode :: ierr
    
    call VecDestroy(self%face_geom, ierr); CHKERRQ(ierr)
    self%dm => null()
    call DMDestroy(self%original_dm, ierr); CHKERRQ(ierr)

    if (allocated(self%bcs)) then
       deallocate(self%bcs)
    end if

    call ISDestroy(self%cell_index, ierr); CHKERRQ(ierr)

    if (allocated(self%ghost_cell)) then
       deallocate(self%ghost_cell)
    end if
    if (allocated(self%ghost_face)) then
       deallocate(self%ghost_face)
    end if

    if (self%has_minc) then
       call DMDestroy(self%minc_dm, ierr); CHKERRQ(ierr)
       do i = 1, size(self%minc)
          call self%minc(i)%destroy()
       end do
       deallocate(self%minc)
    end if

  end subroutine mesh_destroy

!------------------------------------------------------------------------

  subroutine mesh_set_boundary_face_distances(self)
    !! Sets face distances from boundary ghost cells to zero.

    use kinds_module
    use face_module
    use dm_utils_module, only: local_vec_section, section_offset

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscInt :: ibdy, num_faces
    PetscInt :: face_offset, f, iface
    IS :: bdy_IS
    PetscInt, pointer :: bdy_faces(:)
    type(face_type) :: face
    PetscSection :: face_section
    PetscReal, pointer, contiguous :: face_geom_array(:)
    PetscErrorCode :: ierr

    call local_vec_section(self%face_geom, face_section)
    call VecGetArrayF90(self%face_geom, face_geom_array, ierr); CHKERRQ(ierr)
    call face%init()

    if (allocated(self%bcs)) then
       ! Set external boundary face connection distances to zero:
       do ibdy = 1, size(self%bcs, 2)
          call DMGetStratumSize(self%dm, open_boundary_label_name, ibdy, &
               num_faces, ierr); CHKERRQ(ierr)
          if (num_faces > 0) then
             call DMGetStratumIS(self%dm, open_boundary_label_name, ibdy, &
                  bdy_IS, ierr); CHKERRQ(ierr)
             call ISGetIndicesF90(bdy_IS, bdy_faces, ierr); CHKERRQ(ierr)
             do iface = 1, num_faces
                f = bdy_faces(iface)
                call section_offset(face_section, f, face_offset, ierr)
                CHKERRQ(ierr)
                call face%assign_geometry(face_geom_array, face_offset)
                face%distance(2) = 0._dp
             end do
             call ISRestoreIndicesF90(bdy_IS, bdy_faces, ierr); CHKERRQ(ierr)
             call ISDestroy(bdy_IS, ierr); CHKERRQ(ierr)
          end if
       end do
    end if

    call face%destroy()
    call VecRestoreArrayF90(self%face_geom, face_geom_array, ierr); CHKERRQ(ierr)

  end subroutine mesh_set_boundary_face_distances

!------------------------------------------------------------------------

  subroutine mesh_setup_boundaries(self, json, eos, logfile)
    !! Sets up boundary conditions on the mesh.

    use kinds_module
    use eos_module, only: eos_type
    use logfile_module
    use fson
    use fson_value_m, only : TYPE_ARRAY, TYPE_OBJECT, TYPE_INTEGER
    use fson_mpi_module
    use dm_utils_module, only: dm_cell_normal_face

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input file
    class(eos_type), intent(in) :: eos !! EOS object
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    PetscErrorCode :: ierr
    PetscBool :: mesh_has_label
    type(fson_value), pointer :: boundaries, bdy, faces_json, face_json
    type(fson_value), pointer :: cell_normals, cell_normal, item
    PetscInt :: faces_type, face1_type
    PetscInt :: num_boundaries, num_faces, num_cells, ibdy
    PetscInt :: iface, icell, f, np, i, offset
    PetscInt, allocatable :: default_faces(:), default_cells(:)
    PetscInt, allocatable :: faces(:), cells(:)
    PetscInt :: region, cell, normal_len, num_face_items
    PetscReal, allocatable :: primary(:), input_normal(:)
    PetscReal :: normal(3)
    PetscReal, parameter :: default_normal(3) = [0._dp, 0._dp, 1._dp]
    character(len=64) :: bdystr
    character(len=12) :: istr

    default_faces = [PetscInt::] ! empty integer array
    default_cells = [PetscInt::]
    np = eos%num_primary_variables

    call DMHasLabel(self%original_dm, open_boundary_label_name, mesh_has_label, &
         ierr); CHKERRQ(ierr)
    if (.not. mesh_has_label) then
       call DMCreateLabel(self%original_dm, open_boundary_label_name, &
            ierr); CHKERRQ(ierr)
    end if

    if (fson_has_mpi(json, "boundaries")) then

       call fson_get_mpi(json, "boundaries", boundaries)
       num_boundaries = fson_value_count_mpi(boundaries, ".")
       allocate(self%bcs(np + 1, num_boundaries))

       num_faces = 0
       do ibdy = 1, num_boundaries
          write(istr, '(i0)') ibdy - 1
          bdystr = 'boundaries[' // trim(istr) // ']'
          bdy => fson_value_get_mpi(boundaries, ibdy)

          if (fson_has_mpi(bdy, "faces")) then
             call fson_get_mpi(bdy, "faces", faces_json)
             faces_type = fson_type_mpi(faces_json, ".")
             select case (faces_type)
             case (TYPE_ARRAY)
                num_face_items = fson_value_count_mpi(faces_json, ".")
                if (num_face_items > 0) then
                   face_json => fson_value_get_mpi(faces_json, 1)
                   face1_type = fson_type_mpi(face_json, ".")
                   select case (face1_type)
                   case (TYPE_INTEGER)
                      call fson_get_mpi(faces_json, ".", default_faces, faces, &
                           logfile, log_key = trim(bdystr) // ".faces")
                      num_faces = size(faces)
                   case (TYPE_OBJECT)
                      num_faces = 0
                      do i = 1, num_face_items
                         face_json => fson_value_get_mpi(faces_json, i)
                         num_cells = fson_value_count_mpi(face_json, "cells")
                         num_faces = num_faces + num_cells
                      end do
                      allocate(faces(num_faces))
                      offset = 0
                      do i = 1, num_face_items
                         face_json => fson_value_get_mpi(faces_json, i)
                         call fson_get_mpi(face_json, "cells", default_cells, cells, &
                              logfile, log_key = trim(bdystr) // "faces.cells")
                         call fson_get_mpi(face_json, "normal", default_normal, &
                              input_normal, logfile, log_key = trim(bdystr) // "faces.normal")
                         num_cells = size(cells)
                         call get_cell_faces(cells, num_cells, input_normal, offset)
                         offset = offset + num_cells
                      end do
                   case default
                      if (present(logfile)) then
                         call logfile%write(LOG_LEVEL_WARN, "input", &
                              "unrecognised_face_type")
                      end if
                   end select
                end if
             case (TYPE_OBJECT)
                call fson_get_mpi(faces_json, "cells", default_cells, cells, &
                     logfile, log_key = trim(bdystr) // "faces.cells")
                call fson_get_mpi(faces_json, "normal", default_normal, &
                     input_normal, logfile, log_key = trim(bdystr) // "faces.normal")
                num_cells = size(cells)
                num_faces = num_cells
                allocate(faces(num_faces))
                call get_cell_faces(cells, num_cells, input_normal, 0)
             case default
                if (present(logfile)) then
                   call logfile%write(LOG_LEVEL_WARN, "input", &
                        "unrecognised_faces_type")
                end if
             end select

          else if (fson_has_mpi(bdy, "cell_normals")) then
             call fson_get_mpi(bdy, "cell_normals", cell_normals)
             num_faces = fson_value_count_mpi(cell_normals, ".")
             allocate(faces(num_faces))
             do iface = 1, num_faces
                cell_normal => fson_value_get_mpi(cell_normals, iface)
                item => fson_value_get_mpi(cell_normal, 1)
                call fson_get_mpi(item, ".", val = cell)
                item => fson_value_get_mpi(cell_normal, 2)
                call fson_get_mpi(item, ".", val = input_normal)
                normal_len = size(input_normal)
                normal = 0._dp
                normal(1: normal_len) = input_normal
                call dm_cell_normal_face(self%original_dm, cell, normal, f)
                if (f >= 0) then
                   faces(iface) = f
                else
                   if (present(logfile)) then
                      call logfile%write(LOG_LEVEL_WARN, "input", &
                           "faces_not_found", int_keys = ["boundary"], &
                           int_values = [ibdy - 1])
                   end if
                   faces(iface) = -1
                end if
             end do
          end if

          do iface = 1, num_faces
             f = faces(iface)
             if (f >= 0) then
                call DMSetLabelValue(self%original_dm, open_boundary_label_name, &
                     f, ibdy, ierr); CHKERRQ(ierr)
             end if
          end do
          if (allocated(faces)) then
             deallocate(faces)
          end if

          call fson_get_mpi(bdy, "primary", eos%default_primary, &
               primary, logfile, log_key = trim(bdystr) // ".primary")
          call fson_get_mpi(bdy, "region", eos%default_region, &
               region, logfile, log_key = trim(bdystr) // ".region")
          self%bcs(1, ibdy) = dble(region)
          self%bcs(2 : np + 1, ibdy) = primary(1 : np)
       end do
    else if (present(logfile)) then
       call logfile%write(LOG_LEVEL_WARN, "input", "no_boundary_conditions")
    end if

  contains

    subroutine get_cell_faces(cells, num_cells, input_normal, offset)
      ! Get faces for normal vector and specified cells.

      PetscInt, intent(in) :: cells(:), num_cells
      PetscReal, intent(in) :: input_normal(:)
      PetscInt, intent(in) :: offset

      normal_len = size(input_normal)
      normal = 0._dp
      normal(1: normal_len) = input_normal

      do icell = 1, num_cells
         iface = offset + icell
         cell = cells(icell)
         call dm_cell_normal_face(self%original_dm, cell, normal, f)
         if (f >= 0) then
            faces(iface) = f
         else
            if (present(logfile)) then
               call logfile%write(LOG_LEVEL_WARN, "input", &
                    "faces_not_found", int_keys = ["boundary"], &
                    int_values = [ibdy - 1])
            end if
            faces(iface) = -1
         end if
      end do

    end subroutine get_cell_faces

  end subroutine mesh_setup_boundaries

!------------------------------------------------------------------------

  subroutine mesh_set_boundary_values(self, y, fluid_vector, rock_vector, &
       eos, y_range_start, fluid_range_start, rock_range_start)
    !! Sets primary variables (and rock properties) in boundary ghost cells.

    use dm_utils_module, only: global_vec_section, global_section_offset
    use eos_module, only: eos_type
    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(mesh_type), intent(in) :: self
    Vec, intent(in out) :: y !! Primary variables vector
    Vec, intent(in out) :: fluid_vector !! Fluid properties vector
    Vec, intent(in out) :: rock_vector !! Rock properties vector
    class(eos_type), intent(in) :: eos !! EOS module
    PetscInt, intent(in) :: y_range_start !! Start of range for global primary variables vector
    PetscInt, intent(in) :: fluid_range_start !! Start of range for global fluid vector
    PetscInt, intent(in) :: rock_range_start !! Start of range for global rock vector
    ! Locals:
    PetscInt :: ibdy, f, i, num_faces, iface, np, n
    PetscReal, pointer, contiguous :: y_array(:), fluid_array(:), rock_array(:)
    PetscReal, pointer, contiguous :: cell_primary(:), rock1(:), rock2(:)
    PetscSection :: y_section, fluid_section, rock_section
    IS :: bdy_IS
    DMLabel :: ghost_label
    type(fluid_type):: fluid
    type(rock_type) :: rock
    PetscInt :: y_offset, fluid_offset, rock_offsets(2), ghost, num_boundaries
    PetscInt, pointer :: bdy_faces(:), cells(:)
    PetscErrorCode :: ierr

    call global_vec_section(y, y_section)
    call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call global_vec_section(rock_vector, rock_section)
    call VecGetArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
    call DMGetLabel(self%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    np = eos%num_primary_variables
    call fluid%init(eos%num_components, eos%num_phases)
    call rock%init()

    if (allocated(self%bcs)) then
       num_boundaries = size(self%bcs, 2)
       do ibdy = 1, num_boundaries
          call DMGetStratumSize(self%dm, open_boundary_label_name, &
               ibdy, num_faces, ierr); CHKERRQ(ierr)
          if (num_faces > 0) then
             call DMGetStratumIS(self%dm, open_boundary_label_name, &
                  ibdy, bdy_IS, ierr); CHKERRQ(ierr)
             call ISGetIndicesF90(bdy_IS, bdy_faces, ierr); CHKERRQ(ierr)
             do iface = 1, num_faces
                f = bdy_faces(iface)
                call DMPlexGetSupport(self%dm, f, cells, ierr); CHKERRQ(ierr)
                if (size(cells) == 2) then
                   call DMLabelGetValue(ghost_label, cells(1), ghost, ierr)
                   CHKERRQ(ierr)
                   if (ghost < 0) then
                      call global_section_offset(y_section, cells(2), &
                           y_range_start, y_offset, ierr); CHKERRQ(ierr)
                      call global_section_offset(fluid_section, cells(2), &
                           fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)
                      do i = 1, 2
                         call global_section_offset(rock_section, cells(i), &
                              rock_range_start, rock_offsets(i), ierr)
                         CHKERRQ(ierr)
                      end do
                      ! Set primary variables and region:
                      cell_primary => y_array(y_offset : y_offset + np - 1)
                      call fluid%assign(fluid_array, fluid_offset)
                      cell_primary = self%bcs(2: np + 1, ibdy)
                      fluid%region = self%bcs(1, ibdy)
                      ! Copy rock type data from interior cell to boundary ghost cell:
                      n = rock%dof - 1
                      rock1 => rock_array(rock_offsets(1) : rock_offsets(1) + n)
                      rock2 => rock_array(rock_offsets(2) : rock_offsets(2) + n)
                      rock2 = rock1
                   end if
                end if
             end do
             call ISRestoreIndicesF90(bdy_IS, bdy_faces, ierr); CHKERRQ(ierr)
             call ISDestroy(bdy_IS, ierr); CHKERRQ(ierr)
          end if
       end do
    end if

    call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
    call fluid%destroy()
    call rock%destroy()

  end subroutine mesh_set_boundary_values

!------------------------------------------------------------------------

  subroutine mesh_order_vector(self, v, index)
    !! Reorders vector v to correspond to the cell order index of the mesh
    !! DM, rather than that of the given order index set.

    use dm_utils_module, only: vec_reorder

    class(mesh_type), intent(in) :: self
    Vec, intent(in out) :: v !! Global vector to re-order
    IS, intent(in), optional :: index !! Cell order index set for original ordering of v

    call vec_reorder(v, index, self%cell_index)

  end subroutine mesh_order_vector

!------------------------------------------------------------------------

  subroutine mesh_override_face_properties(self, json, logfile)
    !! Sets any face properties overridden in the JSON input-
    !! currently just face permeability directions.

    use fson_mpi_module
    use logfile_module
    use face_module
    use dm_utils_module, only: local_vec_section, section_offset

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    type(fson_value), pointer :: faces_json, face_json
    PetscInt :: num_cells, num_faces, iface, f, i, num_cell_faces
    PetscInt, allocatable :: global_cell_indices(:)
    PetscInt, allocatable :: default_cells(:)
    PetscInt, target :: points(2)
    PetscInt :: permeability_direction, face_offset, num_matching
    PetscSection :: face_section
    PetscInt, pointer :: ppoints(:), cells(:)
    PetscReal, pointer, contiguous :: face_geom_array(:)
    PetscInt, pointer :: cell_faces(:)
    type(face_type) :: face
    character(len=64) :: facestr
    character(len=12) :: istr
    IS :: cell_IS
    PetscErrorCode :: ierr
    PetscInt, parameter :: default_permeability_direction = 1

    default_cells = [PetscInt::] ! empty integer array
    call local_vec_section(self%face_geom, face_section)
    call VecGetArrayF90(self%face_geom, face_geom_array, ierr); CHKERRQ(ierr)
    call face%init()

    if (fson_has_mpi(json, "mesh.faces")) then
       call fson_get_mpi(json, "mesh.faces", faces_json)
       num_faces = fson_value_count_mpi(faces_json, ".")

       do iface = 1, num_faces

          write(istr, '(i0)') iface - 1
          facestr = 'mesh.faces[' // trim(istr) // ']'
          face_json => fson_value_get_mpi(faces_json, iface)
          call fson_get_mpi(face_json, "cells", default_cells, &
               global_cell_indices, logfile, log_key = trim(facestr) // ".cells")
          call fson_get_mpi(face_json, "permeability_direction", &
               default_permeability_direction, permeability_direction, &
               logfile, log_key = trim(facestr) // ".permeability_direction")

          num_cells = size(global_cell_indices)
          if (num_cells == 2) then

             ! get DM mesh points on local processor for both cells:
             call DMGetStratumSize(self%dm, cell_order_label_name, &
                  global_cell_indices(1), num_matching, ierr); CHKERRQ(ierr)
             if (num_matching == 1) then
                call DMGetStratumIS(self%dm, cell_order_label_name, &
                     global_cell_indices(1), cell_IS, ierr); CHKERRQ(ierr)
                call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
                points(1) = cells(1)
                call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
                call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)
                call DMGetStratumSize(self%dm, cell_order_label_name, &
                     global_cell_indices(2), num_matching, ierr); CHKERRQ(ierr)
                if (num_matching == 1) then
                   call DMGetStratumIS(self%dm, cell_order_label_name, &
                        global_cell_indices(2), cell_IS, ierr); CHKERRQ(ierr)
                   call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
                   points(2) = cells(1)
                   call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
                   call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)

                   ppoints => points
                   call DMPlexGetMeet(self%dm, num_cells, ppoints, cell_faces, ierr)
                   CHKERRQ(ierr)
                   num_cell_faces = size(cell_faces)
                   do i = 1, num_cell_faces
                      f = cell_faces(i)
                      if (self%ghost_face(f) < 0) then
                         call section_offset(face_section, f, face_offset, ierr)
                         CHKERRQ(ierr)
                         call face%assign_geometry(face_geom_array, face_offset)
                         face%permeability_direction = dble(permeability_direction)
                      end if
                   end do
                   call DMPlexRestoreMeet(self%dm, num_cells, ppoints, cell_faces, &
                        ierr); CHKERRQ(ierr)
                end if
             end if

          else
             if (present(logfile)) then
                call logfile%write(LOG_LEVEL_WARN, "input", &
                     "incorrect number of cells", int_keys = ["mesh.faces"], &
                     int_values = [iface - 1])
             end if
          end if
       end do
    end if

    call face%destroy()
    call VecRestoreArrayF90(self%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)

  end subroutine mesh_override_face_properties

!------------------------------------------------------------------------

  subroutine mesh_setup_minc(self, json, logfile, err)
    !! Sets up MINC for fractured zones.

    use fson_mpi_module
    use logfile_module
    use fson_value_m, only : TYPE_ARRAY, TYPE_OBJECT

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: num_minc_zones, iminc, num_minc_cells
    type(fson_value), pointer :: minc_json, minci_json
    PetscInt :: minc_type
    character(32) :: imincstr, mincstr
    PetscErrorCode :: ierr

    err = 0

    if (fson_has_mpi(json, "mesh.minc")) then

       call DMCreateLabel(self%original_dm, minc_label_name, ierr)
       CHKERRQ(ierr)

       call fson_get_mpi(json, "mesh.minc", minc_json)
       minc_type = fson_type_mpi(minc_json, ".")

       select case (minc_type)
       case (TYPE_OBJECT)
          num_minc_zones = 1
          allocate(self%minc(num_minc_zones))
          iminc = 1
          mincstr = "minc."
          call self%minc(num_minc_zones)%init(minc_json, self%original_dm, &
               iminc, mincstr, logfile, err)
       case (TYPE_ARRAY)
          num_minc_zones = fson_value_count_mpi(minc_json, ".")
          allocate(self%minc(num_minc_zones))
          do iminc = 1, num_minc_zones
             minci_json => fson_value_get_mpi(minc_json, iminc)
             write(imincstr, '(i0)') iminc - 1
             mincstr = 'minc[' // trim(imincstr) // '].'
             call self%minc(iminc)%init(minci_json, self%original_dm, iminc, &
                  mincstr, logfile, err)
             if (err > 0) exit
          end do
       end select

       call DMGetLabelSize(self%original_dm, minc_label_name, &
            num_minc_cells, ierr); CHKERRQ(ierr)
       self%has_minc = (num_minc_cells > 0)

    else
       self%has_minc = PETSC_FALSE
    end if

  end subroutine mesh_setup_minc

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_dm(self, dof)
    !! Sets up augmented DM for MINC mesh, including MINC cells and
    !! faces. Although they are not used, edges and vertices are also
    !! created for the MINC faces, so that the depth of the DM is
    !! consistent everywhere.

    use dm_utils_module, only: dm_copy_cone_sizes, dm_copy_cones, &
         set_dm_data_layout, dm_set_fv_adjacency, &
         dm_setup_fv_discretization, set_dm_default_data_layout

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: dof !! Degrees of freedom for discretization
    ! Locals:
    PetscInt :: start_chart, end_chart, c, i, iminc, h
    PetscInt :: dim, depth
    PetscInt :: num_minc_zone_cells, num_minc_cells
    PetscInt :: num_minc_zones, num_cells, num_new_points, max_num_levels
    PetscInt, allocatable :: start(:), end(:)
    PetscInt, allocatable :: frac_shift(:), minc_shift(:,:)
    PetscInt, allocatable :: minc_zone(:), num_minc_level_cells(:)
    IS :: minc_IS
    PetscInt, pointer :: minc_cells(:)
    PetscErrorCode :: ierr

    call DMPlexCreate(PETSC_COMM_WORLD, self%minc_dm, ierr); CHKERRQ(ierr)

    call DMGetDimension(self%original_dm, dim, ierr); CHKERRQ(ierr)
    call DMSetDimension(self%minc_dm, dim, ierr); CHKERRQ(ierr)
    call DMPlexGetDepth(self%original_dm, depth, ierr); CHKERRQ(ierr)
    allocate(start(0: depth), end(0: depth))
    do h = 0, depth
       call DMPlexGetHeightStratum(self%original_dm, h, start(h), end(h), ierr)
       CHKERRQ(ierr)
    end do

    num_cells = end(0) - start(0)
    num_minc_zones = size(self%minc)

    max_num_levels = 0
    do iminc = 1, num_minc_zones
       associate(num_levels => self%minc(iminc)%num_levels)
         max_num_levels = max(num_levels, max_num_levels)
       end associate
    end do
    allocate(minc_zone(0: num_cells - 1), &
         num_minc_level_cells(0: max_num_levels))
    minc_zone = -1
    num_minc_level_cells = 0
    do iminc = 1, num_minc_zones
       call DMGetStratumSize(self%original_dm, minc_label_name, iminc, &
            num_minc_zone_cells, ierr); CHKERRQ(ierr)
       if (num_minc_zone_cells > 0) then
          associate(num_levels => self%minc(iminc)%num_levels)
            do i = 0, num_levels
               num_minc_level_cells(i) = num_minc_level_cells(i) + &
                    num_minc_zone_cells
            end do
            call DMGetStratumIS(self%original_dm, minc_label_name, &
                 iminc, minc_IS, ierr); CHKERRQ(ierr)
            call ISGetIndicesF90(minc_IS, minc_cells, ierr); CHKERRQ(ierr)
            do i = 1, num_minc_zone_cells
               c = minc_cells(i)
               minc_zone(c) = iminc
            end do
            call ISRestoreIndicesF90(minc_IS, minc_cells, ierr); CHKERRQ(ierr)
            call ISDestroy(minc_IS, ierr); CHKERRQ(ierr)
          end associate
       end if
    end do

    call DMPlexGetChart(self%original_dm, start_chart, end_chart, ierr)
    CHKERRQ(ierr)
    num_minc_cells = sum(num_minc_level_cells(1:))
    num_new_points = num_minc_cells * (depth + 1)
    call DMPlexSetChart(self%minc_dm, start_chart, &
         end_chart + num_new_points, ierr); CHKERRQ(ierr)

    allocate(frac_shift(0: depth), minc_shift(0: depth, 1: max_num_levels))
    call setup_shifts(frac_shift, minc_shift)

    call set_cone_sizes()
    call DMSetUp(self%minc_dm, ierr); CHKERRQ(ierr)
    call set_cones()

    call DMPlexSymmetrize(self%minc_dm, ierr); CHKERRQ(ierr)
    call DMPlexStratify(self%minc_dm, ierr); CHKERRQ(ierr)
    call transfer_labels(self%dm, self%minc_dm, 1)
    call dm_set_fv_adjacency(self%minc_dm)
    call dm_setup_fv_discretization(self%minc_dm, dof)
    call set_dm_default_data_layout(self%minc_dm, dof)

    call self%assign_dm(self%minc_dm)

    deallocate(start, end, frac_shift, minc_shift, minc_zone, &
         num_minc_level_cells)

  contains

!........................................................................

    subroutine setup_shifts(frac_shift, minc_shift)
      !! Set up shift arrays to determine index shift from original DM
      !! point to corresponding fracture and MINC points in new DM.

      use utils_module, only: array_cumulative_sum

      PetscInt, intent(out) :: frac_shift(0: depth)
      PetscInt, intent(out) :: minc_shift(0: depth, 1: max_num_levels)
      ! Locals:
      PetscInt :: ishift(0: depth)
      PetscInt :: i, s, h, m
      PetscInt :: minc_offset(0: max_num_levels)

      !! Set up ishift array, to take account of the fact that DMPlex
      !! points have the order cells, vertices, faces, edges-
      !! i.e. they are not in depth order.
      ishift(0) = 0
      ishift(depth) = 1
      s = ishift(depth) + 1
      do i = 1, depth - 1
         ishift(i) = s
         s = s + 1
      end do

      minc_offset = 0
      minc_offset(1:) = array_cumulative_sum(num_minc_level_cells(1:))

      frac_shift = ishift * minc_offset(max_num_levels)
      do h = 0, depth
         do m = 1, max_num_levels
            minc_shift(h, m) = end(h) + frac_shift(h) + minc_offset(m - 1)
         end do
      end do

    end subroutine setup_shifts

!........................................................................

    subroutine set_cone_sizes
      !! Sets cone sizes for MINC DM.

      ! Locals:
      PetscInt :: p, cone_size, iminc, m, minc_p, h
      PetscErrorCode :: ierr

      ! Cells:
      do p = start(0), end(0) - 1
         call DMPlexGetConeSize(self%original_dm, p, cone_size, ierr)
         CHKERRQ(ierr)
         call DMPlexSetConeSize(self%minc_dm, p, cone_size + 1, ierr)
         CHKERRQ(ierr)
         iminc = minc_zone(p)
         if (iminc > 0) then
            associate(num_levels => self%minc(iminc)%num_levels)
              do m = 1, num_levels
                 minc_p = p + minc_shift(0, m)
                 if (m < num_levels) then
                    cone_size = 2
                 else
                    cone_size = 1
                 end if
                 call DMPlexSetConeSize(self%minc_dm, minc_p, &
                      cone_size, ierr); CHKERRQ(ierr)
              end do
            end associate
         end if
      end do

      ! Higher level points:
      do h = 1, depth
         call dm_copy_cone_sizes(self%original_dm, self%minc_dm, &
              start(h), end(h) - 1, frac_shift(h))
         if (h < depth) then
            cone_size = 1
         else
            cone_size = 0
         end if
         do p = start(0), end(0) - 1
            iminc = minc_zone(p)
            if (iminc > 0) then
               associate(num_levels => self%minc(iminc)%num_levels)
                 do m = 1, num_levels
                    minc_p = p + minc_shift(h, m)
                    call DMPlexSetConeSize(self%minc_dm, minc_p, &
                         cone_size, ierr); CHKERRQ(ierr)
                 end do
               end associate
            end if
         end do
      end do

    end subroutine set_cone_sizes

!........................................................................

    subroutine set_cones
      !! Sets cones for MINC DM.

      ! Locals:
      PetscInt :: p, m, h, iminc
      PetscInt :: face_p, inner_face_p, minc_p
      PetscInt :: above_p, cone_shift
      PetscInt, pointer :: points(:)
      PetscInt, allocatable :: frac_cone(:), minc_cone(:)
      PetscErrorCode :: ierr

      ! Cells:
      do p = start(0), end(0) - 1
         call DMPlexGetCone(self%original_dm, p, points, ierr)
         CHKERRQ(ierr)
         iminc = minc_zone(p)
         if (iminc > 0) then
            face_p = p + minc_shift(1, 1)
            frac_cone = [points + frac_shift(1), [face_p]]
            associate(num_levels => self%minc(iminc)%num_levels)
              do m = 1, num_levels
                 minc_p = p + minc_shift(0, m)
                 face_p = p + minc_shift(1, m)
                 if (m < num_levels) then
                    inner_face_p = p + minc_shift(1, m + 1)
                    minc_cone = [face_p, inner_face_p]
                 else
                    minc_cone = [face_p]
                 end if
                 call DMPlexSetCone(self%minc_dm, minc_p, &
                      minc_cone, ierr); CHKERRQ(ierr)
              end do
            end associate
         else
            frac_cone = points + frac_shift(1)
         end if
         call DMPlexSetCone(self%minc_dm, p, frac_cone, ierr); CHKERRQ(ierr)
         call DMPlexRestoreCone(self%original_dm, p, points, ierr)
         CHKERRQ(ierr)
      end do

      deallocate(frac_cone, minc_cone)

      ! Higher level points:
      do h = 1, depth
         if (h < depth) then
            cone_shift = frac_shift(h + 1)
         else
            cone_shift = 0
         end if
         call dm_copy_cones(self%original_dm, self%minc_dm, start(h), &
              end(h) - 1, frac_shift(h), cone_shift)
         do p = start(0), end(0) - 1
            iminc = minc_zone(p)
            if (iminc > 0) then
               associate(num_levels => self%minc(iminc)%num_levels)
                 do m = 1, num_levels
                    minc_p = p + minc_shift(h, m)
                    if (h < depth) then
                       above_p = p + minc_shift(h + 1, m)
                       call DMPlexSetCone(self%minc_dm, minc_p, &
                            [above_p], ierr)
                       CHKERRQ(ierr)
                    end if
                 end do
               end associate
            end if
         end do
      end do

    end subroutine set_cones

!........................................................................

    subroutine transfer_labels(dm, minc_dm, max_height)
      !! Transfers relevant labels from original DM to MINC DM
      !! fracture points, applying appropriate shifts to the point
      !! indices.

      use rock_module, only: rocktype_label_name

      DM, intent(in) :: dm
      DM, intent(in out) :: minc_dm
      PetscInt, intent(in) :: max_height
      ! Locals:
      PetscInt :: p, h, l, iid, ip, label_value
      PetscInt :: num_ids, num_points
      PetscInt, parameter :: max_label_name_length = 80
      character(max_label_name_length) :: label_name
      PetscInt, parameter :: num_labels = 6
      character(max_label_name_length) :: label_names(0: num_labels - 1)
      IS :: id_IS, point_IS
      PetscInt, pointer :: ids(:), points(:)
      PetscBool :: has_label
      PetscErrorCode :: ierr

      ! call DMGetNumLabels(dm, num_labels, ierr); CHKERRQ(ierr)
      label_names = [character(max_label_name_length):: &
           "ghost", "Cell Sets", &
           rocktype_label_name, minc_label_name, &
           open_boundary_label_name, cell_order_label_name]

      do l = 0, num_labels - 1
         ! call DMGetLabelName(dm, l, label_name, ierr); CHKERRQ(ierr) ! missing Fortran interface
         label_name = label_names(l)
         call DMHasLabel(dm, label_name, has_label, ierr); CHKERRQ(ierr)
         if (has_label .and. (label_name /= 'depth')) then
            call DMCreateLabel(minc_dm, label_name, ierr); CHKERRQ(ierr)
            call DMGetLabelIdIS(dm, label_name, id_IS, ierr); CHKERRQ(ierr)
            call ISGetLocalSize(id_IS, num_ids, ierr); CHKERRQ(ierr)
            if (num_ids > 0) then
               call ISGetIndicesF90(id_IS, ids, ierr); CHKERRQ(ierr)
               do iid = 1, num_ids
                  call DMGetStratumIS(dm, label_name, ids(iid), point_IS, &
                       ierr); CHKERRQ(ierr)
                  call ISGetLocalSize(point_IS, num_points, ierr)
                  CHKERRQ(ierr)
                  if (num_points > 0) then
                     call ISGetIndicesF90(point_IS, points, ierr)
                     CHKERRQ(ierr)
                     do ip = 1, num_points
                        p = points(ip)
                        h = dm_point_height(dm, p)
                        call DMSetLabelValue(minc_dm, label_name, &
                             p + frac_shift(h), label_value, ierr)
                        CHKERRQ(ierr)
                     end do
                     call ISRestoreIndicesF90(point_IS, points, ierr)
                     CHKERRQ(ierr)
                  end if
                  call ISDestroy(point_IS, ierr); CHKERRQ(ierr)
               end do
               call ISRestoreIndicesF90(id_IS, ids, ierr); CHKERRQ(ierr)
            end if
            call ISDestroy(id_IS, ierr); CHKERRQ(ierr)
         end if
      end do

    end subroutine transfer_labels

!........................................................................

    PetscInt function dm_point_height(dm, p) result(height)
      !! Returns height of point in DM.

      DM, intent(in) :: dm
      PetscInt, intent(in) :: p
      ! Locals:
      PetscInt :: h

      height = -1
      do h = 0, depth
         if ((start(h) <= p) .and. (p < end(h))) then
            height = h
            exit
         end if
      end do

    end function dm_point_height

!........................................................................

  end subroutine mesh_setup_minc_dm

!------------------------------------------------------------------------

end module mesh_module
