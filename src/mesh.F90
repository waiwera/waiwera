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
  use fson
  use kinds_module
  use cell_module
  use face_module
  use zone_module
  use list_module
  use mpi_utils_module
  use minc_module
  use dm_utils_module, only: dm_stratum_type
  use dictionary_module

  implicit none

  private

  PetscInt, parameter, public :: max_mesh_filename_length = 200
  PetscInt, parameter :: partition_overlap = 1 !! Cell overlap for parallel mesh distribution
  character(len = 16), public :: open_boundary_label_name = "open_boundary" !! Name of DMLabel for identifying open boundary faces
  character(len = 16), public :: boundary_ghost_label_name = "boundary_ghost" !! Name of DMLabel for identifying boundary ghost cells
  character(len = 26) :: face_permeability_override_label_name = "face_permeability_override" !! Name of DMLabel for overriding face permeabilities

  type, public :: mesh_type
     !! Mesh type.
     private
     character(max_mesh_filename_length), public :: filename !! Mesh file name
     DM, public :: serial_dm !! Original DM read from file (not distributed)
     DM, public :: original_dm !! Original DM read from file (and distributed)
     DM, public :: dm !! DM representing the mesh topology (may be modified from original_dm)
     Vec, public :: cell_geom !! Vector containing cell geometry data
     Vec, public :: face_geom !! Vector containing face geometry data
     PetscInt, public :: dim !! DM dimension
     PetscInt :: depth !! DM depth
     type(dm_stratum_type), allocatable, public :: strata(:) !! Mesh strata (used for MINC point calculations)
     IS, public :: cell_index !! Index set defining natural to global cell ordering (without boundary cells)
     AO, public :: cell_natural_global !! Application ordering to convert between natural and global cell indices
     AO, public :: original_cell_natural_global !! Natural-to-global AO for original DM
     IS, public :: cell_natural !! Natural indices of local cells
     IS, public :: cell_parent_natural !! Natural indices of parent cells (e.g. MINC fracture cells)
     PetscSF, public :: dist_sf !! Distribution star forest
     PetscInt, public, allocatable :: ghost_cell(:), ghost_face(:) !! Ghost label values for cells and faces
     type(minc_type), allocatable, public :: minc(:) !! Array of MINC zones, with parameters
     PetscReal, public :: permeability_rotation(3, 3) !! Rotation matrix of permeability axes
     PetscReal, public :: thickness !! Mesh thickness (for dimension < 3)
     type(list_type), public :: zones !! Mesh zones
     type(dictionary_type), public :: rock_types !! Dictionary of rock types by name
     PetscBool, public :: radial !! If mesh coordinate system is radial or Cartesian
     PetscBool, public :: has_minc !! If mesh has any MINC cells
     PetscBool, public :: rebalance !! Set false to disable MINC mesh rebalancing
     PetscInt, public :: dof !! Degrees of freedom for default section
   contains
     procedure :: distribute => mesh_distribute
     procedure :: setup_geometry => mesh_setup_geometry
     procedure :: setup_ghost_arrays => mesh_setup_ghost_arrays
     procedure :: destroy_minc => mesh_destroy_minc
     procedure :: destroy_strata => mesh_destroy_strata
     procedure :: setup_coordinate_parameters => mesh_setup_coordinate_parameters
     procedure :: set_permeability_rotation => mesh_set_permeability_rotation
     procedure :: modify_cell_geometry => mesh_modify_cell_geometry
     procedure :: modify_face_geometry => mesh_modify_face_geometry
     procedure :: modify_geometry => mesh_modify_geometry
     procedure :: read_overridden_face_properties => mesh_read_overridden_face_properties
     procedure :: override_face_properties => mesh_override_face_properties
     procedure :: label_cell_array_rock_types => mesh_label_cell_array_rock_types
     procedure :: label_cell_array_zones => mesh_label_cell_array_zones
     procedure :: label_cell_array_minc_zones => mesh_label_cell_array_minc_zones
     procedure :: label_boundaries => mesh_label_boundaries
     procedure :: label_sources => mesh_label_sources
     procedure :: setup_zones => mesh_setup_zones
     procedure :: setup_minc => mesh_setup_minc
     procedure :: setup_minc_dm => mesh_setup_minc_dm
     procedure :: setup_minc_level_cells => mesh_setup_minc_level_cells
     procedure :: setup_minc_dm_strata_shifts => mesh_setup_minc_dm_strata_shifts
     procedure :: set_minc_dm_cone_sizes => mesh_set_minc_dm_cone_sizes
     procedure :: set_minc_dm_cones => mesh_set_minc_dm_cones
     procedure :: setup_minc_dm_depth_label => mesh_setup_minc_dm_depth_label
     procedure :: transfer_labels_to_minc_dm => mesh_transfer_labels_to_minc_dm
     procedure :: setup_minc_output_data => mesh_setup_minc_output_data
     procedure :: setup_minc_dm_cell_natural_global => mesh_setup_minc_dm_cell_natural_global
     procedure :: setup_minc_geometry => mesh_setup_minc_geometry
     procedure :: setup_minc_rock_properties => mesh_setup_minc_rock_properties
     procedure :: setup_minc_point_sf => mesh_setup_minc_point_sf
     procedure :: setup_minc_coordinates => mesh_setup_minc_coordinates
     procedure :: redistribute_dm => mesh_redistribute_dm
     procedure :: redistribute_geometry => mesh_redistribute_geometry
     procedure :: check_face_orientations => mesh_check_face_orientations
     procedure :: geometry_add_boundary => mesh_geometry_add_boundary
     procedure :: boundary_face_geometry => mesh_boundary_face_geometry
     procedure :: setup_cell_natural => mesh_setup_cell_natural
     procedure, public :: init => mesh_init
     procedure, public :: configure => mesh_configure
     procedure, public :: construct_ghost_cells => mesh_construct_ghost_cells
     procedure, public :: set_boundary_conditions => mesh_set_boundary_conditions
     procedure, public :: destroy => mesh_destroy
     procedure, public :: local_to_parent_natural => mesh_local_to_parent_natural
     procedure, public :: global_to_parent_natural => mesh_global_to_parent_natural
     procedure, public :: natural_cell_output_arrays =>  mesh_natural_cell_output_arrays
     procedure, public :: local_cell_minc_level => mesh_local_cell_minc_level
     procedure, public :: destroy_distribution_data => mesh_destroy_distribution_data
     procedure, public :: redistribute => mesh_redistribute
  end type mesh_type

contains

!------------------------------------------------------------------------

  subroutine mesh_distribute(self)
    !! Distributes mesh over processors, and returns star forest from
    !! mesh distribution.

    use dm_utils_module, only: dm_set_default_data_layout, &
         dm_label_partition_ghosts, dm_distribute_index_set
    
    class(mesh_type), intent(in out) :: self
    ! Locals:
    DM :: dm_is
    PetscSection :: section
    PetscErrorCode :: ierr

    call DMClone(self%serial_dm, dm_is, ierr); CHKERRQ(ierr)
    call dm_set_default_data_layout(dm_is, 1)
    call DMGetSection(dm_is, section, ierr); CHKERRQ(ierr)

    call DMPlexDistribute(self%serial_dm, partition_overlap, self%dist_sf, &
         self%original_dm, ierr); CHKERRQ(ierr)
    if (self%original_dm .eq. PETSC_NULL_DM) then
       self%original_dm = self%serial_dm
    else
       call dm_distribute_index_set(self%original_dm, self%dist_sf, &
            section, self%cell_natural)
    end if
    call dm_label_partition_ghosts(self%original_dm)
    call DMDestroy(dm_is, ierr); CHKERRQ(ierr)

  end subroutine mesh_distribute

!------------------------------------------------------------------------

  subroutine mesh_destroy_distribution_data(self)
    !! Destroys distribution start forest and serial DM.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscMPIInt :: np
    PetscErrorCode :: ierr

    call MPI_comm_size(PETSC_COMM_WORLD, np, ierr)
    if (np > 1) then
       call DMDestroy(self%serial_dm, ierr); CHKERRQ(ierr)
    end if
    call PetscSFDestroy(self%dist_sf, ierr); CHKERRQ(ierr)

  end subroutine mesh_destroy_distribution_data

!------------------------------------------------------------------------

  subroutine mesh_construct_ghost_cells(self, gravity)
    !! Constructs ghost cells on open boundary faces.

    use dm_utils_module, only: dm_set_fv_adjacency, &
         dm_set_default_data_layout, dm_label_boundary_ghosts, &
         dm_get_natural_to_global_ao

    class(mesh_type), intent(in out) :: self
    PetscReal, intent(in) :: gravity(:) !! Gravity vector
    ! Locals:
    DM :: ghost_dm
    DMLabel :: label
    PetscErrorCode :: ierr

    call DMRemoveLabel(self%dm, "ghost", label, ierr); CHKERRQ(ierr)
    call DMRemoveLabel(self%dm, "vtk", label, ierr); CHKERRQ(ierr)

    call DMPlexConstructGhostCells(self%dm, open_boundary_label_name, &
         PETSC_NULL_INTEGER, ghost_dm, ierr); CHKERRQ(ierr)
    if (ghost_dm .ne. PETSC_NULL_DM) then
       call DMDestroy(self%dm, ierr); CHKERRQ(ierr);
       self%dm = ghost_dm
       call dm_set_fv_adjacency(self%dm)
       call dm_set_default_data_layout(self%dm, self%dof)
       call dm_label_boundary_ghosts(self%dm, boundary_ghost_label_name)
       self%cell_natural_global = dm_get_natural_to_global_ao(self%dm, self%cell_natural)
       call self%geometry_add_boundary(gravity)
       call self%setup_ghost_arrays()
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
    PetscBool, parameter :: default_radial = PETSC_FALSE
    PetscReal, parameter :: default_thickness = 1._dp

    self%thickness = default_thickness
    select case (self%dim)
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

  subroutine mesh_modify_cell_geometry(self, cell)
    ! Modifies cell geometry parameters for 2-D and radial meshes.

    class(mesh_type), intent(in) :: self
    type(cell_type), intent(in out) :: cell

    if (self%dim == 2) then
       if (self%radial) then
          call modify_cell_geometry_2d_radial(cell)
       else
          call modify_cell_geometry_2d_cartesian(cell)
       end if
    end if

  contains

!........................................................................

    subroutine modify_cell_geometry_2d_cartesian(cell)
      ! Geometry modification for 2D Cartesian cells.

      type(cell_type), intent(in out) :: cell

      cell%volume = cell%volume * self%thickness

    end subroutine modify_cell_geometry_2d_cartesian

!........................................................................

    subroutine modify_cell_geometry_2d_radial(cell)
      ! Geometry modification for 2D radial cells- via Pappus' centroid
      ! theorem.

      use kinds_module, only: dp
      use utils_module, only: pi

      type(cell_type), intent(in out) :: cell
      ! Locals:
      PetscReal :: r

      r = cell%centroid(1)
      cell%volume = cell%volume * 2._dp * pi * r

    end subroutine modify_cell_geometry_2d_radial

  end subroutine mesh_modify_cell_geometry

!------------------------------------------------------------------------

  subroutine mesh_modify_face_geometry(self, face)
    ! Modifies face geometry parameters for 2-D and radial meshes.

    class(mesh_type), intent(in) :: self
    type(face_type), intent(in out) :: face

    if (self%dim == 2) then
       if (self%radial) then
          call modify_face_geometry_2d_radial(face)
       else
          call modify_face_geometry_2d_cartesian(face)
       end if
    end if

  contains

!........................................................................

    subroutine modify_face_geometry_2d_cartesian(face)
      ! Geometry modification for 2D Cartesian faces.

      type(face_type), intent(in out) :: face

      face%area = face%area * self%thickness

    end subroutine modify_face_geometry_2d_cartesian

!........................................................................

    subroutine modify_face_geometry_2d_radial(face)
      ! Geometry modification for 2D radial faces- via Pappus' centroid
      ! theorem.

      use kinds_module, only: dp
      use utils_module, only: pi

      type(face_type), intent(in out) :: face
      ! Locals:
      PetscReal :: r

      r = face%centroid(1)
      face%area = face%area * 2._dp * pi * r

    end subroutine modify_face_geometry_2d_radial

  end subroutine mesh_modify_face_geometry

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
    use utils_module, only: rotation_matrix_2d
    use dm_utils_module, only: section_offset, local_vec_section, dm_set_data_layout

    class(mesh_type), intent(in out) :: self
    Vec, intent(in) :: petsc_face_geom
    PetscReal, intent(in) :: gravity(:)
    ! Locals:
    DM :: dm_cell, dm_face
    PetscSection :: face_section, petsc_face_section, cell_section
    PetscInt :: c, f, i, ghost_cell, ghost_face
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
    PetscInt :: face_variable_dim(num_face_variables)
    PetscInt :: cell_variable_dim(num_cell_variables)
    PetscErrorCode :: ierr

    call DMGetLabel(self%original_dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(self%original_dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    call DMPlexGetHeightStratum(self%original_dm, 1, start_face, end_face, ierr)
    CHKERRQ(ierr)

    ! Set up cell geometry vector:
    call VecGetDM(self%cell_geom, dm_cell, ierr); CHKERRQ(ierr)
    cell_variable_dim = self%dim
    call dm_set_data_layout(dm_cell, cell_variable_num_components, &
         cell_variable_dim, cell_variable_names)

    call local_vec_section(self%cell_geom, cell_section)
    call VecGetArrayF90(self%cell_geom, cell_geom_array, ierr); CHKERRQ(ierr)

    if (self%dim == 2) then
       ! Adjust cell volumes:
       call cell%init(1, 1)
       do c = start_cell, end_cell - 1
          call DMLabelGetValue(ghost_label, c, ghost_cell, ierr); CHKERRQ(ierr)
          if (ghost_cell < 0) then
             offset = section_offset(cell_section, c)
             call cell%assign_geometry(cell_geom_array, offset)
             call self%modify_cell_geometry(cell)
          end if
       end do
       call cell%destroy()
    end if

    ! Set up face geometry vector:
    call DMClone(self%original_dm, dm_face, ierr); CHKERRQ(ierr)
    face_variable_dim = self%dim - 1
    call dm_set_data_layout(dm_face, face_variable_num_components, &
         face_variable_dim, face_variable_names)
    call DMCreateLocalVector(dm_face, self%face_geom, ierr); CHKERRQ(ierr)
    call local_vec_section(self%face_geom, face_section)
    call VecGetArrayF90(self%face_geom, face_geom_array, ierr); CHKERRQ(ierr)
    call local_vec_section(petsc_face_geom, petsc_face_section)
    call VecGetArrayF90(petsc_face_geom, petsc_face_geom_array, ierr)
    CHKERRQ(ierr)
    call face%init()

    do f = start_face, end_face - 1

       call DMLabelGetValue(ghost_label, f, ghost_face, ierr); CHKERRQ(ierr)

       if (ghost_face < 0) then

          face_offset = section_offset(face_section, f)
          petsc_face_offset = section_offset(petsc_face_section, f)

          call DMPlexGetSupport(self%original_dm, f, cells, ierr); CHKERRQ(ierr)
          do i = 1, 2
             cell_offset(i) = section_offset(cell_section, cells(i))
          end do

          call petsc_face%assign_geometry(petsc_face_geom_array, petsc_face_offset)
          call face%assign_geometry(face_geom_array, face_offset)
          call face%assign_cell_geometry(cell_geom_array, cell_offset)

          face%centroid = petsc_face%centroid
          face%area = norm2(petsc_face%area_normal)
          face%normal = petsc_face%area_normal / face%area
          face%gravity_normal = dot_product(gravity, face%normal)
          call self%modify_face_geometry(face)
          call face%calculate_distances()
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

  end subroutine mesh_modify_geometry

!------------------------------------------------------------------------

  subroutine mesh_boundary_face_geometry(self, gravity)
    !! Computes face geometry parameters for open boundary faces. Also
    !! sets the volumes of boundary ghost cells to zero, and their
    !! centroids to the face centroids (these are not used).

    use cell_module, only: cell_type
    use face_module, only: face_type
    use dm_utils_module, only: local_vec_section, section_offset

    class(mesh_type), intent(in out) :: self
    PetscReal, intent(in) :: gravity(:) !! Gravity vector
    ! Locals:
    PetscSection :: cell_section, face_section
    PetscReal, pointer, contiguous :: cell_geom_array(:), face_geom_array(:)
    DMLabel :: bdy_label, ghost_label
    IS :: label_IS, bdy_IS
    PetscInt, pointer :: label_values(:), bdy_faces(:)
    PetscInt :: i, j, num_values, ibdy, num_faces, iface, ghost, f
    PetscInt :: cell_offsets(2), face_offset
    PetscInt, pointer :: cells(:)
    type(cell_type) :: cell
    type(face_type) :: face
    PetscErrorCode :: ierr

    call local_vec_section(self%cell_geom, cell_section)
    call VecGetArrayF90(self%cell_geom, cell_geom_array, ierr); CHKERRQ(ierr)
    call local_vec_section(self%face_geom, face_section)
    call VecGetArrayF90(self%face_geom, face_geom_array, ierr); CHKERRQ(ierr)

    call DMGetLabel(self%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    call DMGetLabel(self%dm, open_boundary_label_name, bdy_label, ierr)
    CHKERRQ(ierr)
    call DMLabelGetValueIS(bdy_label, label_IS, ierr); CHKERRQ(ierr)
    call ISGetIndicesF90(label_IS, label_values, ierr); CHKERRQ(ierr)

    call cell%init(1, 1)
    call face%init()

    num_values = size(label_values)
    do i = 1, num_values
       ibdy = label_values(i)
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
                call DMLabelGetValue(ghost_label, cells(1), ghost, ierr); CHKERRQ(ierr)
                if (ghost < 0) then
                   face_offset = section_offset(face_section, f)
                   do j = 1, 2
                      cell_offsets(j) = section_offset(cell_section, cells(j))
                   end do
                   call face%assign_geometry(face_geom_array, face_offset)
                   call face%assign_cell_geometry(cell_geom_array, cell_offsets)
                   call DMPlexComputeCellGeometryFVM(self%dm, f, face%area, &
                        face%centroid, face%normal, ierr); CHKERRQ(ierr)
                   face%gravity_normal = dot_product(gravity, face%normal)
                   call face%calculate_permeability_direction(self%permeability_rotation)
                   call self%modify_face_geometry(face)
                   face%distance = [dot_product(face%centroid - face%cell(1)%centroid, &
                        face%normal), 0._dp]
                   face%distance12 = face%distance(1)
                   face%cell(2)%volume = 0._dp
                   face%cell(2)%centroid = face%centroid
                end if
             end if
          end do
       end if
    end do

    call ISRestoreIndicesF90(label_IS, label_values, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(self%face_geom, face_geom_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(self%cell_geom, cell_geom_array, ierr); CHKERRQ(ierr)
    call cell%destroy()
    call face%destroy()

  end subroutine mesh_boundary_face_geometry

!------------------------------------------------------------------------

  subroutine mesh_geometry_add_boundary(self, gravity)
    !! Adds space for Dirichlet boundary condition ghost cells to mesh
    !! geometry vectors.

    use dm_utils_module, only: dm_set_data_layout, vec_copy_common_local

    class(mesh_type), intent(in out) :: self
    PetscReal, intent(in) :: gravity(:) !! Gravity vector
    ! Locals:
    Vec :: cell_geom, face_geom
    DM :: dm_cell, dm_face
    PetscInt :: cell_variable_dim(num_cell_variables)
    PetscInt :: face_variable_dim(num_face_variables)
    PetscErrorCode :: ierr

    call DMClone(self%dm, dm_cell, ierr); CHKERRQ(ierr)
    cell_variable_dim = self%dim
    call dm_set_data_layout(dm_cell, cell_variable_num_components, &
         cell_variable_dim, cell_variable_names)
    call DMCreateLocalVector(dm_cell, cell_geom, ierr); CHKERRQ(ierr)
    call vec_copy_common_local(self%cell_geom, cell_geom)
    call VecDestroy(self%cell_geom, ierr); CHKERRQ(ierr)
    self%cell_geom = cell_geom

    call DMClone(self%dm, dm_face, ierr); CHKERRQ(ierr)
    face_variable_dim = self%dim - 1
    call dm_set_data_layout(dm_face, face_variable_num_components, &
         face_variable_dim, face_variable_names)
    call DMCreateLocalVector(dm_face, face_geom, ierr); CHKERRQ(ierr)
    call vec_copy_common_local(self%face_geom, face_geom)
    call VecDestroy(self%face_geom, ierr); CHKERRQ(ierr)
    self%face_geom = face_geom

    call self%boundary_face_geometry(gravity)

  end subroutine mesh_geometry_add_boundary

!------------------------------------------------------------------------

  subroutine mesh_setup_cell_natural(self)
    !! Sets up cell_natural IS on serial DM.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscInt :: c, start_cell, end_cell
    PetscInt, allocatable :: natural(:)
    PetscErrorCode :: ierr

    call DMPlexGetHeightStratum(self%serial_dm, 0, start_cell, &
         end_cell, ierr); CHKERRQ(ierr)
    natural = [(c, c = start_cell, end_cell -1)]
    call ISCreateGeneral(PETSC_COMM_WORLD, end_cell - start_cell, &
            natural, PETSC_COPY_VALUES, self%cell_natural, ierr)
    deallocate(natural)

  end subroutine mesh_setup_cell_natural

!------------------------------------------------------------------------

  subroutine mesh_setup_ghost_arrays(self)
    !! Sets up arrays of ghost label values on cells and faces. This
    !! is just for faster access to these values in the rest of the
    !! code.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscInt :: c, f
    PetscInt :: start_cell, end_cell, start_face, end_face
    DMLabel :: ghost_label
    PetscErrorCode :: ierr

    call DMPlexGetHeightStratum(self%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    call DMPlexGetHeightStratum(self%dm, 1, start_face, end_face, ierr)
    CHKERRQ(ierr)

    if (allocated(self%ghost_cell)) deallocate(self%ghost_cell)
    if (allocated(self%ghost_face)) deallocate(self%ghost_face)

    allocate(self%ghost_cell(start_cell: end_cell - 1))
    allocate(self%ghost_face(start_face: end_face - 1))
    self%ghost_cell = -1
    self%ghost_face = -1

    call DMGetLabel(self%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    do c = start_cell, end_cell - 1
       call DMLabelGetValue(ghost_label, c, self%ghost_cell(c), ierr)
       CHKERRQ(ierr)
    end do

    do f = start_face, end_face - 1
       call DMLabelGetValue(ghost_label, f, self%ghost_face(f), ierr)
       CHKERRQ(ierr)
    end do

  end subroutine mesh_setup_ghost_arrays

!------------------------------------------------------------------------

  subroutine mesh_destroy_minc(self)
    !! Destroys MINC objects.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscInt :: i

    do i = 1, size(self%minc)
       call self%minc(i)%destroy()
    end do
    deallocate(self%minc)

  end subroutine mesh_destroy_minc

!------------------------------------------------------------------------

  subroutine mesh_destroy_strata(self)
    !! Destroys mesh strata.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscInt :: h

    do h = 0, self%depth
       call self%strata(h)%destroy()
    end do
    deallocate(self%strata)

  end subroutine mesh_destroy_strata

!------------------------------------------------------------------------

  subroutine mesh_init(self, eos, json, logfile)
    !! Initializes mesh, reading filename from JSON input file.
    !! If the filename is not present, an error is raised.
    !! Otherwise, the PETSc DM is read in.

    use logfile_module
    use fson_mpi_module
    use fson_value_m, only: TYPE_STRING, TYPE_OBJECT
    use eos_module
    use dm_utils_module

    class(mesh_type), intent(in out) :: self
    class(eos_type), intent(in) :: eos !! EOS object
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
       ! Read in DM:
       call DMPlexCreateFromFile(PETSC_COMM_WORLD, self%filename, PETSC_TRUE, &
            self%serial_dm, ierr); CHKERRQ(ierr)
       call dm_set_fv_adjacency(self%serial_dm)
       self%dof = eos%num_primary_variables
       call DMGetDimension(self%serial_dm, self%dim, ierr); CHKERRQ(ierr)
       call dm_set_default_data_layout(self%serial_dm, self%dof)
       call self%setup_coordinate_parameters(json, logfile)
       call self%set_permeability_rotation(json, logfile)
       call self%read_overridden_face_properties(json, logfile)
       call self%label_cell_array_zones(json)
       call self%label_cell_array_rock_types(json)
       call self%label_cell_array_minc_zones(json)
       call self%label_boundaries(json, logfile)
       call self%label_sources(json)
       call self%setup_cell_natural()
       self%has_minc = PETSC_FALSE
    end if

  end subroutine mesh_init

!------------------------------------------------------------------------

  subroutine mesh_configure(self, gravity, json, logfile, err)
    !! Configures mesh, including distribution over processes and
    !! construction of ghost cells, setup of data layout, geometry and
    !! cell index set.

    use dm_utils_module, only: dm_set_default_data_layout, &
         dm_set_fv_adjacency, dm_get_natural_to_global_ao, &
         dm_get_cell_index
    use logfile_module
    use rock_module, only: setup_rock_types

    class(mesh_type), intent(in out) :: self
    PetscReal, intent(in) :: gravity(:) !! Gravity vector
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    type(logfile_type), intent(in out), optional :: logfile !! Log file
    PetscErrorCode, intent(out) :: err !! Error flag

    err = 0

    call self%distribute()
    call dm_set_default_data_layout(self%original_dm, self%dof)
    call dm_set_fv_adjacency(self%original_dm)
    call self%setup_geometry(gravity)
    self%dm = self%original_dm
    self%original_cell_natural_global = &
         dm_get_natural_to_global_ao(self%original_dm, self%cell_natural)
    self%cell_natural_global = self%original_cell_natural_global

    call self%setup_zones(json, logfile, err)
    if (err == 0) then
       call setup_rock_types(json, self%dm, self%rock_types, logfile, err)
       if (err == 0) then
          call self%setup_minc(json, logfile, err)
          if (err == 0) then
             if (self%has_minc) call self%setup_minc_dm()
             call dm_get_cell_index(self%dm, self%cell_natural_global, &
                  self%cell_index)
          end if
       end if
    end if

    call self%setup_ghost_arrays()
    
  end subroutine mesh_configure

!------------------------------------------------------------------------

  subroutine mesh_destroy(self)
    !! Destroys mesh object.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call VecDestroy(self%cell_geom, ierr); CHKERRQ(ierr)
    call VecDestroy(self%face_geom, ierr); CHKERRQ(ierr)
    call DMDestroy(self%dm, ierr); CHKERRQ(ierr)
    call ISDestroy(self%cell_index, ierr); CHKERRQ(ierr)

    call ISDestroy(self%cell_natural, ierr); CHKERRQ(ierr)
    call AODestroy(self%cell_natural_global, ierr); CHKERRQ(ierr)

    if (allocated(self%ghost_cell)) then
       deallocate(self%ghost_cell)
    end if
    if (allocated(self%ghost_face)) then
       deallocate(self%ghost_face)
    end if

    if (self%has_minc) then
       call DMDestroy(self%original_dm, ierr); CHKERRQ(ierr)
       call AODestroy(self%original_cell_natural_global, ierr); CHKERRQ(ierr)
       call ISDestroy(self%cell_parent_natural, ierr); CHKERRQ(ierr)
       call self%destroy_minc()
       call self%destroy_strata()
    end if
    call self%zones%destroy(mesh_zones_node_data_destroy)

    call self%rock_types%destroy()

  contains

    subroutine mesh_zones_node_data_destroy(node)
      ! Destroys zone data.

      type(list_node_type), pointer, intent(in out) :: node

      select type (zone => node%data)
      class is (zone_type)
         call zone%destroy()
      end select

    end subroutine mesh_zones_node_data_destroy

  end subroutine mesh_destroy

!------------------------------------------------------------------------

  subroutine mesh_set_boundary_conditions(self, json, y, fluid_vector, rock_vector, &
       eos, y_range_start, fluid_range_start, rock_range_start, logfile)
    !! Sets primary variables (and rock properties) in boundary ghost
    !! cells.

    use fson
    use fson_mpi_module
    use kinds_module
    use dm_utils_module, only: global_vec_section, global_section_offset
    use eos_module, only: eos_type
    use fluid_module, only: fluid_type
    use rock_module, only: rock_type
    use logfile_module

    class(mesh_type), intent(in) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input file
    Vec, intent(in out) :: y !! Primary variables vector
    Vec, intent(in out) :: fluid_vector !! Fluid properties vector
    Vec, intent(in out) :: rock_vector !! Rock properties vector
    class(eos_type), intent(in) :: eos !! EOS module
    PetscInt, intent(in) :: y_range_start !! Start of range for global primary variables vector
    PetscInt, intent(in) :: fluid_range_start !! Start of range for global fluid vector
    PetscInt, intent(in) :: rock_range_start !! Start of range for global rock vector
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    type(fson_value), pointer :: boundaries, bdy
    PetscInt :: num_boundaries, ibdy, f, num_faces, iface, np, n
    PetscReal, pointer, contiguous :: y_array(:), fluid_array(:), rock_array(:)
    PetscReal, pointer, contiguous :: cell_primary(:), rock1(:), rock2(:)
    PetscSection :: y_section, fluid_section, rock_section
    IS :: bdy_IS
    DMLabel :: ghost_label
    type(fluid_type):: fluid
    type(rock_type) :: rock
    PetscInt :: y_offset, fluid_offset, rock_offsets(2)
    PetscInt :: ghost, region, i
    PetscInt, pointer :: bdy_faces(:), cells(:)
    PetscReal, allocatable :: primary(:)
    character(len=64) :: bdystr
    character(len=12) :: istr
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

    if (fson_has_mpi(json, "boundaries")) then

       call fson_get_mpi(json, "boundaries", boundaries)
       num_boundaries = fson_value_count_mpi(boundaries, ".")
       bdy => fson_value_children_mpi(boundaries)

       do ibdy = 1, num_boundaries
          write(istr, '(i0)') ibdy - 1
          bdystr = 'boundaries[' // trim(istr) // ']'
          call fson_get_mpi(bdy, "primary", eos%default_primary, &
               primary, logfile, log_key = trim(bdystr) // ".primary")
          call fson_get_mpi(bdy, "region", eos%default_region, &
               region, logfile, log_key = trim(bdystr) // ".region")
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
                      y_offset = global_section_offset(y_section, cells(2), &
                           y_range_start)
                      fluid_offset = global_section_offset(fluid_section, cells(2), &
                           fluid_range_start)
                      do i = 1, 2
                         rock_offsets(i) = global_section_offset(rock_section, cells(i), &
                              rock_range_start)
                      end do
                      ! Set primary variables and region:
                      cell_primary => y_array(y_offset : y_offset + np - 1)
                      call fluid%assign(fluid_array, fluid_offset)
                      cell_primary = primary
                      fluid%region = dble(region)
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
          bdy => fson_value_next_mpi(bdy)
          deallocate(primary)
       end do
    end if

    call fluid%destroy()
    call rock%destroy()
    call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)

  end subroutine mesh_set_boundary_conditions

!------------------------------------------------------------------------

  subroutine mesh_read_overridden_face_properties(self, json, logfile)

    !! Reads in serial data structures for overridden face properties
    !! from JSON input - currently just face permeability directions.

    use fson_utils_module, only: fson_get_default
    use fson_value_m, only : fson_value_count
    use logfile_module

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    PetscMPIInt :: rank
    PetscInt, allocatable, target :: cell_indices(:)
    PetscInt, allocatable :: default_cells(:)
    type(fson_value), pointer :: faces_json, face_json
    PetscInt :: num_faces, iface, f, i, num_cell_faces
    PetscInt :: permeability_direction
    PetscInt :: start_face, end_face
    PetscInt, pointer :: cell_faces(:)
    PetscInt, pointer :: pcells(:)
    character(len=64) :: facestr
    character(len=12) :: istr
    PetscErrorCode :: ierr
    PetscInt, parameter :: default_permeability_direction = 1

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call DMCreateLabel(self%serial_dm, face_permeability_override_label_name, &
         ierr); CHKERRQ(ierr)
    allocate(default_cells(0))

    if (rank == 0) then

       call DMPlexGetHeightStratum(self%serial_dm, 1, start_face, end_face, ierr)
       CHKERRQ(ierr)

       call fson_get(json, "mesh.faces", faces_json)
       if (associated(faces_json)) then
          face_json => faces_json%children
          num_faces = fson_value_count(faces_json)

          do iface = 1, num_faces

             write(istr, '(i0)') iface - 1
             facestr = 'mesh.faces[' // trim(istr) // ']'
             call fson_get_default(face_json, "cells", default_cells, &
                  cell_indices, logfile, log_key = trim(facestr) // ".cells")
             call fson_get_default(face_json, "permeability_direction", &
                  default_permeability_direction, permeability_direction, &
                  logfile, log_key = trim(facestr) // ".permeability_direction")

             associate(num_cells => size(cell_indices))
               if (num_cells == 2) then
                  pcells => cell_indices
                  call DMPlexGetMeet(self%serial_dm, num_cells, pcells, cell_faces, ierr)
                  CHKERRQ(ierr)
                  num_cell_faces = size(cell_faces)
                  if (num_cell_faces == 1) then
                     do i = 1, num_cell_faces
                        f = cell_faces(i)
                        call DMSetLabelValue(self%serial_dm, &
                             face_permeability_override_label_name, f, &
                             permeability_direction, ierr); CHKERRQ(ierr)
                     end do
                  else
                     if (present(logfile)) then
                        call logfile%write(LOG_LEVEL_WARN, "input", &
                             "unrecognised face", int_keys = ["mesh.faces"], &
                             int_values = [iface - 1])
                     end if
                  end if
                  call DMPlexRestoreMeet(self%serial_dm, num_cells, pcells, cell_faces, &
                       ierr); CHKERRQ(ierr)
               else
                  if (present(logfile)) then
                     call logfile%write(LOG_LEVEL_WARN, "input", &
                          "incorrect number of cells", int_keys = ["mesh.faces"], &
                          int_values = [iface - 1])
                  end if
               end if
             end associate

             face_json => face_json%next

          end do
       end if

    end if

  end subroutine mesh_read_overridden_face_properties

!------------------------------------------------------------------------

  subroutine mesh_override_face_properties(self)
    !! Sets any face properties overridden in the JSON input-
    !! currently just face permeability directions.

    use logfile_module
    use dm_utils_module, only: local_vec_section, section_offset
    class(mesh_type), intent(in out) :: self
    ! Locals:
    DMLabel :: label
    IS :: faceIS
    PetscInt :: num_faces
    PetscInt, pointer, contiguous :: faces(:)
    PetscSection :: face_section
    type(face_type) :: face
    PetscReal, pointer, contiguous :: face_geom_array(:)
    PetscInt :: i, face_offset, dirn
    PetscErrorCode :: ierr

    call DMGetLabel(self%dm, face_permeability_override_label_name, &
         label, ierr); CHKERRQ(ierr)

    call local_vec_section(self%face_geom, face_section)
    call VecGetArrayF90(self%face_geom, face_geom_array, ierr); CHKERRQ(ierr)
    call face%init()

    do dirn = 1, self%dim
       call DMLabelGetStratumSize(label, dirn, num_faces, ierr); CHKERRQ(ierr)
       if (num_faces > 0) then
          call DMLabelGetStratumIS(label, dirn, faceIS, ierr); CHKERRQ(ierr)
          call ISGetIndicesF90(faceIS, faces, ierr); CHKERRQ(ierr)
          do i = 1, num_faces
             associate(f => faces(i))
               if (self%ghost_face(f) < 0) then
                  face_offset = section_offset(face_section, f)
                  call face%assign_geometry(face_geom_array, face_offset)
                  face%permeability_direction = dble(dirn)
               end if
             end associate
          end do
          call ISRestoreIndicesF90(faceIS, faces, ierr); CHKERRQ(ierr)
       end if
    end do

    call face%destroy()
    call VecRestoreArrayF90(self%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)

  end subroutine mesh_override_face_properties

!------------------------------------------------------------------------

  subroutine mesh_label_cell_array_zones(self, json)
    !! Labels serial DM for cell array zones, referring to natural
    !! cell indices.

    use fson
    use fson_string_m, only: fson_string_length, fson_string_copy
    use fson_value_m, only : fson_value_count

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscMPIInt :: rank
    type(fson_value), pointer :: zones_json, zone_json
    PetscInt :: num_zones, i, ztype, name_len
    type(zone_cell_array_type) :: zone
    character(:), allocatable :: name
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       call fson_get(json, "mesh.zones", zones_json)
       if (associated(zones_json)) then
          num_zones = fson_value_count(zones_json)
          zone_json => zones_json%children
          do i = 1, num_zones
             ztype = get_zone_type(zone_json)
             if (ztype == ZONE_TYPE_CELL_ARRAY) then
                name_len = fson_string_length(zone_json%name)
                allocate(character(name_len) :: name)
                call fson_string_copy(zone_json%name, name)
                call zone%init_serial(i - 1, name, zone_json)
                call zone%label_serial_dm(self%serial_dm)
                call zone%destroy()
                deallocate(name)
                zone_json => zone_json%next
             end if
          end do
       end if
    end if

  end subroutine mesh_label_cell_array_zones

!------------------------------------------------------------------------

  subroutine mesh_label_cell_array_rock_types(self, json)
    !! Labels serial DM for cell array rock types, referring to
    !! natural cell indices.

    use fson
    use fson_value_m, only : fson_value_count
    use rock_module, only: rock_type_label_name, label_rock_cell

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscMPIInt :: rank
    type(fson_value), pointer :: rocktypes, r
    PetscInt :: start_cell, end_cell, ir, num_rocktypes
    DMLabel :: ghost_label
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call fson_get(json, "rock.types", rocktypes)
       if (associated(rocktypes)) then

          call DMCreateLabel(self%serial_dm, rock_type_label_name, ierr)
          CHKERRQ(ierr)
          call DMPlexGetHeightStratum(self%serial_dm, 0, start_cell, end_cell, &
               ierr); CHKERRQ(ierr)
          call DMGetLabel(self%serial_dm, "ghost", ghost_label, ierr)
          CHKERRQ(ierr)

          num_rocktypes = fson_value_count(rocktypes)
          r => rocktypes%children
          do ir = 1, num_rocktypes
             call label_rock_cells(r, ir)
             r => r%next
          end do

       end if
    end if

  contains

    subroutine label_rock_cells(r, ir)
      !! Sets serial DM rocktype label on specified cells.

      type(fson_value), pointer, intent(in out) :: r
      PetscInt, intent(in) :: ir !! Rock type index
      ! Locals:
      type(fson_value), pointer :: cell_indices_json
      PetscInt, allocatable :: cell_indices(:)
      PetscInt :: ic

      call fson_get(r, "cells", cell_indices_json)
      if (associated(cell_indices_json)) then
         call fson_get(cell_indices_json, ".", cell_indices)
         if (allocated(cell_indices)) then
            associate(num_cells => size(cell_indices))
              do ic = 1, num_cells
                 associate(c => cell_indices(ic))
                   call label_rock_cell(self%serial_dm, &
                        start_cell, end_cell, c, ir)
                 end associate
              end do
            end associate
            deallocate(cell_indices)
         end if
      end if

    end subroutine label_rock_cells

  end subroutine mesh_label_cell_array_rock_types

!------------------------------------------------------------------------

  subroutine mesh_label_cell_array_minc_zones(self, json)
    !! Labels serial DM for cell array MINC zones, referring to
    !! natural cell indices.

    use fson
    use fson_value_m, only : fson_value_count, TYPE_ARRAY, TYPE_OBJECT

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscMPIInt :: rank
    type(fson_value), pointer :: minc_json, minci_json
    PetscInt :: minc_type, iminc, num_minc_zones
    PetscInt :: start_cell, end_cell, minc_rocktype_zone_index
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       call fson_get(json, "mesh.minc", minc_json)
       if (associated(minc_json)) then

          call DMPlexGetHeightStratum(self%serial_dm, 0, start_cell, end_cell, &
               ierr); CHKERRQ(ierr)
          call DMCreateLabel(self%serial_dm, minc_zone_label_name, ierr)
          CHKERRQ(ierr)
          call DMCreateLabel(self%serial_dm, minc_rocktype_zone_label_name, &
               ierr); CHKERRQ(ierr)
          minc_rocktype_zone_index = 1

          minc_type = minc_json%value_type
          select case (minc_type)
          case (TYPE_OBJECT)
             call minc_label_cell_array_zones(minc_json, 1, minc_rocktype_zone_index)
          case (TYPE_ARRAY)
             num_minc_zones = fson_value_count(minc_json)
             minci_json => minc_json%children
             do iminc = 1, num_minc_zones
                call minc_label_cell_array_zones(minci_json, iminc, &
                     minc_rocktype_zone_index)
                minci_json => minci_json%next
             end do
          end select

       end if
    end if

  contains

    subroutine minc_label_cell_array_zones(json, iminc, minc_rocktype_zone_index)
      !! Labels serial DM for cell arrays in a particular MINC zone.

      type(fson_value), pointer, intent(in) :: json
      PetscInt, intent(in) :: iminc
      PetscInt, intent(in out) :: minc_rocktype_zone_index
      ! Locals:
      type(fson_value), pointer :: rock_json, rocki_json, cells_json
      PetscInt :: rock_type, num_rocks, irock, ic
      PetscInt, allocatable :: cells(:)

      call fson_get(json, "rock", rock_json)
      if (associated(rock_json)) then

         rock_type = rock_json%value_type
         select case (rock_type)
         case (TYPE_OBJECT)
            num_rocks = 1
            rocki_json => rock_json
         case (TYPE_ARRAY)
            num_rocks = fson_value_count(rock_json)
            rocki_json => rock_json%children
         end select

         do irock = 1, num_rocks
            call fson_get(rocki_json, "cells", cells_json)
            if (associated(cells_json)) then
               call fson_get(cells_json, ".", cells)
               associate(num_cells => size(cells))
                 do ic = 1, num_cells
                    associate(c => cells(ic))
                      if ((c >= start_cell) .and. (c < end_cell)) then
                         call DMSetLabelValue(self%serial_dm, minc_zone_label_name, &
                              c, iminc, ierr); CHKERRQ(ierr)
                         call DMSetLabelValue(self%serial_dm, &
                              minc_rocktype_zone_label_name, c, &
                              minc_rocktype_zone_index, ierr); CHKERRQ(ierr)
                      end if
                    end associate
                 end do
               end associate
               deallocate(cells)
            end if
            minc_rocktype_zone_index = minc_rocktype_zone_index + 1
            rocki_json => rocki_json%next
         end do

      end if

    end subroutine minc_label_cell_array_zones

  end subroutine mesh_label_cell_array_minc_zones

!------------------------------------------------------------------------

  subroutine mesh_label_boundaries(self, json, logfile)
    !! Labels serial DM for boundary conditions.

    use kinds_module
    use fson
    use fson_value_m, only : fson_value_count, TYPE_ARRAY, TYPE_OBJECT
    use fson_utils_module, only: fson_get_default
    use logfile_module
    use dm_utils_module, only: dm_cell_normal_face, dm_check_create_label

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    PetscMPIInt :: rank
    type(fson_value), pointer :: boundaries_json, bdy_json
    type(fson_value), pointer :: faces_json, face_json, cells_json
    PetscInt :: num_boundaries, num_faces, num_cells, ibdy, offset, i
    PetscInt :: faces_type, num_face_items, face1_type
    PetscInt :: start_cell, end_cell
    PetscInt, allocatable :: faces(:), cells(:)
    PetscInt, allocatable :: default_faces(:), default_cells(:)
    PetscReal, parameter :: default_normal(3) = [0._dp, 0._dp, 1._dp]
    PetscReal, allocatable :: input_normal(:)
    PetscErrorCode :: ierr
    character(len=64) :: bdystr
    character(len=12) :: istr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       default_faces = [PetscInt::] ! empty integer array
       default_cells = [PetscInt::]

       call dm_check_create_label(self%serial_dm, open_boundary_label_name)
       call DMPlexGetHeightStratum(self%serial_dm, 0, start_cell, end_cell, ierr)
       CHKERRQ(ierr)

       call fson_get(json, "boundaries", boundaries_json)
       if (associated(boundaries_json)) then
          num_boundaries = fson_value_count(boundaries_json)
          bdy_json => boundaries_json%children
          num_faces = 0
          do ibdy = 1, num_boundaries
             write(istr, '(i0)') ibdy - 1
             bdystr = 'boundaries[' // trim(istr) // ']'
             call fson_get(bdy_json, "faces", faces_json)
             if (associated(faces_json)) then
                faces_type = faces_json%value_type
                select case (faces_type)
                case (TYPE_ARRAY)
                   num_face_items = fson_value_count(faces_json)
                   if (num_face_items > 0) then
                      face_json => faces_json%children
                      face1_type = face_json%value_type
                      select case (face1_type)
                      case (TYPE_OBJECT)
                         num_faces = 0
                         do i = 1, num_face_items
                            call fson_get(face_json, "cells", cells_json)
                            if (associated(cells_json)) then
                               num_cells = fson_value_count(cells_json)
                               num_faces = num_faces + num_cells
                               face_json => face_json%next
                            end if
                         end do
                         allocate(faces(num_faces))
                         face_json => faces_json%children
                         offset = 0
                         do i = 1, num_face_items
                            call fson_get_default(face_json, "cells", default_cells, cells, &
                                 logfile, log_key = trim(bdystr) // "faces.cells")
                            num_cells = size(cells)
                            call fson_get_default(face_json, "normal", default_normal, &
                                 input_normal, logfile, log_key = trim(bdystr) // "faces.normal")
                            call get_cell_faces(cells, num_cells, input_normal, offset)
                            offset = offset + num_cells
                            deallocate(cells)
                            face_json => face_json%next
                         end do
                      case default
                         if (present(logfile)) then
                            call logfile%write(LOG_LEVEL_WARN, "input", &
                                 "unrecognised_face_type")
                         end if
                      end select
                   end if
                case (TYPE_OBJECT)
                   call fson_get_default(faces_json, "cells", default_cells, cells, &
                        logfile, log_key = trim(bdystr) // "faces.cells")
                   num_cells = size(cells)
                   call fson_get_default(faces_json, "normal", default_normal, &
                        input_normal, logfile, log_key = trim(bdystr) // "faces.normal")
                   num_faces = num_cells
                   allocate(faces(num_faces))
                   call get_cell_faces(cells, num_cells, input_normal, 0)
                   deallocate(cells)
                case default
                   if (present(logfile)) then
                      call logfile%write(LOG_LEVEL_WARN, "input", &
                           "unrecognised_faces_type")
                   end if
                end select
             end if

             do i = 1, num_faces
                associate(f => faces(i))
                  if (f >= 0) then
                     call DMSetLabelValue(self%serial_dm, open_boundary_label_name, &
                          f, ibdy, ierr); CHKERRQ(ierr)
                  end if
                end associate
             end do
             if (allocated(faces)) deallocate(faces)
             bdy_json => bdy_json%next

          end do

       else if (present(logfile)) then
          call logfile%write(LOG_LEVEL_WARN, "input", "no_boundary_conditions")
       end if

    end if

  contains

    subroutine get_cell_faces(cells, num_cells, input_normal, offset)
      ! Get faces for normal vector and specified cells.

      PetscInt, intent(in) :: cells(:), num_cells
      PetscReal, intent(in) :: input_normal(:)
      PetscInt, intent(in) :: offset
      ! Locals:
      PetscInt :: f, normal_len, icell, iface
      PetscReal :: normal(3)

      normal_len = size(input_normal)
      normal = 0._dp
      normal(1: normal_len) = input_normal

      do icell = 1, num_cells
         iface = offset + icell
         associate(c => cells(icell))
           if (c >= 0) then
              if ((start_cell <= c) .and. (c < end_cell)) then
                 call dm_cell_normal_face(self%serial_dm, c, normal, f)
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
              else
                 faces(iface) = -1
              end if
           else
              faces(iface) = -1
           end if
         end associate
      end do

    end subroutine get_cell_faces

  end subroutine mesh_label_boundaries

!------------------------------------------------------------------------

  subroutine mesh_label_sources(self, json)
    !! Labels serial DM for source locations defined by cell indices
    !! (not zones).

    use source_module
    use fson_value_m, only : fson_value_count, TYPE_INTEGER, TYPE_ARRAY
    use dm_utils_module, only: dm_check_create_label

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    ! Locals:
    PetscMPIInt :: rank
    type(fson_value), pointer :: sources_json, source_json
    PetscInt :: source_index, num_source_specs, start_cell, end_cell
    DMLabel :: source_label
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call fson_get(json, "source", sources_json)
       if (associated(sources_json)) then

          call dm_check_create_label(self%serial_dm, source_label_name)
          call DMGetLabel(self%serial_dm, source_label_name, &
               source_label, ierr); CHKERRQ(ierr)
          call DMPlexGetHeightStratum(self%serial_dm, 0, start_cell, &
               end_cell, ierr); CHKERRQ(ierr)

          num_source_specs = fson_value_count(sources_json)
          source_json => sources_json%children
          do source_index = 0, num_source_specs - 1
             call label_source_cells(source_json, source_index)
             source_json => source_json%next
          end do

       end if

    end if

  contains

    subroutine label_source_cells(json, source_index)

      type(fson_value), pointer, intent(in) :: json
      PetscInt, intent(in) :: source_index
      ! Locals:
      type(fson_value), pointer :: cell_json, cells_json
      PetscInt :: cell, i
      PetscInt, allocatable :: cells(:)
      PetscErrorCode :: ierr

      call fson_get(json, "cell", cell_json)
      if (associated(cell_json)) then
         call fson_get(cell_json, ".", cell)
         cells = [cell]
      end if

      call fson_get(json, "cells", cells_json)
      if (associated(cells_json)) then
         select case (cells_json%value_type)
         case (TYPE_INTEGER)
            call fson_get(cells_json, ".", cell)
            cells = [cell]
          case (TYPE_ARRAY)
            call fson_get(cells_json, ".", cells)
         end select
      end if

      if (allocated(cells)) then
         do i = 1, size(cells)
            associate(c => cells(i))
              if ((start_cell <= c) .and. (c < end_cell)) then
                 call DMSetLabelValue(self%serial_dm, source_label_name, &
                      c, source_index, ierr); CHKERRQ(ierr)
              end if
            end associate
         end do
         deallocate(cells)
      end if

    end subroutine label_source_cells

  end subroutine mesh_label_sources

!------------------------------------------------------------------------

  subroutine mesh_setup_zones(self, json, logfile, err)
    !! Sets up zones (for defining e.g. rock types, MINC etc.) in the
    !! mesh, from JSON input.

    use fson_mpi_module
    use logfile_module
    use zone_label_module, only: max_zone_name_length
    use dictionary_module
    use dag_module

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: num_zones, i, ztype
    type(fson_value), pointer :: zones_json, zone_json
    character(:), allocatable :: name
    class(zone_type), pointer :: zone
    type(list_node_type), pointer :: node
    type(dictionary_type) :: zone_dict
    type(dag_type) :: dag
    PetscInt, allocatable :: order(:)
    character(max_zone_name_length), allocatable :: zone_names(:)
    character(max_zone_name_length) :: zone_name

    err = 0
    call self%zones%init(owner = PETSC_TRUE)

    if (fson_has_mpi(json, "mesh.zones")) then

       call fson_get_mpi(json, "mesh.zones", zones_json)

       num_zones = fson_value_count_mpi(zones_json, ".")
       zone_json => fson_value_children_mpi(zones_json)

       do i = 1, num_zones
          zone => null()
          ztype = get_zone_type_mpi(zone_json)
          name = fson_get_name_mpi(zone_json)
          select case (ztype)
          case (ZONE_TYPE_CELL_ARRAY)
             allocate(zone_cell_array_type :: zone)
          case (ZONE_TYPE_BOX)
             allocate(zone_box_type :: zone)
          case (ZONE_TYPE_COMBINE)
             allocate(zone_combine_type :: zone)
          case default
             err = 1
             if (present(logfile)) then
                call logfile%write(LOG_LEVEL_ERR, 'input', &
                     "unrecognised zone type", str_key = "name", &
                     str_value = name)
             end if
          end select
          if (err == 0) then
             call zone%init(i - 1, name, zone_json)
             if (err == 0) then
                call self%zones%append(zone, name)
             else
                if (present(logfile)) then
                   call logfile%write(LOG_LEVEL_ERR, 'input', &
                        "unrecognised zone dependency", &
                        str_key = "name", str_value = name)
                end if
                exit
             end if
          else
             exit
          end if
          zone_json => fson_value_next_mpi(zone_json)
       end do

       call zone_dict%init(self%zones)
       ! Set up zone dependency graph and do topological sort:
       call dag%init(num_zones)
       call self%zones%traverse(zone_dependency_iterator)
       call dag%sort(order, err)

       if (err == 0) then
          call self%zones%tags(zone_names)
          do i = 0, num_zones - 1
             zone_name = zone_names(order(i))
             node => zone_dict%get(zone_name)
             select type(zone => node%data)
                class is (zone_type)
                   call zone%label_dm(self%dm, self%cell_geom, err)
                if (err > 0) then
                   if (present(logfile)) then
                      call logfile%write(LOG_LEVEL_WARN, 'zone', &
                           "can't find cells", str_key = "name", &
                           str_value = zone_name)
                   end if
                   exit
                end if
             end select
          end do

       else
          if (present(logfile)) then
             call logfile%write(LOG_LEVEL_ERR, 'input', &
                  "circular zone dependency")
          end if
       end if

       call zone_dict%destroy()

    end if

  contains

    subroutine zone_dependency_iterator(node, stopped)
      !! Adds dependencies for zone to dependency graph.

      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      character(max_zone_name_length), allocatable :: depends(:)
      PetscInt :: i, num_depends
      PetscInt, allocatable :: edges(:)
      type(list_node_type), pointer :: depend_node

      stopped = PETSC_FALSE
      select type (zone => node%data)
      class is (zone_type)
         call zone%dependencies(depends)
         num_depends = size(depends)
         allocate(edges(num_depends))
         do i = 1, size(depends)
            depend_node => zone_dict%get(depends(i))
            select type (depend_zone => depend_node%data)
            class is (zone_type)
               edges(i) = depend_zone%index
            end select
         end do
         call dag%set_edges(zone%index, edges)
      end select

    end subroutine zone_dependency_iterator
    
  end subroutine mesh_setup_zones

!------------------------------------------------------------------------

  subroutine mesh_setup_minc(self, json, logfile, err)
    !! Sets up MINC for fractured zones.

    use fson_mpi_module
    use logfile_module
    use fson_value_m, only : TYPE_ARRAY, TYPE_OBJECT
    use dictionary_module

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: num_minc_zones, iminc, num_minc_cells, minc_rocktype_zone_index
    type(fson_value), pointer :: minc_json, minci_json
    PetscInt :: minc_type
    character(32) :: imincstr, mincstr
    PetscMPIInt :: rank
    PetscErrorCode :: ierr
    PetscBool :: has_minc_local
    PetscBool :: default_rebalance

    err = 0

    if (fson_has_mpi(json, "mesh.minc")) then

       call fson_get_mpi(json, "mesh.minc", minc_json)
       minc_type = fson_type_mpi(minc_json, ".")
       minc_rocktype_zone_index = 1

       select case (minc_type)
       case (TYPE_OBJECT)
          num_minc_zones = 1
          allocate(self%minc(num_minc_zones))
          iminc = 1
          mincstr = "minc."
          call self%minc(num_minc_zones)%init(minc_json, self%dm, &
               iminc, mincstr, self%rock_types, &
               minc_rocktype_zone_index, logfile, err)
       case (TYPE_ARRAY)
          num_minc_zones = fson_value_count_mpi(minc_json, ".")
          minci_json => fson_value_children_mpi(minc_json)
          allocate(self%minc(num_minc_zones))
          do iminc = 1, num_minc_zones
             write(imincstr, '(i0)') iminc - 1
             mincstr = 'minc[' // trim(imincstr) // '].'
             call self%minc(iminc)%init(minci_json, self%dm, iminc, &
                  mincstr, self%rock_types, minc_rocktype_zone_index, logfile, err)
             if (err > 0) exit
             minci_json => fson_value_next_mpi(minci_json)
          end do
       end select

       call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
       call DMGetLabelSize(self%dm, minc_zone_label_name, &
            num_minc_cells, ierr); CHKERRQ(ierr)
       has_minc_local = (num_minc_cells > 0)
       call MPI_allreduce(has_minc_local, self%has_minc, 1, MPI_LOGICAL, MPI_LOR, &
            PETSC_COMM_WORLD, ierr)

       default_rebalance = PETSC_TRUE

    else
       self%has_minc = PETSC_FALSE
       default_rebalance = PETSC_FALSE
    end if

    call fson_get_mpi(json, "mesh.rebalance", default_rebalance, &
         self%rebalance, logfile)

  end subroutine mesh_setup_minc

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_dm(self)
    !! Sets up augmented DM for MINC mesh, including MINC cells and
    !! faces. Although they are not used, edges and vertices are also
    !! created for the MINC faces, so that the depth of the DM is
    !! consistent everywhere.

    use dm_utils_module
    use list_module

    class(mesh_type), intent(in out) :: self
    ! Locals:
    DM :: minc_dm
    PetscInt :: coord_dim, start_chart, end_chart
    PetscBool :: balance
    PetscInt :: num_minc_cells, minc_end_interior_cell
    PetscInt :: num_new_points, max_num_levels
    PetscInt, allocatable :: stratum_shift(:)
    PetscInt, allocatable :: minc_level_cells(:,:)
    PetscInt, allocatable :: minc_level_cell_count(:)
    PetscErrorCode :: ierr

    call dm_get_strata(self%dm, self%depth, self%strata)
    allocate(stratum_shift(0: self%depth))

    call DMPlexCreate(PETSC_COMM_WORLD, minc_dm, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(minc_dm, 'MINC mesh', ierr); CHKERRQ(ierr)

    call DMSetDimension(minc_dm, self%dim, ierr); CHKERRQ(ierr)
    call DMGetCoordinateDim(self%dm, coord_dim, ierr); CHKERRQ(ierr)
    call DMSetCoordinateDim(minc_dm, coord_dim, ierr); CHKERRQ(ierr)
    call DMPlexGetPartitionBalance(self%dm, balance, ierr); CHKERRQ(ierr)
    call DMPlexSetPartitionBalance(minc_dm, balance, ierr); CHKERRQ(ierr)

    call DMPlexGetChart(self%dm, start_chart, end_chart, ierr)
    CHKERRQ(ierr)

    max_num_levels = maxval(self%minc%num_levels)

    call self%setup_minc_level_cells(max_num_levels, &
         minc_level_cells, minc_level_cell_count, num_minc_cells)

    num_new_points = num_minc_cells * (self%depth + 1)
    call DMPlexSetChart(minc_dm, start_chart, &
         end_chart + num_new_points, ierr); CHKERRQ(ierr)
    call self%setup_minc_dm_strata_shifts(num_minc_cells, &
         max_num_levels, minc_level_cell_count, stratum_shift)

    minc_end_interior_cell = self%strata(0)%end_interior + &
         (stratum_shift(0) + 1) * num_minc_cells
    call DMPlexSetGhostCellStratum(minc_dm, minc_end_interior_cell, -1, &
         ierr); CHKERRQ(ierr)
    self%strata%num_minc_points = num_minc_cells

    call self%set_minc_dm_cone_sizes(minc_dm, num_minc_cells, &
         max_num_levels, minc_level_cells)
    call DMSetUp(minc_dm, ierr); CHKERRQ(ierr)
    call self%set_minc_dm_cones(minc_dm, max_num_levels, &
         minc_level_cells)

    call DMPlexSymmetrize(minc_dm, ierr); CHKERRQ(ierr)
    call self%transfer_labels_to_minc_dm(minc_dm, max_num_levels)
    call self%setup_minc_dm_depth_label(minc_dm, max_num_levels, &
         minc_level_cells)
    call self%setup_minc_output_data(minc_dm, max_num_levels, minc_level_cells)

    call dm_set_fv_adjacency(minc_dm)
    call dm_set_default_data_layout(minc_dm, self%dof)
    call self%setup_minc_point_sf(minc_dm)
    call dm_setup_global_section(minc_dm)

    call self%setup_minc_dm_cell_natural_global(minc_dm, max_num_levels, &
         minc_level_cells, minc_level_cell_count, num_minc_cells)
    self%cell_natural = dm_natural_order_IS(minc_dm, self%cell_natural_global)
    call self%setup_minc_geometry(minc_dm, max_num_levels, &
         minc_level_cells)
    call self%setup_minc_coordinates(minc_dm)

    deallocate(minc_level_cells, minc_level_cell_count, &
         stratum_shift)

    self%dm = minc_dm

  end subroutine mesh_setup_minc_dm

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_level_cells(self, &
       max_num_levels, minc_level_cells, &
       minc_level_cell_count, num_minc_cells)
    !! Set up minc_level_cells and minc_level_cell_count
    !! arrays.

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: max_num_levels
    PetscInt, allocatable, intent(out) :: minc_level_cells(:,:)
    PetscInt, allocatable, intent(out) :: minc_level_cell_count(:)
    PetscInt, intent(out) :: num_minc_cells
    ! Locals:
    PetscInt :: iminc, i, c, m, ghost, num_minc_zone_cells
    IS :: minc_IS
    PetscInt, pointer :: minc_cells(:)
    DMLabel :: ghost_label
    PetscErrorCode :: ierr

    call DMGetLabel(self%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)

    allocate(minc_level_cells(max_num_levels, &
         self%strata(0)%start: self%strata(0)%end - 1))
    minc_level_cells = -1
    allocate(minc_level_cell_count(max_num_levels))
    minc_level_cell_count = 0

    do iminc = 1, size(self%minc)
       call DMGetStratumSize(self%dm, minc_zone_label_name, iminc, &
            num_minc_zone_cells, ierr); CHKERRQ(ierr)
       if (num_minc_zone_cells > 0) then
          call DMGetStratumIS(self%dm, minc_zone_label_name, &
               iminc, minc_IS, ierr); CHKERRQ(ierr)
          call ISGetIndicesF90(minc_IS, minc_cells, ierr); CHKERRQ(ierr)
          associate(num_levels => self%minc(iminc)%num_levels)
            do i = 1, num_minc_zone_cells
               c = minc_cells(i)
               call DMLabelGetValue(ghost_label, c, ghost, ierr)
               if (ghost < 0) then
                  do m = 1, self%minc(iminc)%num_levels
                     minc_level_cells(m, c) = minc_level_cell_count(m)
                     minc_level_cell_count(m) = minc_level_cell_count(m) + 1
                  end do
               end if
            end do
          end associate
          call ISRestoreIndicesF90(minc_IS, minc_cells, ierr)
          CHKERRQ(ierr)
          call ISDestroy(minc_IS, ierr); CHKERRQ(ierr)
       end if
    end do

    num_minc_cells = sum(minc_level_cell_count)

  end subroutine mesh_setup_minc_level_cells

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_dm_strata_shifts(self, num_minc_cells, &
       max_num_levels, minc_level_cell_count, stratum_shift)
    !! Set up minc_shift array in DM strata, to determine index shift
    !! from original DM point to corresponding point in MINC DM.

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: num_minc_cells, max_num_levels
    PetscInt, intent(in) :: minc_level_cell_count(max_num_levels)
    PetscInt, intent(out) :: stratum_shift(0: self%depth)
    ! Locals:
    PetscInt :: h, m
    PetscInt :: minc_offset(0: max_num_levels)

    ! Set up stratum_shift array (taking account of the fact that
    ! DMPlex points have the order cells, vertices, faces, edges-
    ! i.e. they are not in depth order.) This represents cumulative shift
    ! resulting from points added to lower-index strata in the DM.
    stratum_shift(0) = 0
    stratum_shift(self%depth) = 1
    stratum_shift(1: self%depth - 1) = [(h + 1, h = 1, self%depth - 1)]

    ! Offset of the start of each MINC level within the stratum (note
    ! this is the same for all strata):
    minc_offset = 0
    do m = 1, max_num_levels
       minc_offset(m) = minc_offset(m - 1) + minc_level_cell_count(m)
    end do

    do h = 0, self%depth
       allocate(self%strata(h)%minc_shift(0: max_num_levels))
       self%strata(h)%minc_shift(0) = stratum_shift(h) * num_minc_cells
       do m = 1, max_num_levels
          self%strata(h)%minc_shift(m) = stratum_shift(h) * num_minc_cells + &
               self%strata(h)%end_non_ghost + minc_offset(m - 1)
       end do
    end do

  end subroutine mesh_setup_minc_dm_strata_shifts

!------------------------------------------------------------------------

  subroutine mesh_set_minc_dm_cone_sizes(self, minc_dm, &
       num_minc_cells, max_num_levels, &
       minc_level_cells)
    !! Sets cone sizes for MINC DM.

    class(mesh_type), intent(in out) :: self
    DM, intent(in out) :: minc_dm
    PetscInt, intent(in) :: num_minc_cells, max_num_levels
    PetscInt, intent(in) :: minc_level_cells(max_num_levels, &
         self%strata(0)%start: self%strata(0)%end - 1)
    ! Locals:
    PetscInt :: ic, p, minc_p, orig_cone_size, cone_size
    PetscInt :: iminc, m, h, c, ghost
    DMLabel :: minc_zone_label, ghost_label
    PetscErrorCode :: ierr

    call DMGetLabel(self%dm, minc_zone_label_name, minc_zone_label, &
         ierr); CHKERRQ(ierr)
    call DMGetLabel(self%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)

    do c = self%strata(0)%start, self%strata(0)%end - 1

       minc_p = self%strata(0)%minc_point(c, 0)
       call DMLabelGetValue(minc_zone_label, c, iminc, ierr); CHKERRQ(ierr)
       call DMLabelGetValue(ghost_label, c, ghost, ierr)
       call DMPlexGetConeSize(self%dm, c, orig_cone_size, ierr); CHKERRQ(ierr)
       if ((iminc > 0) .and. (ghost < 0)) then

          ! Fracture cells:
          cone_size = orig_cone_size + 1
          call DMPlexSetConeSize(minc_dm, minc_p, cone_size, ierr)
          CHKERRQ(ierr)

          do m = 1, self%minc(iminc)%num_levels
             ! MINC cells:
             ic = minc_level_cells(m, c)
             minc_p = self%strata(0)%minc_point(ic, m)
             if (m < self%minc(iminc)%num_levels) then
                cone_size = 2
             else
                cone_size = 1
             end if
             call DMPlexSetConeSize(minc_dm, minc_p, cone_size, &
                  ierr); CHKERRQ(ierr)
             ! MINC DAG points (height h > 0):
             do h = 1, self%depth -1
                minc_p = self%strata(h)%minc_point(ic, m)
                cone_size = 1
                call DMPlexSetConeSize(minc_dm, minc_p, &
                     cone_size, ierr); CHKERRQ(ierr)
             end do
          end do

       else
          ! Non-MINC and ghost cells:
          cone_size = orig_cone_size
          call DMPlexSetConeSize(minc_dm, minc_p, cone_size, ierr)
          CHKERRQ(ierr)
       end if

    end do

    ! Non-MINC and fracture DAG points for height h > 0:
    do h = 1, self%depth
       do p = self%strata(h)%start, self%strata(h)%end - 1
          minc_p = self%strata(h)%minc_point(p, 0)
          call DMPlexGetConeSize(self%dm, p, orig_cone_size, ierr)
          CHKERRQ(ierr)
          cone_size = orig_cone_size
          call DMPlexSetConeSize(minc_dm, minc_p, cone_size, ierr)
          CHKERRQ(ierr)
       end do
    end do

  end subroutine mesh_set_minc_dm_cone_sizes

!------------------------------------------------------------------------

  subroutine mesh_set_minc_dm_cones(self, minc_dm, max_num_levels, &
       minc_level_cells)
    !! Sets cones and cone orientations for MINC DM.

    use dm_utils_module, only: dm_stratum_type, dm_copy_cone_orientation

    class(mesh_type), intent(in out) :: self
    DM, intent(in out) :: minc_dm
    PetscInt, intent(in) :: max_num_levels
    PetscInt, intent(in) :: minc_level_cells(max_num_levels, &
         self%strata(0)%start: self%strata(0)%end - 1)
    ! Locals:
    PetscInt :: c, p, m, h, iminc, ic, ic_m1, ghost
    DMLabel :: minc_zone_label, ghost_label
    PetscInt :: minc_p, above_p, face_p, inner_face_p
    PetscInt, pointer :: points(:)
    PetscInt, allocatable :: cell_cone(:), minc_cone(:), minc_orientation(:)
    PetscInt, pointer :: orientation(:)
    PetscErrorCode :: ierr

    call DMGetLabel(self%dm, minc_zone_label_name, minc_zone_label, &
         ierr); CHKERRQ(ierr)
    call DMGetLabel(self%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)

    do c = self%strata(0)%start, self%strata(0)%end - 1

       minc_p = self%strata(0)%minc_point(c, 0)
       call DMPlexGetCone(self%dm, c, points, ierr); CHKERRQ(ierr)
       call DMLabelGetValue(minc_zone_label, c, iminc, ierr); CHKERRQ(ierr)
       call DMLabelGetValue(ghost_label, c, ghost, ierr)
       if ((iminc > 0) .and. (ghost < 0)) then

          ! Fracture cells:
          ic_m1 = minc_level_cells(1, c)
          face_p = self%strata(1)%minc_point(ic_m1, 1)
          cell_cone = [self%strata(1)%minc_point(points, 0), [face_p]]
          call DMPlexSetCone(minc_dm, minc_p, cell_cone, ierr); CHKERRQ(ierr)
          deallocate(cell_cone)
          call DMPlexGetConeOrientation(self%dm, c, orientation, ierr); CHKERRQ(ierr)
          minc_orientation = [orientation, [0]]
          call DMPlexSetConeOrientation(minc_dm, minc_p, minc_orientation, &
               ierr); CHKERRQ(ierr)
          call DMPlexRestoreConeOrientation(self%dm, c, orientation, &
               ierr); CHKERRQ(ierr)
          deallocate(minc_orientation)

          do m = 1, self%minc(iminc)%num_levels
             ! MINC cells:
             ic = minc_level_cells(m, c)
             minc_p = self%strata(0)%minc_point(ic, m)
             face_p = self%strata(1)%minc_point(ic, m)
             if (m < self%minc(iminc)%num_levels) then
                ic_m1 = minc_level_cells(m + 1, c)
                inner_face_p = self%strata(1)%minc_point(ic_m1, m + 1)
                minc_cone = [face_p, inner_face_p]
             else
                minc_cone = [face_p]
             end if
             call DMPlexSetCone(minc_dm, minc_p, minc_cone, ierr); CHKERRQ(ierr)
             deallocate(minc_cone)
             ! MINC DAG points for height h > 0:
             do h = 1, self%depth - 1
                minc_p = self%strata(h)%minc_point(ic, m)
                above_p = self%strata(h + 1)%minc_point(ic, m)
                call DMPlexSetCone(minc_dm, minc_p, [above_p], ierr)
                CHKERRQ(ierr)
             end do
          end do

       else
          ! Non-MINC and ghost cells:
          cell_cone = self%strata(1)%minc_point(points, 0)
          call DMPlexSetCone(minc_dm, minc_p, cell_cone, ierr); CHKERRQ(ierr)
          deallocate(cell_cone)
          call dm_copy_cone_orientation(self%dm, c, minc_dm, minc_p)
          CHKERRQ(ierr)
       end if

       call DMPlexRestoreCone(self%dm, c, points, ierr); CHKERRQ(ierr)

    end do

    ! Non-MINC and fracture DAG points for height h > 0:
    do h = 1, self%depth - 1
       do p = self%strata(h)%start, self%strata(h)%end - 1
          minc_p = self%strata(h)%minc_point(p, 0)
          call DMPlexGetCone(self%dm, p, points, ierr); CHKERRQ(ierr)
          minc_cone = self%strata(h + 1)%minc_point(points, 0)
          call DMPlexSetCone(minc_dm, minc_p, minc_cone, ierr); CHKERRQ(ierr)
          deallocate(minc_cone)
          call DMPlexRestoreCone(self%dm, p, points, ierr); CHKERRQ(ierr)
          call dm_copy_cone_orientation(self%dm, p, minc_dm, minc_p)
       end do
    end do

  end subroutine mesh_set_minc_dm_cones

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_dm_depth_label(self, minc_dm, max_num_levels, &
       minc_level_cells)
    !! Set DAG depth label for MINC points added to the DM.

    class(mesh_type), intent(in out) :: self
    DM, intent(in out) :: minc_dm
    PetscInt, intent(in) :: max_num_levels
    PetscInt, intent(in) :: minc_level_cells(max_num_levels, &
         self%strata(0)%start: self%strata(0)%end - 1)
    ! Locals:
    PetscInt :: iminc, num_minc_zone_cells
    PetscInt :: i, c, m, ic, minc_p, h
    IS :: minc_IS
    PetscInt, pointer :: minc_cells(:)
    PetscErrorCode :: ierr

    do iminc = 1, size(self%minc)
       call DMGetStratumSize(self%dm, minc_zone_label_name, iminc, &
            num_minc_zone_cells, ierr); CHKERRQ(ierr)
       if (num_minc_zone_cells > 0) then
          call DMGetStratumIS(self%dm, minc_zone_label_name, &
               iminc, minc_IS, ierr); CHKERRQ(ierr)
          call ISGetIndicesF90(minc_IS, minc_cells, ierr); CHKERRQ(ierr)
          associate(num_levels => self%minc(iminc)%num_levels)
            do i = 1, num_minc_zone_cells
               c = minc_cells(i)
               do m = 1, self%minc(iminc)%num_levels
                  ic = minc_level_cells(m, c)
                  do h = 0, self%depth
                     minc_p = self%strata(h)%minc_point(ic, m)
                     associate(p_depth => self%depth - h)
                       call DMSetLabelValue(minc_dm, "depth", &
                            minc_p, p_depth, ierr); CHKERRQ(ierr)
                     end associate
                  end do
               end do
            end do
          end associate
       end if
    end do

  end subroutine mesh_setup_minc_dm_depth_label

!------------------------------------------------------------------------

  subroutine mesh_transfer_labels_to_minc_dm(self, minc_dm, max_num_levels)
    !! Transfers labels from original DM to MINC DM fracture points,
    !! applying appropriate shifts to the point indices.

    use dm_utils_module, only: dm_point_stratum_height

    class(mesh_type), intent(in out) :: self
    DM, intent(in out) :: minc_dm
    PetscInt, intent(in) :: max_num_levels

    ! Locals:
    PetscInt :: p, minc_p, h, l, iid, ip, label_value
    PetscInt :: num_ids, num_points
    PetscInt, parameter :: max_label_name_length = 80
    character(max_label_name_length) :: label_name
    PetscInt :: num_labels
    IS :: id_IS, point_IS
    PetscInt, pointer :: ids(:), points(:)
    PetscErrorCode :: ierr

    call DMGetNumLabels(self%dm, num_labels, ierr); CHKERRQ(ierr)

    do l = 0, num_labels - 1
       call DMGetLabelName(self%dm, l, label_name, ierr); CHKERRQ(ierr)
       call DMCreateLabel(minc_dm, label_name, ierr); CHKERRQ(ierr)
       call DMGetLabelIdIS(self%dm, label_name, id_IS, ierr); CHKERRQ(ierr)
       call ISGetLocalSize(id_IS, num_ids, ierr); CHKERRQ(ierr)
       if (num_ids > 0) then
          call ISGetIndicesF90(id_IS, ids, ierr); CHKERRQ(ierr)
          do iid = 1, num_ids
             label_value = ids(iid)
             call DMGetStratumIS(self%dm, label_name, ids(iid), point_IS, &
                  ierr); CHKERRQ(ierr)
             call ISGetLocalSize(point_IS, num_points, ierr)
             CHKERRQ(ierr)
             if (num_points > 0) then
                call ISGetIndicesF90(point_IS, points, ierr)
                CHKERRQ(ierr)
                do ip = 1, num_points
                   p = points(ip)
                   h = dm_point_stratum_height(self%strata, p)
                   minc_p = self%strata(h)%minc_point(p, 0)
                   call DMSetLabelValue(minc_dm, label_name, &
                        minc_p, label_value, ierr)
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
    end do

  end subroutine mesh_transfer_labels_to_minc_dm

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_output_data(self, minc_dm, max_num_levels, &
       minc_level_cells)
    !! Sets up minc_level DM label on MINC DM, and
    !! cell_parent_natural IS, both used for output purposes. The
    !! minc_level label contains the MINC level assigned to each
    !! cell. Non-MINC cells are assigned level 0 (as are fracture
    !! cells in MINC zones). The cell_parent_natural IS contains,
    !! for each cell, the natural index of the corresponding original
    !! single-porosity cell.

    use dm_utils_module, only: local_to_natural_cell_index

    class(mesh_type), intent(in out) :: self
    DM, intent(in out) :: minc_dm
    PetscInt, intent(in) :: max_num_levels
    PetscInt, intent(in) :: minc_level_cells(max_num_levels, &
         self%strata(0)%start: self%strata(0)%end - 1)
    ! Locals:
    PetscInt :: m, ic, i, c, ghost, minc_p, start_cell, end_cell
    PetscInt, allocatable :: natural(:)
    PetscInt :: iminc, num_minc_zone_cells
    IS :: minc_IS
    PetscInt, pointer :: minc_cells(:)
    DMLabel :: ghost_label
    ISLocalToGlobalMapping :: l2g
    PetscErrorCode :: ierr

    call DMGetLabel(minc_dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)
    call DMCreateLabel(minc_dm, minc_level_label_name, ierr)
    CHKERRQ(ierr)
    call DMGetLocalToGlobalMapping(self%dm, l2g, ierr); CHKERRQ(ierr)

    call DMPlexGetHeightStratum(minc_dm, 0, start_cell, end_cell, &
         ierr); CHKERRQ(ierr)

    allocate(natural(start_cell: end_cell - 1))
    natural = -1

    ! Non-MINC and fracture cells:
    do c = self%strata(0)%start, self%strata(0)%end - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr)
       if (ghost < 0) then
          minc_p = self%strata(0)%minc_point(c, 0)
          call DMSetLabelValue(minc_dm, minc_level_label_name, &
               minc_p, 0, ierr); CHKERRQ(ierr)
          natural(c) = local_to_natural_cell_index( &
               self%cell_natural_global, l2g, c)
       end if
    end do

    ! MINC cells:
    do iminc = 1, size(self%minc)
       call DMGetStratumSize(self%dm, minc_zone_label_name, iminc, &
            num_minc_zone_cells, ierr); CHKERRQ(ierr)
       if (num_minc_zone_cells > 0) then
          call DMGetStratumIS(self%dm, minc_zone_label_name, &
               iminc, minc_IS, ierr); CHKERRQ(ierr)
          call ISGetIndicesF90(minc_IS, minc_cells, ierr); CHKERRQ(ierr)
          associate(num_levels => self%minc(iminc)%num_levels)
            do i = 1, num_minc_zone_cells
               c = minc_cells(i)
               call DMLabelGetValue(ghost_label, c, ghost, ierr)
               if (ghost < 0) then
                  do m = 1, self%minc(iminc)%num_levels
                     ic = minc_level_cells(m, c)
                     minc_p = self%strata(0)%minc_point(ic, m)
                     call DMSetLabelValue(minc_dm, minc_level_label_name, &
                          minc_p, m, ierr); CHKERRQ(ierr)
                     natural(minc_p) = local_to_natural_cell_index( &
                          self%cell_natural_global, l2g, c)
                  end do
               end if
            end do
          end associate
       end if
    end do

    call ISCreateGeneral(PETSC_COMM_WORLD, end_cell - start_cell, &
         natural, PETSC_COPY_VALUES, self%cell_parent_natural, ierr); CHKERRQ(ierr)
    deallocate(natural)

  end subroutine mesh_setup_minc_output_data

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_dm_cell_natural_global(self, minc_dm, &
       max_num_levels, minc_level_cells, minc_level_cell_count, &
       num_minc_cells)

    !! Sets up natural-to-global AO for MINC DM, and overwrites
    !! original AO.

    use dm_utils_module
    use utils_module, only: array_cumulative_sum, get_mpi_int_gather_array

    class(mesh_type), intent(in out) :: self
    DM, intent(in) :: minc_dm
    PetscInt, intent(in) :: max_num_levels
    PetscInt, intent(in) :: minc_level_cells(max_num_levels, &
         self%strata(0)%start: self%strata(0)%end - 1)
    PetscInt, intent(in) :: minc_level_cell_count(max_num_levels)
    PetscInt, intent(in) :: num_minc_cells
    ! Locals:
    AO :: minc_ao
    PetscMPIInt :: rank, num_procs
    PetscInt, allocatable :: natural(:), global(:)
    PetscInt :: start_cell, end_interior_cell, offset, n_all
    PetscInt :: mapping_count, num_ghost_cells, num_non_ghost_cells
    ISLocalToGlobalMapping :: l2g, minc_l2g
    PetscInt :: total_minc_cell_count
    PetscInt :: ic, c, m, inatural, h, p(1)
    PetscInt, allocatable :: minc_global(:), minc_frac_natural(:)
    PetscInt, allocatable :: minc_global_all(:), minc_frac_natural_all(:)
    PetscInt, allocatable :: minc_counts(:), minc_displacements(:)
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD, num_procs, ierr)

    start_cell = self%strata(0)%start
    end_interior_cell = self%strata(0)%end_interior
    num_ghost_cells = dm_get_num_partition_ghost_points(self%dm, 0)
    num_non_ghost_cells = end_interior_cell - start_cell - num_ghost_cells
    ! Each rank has its own array elements for MINC level 0:
    mapping_count = num_non_ghost_cells

    ! Rank 0 has array elements for MINC cells from all ranks:
    call MPI_reduce(num_minc_cells, total_minc_cell_count, 1, &
         MPI_INTEGER, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
       mapping_count = mapping_count + total_minc_cell_count
    end if
    allocate(natural(mapping_count), global(mapping_count))

    call DMGetLocalToGlobalMapping(self%dm, l2g, ierr); CHKERRQ(ierr)
    call DMGetLocalToGlobalMapping(minc_dm, minc_l2g, ierr); CHKERRQ(ierr)

    call get_original_cell_natural_global_indices()
    inatural = get_new_natural_index()
    offset = num_non_ghost_cells + 1

    ! MINC cells (level > 0):
    minc_counts = get_mpi_int_gather_array()
    minc_displacements = get_mpi_int_gather_array()
    h = 0

    do m = 1, max_num_levels
       associate(n => minc_level_cell_count(m))

         allocate(minc_frac_natural(n), minc_global(n))
         call MPI_gather(n, 1, MPI_INTEGER, minc_counts, 1, &
              MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
         if (rank == 0) then
            minc_displacements = [[0], &
                 array_cumulative_sum(minc_counts(1: num_procs - 1))]
            n_all = sum(minc_counts)
         else
            n_all = 1
         end if

         do c = self%strata(0)%start, self%strata(0)%end - 1
            ic = minc_level_cells(m, c)
            if (ic >= 0) then
               p = self%strata(0)%minc_point(ic, m)
               ic = ic + 1
               minc_frac_natural(ic) = local_to_natural_cell_index( &
                    self%cell_natural_global, l2g, c)
               call ISLocalToGlobalMappingApplyBlock(minc_l2g, 1, p, &
                    minc_global(ic:ic), ierr); CHKERRQ(ierr)
            end if
         end do

         allocate(minc_global_all(n_all), minc_frac_natural_all(n_all))
         call MPI_gatherv(minc_global, n, MPI_INTEGER, &
              minc_global_all, minc_counts, minc_displacements, &
              MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
         call MPI_gatherv(minc_frac_natural, n, MPI_INTEGER, &
              minc_frac_natural_all, minc_counts, minc_displacements, &
              MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)

         call assign_minc_natural_indices(n_all, minc_frac_natural_all, &
              minc_global_all, natural, global, inatural, offset)

         deallocate(minc_frac_natural, minc_global, minc_global_all, &
              minc_frac_natural_all)

       end associate
    end do
    deallocate(minc_counts, minc_displacements)

    call AOCreateMapping(PETSC_COMM_WORLD, mapping_count, &
         natural, global, minc_ao, ierr); CHKERRQ(ierr)
    deallocate(natural, global)

    self%cell_natural_global = minc_ao

  contains

!........................................................................

    subroutine get_original_cell_natural_global_indices()
      !! Gets natural and global indices for original cells on each
      !! process.

      ! Locals:
      PetscInt :: ic, ghost, carray(1)
      DMLabel :: ghost_label
      PetscErrorCode :: ierr

      call DMGetLabel(self%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
      ic = 0
      do c = start_cell, end_interior_cell - 1
         call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
         if (ghost < 0) then
            ic = ic + 1
            natural(ic) = local_to_natural_cell_index(self%cell_natural_global, &
                 l2g, c)
            carray = c
            call ISLocalToGlobalMappingApplyBlock(minc_l2g, 1, carray, &
                 global(ic:ic), ierr); CHKERRQ(ierr)
         end if
      end do

    end subroutine get_original_cell_natural_global_indices

!........................................................................

    PetscInt function get_new_natural_index() result(inatural)
      !! Initialises new natural index for MINC cells.

      ! Locals:
      PetscInt :: local_max_natural, max_natural

      local_max_natural = maxval(natural(1: num_non_ghost_cells))
      call MPI_reduce(local_max_natural, max_natural, 1, MPI_INTEGER, &
           MPI_MAX, 0, PETSC_COMM_WORLD, ierr)
      if (rank == 0) then
         inatural = max_natural + 1
      else
         inatural = 0
      end if

    end function get_new_natural_index

!........................................................................

    subroutine assign_minc_natural_indices(n_all, minc_frac_natural_all, &
         minc_global_all, natural, global, inatural, offset)
      !! Assigns natural indices to MINC level cells on root
      !! process. The 'natural' ordering assigned to MINC cells must
      !! be independent of the mesh partitioning. Here the MINC cells
      !! are ordered first by MINC level and then by the natural order
      !! of the corresponding fracture cells.

      PetscInt, intent(in) :: n_all
      PetscInt, intent(in) :: minc_frac_natural_all(:), minc_global_all(:)
      PetscInt, intent(in out), target :: natural(:), global(:)
      PetscInt, intent(in out) :: inatural, offset
      ! Locals:
      PetscInt, allocatable :: isort(:)
      PetscInt, pointer :: minc_natural_sorted(:), minc_global_sorted(:)
      PetscInt :: i
      PetscErrorCode :: ierr

      if (rank == 0) then

         minc_natural_sorted => natural(offset: offset + n_all - 1)
         minc_global_sorted => global(offset: offset + n_all - 1)

         isort = [(i - 1, i = 1, n_all)]
         call PetscSortIntWithPermutation(n_all, &
              minc_frac_natural_all, isort, ierr); CHKERRQ(ierr)
         isort = isort + 1 ! convert to 1-based

         do i = 1, n_all
            minc_natural_sorted(i) = inatural
            minc_global_sorted(i) = minc_global_all(isort(i))
            inatural = inatural + 1
         end do

         offset = offset + n_all

      end if

    end subroutine assign_minc_natural_indices

!........................................................................

  end subroutine mesh_setup_minc_dm_cell_natural_global

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_geometry(self, minc_dm, max_num_levels, &
       minc_level_cells)
    !! Sets up cell and face geometry vectors for MINC mesh. These
    !! overwrite the original geometry vectors.

    use kinds_module
    use dm_utils_module, only: section_offset, local_vec_section, &
         dm_set_data_layout

    class(mesh_type), intent(in out) :: self
    DM, intent(in out) :: minc_dm
    PetscInt, intent(in) :: max_num_levels
    PetscInt, intent(in) :: minc_level_cells(max_num_levels, &
         self%strata(0)%start: self%strata(0)%end - 1)
    ! Locals:
    PetscInt :: cell_variable_dim(num_cell_variables), &
         face_variable_dim(num_face_variables)
    PetscInt :: iminc, c, f, i, m, minc_p, face_p, ghost
    PetscInt :: cell_offset, minc_cell_offset, face_offset, minc_face_offset
    PetscInt :: ic, num_minc_zone_cells
    DM :: dm_cell, dm_face
    DMLabel :: ghost_label
    IS :: minc_IS
    PetscInt, pointer :: minc_cells(:)
    Vec :: minc_cell_geom, minc_face_geom
    PetscSection :: cell_section, minc_cell_section
    PetscSection :: face_section, minc_face_section
    PetscReal, pointer, contiguous :: cell_geom_array(:), minc_cell_geom_array(:)
    PetscReal, pointer, contiguous :: face_geom_array(:), minc_face_geom_array(:)
    type(cell_type) :: cell
    type(face_type) :: face
    PetscReal :: orig_volume, orig_centroid(3)
    character(80) :: name
    PetscErrorCode :: ierr
    PetscInt, parameter :: nc = 1, np = 1 ! dummy values for cell & face init

    call DMGetLabel(self%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    call DMClone(minc_dm, dm_cell, ierr); CHKERRQ(ierr)
    cell_variable_dim = self%dim
    call dm_set_data_layout(dm_cell, cell_variable_num_components, &
         cell_variable_dim, cell_variable_names)
    call DMCreateLocalVector(dm_cell, minc_cell_geom, ierr); CHKERRQ(ierr)
    call PetscObjectGetName(self%cell_geom, name, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(minc_cell_geom, name, ierr); CHKERRQ(ierr)
    call local_vec_section(minc_cell_geom, minc_cell_section)
    call VecGetArrayF90(minc_cell_geom, minc_cell_geom_array, ierr)
    CHKERRQ(ierr)

    call local_vec_section(self%cell_geom, cell_section)
    call VecGetArrayReadF90(self%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)

    call DMClone(minc_dm, dm_face, ierr); CHKERRQ(ierr)
    face_variable_dim = self%dim - 1
    call dm_set_data_layout(dm_face, face_variable_num_components, &
         face_variable_dim, face_variable_names)
    call DMCreateLocalVector(dm_face, minc_face_geom, ierr); CHKERRQ(ierr)
    call PetscObjectGetName(self%face_geom, name, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(minc_face_geom, name, ierr); CHKERRQ(ierr)
    call local_vec_section(minc_face_geom, minc_face_section)
    call VecGetArrayF90(minc_face_geom, minc_face_geom_array, ierr)
    CHKERRQ(ierr)

    call local_vec_section(self%face_geom, face_section)
    call VecGetArrayReadF90(self%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)

    ! Copy original cell geometry:
    call cell%init(nc, np)
    do c = self%strata(0)%start, self%strata(0)%end - 1
       cell_offset = section_offset(cell_section, c)
       minc_p = self%strata(0)%minc_point(c, 0)
       minc_cell_offset = section_offset(minc_cell_section, minc_p)
       minc_cell_geom_array(minc_cell_offset: &
            minc_cell_offset + cell%dof - 1) = cell_geom_array(cell_offset: &
            cell_offset + cell%dof - 1)
    end do

    ! Copy original face geometry:
    call face%init(nc, np)
    do f = self%strata(1)%start, self%strata(1)%end - 1
       call DMLabelGetValue(ghost_label, f, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          face_offset = section_offset(face_section, f)
          minc_p = self%strata(1)%minc_point(f, 0)
          minc_face_offset = section_offset(minc_face_section, minc_p)
          minc_face_geom_array(minc_face_offset: &
               minc_face_offset + face%dof - 1) = face_geom_array(face_offset: &
               face_offset + face%dof - 1)
       end if
    end do

    do iminc = 1, size(self%minc)
       associate(minc => self%minc(iminc))
         call DMGetStratumSize(self%dm, minc_zone_label_name, iminc, &
              num_minc_zone_cells, ierr); CHKERRQ(ierr)
         if (num_minc_zone_cells > 0) then
            call DMGetStratumIS(self%dm, minc_zone_label_name, &
                 iminc, minc_IS, ierr); CHKERRQ(ierr)
            call ISGetIndicesF90(minc_IS, minc_cells, ierr); CHKERRQ(ierr)
            do i = 1, num_minc_zone_cells
               c = minc_cells(i)
               call DMLabelGetValue(ghost_label, c, ghost, ierr)
               if (ghost < 0) then
                  minc_p = self%strata(0)%minc_point(c, 0)
                  minc_cell_offset = section_offset(minc_cell_section, minc_p)
                  call cell%assign_geometry(minc_cell_geom_array, minc_cell_offset)
                  orig_volume = cell%volume
                  orig_centroid = cell%centroid
                  ! Modify fracture cell volume:
                  cell%volume = orig_volume * minc%volume(1)
                  do m = 1, minc%num_levels
                     ! Assign MINC cell geometry:
                     ic = minc_level_cells(m, c)
                     minc_p = self%strata(0)%minc_point(ic, m)
                     minc_cell_offset = section_offset(minc_cell_section, minc_p)
                     call cell%assign_geometry(minc_cell_geom_array, minc_cell_offset)
                     cell%volume = orig_volume * minc%volume(m + 1)
                     cell%centroid = orig_centroid
                     ! Assign MINC face geometry:
                     face_p = self%strata(1)%minc_point(ic, m)
                     minc_face_offset = section_offset(minc_face_section, face_p)
                     call face%assign_geometry(minc_face_geom_array, minc_face_offset)
                     face%area = orig_volume * minc%connection_area(m)
                     face%distance = minc%connection_distance(m: m + 1)
                     face%distance12 = sum(face%distance)
                     face%normal = 0._dp
                     face%gravity_normal = 0._dp
                     face%centroid = orig_centroid
                     face%permeability_direction = dble(1)
                  end do
               end if
            end do
         end if
       end associate
    end do

    call cell%destroy()
    call VecRestoreArrayReadF90(self%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(minc_cell_geom, minc_cell_geom_array, ierr)
    CHKERRQ(ierr)

    call face%destroy()
    call VecRestoreArrayReadF90(self%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(minc_face_geom, minc_face_geom_array, ierr)
    CHKERRQ(ierr)

    call VecDestroy(self%cell_geom, ierr); CHKERRQ(ierr)
    self%cell_geom = minc_cell_geom
    call VecDestroy(self%face_geom, ierr); CHKERRQ(ierr)
    self%face_geom = minc_face_geom

  end subroutine mesh_setup_minc_geometry

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_rock_properties(self, json, rock_vector, &
       rock_range_start, logfile, err)
    !! Sets rock properties in MINC cells. Rock types can be specified
    !! for fracture and matrix cells in each MINC zone. If fracture
    !! cell rock properties are not specified in the fracture rock
    !! type, their values are left unchanged from the original cell
    !! rock properties. If matrix cell rock properties are not
    !! specified in the matrix rock type, their values are also set
    !! equal to the properties from the original rock type, except the
    !! porosity, which is calculated to preserve the total void
    !! fraction of the rock in the original cell.

    use kinds_module
    use fson_mpi_module
    use logfile_module
    use dictionary_module
    use rock_module
    use dm_utils_module, only: global_vec_section, global_section_offset
    use fson_value_m, only : TYPE_ARRAY, TYPE_OBJECT

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    Vec, intent(in out) :: rock_vector
    PetscInt, intent(in) :: rock_range_start
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: iminc, h, i, c, m, cell_p
    PetscInt :: rock_json_type, irock, num_rocks
    PetscInt :: num_minc_zone_cells, orig_offset, offset
    IS :: minc_IS
    PetscSection :: section
    PetscReal, pointer, contiguous :: rock_array(:)
    PetscInt, pointer :: minc_cells(:)
    type(rock_type) :: orig_rock, rock
    type(fson_value), pointer :: minc_json, minci_json
    type(fson_value), pointer :: rock_json, rocki_json
    PetscInt :: minc_type, max_num_levels, dof, minc_rocktype_zone_index
    PetscInt, allocatable :: ic(:)
    type(rock_type) :: fracture_rock, matrix_rock
    PetscReal :: fracture_porosity, matrix_porosity
    PetscReal, pointer, contiguous :: fracture_rock_array(:), &
         matrix_rock_array(:)
    PetscErrorCode :: ierr
    character(len=64) :: minc_str, minc_rock_str
    character(len=12) :: imstr, irockstr

    err = 0

    call fson_get_mpi(json, "mesh.minc", minc_json)
    minc_type = fson_type_mpi(minc_json, ".")
    select case (minc_type)
    case (TYPE_OBJECT)
       minci_json => minc_json
    case (TYPE_ARRAY)
       minci_json => fson_value_children_mpi(minc_json)
    end select

    call VecGetArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
    call global_vec_section(rock_vector, section)

    call orig_rock%init()
    call rock%init()
    dof = rock%dof

    allocate(fracture_rock_array(dof), matrix_rock_array(dof))
    call fracture_rock%assign(fracture_rock_array, 1)
    call matrix_rock%assign(matrix_rock_array, 1)

    max_num_levels = maxval(self%minc%num_levels)
    allocate(ic(max_num_levels))
    ic = 0
    h = 0
    minc_rocktype_zone_index = 1
    
    do iminc = 1, size(self%minc)
       associate(minc => self%minc(iminc))

         write(imstr, '(i0)') iminc - 1
         minc_str = 'mesh.minc[' // trim(imstr) // ']'

         if (fson_has_mpi(minci_json, "rock")) then

            call fson_get_mpi(minci_json, "rock", rock_json)
            rock_json_type = fson_type_mpi(rock_json, ".")
            select case (rock_json_type)
            case (TYPE_OBJECT)
               num_rocks = 1
               rocki_json => rock_json
            case (TYPE_ARRAY)
               num_rocks = fson_value_count_mpi(rock_json, ".")
               rocki_json => fson_value_children_mpi(rock_json)
            end select

            do irock = 1, num_rocks

               write(irockstr, '(i0)') irock - 1
               minc_rock_str = 'rock[' // trim(irockstr) // ']'

               call read_rock_type("fracture", rocki_json, json, fracture_rock, err)
               if (err == 0) then
                  call read_rock_type("matrix", rocki_json, json, matrix_rock, err)
                  if (err == 0) then

                     call DMGetStratumSize(self%dm, minc_rocktype_zone_label_name, &
                          minc_rocktype_zone_index, num_minc_zone_cells, ierr)
                     CHKERRQ(ierr)
                     if (num_minc_zone_cells > 0) then
                        call DMGetStratumIS(self%dm, minc_rocktype_zone_label_name, &
                             minc_rocktype_zone_index, minc_IS, ierr); CHKERRQ(ierr)
                        call ISGetIndicesF90(minc_IS, minc_cells, ierr); CHKERRQ(ierr)
                        do i = 1, num_minc_zone_cells
                           c = minc_cells(i)
                           if (self%ghost_cell(c) < 0) then
                              orig_offset = global_section_offset(section, c, rock_range_start)
                              call orig_rock%assign(rock_array, orig_offset)

                              if (fracture_rock%porosity < 0._dp) then
                                 fracture_porosity = orig_rock%porosity
                              else
                                 fracture_porosity = fracture_rock%porosity
                              end if

                              if (matrix_rock%porosity < 0._dp) then
                                 ! Default matrix porosity preserves void fraction of original rock:
                                 associate(fracture_volume => minc%volume(1))
                                   matrix_porosity = (orig_rock%porosity - &
                                        fracture_porosity * fracture_volume) / &
                                        (1._dp - fracture_volume)
                                 end associate
                              else
                                 matrix_porosity = matrix_rock%porosity
                              end if

                              associate(orig_properties => rock_array(orig_offset: &
                                   orig_offset + dof - 1))

                                do m = 1, minc%num_levels
                                   cell_p = self%strata(h)%minc_point(ic(m), m)
                                   offset = global_section_offset(section, cell_p, rock_range_start)
                                   ! Update specified matrix properties:
                                   associate(matrix_properties => rock_array(offset: &
                                        offset + dof - 1))
                                     matrix_properties = merge( &
                                          matrix_rock_array, orig_properties, &
                                          matrix_rock_array > 0._dp)
                                   end associate
                                   call rock%assign(rock_array, offset)
                                   rock%porosity = matrix_porosity
                                   ic(m) = ic(m) + 1
                                end do

                                ! Update specified fracture properties:
                                orig_properties = merge( &
                                     fracture_rock_array, orig_properties, &
                                     fracture_rock_array > 0._dp)
                              end associate
                           end if
                        end do
                     end if

                  else
                     exit
                  end if
               else
                  exit
               end if

               if (rock_json_type == TYPE_ARRAY) then
                  rocki_json => fson_value_next_mpi(rocki_json)
               end if
               minc_rocktype_zone_index = minc_rocktype_zone_index + 1

            end do

         end if

         if (minc_type == TYPE_ARRAY) then
            minci_json => fson_value_next_mpi(minci_json)
         end if

       end associate
    end do

    deallocate(ic)
    call fracture_rock%destroy()
    call matrix_rock%destroy()
    deallocate(fracture_rock_array, matrix_rock_array)
    call rock%destroy()
    call orig_rock%destroy()
    call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)

  contains

 !........................................................................

    subroutine read_rock_type(name, rocki_json, json, rock, err)

      character(*), intent(in) :: name
      type(fson_value), pointer, intent(in) :: rocki_json, json
      type(rock_type), intent(out) :: rock
      PetscErrorCode, intent(out) :: err
      ! Locals:
      character(max_rockname_length) :: rock_name
      character(len=64) :: rockstr
      type(list_node_type), pointer :: node

      err = 0

      if (fson_has_mpi(rocki_json, trim(name) // ".type")) then
         call fson_get_mpi(rocki_json, trim(name) // ".type", val = rock_name)
         node => self%rock_types%get(rock_name)
         if (associated(node)) then
            select type (item => node%data)
            type is (rock_dict_item_type)
               rockstr = trim(minc_str) // "." // trim(minc_rock_str) // "." &
                    // trim(name) // ".type."
               call read_rock_parameters(item%rock, rock, rockstr)
            end select
         else
            if (present(logfile)) then
               call logfile%write(LOG_LEVEL_ERR, "input", "unrecognised", &
                    str_key = trim(minc_str) // "." // trim(minc_rock_str) // "." &
                    // trim(name) // ".type", str_value = rock_name)
               err = 1
            end if
         end if
      else
         if (present(logfile)) then
            call logfile%write(LOG_LEVEL_ERR, "input", "not found", &
                 str_key = trim(minc_str), str_value = trim(name) // ".type")
            err = 1
         end if
      end if

    end subroutine read_rock_type

!........................................................................

    subroutine read_rock_parameters(json, rock, rockstr)

      type(fson_value), pointer, intent(in) :: json
      type(rock_type), intent(out) :: rock
      character(*), intent(in) :: rockstr
      ! Locals:
      PetscReal, parameter :: default = -1._dp ! Flag missing values with -1
      PetscReal, parameter :: default_permeability(3) = [-1._dp, -1._dp, -1._dp]
      PetscReal, allocatable :: permeability(:)

      call fson_get_mpi(json, "permeability", default_permeability, &
           permeability, logfile, trim(rockstr) // "permeability")
      rock%permeability = 0._dp
      rock%permeability(1: size(permeability)) = permeability
      call fson_get_mpi(json, "wet_conductivity", default, &
           rock%wet_conductivity, logfile, trim(rockstr) // "wet_conductivity")
      call fson_get_mpi(json, "dry_conductivity", default, &
           rock%dry_conductivity, logfile, trim(rockstr) // "dry_conductivity")
      call fson_get_mpi(json, "porosity", default, rock%porosity, logfile, &
           trim(rockstr) // "porosity")
      call fson_get_mpi(json, "density", default, rock%density, logfile, &
           trim(rockstr) // "density")
      call fson_get_mpi(json, "specific_heat", default, &
           rock%specific_heat, logfile, trim(rockstr) // "specific_heat")

    end subroutine read_rock_parameters

  end subroutine mesh_setup_minc_rock_properties

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_point_sf(self, minc_dm)
    !! Sets up point SF for MINC DM, for communicating ghost values between
    !! processors.

    use dm_utils_module, only: dm_point_stratum_height

    class(mesh_type), intent(in out) :: self
    DM, intent(in out) :: minc_dm
    ! Locals:
    PetscSF :: sf, minc_sf
    PetscInt :: start_chart, end_chart, start_minc_chart, end_minc_chart
    PetscInt :: h, i, p
    PetscInt :: num_roots, num_leaves, num_minc_roots
    PetscInt, allocatable :: new_locations(:), new_remote_locations(:)
    PetscInt, pointer :: leaves(:)
    type(PetscSFNode), pointer :: roots(:)
    PetscInt, allocatable :: minc_leaves(:)
    type(PetscSFNode), allocatable :: minc_roots(:)
    PetscErrorCode :: ierr
    PetscInt, parameter :: m = 0 ! MINC level

    call DMGetPointSF(minc_dm, minc_sf, ierr); CHKERRQ(ierr)
    call DMGetPointSF(self%dm, sf, ierr); CHKERRQ(ierr)
    call PetscSFGetGraph(sf, num_roots, num_leaves, leaves, roots, ierr)
    CHKERRQ(ierr)

    if (num_roots >= 0) then

       call DMPlexGetChart(self%dm, start_chart, end_chart, ierr)
       CHKERRQ(ierr)
       ! All DM points are considered potential roots for leaves on
       ! another process:
       call DMPlexGetChart(minc_dm, start_minc_chart, end_minc_chart, ierr)
       CHKERRQ(ierr)
       num_minc_roots = end_minc_chart - start_minc_chart

       ! Update root locations and broadcast:
       allocate(new_locations(0: num_roots - 1), &
            new_remote_locations(0: end_chart - start_chart))
       do p = 0, num_roots - 1
          h = dm_point_stratum_height(self%strata, p)
          new_locations(p) = self%strata(h)%minc_point(p, m)
       end do
       call PetscSFBcastBegin(sf, MPI_INTEGER, new_locations, &
            new_remote_locations, ierr); CHKERRQ(ierr)
       call PetscSFBcastEnd(sf, MPI_INTEGER, new_locations, &
            new_remote_locations, ierr); CHKERRQ(ierr)

       ! Determine MINC roots and leaves:
       allocate(minc_leaves(num_leaves), minc_roots(num_leaves))
       do i = 1, num_leaves
          p = leaves(i)
          h = dm_point_stratum_height(self%strata, p)
          minc_leaves(i) = self%strata(h)%minc_point(p, m)
          minc_roots(i)%rank = roots(i)%rank
          minc_roots(i)%index = new_remote_locations(p)
       end do
       deallocate(new_locations, new_remote_locations)

       call PetscSFSetGraph(minc_sf, num_minc_roots, num_leaves, &
            minc_leaves, PETSC_COPY_VALUES, minc_roots, PETSC_COPY_VALUES, &
            ierr); CHKERRQ(ierr)
       deallocate(minc_leaves, minc_roots)

    end if

  end subroutine mesh_setup_minc_point_sf

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_coordinates(self, minc_dm)
    !! Sets coordinate section and vector in MINC DM. Coordinates are
    !! copied from the vertices of the original single-porosity DM. No
    !! coordinates are assigned to the (dummy) vertices added when
    !! the MINC DM is constructed.

    use dm_utils_module, only: section_offset

    class(mesh_type), intent(in out) :: self
    DM, intent(in out) :: minc_dm
    ! Locals:
    PetscInt :: dim, v, coord_size, dof, offset, minc_offset
    PetscInt :: minc_v, depth
    PetscInt :: vstart, vend, minc_vstart, minc_vend
    PetscSection :: section, minc_section
    Vec :: coordinates, minc_coordinates
    PetscReal, pointer :: coord_array(:), minc_coord_array(:)
    PetscReal, pointer :: pos(:), minc_pos(:)
    PetscErrorCode :: ierr

    PetscMPIInt :: rank
    call mpi_comm_rank(PETSC_COMM_WORLD, rank, ierr)

    call DMGetCoordinateDim(self%dm, dim, ierr); CHKERRQ(ierr)
    call DMGetCoordinatesLocal(self%dm, coordinates, ierr); CHKERRQ(ierr)
    call DMGetCoordinateSection(self%dm, section, ierr); CHKERRQ(ierr)
    call DMPlexGetDepthStratum(self%dm, 0, vstart, vend, ierr); CHKERRQ(ierr)

    call DMSetCoordinateDim(minc_dm, dim, ierr); CHKERRQ(ierr)
    call PetscSectionCreate(PETSC_COMM_WORLD, minc_section, ierr); CHKERRQ(ierr)
    call PetscSectionSetNumFields(minc_section, 1, ierr); CHKERRQ(ierr)
    call PetscSectionSetFieldComponents(minc_section, 0, dim, ierr)
    CHKERRQ(ierr)
    call DMPlexGetDepthStratum(minc_dm, 0, minc_vstart, minc_vend, ierr)
    CHKERRQ(ierr)
    call PetscSectionSetChart(minc_section, minc_vstart, minc_vend, ierr)
    CHKERRQ(ierr)
    do v = minc_vstart, minc_vend - 1
       call PetscSectionSetDof(minc_section, v, dim, ierr); CHKERRQ(ierr)
       call PetscSectionSetFieldDof(minc_section, v, 0, dim, ierr)
       CHKERRQ(ierr)
    end do
    call PetscSectionSetUp(minc_section, ierr); CHKERRQ(ierr)
    call DMSetCoordinateSection(minc_dm, PETSC_DETERMINE, minc_section, ierr)
    CHKERRQ(ierr)

    call VecCreate(PETSC_COMM_SELF, minc_coordinates, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(minc_coordinates, "coordinates", ierr); CHKERRQ(ierr)
    call PetscSectionGetStorageSize(minc_section, coord_size, ierr); CHKERRQ(ierr)
    call VecSetSizes(minc_coordinates, coord_size, PETSC_DETERMINE, ierr)
    CHKERRQ(ierr)
    call VecSetBlockSize(minc_coordinates, dim, ierr); CHKERRQ(ierr)
    call VecSetType(minc_coordinates, VECSTANDARD, ierr); CHKERRQ(ierr)
    call DMSetCoordinatesLocal(minc_dm, minc_coordinates, ierr); CHKERRQ(ierr)

    ! Copy vertex coordinates from original DM:
    call VecGetArrayReadF90(coordinates, coord_array, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(minc_coordinates, minc_coord_array, ierr); CHKERRQ(ierr)
    minc_coord_array = 0._dp
    call DMPlexGetDepth(self%dm, depth, ierr); CHKERRQ(ierr)
    do v = vstart, vend - 1
       minc_v = self%strata(depth)%minc_point(v, 0)
       call PetscSectionGetDof(section, v, dof, ierr); CHKERRQ(ierr)
       offset = section_offset(section, v)
       minc_offset = section_offset(minc_section, minc_v)
       pos => coord_array(offset: offset + dof - 1)
       minc_pos => minc_coord_array(minc_offset: minc_offset + dof - 1)
       minc_pos = pos
    end do
    call VecRestoreArrayF90(minc_coordinates, minc_coord_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(coordinates, coord_array, ierr); CHKERRQ(ierr)

    call VecDestroy(minc_coordinates, ierr); CHKERRQ(ierr)
    call PetscSectionDestroy(minc_section, ierr); CHKERRQ(ierr)

  end subroutine mesh_setup_minc_coordinates

!------------------------------------------------------------------------

  subroutine mesh_redistribute(self, sf)
    !! Redistributes DM and other mesh data for load balancing. Should
    !! only be called when running in parallel. The redistribution
    !! star forest is returned.

    use dm_utils_module, only: dm_set_default_data_layout, &
         dm_natural_order_IS, dm_get_natural_to_global_ao, &
         dm_get_cell_index, dm_distribute_index_set

    class(mesh_type), intent(in out) :: self
    PetscSF, intent(out) :: sf !! Redistribution star forest
    ! Locals:
    DM :: dm_is
    PetscSection :: section
    PetscErrorCode :: ierr

    call DMClone(self%dm, dm_is, ierr); CHKERRQ(ierr)
    call dm_set_default_data_layout(dm_is, 1)
    call DMGetSection(dm_is, section, ierr); CHKERRQ(ierr)

    call self%redistribute_dm(sf)
    if (sf .ne. PETSC_NULL_SF) then
       call self%redistribute_geometry(sf)
       call dm_distribute_index_set(self%dm, sf, section, self%cell_natural)
       if (self%has_minc) then
          call dm_distribute_index_set(self%dm, sf, section, self%cell_parent_natural)
       end if
       self%cell_natural_global = &
            dm_get_natural_to_global_ao(self%dm, self%cell_natural)
       call ISDestroy(self%cell_index, ierr); CHKERRQ(ierr)
       call dm_get_cell_index(self%dm, self%cell_natural_global, &
                  self%cell_index)
       call self%setup_ghost_arrays()
    end if

    call DMDestroy(dm_is, ierr); CHKERRQ(ierr)

  end subroutine mesh_redistribute

!------------------------------------------------------------------------

  subroutine mesh_redistribute_dm(self, sf)
    !! Redistributes DM (e.g. to improve load balancing for MINC
    !! simulations), and returns SF for the redistribution.

    use dm_utils_module, only: dm_set_default_data_layout, &
         dm_set_fv_adjacency, dm_label_partition_ghosts

    class(mesh_type), intent(in out) :: self
    PetscSF, intent(out) :: sf
    ! Locals:
    DM :: dm_dist
    PetscErrorCode :: ierr

    call DMPlexDistribute(self%dm, partition_overlap, sf, &
         dm_dist, ierr); CHKERRQ(ierr)
    if (dm_dist .ne. PETSC_NULL_DM) then
       call DMDestroy(self%dm, ierr); CHKERRQ(ierr)
       self%dm = dm_dist
       call dm_set_fv_adjacency(self%dm)
       call dm_set_default_data_layout(self%dm, self%dof)
       call dm_label_partition_ghosts(self%dm)
    end if

  end subroutine mesh_redistribute_dm

!------------------------------------------------------------------------

  subroutine mesh_redistribute_geometry(self, sf)
    !! Redistributes cell and face geometry vectors after
    !! redistributing DM according to the specified SF.

    use dm_utils_module, only: dm_distribute_local_vec

    class(mesh_type), intent(in out) :: self
    PetscSF, intent(in out) :: sf !! Redistribution star forest

    if (sf .ne. PETSC_NULL_SF) then
       call dm_distribute_local_vec(self%dm, sf, self%cell_geom)
       call dm_distribute_local_vec(self%dm, sf, self%face_geom)
       call self%check_face_orientations()
    end if

  end subroutine mesh_redistribute_geometry

!------------------------------------------------------------------------

  subroutine mesh_check_face_orientations(self)
    !! Checks interior faces to make sure they match the orientation
    !! of the cells in their support (sometimes the orientations can
    !! be reversed during redistribution). If the order of cells in
    !! the face has been reversed, it is necessary also to reverse the
    !! normal vector and centroid distance array in the face
    !! geometry. For MINC faces a reversed orientation is detected by
    !! MINC levels that decrease from the first cell to the second.

    use dm_utils_module, only: local_vec_section, section_offset, &
         dm_get_end_interior_cell
    use face_module, only: face_type
    use minc_module, only: minc_level_label_name

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscSection :: cell_geom_section, face_geom_section
    PetscReal, pointer, contiguous :: cell_geom_array(:), face_geom_array(:)
    type(face_type) :: face
    PetscInt :: start_cell, end_cell, end_interior_cell
    PetscInt :: start_face, end_face, f, offset, i, cell_offset(2)
    PetscInt :: ghost, minc_levels(2), num_cells
    DMLabel :: ghost_label, minc_level_label
    PetscInt, pointer :: cells(:)
    PetscBool :: reversed
    PetscErrorCode :: ierr

    call local_vec_section(self%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(self%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)
    call local_vec_section(self%face_geom, face_geom_section)
    call VecGetArrayF90(self%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)

    call face%init(1, 1)

    call DMGetLabel(self%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    if (self%has_minc) then
       call DMGetLabel(self%dm, minc_level_label_name, minc_level_label, ierr)
       CHKERRQ(ierr)
    end if
    call DMPlexGetHeightStratum(self%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    end_interior_cell = dm_get_end_interior_cell(self%dm, end_cell)
    call DMPlexGetHeightStratum(self%dm, 1, start_face, end_face, ierr)
    CHKERRQ(ierr)
    minc_levels = 0

    do f = start_face, end_face - 1
       call DMLabelGetValue(ghost_label, f, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          call DMPlexGetSupportSize(self%dm, f, num_cells, ierr); CHKERRQ(ierr)
          call DMPlexGetSupport(self%dm, f, cells, ierr); CHKERRQ(ierr)
          do i = 1, 2
             if (self%has_minc) then
                call DMLabelGetValue(minc_level_label, cells(i), &
                     minc_levels(i), ierr); CHKERRQ(ierr)
             end if
             cell_offset(i) = section_offset(cell_geom_section, cells(i))
          end do
          offset = section_offset(face_geom_section, f)
          call face%assign_geometry(face_geom_array, offset)
          call face%assign_cell_geometry(cell_geom_array, cell_offset)
          if (all(cells < end_interior_cell)) then
             if (all(minc_levels <= 0)) then
                reversed = face%reversed_orientation()
             else
                reversed = (minc_levels(1) > minc_levels(2))
             end if
             if (reversed) call face%reverse_geometry()
          end if
       end if
    end do

    call face%destroy()
    call VecRestoreArrayReadF90(self%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(self%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)

  end subroutine mesh_check_face_orientations

!------------------------------------------------------------------------

  PetscInt function mesh_local_cell_minc_level(self, local) result(minc_level)
    !! Takes local cell index and returns MINC level. If the mesh is
    !! not MINC, the returned value is zero.

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: local !! Local cell index
    ! Locals:
    DMLabel :: minc_level_label
    PetscErrorCode :: ierr

    if (self%has_minc) then
       call DMGetLabel(self%dm, minc_level_label_name, &
            minc_level_label, ierr); CHKERRQ(ierr)
       call DMLabelGetValue(minc_level_label, local, minc_level, ierr)
    else
       minc_level = 0
    end if

  end function mesh_local_cell_minc_level

!------------------------------------------------------------------------

  PetscInt function mesh_local_to_parent_natural(self, local) &
       result(natural)
    !! Takes a local cell index and returns natural index of the
    !! corresponding parent cell. (For a non-MINC mesh, the
    !! 'parent' cell is just the original cell itself.)

    class(mesh_type), intent(in) :: self
    PetscInt, intent(in) :: local !! Local cell index
    ! Locals:
    PetscInt, pointer :: natural_array(:)
    PetscErrorCode :: ierr

    if (self%has_minc) then
       call ISGetIndicesF90(self%cell_parent_natural, natural_array, &
            ierr); CHKERRQ(ierr)
       natural = natural_array(local + 1)
       call ISRestoreIndicesF90(self%cell_parent_natural, natural_array, &
            ierr); CHKERRQ(ierr)
    else
       call ISGetIndicesF90(self%cell_natural, natural_array, &
            ierr); CHKERRQ(ierr)
       natural = natural_array(local + 1)
       call ISRestoreIndicesF90(self%cell_natural, natural_array, &
            ierr); CHKERRQ(ierr)
    end if

  end function mesh_local_to_parent_natural

!------------------------------------------------------------------------

  subroutine mesh_global_to_parent_natural(self, global, &
       parent_natural, minc_level)
    !! Takes a global cell index and returns natural index of
    !! corresponding parent cell, together with the MINC level of
    !! the cell. (For a non-MINC mesh, the 'parent' cell is just the
    !! original cell itself, and the MINC level is zero.)

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: global !! Global cell index
    PetscInt, intent(out) :: parent_natural !! Natural index of parent cell
    PetscInt, intent(out) :: minc_level !! MINC level of cell
    ! Locals:
    PetscInt :: idx(1), local_array(1), n
    ISLocalToGlobalMapping :: l2g
    PetscMPIInt :: rank, found_rank, process_found_rank
    PetscErrorCode :: ierr

    idx(1) = global
    if (self%has_minc) then
       call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
       process_found_rank = 0
       call DMGetLocalToGlobalMapping(self%dm, l2g, ierr); CHKERRQ(ierr)
       call ISGlobalToLocalMappingApplyBlock(l2g, IS_GTOLM_MASK, 1, &
            idx, n, local_array, ierr); CHKERRQ(ierr)
       associate(c => local_array(1))
         if (c >= 0) then
            if (self%ghost_cell(c) < 0) then
               parent_natural = self%local_to_parent_natural(c)
               minc_level = self%local_cell_minc_level(c)
               CHKERRQ(ierr)
               process_found_rank = rank
            end if
         end if
         call MPI_Allreduce(process_found_rank, found_rank, 1, MPI_INT, &
              MPI_MAX, PETSC_COMM_WORLD, ierr)
         call MPI_bcast(parent_natural, 1, MPI_LOGICAL, found_rank, &
            PETSC_COMM_WORLD, ierr)
         call MPI_bcast(minc_level, 1, MPI_LOGICAL, found_rank, &
            PETSC_COMM_WORLD, ierr)
       end associate
    else
       call AOPetscToApplication(self%cell_natural_global, 1, idx, ierr)
       CHKERRQ(ierr)
       parent_natural = idx(1)
       minc_level = 0
    end if

  end subroutine mesh_global_to_parent_natural

!------------------------------------------------------------------------

  subroutine mesh_natural_cell_output_arrays(self, natural, minc_level, &
       keys, values)
    !! Takes a natural cell index and MINC level, and returns arrays
    !! of keys and values for output. If the mesh is MINC, these
    !! arrays have two entries, one for natural cell index
    !! and one for MINC level; otherwise, the arrays have only one
    !! entry each, for natural cell index. These arrays are intended
    !! for output to logfile.

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: natural !! Natural cell index
    PetscInt, intent(in) :: minc_level !! MINC level
    character(len = *), allocatable :: keys(:) !! Key array
    PetscInt, allocatable :: values(:)

    if (self%has_minc) then
       keys = ['cell', 'minc']
       values = [natural, minc_level]
    else
       keys = ['cell']
       values = [natural]
    end if

  end subroutine mesh_natural_cell_output_arrays

!------------------------------------------------------------------------

end module mesh_module
