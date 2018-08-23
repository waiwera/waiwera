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
  use zone_module
  use list_module
  use mpi_utils_module
  use fson
  use minc_module
  use dm_utils_module, only: dm_stratum_type
  use dictionary_module

  implicit none

  private

  PetscInt, parameter, public :: max_mesh_filename_length = 200
  character(len = 16), public :: open_boundary_label_name = "open_boundary" !! Name of DMLabel for identifying open boundaries
  character(len = 22) :: face_property_override_label_name = "face_property_override"

  type, public :: mesh_type
     !! Mesh type.
     private
     character(max_mesh_filename_length), public :: filename !! Mesh file name
     DM, public :: serial_dm !! Original DM read from file (not distributed)
     DM, public :: original_dm !! Original DM read from file (and distributed)
     DM, public :: dm !! DM representing the mesh topology (may be modified from original_dm)
     Vec, public :: cell_geom !! Vector containing cell geometry data
     Vec, public :: face_geom !! Vector containing face geometry data
     PetscInt :: depth !! DM depth
     type(dm_stratum_type), allocatable, public :: strata(:) !! Mesh strata (used for MINC point calculations)
     PetscReal, allocatable, public :: bcs(:,:) !! Array containing boundary conditions
     IS, public :: cell_index !! Index set defining natural to global cell ordering (without boundary cells)
     AO, public :: cell_order !! Application ordering to convert between global and natural cell indices
     AO, public :: original_cell_order !! Global-to-natural AO for original DM
     PetscSF, public :: dist_sf !! Distribution star forest
     PetscInt, allocatable, public :: minc_cell_map(:) !! Mapping from MINC cell local indices to original single-porosity cell local indices
     PetscInt, public, allocatable :: ghost_cell(:), ghost_face(:) !! Ghost label values for cells and faces
     type(minc_type), allocatable, public :: minc(:) !! Array of MINC zones, with parameters
     PetscReal, public :: permeability_rotation(3, 3) !! Rotation matrix of permeability axes
     PetscReal, public :: thickness !! Mesh thickness (for dimension < 3)
     type(list_type), public :: zones !! Mesh zones
     type(dictionary_type), public :: rock_types !! Dictionary of rock types by name
     PetscBool, public :: radial !! If mesh coordinate system is radial or Cartesian
     PetscBool, public :: has_minc !! If mesh has any MINC cells
   contains
     procedure :: distribute => mesh_distribute
     procedure :: construct_ghost_cells => mesh_construct_ghost_cells
     procedure :: setup_geometry => mesh_setup_geometry
     procedure :: setup_ghost_arrays => mesh_setup_ghost_arrays
     procedure :: destroy_minc => mesh_destroy_minc
     procedure :: destroy_strata => mesh_destroy_strata
     procedure :: setup_coordinate_parameters => mesh_setup_coordinate_parameters
     procedure :: set_boundary_face_distances => mesh_set_boundary_face_distances
     procedure :: set_permeability_rotation => mesh_set_permeability_rotation
     procedure :: modify_geometry => mesh_modify_geometry
     procedure :: read_overridden_face_properties => mesh_read_overridden_face_properties
     procedure :: override_face_properties => mesh_override_face_properties
     procedure :: setup_zones => mesh_setup_zones
     procedure :: setup_minc => mesh_setup_minc
     procedure :: setup_minc_dm => mesh_setup_minc_dm
     procedure :: setup_minc_dm_zone_and_cells => mesh_setup_minc_dm_zone_and_cells
     procedure :: setup_minc_dm_strata_shifts => mesh_setup_minc_dm_strata_shifts
     procedure :: set_minc_dm_cone_sizes => mesh_set_minc_dm_cone_sizes
     procedure :: set_minc_dm_cones => mesh_set_minc_dm_cones
     procedure :: setup_minc_dm_depth_label => mesh_setup_minc_dm_depth_label
     procedure :: transfer_labels_to_minc_dm => mesh_transfer_labels_to_minc_dm
     procedure :: setup_minc_dm_level_label_and_cell_map => mesh_setup_minc_dm_level_label_and_cell_map
     procedure :: setup_minc_dm_cell_order => mesh_setup_minc_dm_cell_order
     procedure :: setup_minc_geometry => mesh_setup_minc_geometry
     procedure :: setup_minc_rock_properties => mesh_setup_minc_rock_properties
     procedure :: setup_minc_point_sf => mesh_setup_minc_point_sf
     procedure, public :: init => mesh_init
     procedure, public :: configure => mesh_configure
     procedure, public :: setup_boundaries => mesh_setup_boundaries
     procedure, public :: set_boundary_values => mesh_set_boundary_values
     procedure, public :: destroy => mesh_destroy
     procedure, public :: local_to_fracture_natural => mesh_local_to_fracture_natural
     procedure, public :: global_to_fracture_natural => mesh_global_to_fracture_natural
     procedure, public :: natural_cell_output_arrays =>  mesh_natural_cell_output_arrays
     procedure, public :: local_cell_minc_level => mesh_local_cell_minc_level
     procedure, public :: destroy_distribution_data => mesh_destroy_distribution_data
  end type mesh_type

contains

!------------------------------------------------------------------------

  subroutine mesh_distribute(self)
    !! Distributes mesh over processors, and returns star forest from
    !! mesh distribution.
    
    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    PetscInt, parameter :: overlap = 1

    call DMPlexDistribute(self%serial_dm, overlap, self%dist_sf, &
         self%original_dm, ierr); CHKERRQ(ierr)
    if (self%original_dm .eq. PETSC_NULL_DM) then
       self%original_dm = self%serial_dm
    end if

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

  subroutine mesh_construct_ghost_cells(self)
    !! Constructs ghost cells on open boundary faces.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    DM :: ghost_dm
    PetscErrorCode :: ierr

    call DMPlexConstructGhostCells(self%original_dm, open_boundary_label_name, &
         PETSC_NULL_INTEGER, ghost_dm, ierr); CHKERRQ(ierr)
    if (ghost_dm .ne. PETSC_NULL_DM) then
       call DMDestroy(self%original_dm, ierr); CHKERRQ(ierr);
       self%original_dm = ghost_dm
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
    call DMGetDimension(self%serial_dm, dim, ierr); CHKERRQ(ierr)
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
    DM :: dm_cell, dm_face
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
    PetscInt :: cell_variable_dim(num_cell_variables)
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
    call VecGetDM(self%cell_geom, dm_cell, ierr); CHKERRQ(ierr)
    cell_variable_dim = dim
    call set_dm_data_layout(dm_cell, cell_variable_num_components, &
         cell_variable_dim, cell_variable_names)

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
    PetscInt :: start_cell, end_cell, start_face, end_face
    DMLabel :: ghost_label
    PetscErrorCode :: ierr

    call DMPlexGetHeightStratum(self%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    call DMPlexGetHeightStratum(self%dm, 1, start_face, end_face, ierr)
    CHKERRQ(ierr)

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

  subroutine mesh_init(self, json, logfile)
    !! Initializes mesh, reading filename from JSON input file.
    !! If the filename is not present, an error is raised.
    !! Otherwise, the PETSc DM is read in.

    use logfile_module
    use fson_mpi_module
    use fson_value_m, only: TYPE_STRING, TYPE_OBJECT

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
       ! Read in DM:
       call DMPlexCreateFromFile(PETSC_COMM_WORLD, self%filename, PETSC_TRUE, &
            self%serial_dm, ierr); CHKERRQ(ierr)
       call self%setup_coordinate_parameters(json, logfile)
       call self%set_permeability_rotation(json, logfile)
       call self%read_overridden_face_properties(json, logfile)
       call self%rock_types%init(owner = PETSC_TRUE)
       self%has_minc = PETSC_FALSE
    end if

  end subroutine mesh_init

!------------------------------------------------------------------------

  subroutine mesh_configure(self, eos, gravity, json, logfile, viewer, err)
    !! Configures mesh, including distribution over processes and
    !! construction of ghost cells, setup of data layout, geometry and
    !! cell index set.

    use eos_module, only: eos_type
    use dm_utils_module
    use logfile_module
    use rock_module, only: setup_rock_types

    class(mesh_type), intent(in out) :: self
    class(eos_type), intent(in) :: eos !! EOS object
    PetscReal, intent(in) :: gravity(:) !! Gravity vector
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    type(logfile_type), intent(in out), optional :: logfile !! Log file
    PetscViewer, intent(in out) :: viewer !! PetscViewer for output of cell index sets to HDF5 file
    PetscErrorCode, intent(out) :: err !! Error flag
    ! Locals:
    PetscMPIInt :: np
    PetscErrorCode :: ierr

    err = 0
    call MPI_comm_size(PETSC_COMM_WORLD, np, ierr)

    call self%distribute()

    associate(dof => eos%num_primary_variables)

      call dm_setup_fv_discretization(self%original_dm, dof)
      call self%setup_boundaries(json, eos, logfile)
      call self%construct_ghost_cells()
      call set_dm_default_data_layout(self%original_dm, dof)
      call dm_set_fv_adjacency(self%original_dm)

      call self%setup_geometry(gravity)
      self%dm = self%original_dm
      self%original_cell_order = dm_get_natural_to_global_ao(self%original_dm, &
           self%dist_sf)
      self%cell_order = self%original_cell_order

      call self%setup_zones(json, logfile, err)
      if (err == 0) then
         call setup_rock_types(json, self%dm, self%cell_order, &
              self%rock_types, logfile, err)
         if (err == 0) then
            call self%setup_minc(json, logfile, err)
            if (err == 0) then
               if (self%has_minc) call self%setup_minc_dm(dof)
               call dm_get_cell_index(self%dm, self%cell_order, &
                    self%cell_index)
               if (viewer /= PETSC_NULL_VIEWER) then
                  call ISView(self%cell_index, viewer, ierr); CHKERRQ(ierr)
               end if
            end if
         end if
      end if

    end associate

    call self%setup_ghost_arrays()
    
  end subroutine mesh_configure

!------------------------------------------------------------------------

  subroutine mesh_destroy(self)
    !! Destroys mesh object.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call VecDestroy(self%face_geom, ierr); CHKERRQ(ierr)
    call DMDestroy(self%dm, ierr); CHKERRQ(ierr)
    call ISDestroy(self%cell_index, ierr); CHKERRQ(ierr)

    if (allocated(self%bcs)) then
       deallocate(self%bcs)
    end if

    call AODestroy(self%cell_order, ierr); CHKERRQ(ierr)

    if (allocated(self%ghost_cell)) then
       deallocate(self%ghost_cell)
    end if
    if (allocated(self%ghost_face)) then
       deallocate(self%ghost_face)
    end if

    if (self%has_minc) then
       call DMDestroy(self%original_dm, ierr); CHKERRQ(ierr)
       call AODestroy(self%original_cell_order, ierr); CHKERRQ(ierr)
       call self%destroy_minc()
       call self%destroy_strata()
       deallocate(self%minc_cell_map)
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
          call DMGetStratumSize(self%original_dm, open_boundary_label_name, ibdy, &
               num_faces, ierr); CHKERRQ(ierr)
          if (num_faces > 0) then
             call DMGetStratumIS(self%original_dm, open_boundary_label_name, ibdy, &
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
    use dm_utils_module, only: dm_get_natural_to_global_ao, &
         dm_cell_normal_face, dm_get_num_partition_ghost_cells, &
         dm_get_end_interior_cell

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input file
    class(eos_type), intent(in) :: eos !! EOS object
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    PetscErrorCode :: ierr
    PetscBool :: mesh_has_label
    PetscInt :: start_cell, end_cell, end_interior_cell
    PetscInt :: num_ghost_cells, num_non_ghost_cells, end_non_ghost_cell
    ISLocalToGlobalMapping :: l2g
    AO :: ao
    type(fson_value), pointer :: boundaries, bdy, faces_json, face_json
    PetscInt :: faces_type, face1_type
    PetscInt :: num_boundaries, num_faces, num_cells, ibdy
    PetscInt :: iface, icell, np, i, offset, nout
    PetscInt, allocatable :: default_faces(:), default_cells(:)
    PetscInt, allocatable :: faces(:), cells(:), local_cells(:)
    PetscInt :: region, normal_len, num_face_items
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
       bdy => fson_value_children_mpi(boundaries)

       call DMPlexGetHeightStratum(self%original_dm, 0, start_cell, end_cell, ierr)
       CHKERRQ(ierr)
       end_interior_cell = dm_get_end_interior_cell(self%original_dm, end_cell)
       num_ghost_cells = dm_get_num_partition_ghost_cells(self%original_dm)
       num_non_ghost_cells = end_interior_cell - start_cell - num_ghost_cells
       end_non_ghost_cell = start_cell + num_non_ghost_cells
       ao = dm_get_natural_to_global_ao(self%original_dm, self%dist_sf)
       call DMGetLocalToGlobalMapping(self%original_dm, l2g, ierr); CHKERRQ(ierr)

       num_faces = 0
       do ibdy = 1, num_boundaries
          write(istr, '(i0)') ibdy - 1
          bdystr = 'boundaries[' // trim(istr) // ']'

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
                   case (TYPE_OBJECT)
                      num_faces = 0
                      face_json => fson_value_children_mpi(faces_json)
                      do i = 1, num_face_items
                         num_cells = fson_value_count_mpi(face_json, "cells")
                         num_faces = num_faces + num_cells
                         face_json => fson_value_next_mpi(face_json)
                      end do
                      allocate(faces(num_faces))
                      face_json => fson_value_children_mpi(faces_json)
                      offset = 0
                      do i = 1, num_face_items
                         call fson_get_mpi(face_json, "cells", default_cells, cells, &
                              logfile, log_key = trim(bdystr) // "faces.cells")
                         num_cells = size(cells)
                         allocate(local_cells(num_cells))
                         call AOApplicationToPetsc(ao, num_cells, cells, ierr); CHKERRQ(ierr)
                         call ISGlobalToLocalMappingApplyBlock(l2g, IS_GTOLM_MASK, num_cells, &
                              cells, nout, local_cells, ierr); CHKERRQ(ierr)
                         call fson_get_mpi(face_json, "normal", default_normal, &
                              input_normal, logfile, log_key = trim(bdystr) // "faces.normal")
                         call get_cell_faces(local_cells, num_cells, input_normal, offset)
                         offset = offset + num_cells
                         deallocate(cells, local_cells)
                         face_json => fson_value_next_mpi(face_json)
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
                num_cells = size(cells)
                allocate(local_cells(num_cells))
                call AOApplicationToPetsc(ao, num_cells, cells, ierr); CHKERRQ(ierr)
                call ISGlobalToLocalMappingApplyBlock(l2g, IS_GTOLM_MASK, num_cells, &
                     cells, nout, local_cells, ierr); CHKERRQ(ierr)
                call fson_get_mpi(faces_json, "normal", default_normal, &
                     input_normal, logfile, log_key = trim(bdystr) // "faces.normal")
                num_faces = num_cells
                allocate(faces(num_faces))
                call get_cell_faces(local_cells, num_cells, input_normal, 0)
                deallocate(cells, local_cells)
             case default
                if (present(logfile)) then
                   call logfile%write(LOG_LEVEL_WARN, "input", &
                        "unrecognised_faces_type")
                end if
             end select

          end if

          do iface = 1, num_faces
             associate(f => faces(iface))
               if (f >= 0) then
                  call DMSetLabelValue(self%original_dm, open_boundary_label_name, &
                       f, ibdy, ierr); CHKERRQ(ierr)
               end if
             end associate
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

          bdy => fson_value_next_mpi(bdy)

       end do
       call AODestroy(ao, ierr); CHKERRQ(ierr)

    else if (present(logfile)) then
       call logfile%write(LOG_LEVEL_WARN, "input", "no_boundary_conditions")
    end if

  contains

    subroutine get_cell_faces(cells, num_cells, input_normal, offset)
      ! Get faces for normal vector and specified cells.

      PetscInt, intent(in) :: cells(:), num_cells
      PetscReal, intent(in) :: input_normal(:)
      PetscInt, intent(in) :: offset
      ! Locals:
      PetscInt :: f

      normal_len = size(input_normal)
      normal = 0._dp
      normal(1: normal_len) = input_normal

      do icell = 1, num_cells
         iface = offset + icell
         associate(c => cells(icell))
           if (c >= 0) then
              if ((start_cell <= c) .and. (c <= end_non_ghost_cell)) then
                 call dm_cell_normal_face(self%original_dm, c, normal, f)
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

  subroutine mesh_read_overridden_face_properties(self, json, logfile)

    !! Reads in serial data structures for overridden face properties
    !! from JSON input - currently just face permeability directions.

    use fson_mpi_module
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
    call DMCreateLabel(self%serial_dm, face_property_override_label_name, &
         ierr); CHKERRQ(ierr)

    if (rank == 0) then

       default_cells = [PetscInt::] ! empty integer array
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
                             face_property_override_label_name, f, &
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
    use face_module
    use dm_utils_module, only: local_vec_section, section_offset
    class(mesh_type), intent(in out) :: self
    ! Locals:
    DMLabel :: label
    IS :: faceIS
    PetscInt :: dim, num_faces
    PetscInt, pointer, contiguous :: faces(:)
    PetscSection :: face_section
    type(face_type) :: face
    PetscReal, pointer, contiguous :: face_geom_array(:)
    PetscInt :: i, face_offset, dirn
    PetscErrorCode :: ierr

    call DMGetDimension(self%original_dm, dim, ierr); CHKERRQ(ierr)
    call DMGetLabel(self%original_dm, face_property_override_label_name, &
         label, ierr); CHKERRQ(ierr)

    call local_vec_section(self%face_geom, face_section)
    call VecGetArrayF90(self%face_geom, face_geom_array, ierr); CHKERRQ(ierr)
    call face%init()

    do dirn = 1, dim
       call DMLabelGetStratumSize(label, dirn, num_faces, ierr); CHKERRQ(ierr)
       if (num_faces > 0) then
          call DMLabelGetStratumIS(label, dirn, faceIS, ierr); CHKERRQ(ierr)
          call ISGetIndicesF90(faceIS, faces, ierr); CHKERRQ(ierr)
          do i = 1, num_faces
             associate(f => faces(i))
               if (self%ghost_face(f) < 0) then
                  call section_offset(face_section, f, face_offset, ierr)
                  CHKERRQ(ierr)
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
          ztype = get_zone_type(zone_json)
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
                   call zone%label_dm(self%dm, self%cell_order, &
                        self%cell_geom, err)
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

    err = 0

    if (fson_has_mpi(json, "mesh.minc")) then

       call DMCreateLabel(self%dm, minc_zone_label_name, ierr)
       CHKERRQ(ierr)
       call DMCreateLabel(self%dm, minc_rocktype_zone_label_name, ierr)
       CHKERRQ(ierr)

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
               self%cell_order, iminc, mincstr, self%rock_types, &
               minc_rocktype_zone_index, logfile, err)
       case (TYPE_ARRAY)
          num_minc_zones = fson_value_count_mpi(minc_json, ".")
          minci_json => fson_value_children_mpi(minc_json)
          allocate(self%minc(num_minc_zones))
          do iminc = 1, num_minc_zones
             write(imincstr, '(i0)') iminc - 1
             mincstr = 'minc[' // trim(imincstr) // '].'
             call self%minc(iminc)%init(minci_json, self%dm, self%cell_order, iminc, &
                  mincstr, self%rock_types, minc_rocktype_zone_index, logfile, err)
             if (err > 0) exit
             minci_json => fson_value_next_mpi(minci_json)
          end do
       end select

       call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
       call DMGetLabelSize(self%dm, minc_zone_label_name, &
            num_minc_cells, ierr); CHKERRQ(ierr)
       has_minc_local = (num_minc_cells > 0)
       call MPI_reduce(has_minc_local, self%has_minc, 1, MPI_LOGICAL, MPI_LOR, &
         0, PETSC_COMM_WORLD, ierr)
       call MPI_bcast(self%has_minc, 1, MPI_LOGICAL, 0, &
            PETSC_COMM_WORLD, ierr)
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

    use dm_utils_module
    use list_module

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: dof !! Degrees of freedom for discretization
    ! Locals:
    DM :: minc_dm
    PetscInt :: start_chart, end_chart, m
    PetscInt :: dim
    PetscInt :: num_minc_cells
    PetscInt :: num_minc_zones, num_cells, num_new_points, max_num_levels
    PetscInt, allocatable :: minc_zone(:), minc_end_interior(:)
    PetscInt, allocatable :: stratum_shift(:)
    type(list_type), allocatable :: minc_level_cells(:)
    PetscErrorCode :: ierr

    call dm_get_strata(self%dm, self%depth, self%strata)
    allocate(stratum_shift(0: self%depth))

    call DMPlexCreate(PETSC_COMM_WORLD, minc_dm, ierr); CHKERRQ(ierr)

    call DMGetDimension(self%dm, dim, ierr); CHKERRQ(ierr)
    call DMSetDimension(minc_dm, dim, ierr); CHKERRQ(ierr)
    call DMPlexGetChart(self%dm, start_chart, end_chart, ierr)
    CHKERRQ(ierr)

    num_cells = self%strata(0)%size()
    num_minc_zones = size(self%minc)
    max_num_levels = maxval(self%minc%num_levels)

    call self%setup_minc_dm_zone_and_cells(num_cells, &
       max_num_levels, num_minc_zones, minc_zone, minc_level_cells, &
       num_minc_cells)

    num_new_points = num_minc_cells * (self%depth + 1)
    call DMPlexSetChart(minc_dm, start_chart, &
         end_chart + num_new_points, ierr); CHKERRQ(ierr)
    call self%setup_minc_dm_strata_shifts(num_minc_cells, &
         max_num_levels, minc_level_cells, stratum_shift)

    allocate(minc_end_interior(0: self%depth))
    minc_end_interior = self%strata%end_interior + &
         (stratum_shift + 1) * num_minc_cells
    call DMPlexSetHybridBounds(minc_dm, minc_end_interior(0), &
         minc_end_interior(1), minc_end_interior(self%depth - 1), &
         minc_end_interior(self%depth), ierr); CHKERRQ(ierr)
    deallocate(minc_end_interior)
    self%strata%num_minc_points = num_minc_cells

    call self%set_minc_dm_cone_sizes(minc_dm, num_cells, num_minc_cells, &
         max_num_levels, minc_zone, minc_level_cells)
    call DMSetUp(minc_dm, ierr); CHKERRQ(ierr)
    call self%set_minc_dm_cones(minc_dm, num_cells, max_num_levels, &
         minc_zone, minc_level_cells)

    call DMPlexSymmetrize(minc_dm, ierr); CHKERRQ(ierr)
    call self%transfer_labels_to_minc_dm(minc_dm, max_num_levels)
    call self%setup_minc_dm_depth_label(minc_dm, max_num_levels, &
         minc_level_cells)
    call self%setup_minc_dm_level_label_and_cell_map(minc_dm, max_num_levels, &
         minc_level_cells)

    call dm_set_fv_adjacency(minc_dm)
    call dm_setup_fv_discretization(minc_dm, dof)
    call set_dm_default_data_layout(minc_dm, dof)
    call self%setup_minc_point_sf(minc_dm)
    call dm_setup_global_section(minc_dm)

    call self%setup_minc_dm_cell_order(minc_dm, max_num_levels, &
         minc_level_cells)
    call self%setup_minc_geometry(minc_dm, num_cells, max_num_levels, &
         num_minc_zones, minc_zone)

    do m = 0, max_num_levels
       call minc_level_cells(m)%destroy()
    end do
    deallocate(minc_zone, minc_level_cells, stratum_shift)

    self%dm = minc_dm

  end subroutine mesh_setup_minc_dm

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_dm_zone_and_cells(self, num_cells, &
       max_num_levels, num_minc_zones, minc_zone, minc_level_cells, &
       num_minc_cells)
    !! Set up minc_zone and minc_level_cells arrays, and minc
    !! level label on DM.

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: num_cells, max_num_levels, num_minc_zones 
    PetscInt, allocatable, intent(out) :: minc_zone(:)
    type(list_type), allocatable, intent(out) :: minc_level_cells(:)
    PetscInt, intent(out) :: num_minc_cells
    ! Locals:
    PetscInt :: iminc, i, c, m, ghost, num_minc_zone_cells
    IS :: minc_IS
    PetscInt, pointer :: minc_cells(:)
    DMLabel :: ghost_label
    PetscErrorCode :: ierr
    PetscInt, pointer :: pc

    allocate(minc_zone(0: num_cells - 1), &
         minc_level_cells(0: max_num_levels))
    minc_zone = -1
    do m = 0, max_num_levels
       call minc_level_cells(m)%init(owner = PETSC_TRUE)
    end do

    call DMGetLabel(self%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    do iminc = 1, num_minc_zones
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
                  minc_zone(c) = iminc
                  do m = 0, num_levels
                     allocate(pc); pc = c
                     call minc_level_cells(m)%append(pc)
                  end do
               end if
            end do
          end associate
          call ISRestoreIndicesF90(minc_IS, minc_cells, ierr)
          CHKERRQ(ierr)
          call ISDestroy(minc_IS, ierr); CHKERRQ(ierr)
       end if
    end do

    num_minc_cells = 0
    do m = 1, max_num_levels
       num_minc_cells = num_minc_cells + minc_level_cells(m)%count
    end do

  end subroutine mesh_setup_minc_dm_zone_and_cells

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_dm_strata_shifts(self, num_minc_cells, &
       max_num_levels, minc_level_cells, stratum_shift)
    !! Set up minc_shift array in DM strata, to determine index shift
    !! from original DM point to corresponding point in MINC DM.

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: num_minc_cells, max_num_levels
    type(list_type), intent(in out) :: minc_level_cells(0: max_num_levels)
    PetscInt, intent(out) :: stratum_shift(0: self%depth)
    ! Locals:
    PetscInt :: i, h, m
    PetscInt :: minc_offset(0: max_num_levels)

    ! Set up stratum_shift array (taking account of the fact that
    ! DMPlex points have the order cells, vertices, faces, edges-
    ! i.e. they are not in depth order.) This represents cumulative shift
    ! resulting from points added to lower-index strata in the DM.
    stratum_shift(0) = 0
    stratum_shift(self%depth) = 1
    stratum_shift(1: self%depth - 1) = [(i + 1, i = 1, self%depth - 1)]

    ! Offset of the start of each MINC level within the stratum (note
    ! this is the same for all strata):
    minc_offset = 0
    do i = 1, max_num_levels
       minc_offset(i) = minc_offset(i - 1) + minc_level_cells(i)%count
    end do

    do h = 0, self%depth
       allocate(self%strata(h)%minc_shift(0: max_num_levels))
       self%strata(h)%minc_shift(0) = stratum_shift(h) * num_minc_cells
       do m = 1, max_num_levels
          self%strata(h)%minc_shift(m) = stratum_shift(h) * num_minc_cells + &
               self%strata(h)%end_interior + minc_offset(m - 1)
       end do
    end do

  end subroutine mesh_setup_minc_dm_strata_shifts

!------------------------------------------------------------------------

  subroutine mesh_set_minc_dm_cone_sizes(self, minc_dm, num_cells, &
       num_minc_cells, max_num_levels, minc_zone, &
       minc_level_cells)
    !! Sets cone sizes for MINC DM.

    class(mesh_type), intent(in out) :: self
    DM, intent(in out) :: minc_dm
    PetscInt, intent(in) :: num_cells, num_minc_cells, max_num_levels
    PetscInt, intent(in) :: minc_zone(0: num_cells - 1)
    type(list_type), intent(in out) :: minc_level_cells(0: max_num_levels)
    ! Locals:
    PetscInt :: ic, p, minc_p, cone_size, new_cone_size
    PetscInt :: iminc, m, h
    PetscErrorCode :: ierr

    ! Non-MINC and fracture cells:
    h = 0
    do p = self%strata(h)%start, self%strata(h)%end - 1
       iminc = minc_zone(p)
       call DMPlexGetConeSize(self%dm, p, cone_size, ierr)
       CHKERRQ(ierr)
       if (iminc > 0) then
          new_cone_size = cone_size + 1
       else
          new_cone_size = cone_size
       end if
       minc_p = self%strata(h)%minc_point(p, 0)
       call DMPlexSetConeSize(minc_dm, minc_p, new_cone_size, ierr)
       CHKERRQ(ierr)
    end do
    ! MINC cells:
    do m = 1, max_num_levels
       ic = 0
       call minc_level_cells(m)%traverse(minc_cell_cone_size_iterator)
    end do

    ! Non-MINC and fracture DAG points for height h > 0:
    do h = 1, self%depth
       do p = self%strata(h)%start, self%strata(h)%end - 1
          minc_p = self%strata(h)%minc_point(p, 0)
          call DMPlexGetConeSize(self%dm, p, cone_size, ierr)
          CHKERRQ(ierr)
          call DMPlexSetConeSize(minc_dm, minc_p, cone_size, ierr)
          CHKERRQ(ierr)
       end do
    end do

    do m = 1, max_num_levels
       ic = 0
       call minc_level_cells(m)%traverse(minc_dag_cone_size_iterator)
    end do

  contains

    subroutine minc_cell_cone_size_iterator(node, stopped)
      !! Sets cone size for MINC cells.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: minc_p, iminc, cone_size
      PetscInt, parameter :: h = 0

      select type (c => node%data)
      type is (PetscInt)
         minc_p = self%strata(h)%minc_point(ic, m)
         iminc = minc_zone(c)
         associate(num_levels => self%minc(iminc)%num_levels)
           if (m < num_levels) then
              cone_size = 2
           else
              cone_size = 1
           end if
           call DMPlexSetConeSize(minc_dm, minc_p, &
                cone_size, ierr); CHKERRQ(ierr)
         end associate
      end select
      ic = ic + 1

    end subroutine minc_cell_cone_size_iterator

    subroutine minc_dag_cone_size_iterator(node, stopped)
      !! Sets cone size for MINC DAG points (height h > 0).

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: minc_p, cone_size, h

      select type (c => node%data)
      type is (PetscInt)
         do h = 1, self%depth
            minc_p = self%strata(h)%minc_point(ic, m)
            if (h < self%depth) then
               cone_size = 1
            else
               cone_size = 0
            end if
            call DMPlexSetConeSize(minc_dm, minc_p, &
                 cone_size, ierr); CHKERRQ(ierr)
         end do
      end select
      ic = ic + 1

    end subroutine minc_dag_cone_size_iterator

  end subroutine mesh_set_minc_dm_cone_sizes

!------------------------------------------------------------------------

  subroutine mesh_set_minc_dm_cones(self, minc_dm, num_cells, max_num_levels, &
       minc_zone, minc_level_cells)
    !! Sets cones for MINC DM.

    use dm_utils_module, only: dm_stratum_type

    class(mesh_type), intent(in out) :: self
    DM, intent(in out) :: minc_dm
    PetscInt, intent(in) :: num_cells, max_num_levels
    PetscInt, intent(in) :: minc_zone(0: num_cells - 1)
    type(list_type), intent(in out) :: minc_level_cells(0: max_num_levels)
    ! Locals:
    PetscInt :: p, m, h, iminc, ic, minc_p
    PetscInt, pointer :: points(:)
    PetscInt, allocatable :: cell_cone(:), minc_cone(:)
    PetscErrorCode :: ierr

    h = 0
    ! Non-MINC cells:
    do p = self%strata(h)%start, self%strata(h)%end - 1
       iminc = minc_zone(p)
       if (iminc <= 0) then
          minc_p = self%strata(h)%minc_point(p, 0)
          call DMPlexGetCone(self%dm, p, points, ierr)
          CHKERRQ(ierr)
          cell_cone = self%strata(h + 1)%minc_point(points, 0)
          call DMPlexSetCone(minc_dm, minc_p, cell_cone, ierr); CHKERRQ(ierr)
          deallocate(cell_cone)
          call DMPlexRestoreCone(self%dm, p, points, ierr)
          CHKERRQ(ierr)
       end if
    end do
    ! Fracture cells:
    m = 0; ic = 0
    call minc_level_cells(m)%traverse(fracture_cell_cone_iterator)
    ! MINC cells:
    do m = 1, max_num_levels
       ic = 0
       call minc_level_cells(m)%traverse(minc_cell_cone_iterator)
    end do

    ! Fracture and non-MINC DAG points for height h > 0:
    do h = 1, self%depth - 1
       do p = self%strata(h)%start, self%strata(h)%end - 1
          minc_p = self%strata(h)%minc_point(p, 0)
          call DMPlexGetCone(self%dm, p, points, ierr); CHKERRQ(ierr)
          minc_cone = self%strata(h + 1)%minc_point(points, 0)
          call DMPlexSetCone(minc_dm, minc_p, minc_cone, ierr)
          CHKERRQ(ierr)
          deallocate(minc_cone)
          call DMPlexRestoreCone(self%dm, p, points, ierr)
          CHKERRQ(ierr)
       end do
    end do
    ! MINC DAG points for height h > 0:
    do m = 1, max_num_levels
       ic = 0
       call minc_level_cells(m)%traverse(minc_dag_cone_iterator)
    end do

  contains

    subroutine fracture_cell_cone_iterator(node, stopped)
      !! Sets cone for fracture cells.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: minc_p, face_p
      PetscInt, pointer :: points(:)
      PetscInt, allocatable :: cell_cone(:)
      PetscInt, parameter :: h = 0, m = 0

      select type (c => node%data)
      type is (PetscInt)
         call DMPlexGetCone(self%dm, c, points, ierr)
         CHKERRQ(ierr)
         minc_p = self%strata(h)%minc_point(c, m)
         face_p = self%strata(h + 1)%minc_point(ic, m + 1)
         cell_cone = [self%strata(h + 1)%minc_point(points, m), [face_p]]
         call DMPlexSetCone(minc_dm, minc_p, cell_cone, ierr); CHKERRQ(ierr)
         call DMPlexRestoreCone(self%dm, c, points, ierr)
         deallocate(cell_cone)
         CHKERRQ(ierr)
      end select
      ic = ic + 1

    end subroutine fracture_cell_cone_iterator

    subroutine minc_cell_cone_iterator(node, stopped)
      !! Sets cone for MINC cells (m > 0).

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: minc_p, iminc, face_p, inner_face_p
      PetscInt, parameter :: h = 0

      select type (c => node%data)
      type is (PetscInt)
         minc_p = self%strata(h)%minc_point(ic, m)
         face_p = self%strata(h + 1)%minc_point(ic, m)
         iminc = minc_zone(c)
         associate(num_levels => self%minc(iminc)%num_levels)
           if (m < num_levels) then
              inner_face_p = self%strata(h + 1)%minc_point(ic, m + 1)
              minc_cone = [face_p, inner_face_p]
           else
              minc_cone = [face_p]
           end if
           call DMPlexSetCone(minc_dm, minc_p, minc_cone, ierr)
           CHKERRQ(ierr)
           deallocate(minc_cone)
         end associate
      end select
      ic = ic + 1

    end subroutine minc_cell_cone_iterator

    subroutine minc_dag_cone_iterator(node, stopped)
      !! Sets cone for MINC DAG points (height h > 0).

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: minc_p, h, above_p

      select type (c => node%data)
      type is (PetscInt)
         do h = 1, self%depth - 1
            minc_p = self%strata(h)%minc_point(ic, m)
            above_p = self%strata(h + 1)%minc_point(ic, m)
            call DMPlexSetCone(minc_dm, minc_p, [above_p], ierr)
            CHKERRQ(ierr)
         end do
      end select
      ic = ic + 1

    end subroutine minc_dag_cone_iterator

  end subroutine mesh_set_minc_dm_cones

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_dm_depth_label(self, minc_dm, max_num_levels, &
       minc_level_cells)
    !! Set DAG depth label for MINC points added to the DM.

    class(mesh_type), intent(in out) :: self
    DM, intent(in out) :: minc_dm
    PetscInt, intent(in) :: max_num_levels
    type(list_type), intent(in out) :: minc_level_cells(0: max_num_levels)
    ! Locals:
    PetscInt :: m, ic

    do m = 1, max_num_levels
       ic = 0
       call minc_level_cells(m)%traverse(minc_depth_label_iterator)
    end do

  contains

    subroutine minc_depth_label_iterator(node, stopped)
      !! Sets depth label for all MINC DAG points.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: minc_p, h
      PetscErrorCode :: ierr

      select type (c => node%data)
      type is (PetscInt)
         do h = 0, self%depth
            minc_p = self%strata(h)%minc_point(ic, m)
            associate(p_depth => self%depth - h)
              call DMSetLabelValue(minc_dm, "depth", &
                   minc_p, p_depth, ierr); CHKERRQ(ierr)
            end associate
         end do
      end select
      ic = ic + 1

    end subroutine minc_depth_label_iterator

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

  subroutine mesh_setup_minc_dm_level_label_and_cell_map(self, minc_dm, &
       max_num_levels, minc_level_cells)
    !! Sets up minc_level label on MINC DM, with MINC level assigned
    !! to all cells. Non-MINC cells are assigned level 0 (as are
    !! fracture cells in MINC zones). Also creates minc_cell_map
    !! array, mapping MINC DM cell local indices to their original
    !! single-porosity cell local indices.

    class(mesh_type), intent(in out) :: self
    DM, intent(in out) :: minc_dm
    PetscInt, intent(in) :: max_num_levels
    type(list_type), intent(in out) :: minc_level_cells(0: max_num_levels)
    ! Locals:
    PetscInt :: m, ic, h, c, ghost, minc_p
    PetscInt :: start_cell, end_cell
    DMLabel :: ghost_label
    PetscErrorCode :: ierr

    call DMGetLabel(minc_dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)
    call DMCreateLabel(minc_dm, minc_level_label_name, ierr)
    CHKERRQ(ierr)

    h = 0
    call DMPlexGetHeightStratum(minc_dm, h, start_cell, end_cell, &
         ierr); CHKERRQ(ierr)
    associate(num_cells => end_cell - start_cell)
      allocate(self%minc_cell_map(0: num_cells - 1))
    end associate

    m = 0
    do c = self%strata(h)%start, self%strata(h)%end - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr)
       if (ghost < 0) then
          minc_p = self%strata(h)%minc_point(c, m)
          call DMSetLabelValue(minc_dm, minc_level_label_name, &
               minc_p, m, ierr); CHKERRQ(ierr)
          self%minc_cell_map(minc_p) = c
       end if
    end do

    do m = 1, max_num_levels
       ic = 0
       call minc_level_cells(m)%traverse(minc_level_label_iterator)
    end do

  contains

    subroutine minc_level_label_iterator(node, stopped)
      !! Sets level label for all MINC cells and faces.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: minc_p
      PetscErrorCode :: ierr

      select type (c => node%data)
      type is (PetscInt)
         minc_p = self%strata(h)%minc_point(ic, m)
         call DMSetLabelValue(minc_dm, minc_level_label_name, &
              minc_p, m, ierr); CHKERRQ(ierr)
         self%minc_cell_map(minc_p) = c
      end select
      ic = ic + 1

    end subroutine minc_level_label_iterator

  end subroutine mesh_setup_minc_dm_level_label_and_cell_map

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_dm_cell_order(self, minc_dm, max_num_levels, &
       minc_level_cells)
    !! Sets up natural-to-global AO for MINC DM, and overwrites
    !! original AO.

    use dm_utils_module
    use utils_module, only: array_cumulative_sum, get_mpi_int_gather_array

    class(mesh_type), intent(in out) :: self
    DM, intent(in) :: minc_dm
    PetscInt, intent(in) :: max_num_levels
    type(list_type), intent(in out) :: minc_level_cells(0: max_num_levels)
    ! Locals:
    AO :: minc_ao
    PetscMPIInt :: rank, num_procs
    PetscInt, allocatable :: natural(:), global(:)
    PetscInt :: start_cell, end_interior_cell, offset, n_all
    PetscInt :: mapping_count, num_ghost_cells, num_non_ghost_cells
    ISLocalToGlobalMapping :: l2g, minc_l2g
    PetscInt :: local_minc_cell_count, total_minc_cell_count
    PetscInt :: ic, c, m, inatural
    PetscInt, allocatable :: minc_global(:), minc_frac_natural(:)
    PetscInt, allocatable :: minc_global_all(:), minc_frac_natural_all(:)
    PetscInt, allocatable :: minc_counts(:), minc_displacements(:)
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD, num_procs, ierr)

    start_cell = self%strata(0)%start
    end_interior_cell = self%strata(0)%end_interior
    num_ghost_cells = dm_get_num_partition_ghost_cells(self%dm)
    num_non_ghost_cells = end_interior_cell - start_cell - num_ghost_cells
    ! Each rank has its own array elements for MINC level 0:
    mapping_count = num_non_ghost_cells

    ! Rank 0 has array elements for MINC cells from all ranks:
    local_minc_cell_count = sum(minc_level_cells(1:)%count)
    call MPI_reduce(local_minc_cell_count, total_minc_cell_count, 1, &
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

    do m = 1, max_num_levels
       associate(n => minc_level_cells(m)%count)

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

         ic = 0
         call minc_level_cells(m)%traverse(minc_indices_iterator)

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

    self%cell_order = minc_ao

  contains

!........................................................................

    subroutine get_original_cell_natural_global_indices()
      !! Gets natural and global indices for original cells on each
      !! process, and returns number of cells.

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
            natural(ic) = local_to_natural_cell_index(self%cell_order, &
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

    subroutine minc_indices_iterator(node, stopped)
      !! Gets global indices and natural indices of fracture cells for
      !! MINC cells on current process in specified MINC level.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: p(1)
      PetscErrorCode :: ierr

      select type (c => node%data)
      type is (PetscInt)
         p = self%strata(0)%minc_point(ic, m)
         ic = ic + 1
         minc_frac_natural(ic) = local_to_natural_cell_index( &
              self%cell_order, l2g, c)
         call ISLocalToGlobalMappingApplyBlock(minc_l2g, 1, p, &
              minc_global(ic:ic), ierr); CHKERRQ(ierr)
      end select

    end subroutine minc_indices_iterator

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

  end subroutine mesh_setup_minc_dm_cell_order

!------------------------------------------------------------------------

  subroutine mesh_setup_minc_geometry(self, minc_dm, num_cells, &
       max_num_levels, num_minc_zones, minc_zone)
    !! Sets up cell and face geometry vectors for MINC mesh. These
    !! overwrite the original geometry vectors.

    use kinds_module
    use dm_utils_module, only: section_offset, local_vec_section, &
         set_dm_data_layout
    use cell_module
    use face_module

    class(mesh_type), intent(in out) :: self
    DM, intent(in out) :: minc_dm
    PetscInt, intent(in) :: num_cells, max_num_levels, num_minc_zones
    PetscInt, intent(in) :: minc_zone(0: num_cells - 1)
    ! Locals:
    PetscInt :: dim, cell_variable_dim(num_cell_variables), &
         face_variable_dim(num_face_variables)
    PetscInt :: iminc, c, f, h, i, m, minc_p, face_p, ghost
    PetscInt :: cell_offset, minc_cell_offset, face_offset, minc_face_offset
    PetscInt :: ic(max_num_levels), num_minc_zone_cells
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
    PetscErrorCode :: ierr
    PetscInt, parameter :: nc = 1, np = 1 ! dummy values for cell & face init

    call DMGetDimension(self%dm, dim, ierr); CHKERRQ(ierr)
    call DMGetLabel(self%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    call DMClone(minc_dm, dm_cell, ierr); CHKERRQ(ierr)
    cell_variable_dim = dim
    call set_dm_data_layout(dm_cell, cell_variable_num_components, &
         cell_variable_dim, cell_variable_names)
    call DMCreateLocalVector(dm_cell, minc_cell_geom, ierr); CHKERRQ(ierr)
    call local_vec_section(minc_cell_geom, minc_cell_section)
    call VecGetArrayF90(minc_cell_geom, minc_cell_geom_array, ierr)
    CHKERRQ(ierr)

    call local_vec_section(self%cell_geom, cell_section)
    call VecGetArrayReadF90(self%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)

    call DMClone(minc_dm, dm_face, ierr); CHKERRQ(ierr)
    face_variable_dim = dim - 1
    call set_dm_data_layout(dm_face, face_variable_num_components, &
         face_variable_dim, face_variable_names)
    call DMCreateLocalVector(dm_face, minc_face_geom, ierr); CHKERRQ(ierr)
    call local_vec_section(minc_face_geom, minc_face_section)
    call VecGetArrayF90(minc_face_geom, minc_face_geom_array, ierr)
    CHKERRQ(ierr)

    call local_vec_section(self%face_geom, face_section)
    call VecGetArrayReadF90(self%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)

    ! Copy original cell geometry:
    h = 0
    call cell%init(nc, np)
    do c = self%strata(h)%start, self%strata(h)%end - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          call section_offset(cell_section, c, cell_offset, ierr)
          CHKERRQ(ierr)
          call section_offset(minc_cell_section, c, minc_cell_offset, ierr)
          CHKERRQ(ierr)
          minc_cell_geom_array(minc_cell_offset: &
               minc_cell_offset + cell%dof - 1) = cell_geom_array(cell_offset: &
               cell_offset + cell%dof - 1)
       end if
    end do

    ! Copy original face geometry:
    h = 1
    call face%init(nc, np)
    do f = self%strata(h)%start, self%strata(h)%end - 1
       call DMLabelGetValue(ghost_label, f, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          call section_offset(face_section, f, face_offset, ierr)
          CHKERRQ(ierr)
          minc_p = self%strata(h)%minc_point(f, 0)
          call section_offset(minc_face_section, minc_p, minc_face_offset, ierr)
          CHKERRQ(ierr)
          minc_face_geom_array(minc_face_offset: &
               minc_face_offset + face%dof - 1) = face_geom_array(face_offset: &
               face_offset + face%dof - 1)
       end if
    end do

    h = 0
    ic = 0
    do iminc = 1, num_minc_zones
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
                  minc_p = self%strata(h)%minc_point(c, 0)
                  call section_offset(minc_cell_section, minc_p, minc_cell_offset, ierr)
                  CHKERRQ(ierr)
                  call cell%assign_geometry(minc_cell_geom_array, minc_cell_offset)
                  orig_volume = cell%volume
                  orig_centroid = cell%centroid
                  ! Modify fracture cell volume:
                  cell%volume = orig_volume * minc%volume(1)
                  do m = 1, minc%num_levels
                     ! Assign MINC cell geometry:
                     minc_p = self%strata(h)%minc_point(ic(m), m)
                     call section_offset(minc_cell_section, minc_p, &
                          minc_cell_offset, ierr); CHKERRQ(ierr)
                     call cell%assign_geometry(minc_cell_geom_array, minc_cell_offset)
                     cell%volume = orig_volume * minc%volume(m + 1)
                     cell%centroid = orig_centroid
                     ! Assign MINC face geometry:
                     face_p = self%strata(h + 1)%minc_point(ic(m), m)
                     call section_offset(minc_face_section, face_p, &
                          minc_face_offset, ierr); CHKERRQ(ierr)
                     call face%assign_geometry(minc_face_geom_array, minc_face_offset)
                     face%area = orig_volume * minc%connection_area(m)
                     face%distance = minc%connection_distance(m: m + 1)
                     face%normal = 0._dp
                     face%gravity_normal = 0._dp
                     face%centroid = 0._dp
                     face%permeability_direction = dble(1)
                     ic(m) = ic(m) + 1
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
    PetscInt, intent(out) :: rock_range_start
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
                              call global_section_offset(section, c, rock_range_start, &
                                   orig_offset, ierr); CHKERRQ(ierr)
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
                                   call global_section_offset(section, cell_p, rock_range_start, &
                                        offset, ierr); CHKERRQ(ierr)
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
                    // trim(name) // ".type"
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

  PetscInt function mesh_local_to_fracture_natural(self, local) &
       result(natural)
    !! Takes a local cell index and returns natural index of the
    !! corresponding fracture cell. (For a non-MINC mesh, the
    !! 'fracture' cell is just the original cell itself.)

    use dm_utils_module, only: local_to_natural_cell_index

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: local !! Local cell index
    ! Locals:
    PetscInt :: fracture_local
    ISLocalToGlobalMapping :: l2g
    PetscErrorCode :: ierr

    if (self%has_minc) then
       fracture_local = self%minc_cell_map(local)
    else
       fracture_local = local
    end if

    call DMGetLocalToGlobalMapping(self%dm, l2g, ierr); CHKERRQ(ierr)
    natural = local_to_natural_cell_index(self%cell_order, &
         l2g, fracture_local)

  end function mesh_local_to_fracture_natural

!------------------------------------------------------------------------

  subroutine mesh_global_to_fracture_natural(self, global, &
       fracture_natural, minc_level)
    !! Takes a global cell index and returns natural index of
    !! corresponding fracture cell, together with the MINC level of
    !! the cell. (For a non-MINC mesh, the 'fracture' cell is just the
    !! original cell itself, and the MINC level is zero.)

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: global !! Global cell index
    PetscInt, intent(out) :: fracture_natural !! Natural index of fracture cell
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
               fracture_natural = self%local_to_fracture_natural(c)
               minc_level = self%local_cell_minc_level(c)
               CHKERRQ(ierr)
               process_found_rank = rank
            end if
         end if
         call MPI_Allreduce(process_found_rank, found_rank, 1, MPI_INT, &
              MPI_MAX, PETSC_COMM_WORLD, ierr)
         call MPI_bcast(fracture_natural, 1, MPI_LOGICAL, found_rank, &
            PETSC_COMM_WORLD, ierr)
         call MPI_bcast(minc_level, 1, MPI_LOGICAL, found_rank, &
            PETSC_COMM_WORLD, ierr)
       end associate
    else
       call AOPetscToApplication(self%cell_order, 1, idx, ierr); CHKERRQ(ierr)
       fracture_natural = idx(1)
       minc_level = 0
    end if

  end subroutine mesh_global_to_fracture_natural

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
