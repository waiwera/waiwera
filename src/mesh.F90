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
  use cell_order_module
  use zone_module
  use list_module
  use mpi_utils_module
  use fson
  use minc_module

  implicit none

  private

  PetscInt, parameter, public :: max_mesh_filename_length = 200
  character(len = 16), public :: open_boundary_label_name = "open_boundary" !! Name of DMLabel for identifying open boundaries

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
     type(list_type), public :: zones !! Mesh zones
     PetscBool, public :: radial !! If mesh coordinate system is radial or Cartesian
     PetscBool, public :: has_minc !! If mesh has any MINC cells
   contains
     procedure :: assign_dm => mesh_assign_dm
     procedure :: setup_cell_order_label => mesh_setup_cell_order_label
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
     procedure :: setup_zones => mesh_setup_zones
     procedure :: setup_minc => mesh_setup_minc
     procedure :: setup_minc_dm => mesh_setup_minc_dm
     procedure :: setup_minc_dm_zone_and_cells => mesh_setup_minc_dm_zone_and_cells
     procedure :: setup_minc_dm_shifts => mesh_setup_minc_dm_shifts
     procedure :: set_minc_dm_cone_sizes => mesh_set_minc_dm_cone_sizes
     procedure :: set_minc_dm_cones => mesh_set_minc_dm_cones
     procedure :: set_minc_dm_depth_labels => mesh_set_minc_dm_depth_labels
     procedure :: transfer_labels_to_minc_dm => mesh_transfer_labels_to_minc_dm
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

    self%ghost_cell = -1
    self%ghost_face = -1

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

  subroutine mesh_configure(self, dof, gravity, json, logfile, viewer, err)
    !! Configures mesh, including distribution over processes and
    !! construction of ghost cells, setup of data layout, geometry and
    !! cell index set.

    use dm_utils_module, only: dm_setup_fv_discretization, &
         set_dm_default_data_layout
    use cell_order_module
    use logfile_module

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: dof !! Number of degrees of freedom (primary thermodynamic variables)
    PetscReal, intent(in) :: gravity(:)
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    type(logfile_type), intent(in out), optional :: logfile !! Log file
    PetscViewer, intent(in out), optional :: viewer !! PetscViewer for output of cell index set to HDF5 file
    PetscErrorCode, intent(out) :: err !! Error flag

    err = 0
    call dm_setup_cell_order_label(self%dm)

    call self%distribute()
    call self%construct_ghost_cells()
    call dm_setup_fv_discretization(self%original_dm, dof)
    call set_dm_default_data_layout(self%original_dm, dof)

    call self%setup_geometry(gravity)
    call self%setup_zones(json, logfile, err)

    if (err == 0) then
       call self%setup_minc(json, logfile, err)
       if ((err == 0) .and. self%has_minc) then
          call self%setup_minc_dm(dof)
       end if
    end if

    call self%get_bounds()
    call self%setup_ghost_arrays()
    call dm_setup_cell_index(self%dm, self%cell_index, viewer)

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

    call self%zones%destroy(mesh_zones_node_data_destroy)

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

       do i = 1, num_zones
          zone => null()
          zone_json => fson_value_get_mpi(zones_json, i)
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

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: num_minc_zones, iminc, num_minc_cells
    type(fson_value), pointer :: minc_json, minci_json
    PetscInt :: minc_type
    character(32) :: imincstr, mincstr
    PetscMPIInt :: rank
    PetscErrorCode :: ierr
    PetscBool :: has_minc_local

    err = 0

    if (fson_has_mpi(json, "mesh.minc")) then

       call DMCreateLabel(self%original_dm, minc_zone_label_name, ierr)
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

       call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
       call DMGetLabelSize(self%original_dm, minc_zone_label_name, &
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
    PetscInt :: start_chart, end_chart, h, m
    PetscInt :: dim, depth
    PetscInt :: num_minc_cells
    PetscInt :: num_minc_zones, num_cells, num_new_points, max_num_levels
    PetscInt, allocatable :: stratum_start(:), stratum_end(:)
    PetscInt, allocatable :: frac_shift(:), minc_shift(:,:)
    PetscInt, allocatable :: minc_zone(:)
    type(list_type), allocatable :: minc_level_cells(:)
    PetscErrorCode :: ierr

    call DMPlexCreate(PETSC_COMM_WORLD, self%minc_dm, ierr); CHKERRQ(ierr)

    call DMGetDimension(self%original_dm, dim, ierr); CHKERRQ(ierr)
    call DMSetDimension(self%minc_dm, dim, ierr); CHKERRQ(ierr)
    call DMPlexGetDepth(self%original_dm, depth, ierr); CHKERRQ(ierr)
    call DMPlexGetChart(self%original_dm, start_chart, end_chart, ierr)
    CHKERRQ(ierr)
    allocate(stratum_start(0: depth), stratum_end(0: depth))
    do h = 0, depth
       call DMPlexGetHeightStratum(self%original_dm, h, stratum_start(h), &
            stratum_end(h), ierr); CHKERRQ(ierr)
    end do

    num_cells = stratum_end(0) - stratum_start(0)
    num_minc_zones = size(self%minc)
    max_num_levels = maxval(self%minc%num_levels)

    call self%setup_minc_dm_zone_and_cells(num_cells, &
       max_num_levels, num_minc_zones, minc_zone, minc_level_cells, &
       num_minc_cells)

    num_new_points = num_minc_cells * (depth + 1)
    call DMPlexSetChart(self%minc_dm, start_chart, &
         end_chart + num_new_points, ierr); CHKERRQ(ierr)

    allocate(frac_shift(0: depth), minc_shift(0: depth, 1: max_num_levels))
    call self%setup_minc_dm_shifts(num_minc_cells, depth, &
         max_num_levels, minc_level_cells, stratum_start, stratum_end, &
         frac_shift, minc_shift)

    call self%set_minc_dm_cone_sizes(depth, num_cells, max_num_levels, &
         minc_zone, frac_shift, minc_shift, stratum_start, stratum_end, &
         minc_level_cells)
    call DMSetUp(self%minc_dm, ierr); CHKERRQ(ierr)
    call self%set_minc_dm_cones(depth, num_cells, max_num_levels, &
         minc_zone, frac_shift, minc_shift, stratum_start, stratum_end, &
         minc_level_cells)

    call DMPlexSymmetrize(self%minc_dm, ierr); CHKERRQ(ierr)
    call self%transfer_labels_to_minc_dm(1, depth, stratum_start, &
         stratum_end, frac_shift)
    call self%set_minc_dm_depth_labels(depth, max_num_levels, &
         minc_shift, minc_level_cells)

    call dm_set_fv_adjacency(self%minc_dm)
    call dm_setup_fv_discretization(self%minc_dm, dof)
    call set_dm_default_data_layout(self%minc_dm, dof)

    call self%assign_dm(self%minc_dm)

    do m = 0, max_num_levels
       call minc_level_cells(m)%destroy()
    end do
    deallocate(stratum_start, stratum_end, frac_shift, &
         minc_shift, minc_zone, minc_level_cells)

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

    call DMGetLabel(self%original_dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    do iminc = 1, num_minc_zones
       call DMGetStratumSize(self%original_dm, minc_zone_label_name, iminc, &
            num_minc_zone_cells, ierr); CHKERRQ(ierr)
       if (num_minc_zone_cells > 0) then
          call DMGetStratumIS(self%original_dm, minc_zone_label_name, &
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

  subroutine mesh_setup_minc_dm_shifts(self, num_minc_cells, depth, &
       max_num_levels, minc_level_cells, stratum_start, stratum_end, &
       frac_shift, minc_shift)
    !! Set up shift arrays to determine index shift from original DM
    !! point to corresponding fracture and MINC points in new DM.

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: num_minc_cells, depth, max_num_levels
    type(list_type), intent(in out) :: minc_level_cells(0: max_num_levels)
    PetscInt, intent(in) :: stratum_start(0: depth), stratum_end(0: depth)
    PetscInt, intent(out) :: frac_shift(0: depth)
    PetscInt, intent(out) :: minc_shift(0: depth, 1: max_num_levels)
    ! Locals:
    PetscInt :: ishift(0: depth)
    PetscInt :: i, h, m
    PetscInt :: minc_offset(0: max_num_levels)

    !! Set up ishift array, to take account of the fact that DMPlex
    !! points have the order cells, vertices, faces, edges-
    !! i.e. they are not in depth order.
    ishift(0) = 0
    ishift(depth) = 1
    ishift(1: depth - 1) = [(i + 1, i = 1, depth - 1)]

    minc_offset = 0
    do i = 1, max_num_levels
       minc_offset(i) = minc_offset(i - 1) + minc_level_cells(i)%count
    end do

    frac_shift = ishift * num_minc_cells
    minc_shift = 0
    do h = 0, depth
       do m = 1, max_num_levels
          minc_shift(h, m) = stratum_end(h) + frac_shift(h) + &
               minc_offset(m - 1)
       end do
    end do

  end subroutine mesh_setup_minc_dm_shifts

!------------------------------------------------------------------------

  subroutine mesh_set_minc_dm_cone_sizes(self, depth, num_cells, &
       max_num_levels, minc_zone, frac_shift, minc_shift, &
       stratum_start, stratum_end, minc_level_cells)
    !! Sets cone sizes for MINC DM.

    use dm_utils_module, only: dm_copy_cone_sizes

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: depth, num_cells, max_num_levels
    PetscInt, intent(in) :: minc_zone(0: num_cells - 1)
    PetscInt, intent(in) :: frac_shift(0: depth)
    PetscInt, intent(in) :: minc_shift(0: depth, 1: max_num_levels)
    PetscInt, intent(in) :: stratum_start(0: depth), stratum_end(0: depth)
    type(list_type), intent(in out) :: minc_level_cells(0: max_num_levels)
    ! Locals:
    PetscInt :: ic, c, cone_size, new_cone_size
    PetscInt :: iminc, m, h
    PetscErrorCode :: ierr

    ! Fracture cells:
    h = 0
    do c = stratum_start(h), stratum_end(h) - 1
       iminc = minc_zone(c)
       call DMPlexGetConeSize(self%original_dm, c, cone_size, ierr)
       CHKERRQ(ierr)
       if (iminc > 0) then
          new_cone_size = cone_size + 1
       else
          new_cone_size = cone_size
       end if
       call DMPlexSetConeSize(self%minc_dm, c, new_cone_size, ierr)
       CHKERRQ(ierr)
    end do
    ! MINC cells:
    do m = 1, max_num_levels
       ic = 0
       call minc_level_cells(m)%traverse(minc_cell_cone_size_iterator)
    end do

    ! Fracture DAG points for height h > 0:
    do h = 1, depth
       call dm_copy_cone_sizes(self%original_dm, self%minc_dm, &
            stratum_start(h), stratum_end(h) - 1, frac_shift(h))
    end do
    ! MINC DAG points for height h > 0:
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
         minc_p = ic + minc_shift(h, m)
         iminc = minc_zone(c)
         associate(num_levels => self%minc(iminc)%num_levels)
           if (m < num_levels) then
              cone_size = 2
           else
              cone_size = 1
           end if
           call DMPlexSetConeSize(self%minc_dm, minc_p, &
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
         do h = 1, depth
            minc_p = ic + minc_shift(h, m)
            if (h < depth) then
               cone_size = 1
            else
               cone_size = 0
            end if
            call DMPlexSetConeSize(self%minc_dm, minc_p, &
                 cone_size, ierr); CHKERRQ(ierr)
         end do
      end select
      ic = ic + 1

    end subroutine minc_dag_cone_size_iterator

  end subroutine mesh_set_minc_dm_cone_sizes

!------------------------------------------------------------------------

  subroutine mesh_set_minc_dm_cones(self, depth, num_cells, &
       max_num_levels, minc_zone, frac_shift, minc_shift, &
       stratum_start, stratum_end, minc_level_cells)
    !! Sets cones for MINC DM.

    use dm_utils_module, only: dm_copy_cones

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: depth, num_cells, max_num_levels
    PetscInt, intent(in) :: minc_zone(0: num_cells - 1)
    PetscInt, intent(in) :: frac_shift(0: depth)
    PetscInt, intent(in) :: minc_shift(0: depth, 1: max_num_levels)
    PetscInt, intent(in) :: stratum_start(0: depth), stratum_end(0: depth)
    type(list_type), intent(in out) :: minc_level_cells(0: max_num_levels)
    ! Locals:
    PetscInt :: c, m, h, iminc, ic
    PetscInt :: cone_shift
    PetscInt, pointer :: points(:)
    PetscInt, allocatable :: cell_cone(:), minc_cone(:)
    PetscErrorCode :: ierr

    h = 0
    ! Non-MINC cells:
    do c = stratum_start(h), stratum_end(h) - 1
       iminc = minc_zone(c)
       if (iminc <= 0) then
          call DMPlexGetCone(self%original_dm, c, points, ierr)
          CHKERRQ(ierr)
          cell_cone = points + frac_shift(h + 1)
          call DMPlexSetCone(self%minc_dm, c, cell_cone, ierr); CHKERRQ(ierr)
          call DMPlexRestoreCone(self%original_dm, c, points, ierr)
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

    ! Fracture DAG points for height h > 0:
    do h = 1, depth
       if (h < depth) then
          cone_shift = frac_shift(h + 1)
       else
          cone_shift = 0
       end if
       call dm_copy_cones(self%original_dm, self%minc_dm, &
            stratum_start(h), stratum_end(h) - 1, &
            frac_shift(h), cone_shift)
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
      PetscInt :: face_p
      PetscInt, pointer :: points(:)
      PetscInt, allocatable :: cell_cone(:)
      PetscInt, parameter :: h = 0, m = 0

      select type (c => node%data)
      type is (PetscInt)
         call DMPlexGetCone(self%original_dm, c, points, ierr)
         CHKERRQ(ierr)
         face_p = ic + minc_shift(h + 1, m + 1)
         cell_cone = [points + frac_shift(h + 1), [face_p]]
         call DMPlexSetCone(self%minc_dm, c, cell_cone, ierr); CHKERRQ(ierr)
         call DMPlexRestoreCone(self%original_dm, c, points, ierr)
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
         minc_p = ic + minc_shift(h, m)
         face_p = ic + minc_shift(h + 1, m)
         iminc = minc_zone(c)
         associate(num_levels => self%minc(iminc)%num_levels)
           if (m < num_levels) then
              inner_face_p = ic + minc_shift(h + 1, m + 1)
              minc_cone = [face_p, inner_face_p]
           else
              minc_cone = [face_p]
           end if
           call DMPlexSetCone(self%minc_dm, minc_p, &
                minc_cone, ierr); CHKERRQ(ierr)
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
         do h = 1, depth - 1
            minc_p = ic + minc_shift(h, m)
            above_p = ic + minc_shift(h + 1, m)
            call DMPlexSetCone(self%minc_dm, minc_p, &
                 [above_p], ierr); CHKERRQ(ierr)
         end do
      end select
      ic = ic + 1

    end subroutine minc_dag_cone_iterator

  end subroutine mesh_set_minc_dm_cones

!------------------------------------------------------------------------

  subroutine mesh_set_minc_dm_depth_labels(self, depth, max_num_levels, &
       minc_shift, minc_level_cells)
    !! Set DAG depth labels for MINC points added to the DM.

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: depth, max_num_levels
    PetscInt, intent(in) :: minc_shift(0: depth, 1: max_num_levels)
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
         do h = 0, depth
            minc_p = ic + minc_shift(h, m)
            associate(p_depth => depth - h)
              call DMSetLabelValue(self%minc_dm, "depth", &
                   minc_p, p_depth, ierr); CHKERRQ(ierr)
            end associate
         end do
      end select
      ic = ic + 1

    end subroutine minc_depth_label_iterator

  end subroutine mesh_set_minc_dm_depth_labels

!------------------------------------------------------------------------

  subroutine mesh_transfer_labels_to_minc_dm(self, max_height, depth, &
       stratum_start, stratum_end, frac_shift)
    !! Transfers labels from original DM to MINC DM fracture points,
    !! applying appropriate shifts to the point indices.

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: max_height, depth
    PetscInt, intent(in) :: stratum_start(0: depth), stratum_end(0: depth)
    PetscInt, intent(in) :: frac_shift(0: depth)
    ! Locals:
    PetscInt :: p, h, l, iid, ip, label_value
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
       call DMCreateLabel(self%minc_dm, label_name, ierr); CHKERRQ(ierr)
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
                   h = dm_point_height(self%dm, p)
                   call DMSetLabelValue(self%minc_dm, label_name, &
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
    end do

  contains

    PetscInt function dm_point_height(dm, p) result(height)
      !! Returns height of point in DM.

      DM, intent(in) :: dm
      PetscInt, intent(in) :: p
      ! Locals:
      PetscInt :: h

      height = -1
      do h = 0, depth
         if ((stratum_start(h) <= p) .and. (p < stratum_end(h))) then
            height = h
            exit
         end if
      end do

    end function dm_point_height

  end subroutine mesh_transfer_labels_to_minc_dm

!------------------------------------------------------------------------

end module mesh_module
