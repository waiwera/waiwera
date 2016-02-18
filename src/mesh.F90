module mesh_module
  !! Module for mesh type.

  use mpi_module
  use fson

  implicit none

  private

#include <petsc/finclude/petsc.h90>

  PetscInt, parameter, public :: max_mesh_filename_length = 200
  character(len = 13), parameter, public :: &
       natural_order_label_name = "natural order"

  type, public :: mesh_type
     !! Mesh type.
     private
     character(max_mesh_filename_length), public :: filename
     DM, public :: dm
     Vec, public :: cell_geom, face_geom
     PetscInt, public :: start_cell, end_cell, end_interior_cell
     PetscInt, public :: start_face, end_face
     PetscReal, allocatable, public :: bcs(:,:)
     IS, public :: natural_order
   contains
     procedure :: setup_natural_order_label => mesh_setup_natural_order_label
     procedure :: setup_natural_order => mesh_setup_natural_order
     procedure :: distribute => mesh_distribute
     procedure :: construct_ghost_cells => mesh_construct_ghost_cells
     procedure :: setup_data_layout => mesh_setup_data_layout
     procedure :: setup_geometry => mesh_setup_geometry
     procedure :: setup_discretization => mesh_setup_discretization
     procedure :: get_bounds => mesh_get_bounds
     procedure, public :: init => mesh_init
     procedure, public :: configure => mesh_configure
     procedure, public :: setup_boundaries => mesh_setup_boundaries
     procedure, public :: set_boundary_values => mesh_set_boundary_values
     procedure, public :: order_vector => mesh_order_vector
     procedure, public :: destroy => mesh_destroy
  end type mesh_type

contains

!------------------------------------------------------------------------

  subroutine mesh_setup_natural_order_label(self)
    !! Sets up natural ordering label on cells. This assumes the mesh
    !! has not yet been distributed.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscBool :: has_label
    PetscInt :: start_cell, end_cell, c
    PetscErrorCode :: ierr

    call DMHasLabel(self%dm, natural_order_label_name, has_label, &
         ierr); CHKERRQ(ierr)

    if (.not. has_label) then
       call DMCreateLabel(self%dm, natural_order_label_name, ierr)
       CHKERRQ(ierr)
       call DMPlexGetHeightStratum(self%dm, 0, start_cell, end_cell, &
            ierr); CHKERRQ(ierr)
       do c = start_cell, end_cell - 1
          call DMSetLabelValue(self%dm, natural_order_label_name, c, &
               c, ierr); CHKERRQ(ierr)
       end do
    end if

  end subroutine mesh_setup_natural_order_label

!------------------------------------------------------------------------

  subroutine mesh_setup_natural_order(self)
    !! Sets up natural order index set from natural order label on DM.
    !! This index set corresponds to a block size of 1.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscInt :: n, c, i, ghost, order
    DMLabel :: ghost_label, order_label
    PetscInt, allocatable :: order_array(:)
    PetscErrorCode :: ierr

    call DMGetLabel(self%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)

    ! Count interior cells:
    n = 0
    do c = self%start_cell, self%end_interior_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) n = n + 1
    end do
    allocate(order_array(n))

    call DMGetLabel(self%dm, natural_order_label_name, order_label, ierr)
    CHKERRQ(ierr)
    i = 1
    do c = self%start_cell, self%end_interior_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          call DMLabelGetValue(order_label, c, order, ierr); CHKERRQ(ierr)
          order_array(i) = order
          i = i + 1
       end if
    end do

    call ISCreateGeneral(mpi%comm, n, order_array, PETSC_COPY_VALUES, &
         self%natural_order, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(self%natural_order, "natural_order", ierr)
    CHKERRQ(ierr)

    deallocate(order_array)

  end subroutine mesh_setup_natural_order

!------------------------------------------------------------------------

  subroutine mesh_distribute(self)
    !! Distributes mesh over processors.
    
    class(mesh_type), intent(in out) :: self
    ! Locals:
    DM :: dist_dm
    PetscErrorCode :: ierr
    PetscInt, parameter :: overlap = 1

    call DMPlexDistribute(self%dm, overlap, PETSC_NULL_OBJECT, &
         dist_dm, ierr); CHKERRQ(ierr)
    if (dist_dm /= 0) then
       call DMDestroy(self%dm, ierr); CHKERRQ(ierr)
       self%dm = dist_dm
    end if
    
  end subroutine mesh_distribute

!------------------------------------------------------------------------

  subroutine mesh_construct_ghost_cells(self)
    !! Constructs ghost cells on open boundary faces.

    use boundary_module, only: open_boundary_label_name

    class(mesh_type), intent(in out) :: self
    ! Locals:
    DM :: ghost_dm
    PetscErrorCode :: ierr

    call DMPlexConstructGhostCells(self%dm, open_boundary_label_name, &
         PETSC_NULL_INTEGER, ghost_dm, ierr); CHKERRQ(ierr)
    if (ghost_dm /= 0) then
       call DMDestroy(self%dm, ierr); CHKERRQ(ierr);
       self%dm = ghost_dm
    end if

  end subroutine mesh_construct_ghost_cells

!------------------------------------------------------------------------

  subroutine mesh_setup_data_layout(self, dof)
    !! Sets up default section data layout for the mesh.

    use dm_utils_module, only: set_dm_data_layout

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: dof !! Degrees of freedom
    ! Locals:
    PetscInt :: num_components(1), field_dim(1)
    character(7) :: field_names(1)

    num_components = dof
    field_dim = 3
    field_names(1) = "Primary"

    call set_dm_data_layout(self%dm, num_components, field_dim, &
         field_names)

  end subroutine mesh_setup_data_layout

!------------------------------------------------------------------------

  subroutine mesh_setup_discretization(self, dof)
    !! Sets up finite volume discretization on the DM.

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: dof
    ! Locals:
    PetscFV :: fvm
    PetscDS :: ds
    PetscInt :: dim
    PetscErrorCode :: ierr

    call PetscFVCreate(mpi%comm, fvm, ierr); CHKERRQ(ierr)
    call PetscFVSetFromOptions(fvm, ierr); CHKERRQ(ierr)
    call PetscFVSetNumComponents(fvm, dof, ierr); CHKERRQ(ierr)
    call DMGetDimension(self%dm, dim, ierr); CHKERRQ(ierr)
    call PetscFVSetSpatialDimension(fvm, dim, ierr); CHKERRQ(ierr)
    call DMGetDS(self%dm, ds, ierr); CHKERRQ(ierr)
    call PetscDSAddDiscretization(ds, fvm, ierr); CHKERRQ(ierr)
    call PetscFVDestroy(fvm, ierr); CHKERRQ(ierr)

  end subroutine mesh_setup_discretization

!------------------------------------------------------------------------

  subroutine mesh_setup_geometry(self)
    !! Sets up global vectors containing geometry data (e.g. cell volumes,
    !! cell centroids, face areas, face-to-centroid distances) for the mesh.

    use kinds_module
    use face_module
    use dm_utils_module, only: section_offset, local_vec_section
    use boundary_module, only: open_boundary_label_name

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    Vec :: petsc_face_geom
    DM :: dm_face
    PetscSection :: face_section, petsc_face_section, cell_section
    PetscInt :: f, face_dof, ghost_face, i, bdy_face, ibdy
    PetscInt :: num_faces, iface
    PetscInt :: face_offset, petsc_face_offset
    PetscInt :: cell_offset(2)
    type(face_type) :: face
    type(petsc_face_type) :: petsc_face
    PetscReal, pointer :: face_geom_array(:), petsc_face_geom_array(:)
    PetscReal, pointer :: cell_geom_array(:)
    DMLabel :: ghost_label, bdy_label
    PetscInt, pointer :: cells(:)
    IS :: bdy_IS
    PetscInt, pointer :: bdy_faces(:)

    ! First call PETSc geometry routine- we use the cell geometry vector but need to 
    ! create our own face geometry vector, containing additional parameters:
    call DMPlexTSGetGeometryFVM(self%dm, petsc_face_geom, self%cell_geom, &
         PETSC_NULL_REAL, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(self%cell_geom, "cell_geometry", ierr); CHKERRQ(ierr)

    ! Set up face geometry vector:
    call DMClone(self%dm, dm_face, ierr); CHKERRQ(ierr)
    ! replace this with set_dm_data_layout()?
    call PetscSectionCreate(mpi%comm, face_section, ierr); CHKERRQ(ierr)
    call PetscSectionSetChart(face_section, self%start_face, &
         self%end_face, ierr); CHKERRQ(ierr)
    call face%init()
    face_dof = face%dof()
    do f = self%start_face, self%end_face - 1
       call PetscSectionSetDof(face_section, f, face_dof, ierr); CHKERRQ(ierr)
    end do
    call PetscSectionSetUp(face_section, ierr); CHKERRQ(ierr)
    call DMSetDefaultSection(dm_face, face_section, ierr); CHKERRQ(ierr)

    call DMCreateLocalVector(dm_face, self%face_geom, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(self%face_geom, face_geom_array, ierr); CHKERRQ(ierr)

    call local_vec_section(self%cell_geom, cell_section)
    call VecGetArrayF90(self%cell_geom, cell_geom_array, ierr); CHKERRQ(ierr)

    call local_vec_section(petsc_face_geom, petsc_face_section)
    call VecGetArrayF90(petsc_face_geom, petsc_face_geom_array, ierr); CHKERRQ(ierr)

    call DMGetLabel(self%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    call DMGetLabel(self%dm, open_boundary_label_name, bdy_label, ierr)
    CHKERRQ(ierr)

    do f = self%start_face, self%end_face - 1

       call DMLabelGetValue(ghost_label, f, ghost_face, ierr); CHKERRQ(ierr)

       if (ghost_face < 0) then

          call DMLabelGetValue(bdy_label, f, bdy_face, ierr); CHKERRQ(ierr)

          call section_offset(face_section, f, face_offset, ierr); CHKERRQ(ierr)
          call section_offset(petsc_face_section, f, petsc_face_offset, ierr)
          CHKERRQ(ierr)

          call DMPlexGetSupport(self%dm, f, cells, ierr); CHKERRQ(ierr)
          do i = 1, 2
             call section_offset(cell_section, cells(i), cell_offset(i), ierr)
             CHKERRQ(ierr)
          end do

          call petsc_face%assign(petsc_face_geom_array, petsc_face_offset)
          call face%assign(face_geom_array, face_offset, cell_geom_array, cell_offset)

          face%centroid = petsc_face%centroid
          face%area = norm2(petsc_face%area_normal)
          face%normal = petsc_face%area_normal / face%area
          do i = 1, 2
             face%distance(i) = norm2(face%centroid - face%cell(i)%centroid)
          end do
          call face%calculate_permeability_direction()

       end if

    end do

    ! Set external boundary face connection distances to zero:
    do ibdy = 1, size(self%bcs, 2)
       call DMGetStratumIS(self%dm, open_boundary_label_name, ibdy, bdy_IS, &
            ierr); CHKERRQ(ierr)
       if (bdy_IS /= 0) then
          call ISGetIndicesF90(bdy_IS, bdy_faces, ierr); CHKERRQ(ierr)
          num_faces = size(bdy_faces)
          do iface = 1, num_faces
             f = bdy_faces(iface)
             call section_offset(face_section, f, face_offset, ierr)
             CHKERRQ(ierr)
             call face%assign(face_geom_array, face_offset)
             face%distance(2) = 0._dp
          end do
       end if
       call ISDestroy(bdy_IS, ierr); CHKERRQ(ierr)
    end do

    call face%destroy()
    call petsc_face%destroy()
    call VecRestoreArrayF90(self%cell_geom, cell_geom_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(self%face_geom, face_geom_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(petsc_face_geom, petsc_face_geom_array, ierr)
    CHKERRQ(ierr)

    call PetscSectionDestroy(face_section, ierr); CHKERRQ(ierr)
    call DMDestroy(dm_face, ierr); CHKERRQ(ierr)

  end subroutine mesh_setup_geometry

!------------------------------------------------------------------------

  subroutine mesh_get_bounds(self)
    !! Gets cell bounds on current processor.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call DMPlexGetHybridBounds(self%dm, self%end_interior_cell, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         ierr)
    CHKERRQ(ierr)

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

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscErrorCode :: ierr
    type(fson_value), pointer :: mesh

    if (mpi%rank == mpi%input_rank) then
       call fson_get(json, "mesh", mesh)
       if (associated(mesh)) then
          call fson_get(mesh, ".", self%filename)
       else
          if (present(logfile)) then
             call logfile%write(LOG_LEVEL_ERR, 'mesh', 'init', &
                  str_key = 'stop            ', &
                  str_value = 'mesh not found in input.')
          end if
          stop
       end if
    end if

    call MPI_bcast(self%filename, max_mesh_filename_length, MPI_CHARACTER, &
         mpi%input_rank, mpi%comm, ierr)

    ! Read in DM:
    call DMPlexCreateFromFile(mpi%comm, self%filename, PETSC_TRUE, self%dm, ierr)
    CHKERRQ(ierr)

  end subroutine mesh_init

!------------------------------------------------------------------------

  subroutine mesh_configure(self, primary_variable_names)
    !! Configures mesh.

    class(mesh_type), intent(in out) :: self
    character(*), intent(in) :: primary_variable_names(:) !! Names of primary thermodynamic variables
    ! Locals:
    PetscInt :: dof

    call self%setup_natural_order_label()
    call self%distribute()
    call self%construct_ghost_cells()

    dof = size(primary_variable_names)
    call self%setup_discretization(dof)

    call self%get_bounds()

    call self%setup_data_layout(dof)

    call self%setup_geometry()

    call self%setup_natural_order()

  end subroutine mesh_configure

!------------------------------------------------------------------------

  subroutine mesh_destroy(self)
    !! Destroys mesh object.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    
    call VecDestroy(self%face_geom, ierr); CHKERRQ(ierr)
    call DMDestroy(self%dm, ierr); CHKERRQ(ierr)

    if (allocated(self%bcs)) then
       deallocate(self%bcs)
    end if

    call ISDestroy(self%natural_order, ierr); CHKERRQ(ierr)

  end subroutine mesh_destroy

!------------------------------------------------------------------------

  subroutine mesh_setup_boundaries(self, json, eos, logfile)
    !! Sets up boundary conditions on the mesh.

    use eos_module, only: eos_type
    use boundary_module, only: setup_boundaries
    use logfile_module

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    class(eos_type), intent(in) :: eos
    type(logfile_type), intent(in out), optional :: logfile

    call setup_boundaries(json, eos, self%dm, self%bcs, logfile)

  end subroutine mesh_setup_boundaries

!------------------------------------------------------------------------

  subroutine mesh_set_boundary_values(self, y, rock_vector, eos, &
       y_range_start, rock_range_start)
    !! Sets primary variables (and rock properties) in boundary ghost cells.

    use dm_utils_module, only: global_vec_section, global_section_offset
    use boundary_module, only: open_boundary_label_name
    use eos_module, only: eos_type
    use rock_module, only: rock_type

    class(mesh_type), intent(in) :: self
    Vec, intent(in out) :: y, rock_vector
    class(eos_type), intent(in) :: eos
    PetscInt, intent(in) :: y_range_start, rock_range_start
    ! Locals:
    PetscInt :: ibdy, f, i, num_faces, iface, rock_dof, np, n
    PetscReal, pointer :: y_array(:), rock_array(:)
    PetscReal, pointer :: cell_primary(:), rock1(:), rock2(:)
    PetscSection :: y_section, rock_section
    IS :: bdy_IS
    DMLabel :: bdy_label
    type(rock_type) :: rock
    PetscInt :: y_offset, rock_offsets(2)
    PetscInt, pointer :: bdy_faces(:), cells(:)
    PetscErrorCode :: ierr

    call global_vec_section(y, y_section)
    call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call global_vec_section(rock_vector, rock_section)
    call VecGetArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
    call DMGetLabel(self%dm, open_boundary_label_name, &
         bdy_label, ierr); CHKERRQ(ierr)
    rock_dof = rock%dof()
    np = eos%num_primary_variables

    do ibdy = 1, size(self%bcs, 1)
       call DMGetStratumIS(self%dm, open_boundary_label_name, &
            ibdy, bdy_IS, ierr); CHKERRQ(ierr)
       if (bdy_IS /= 0) then
          call ISGetIndicesF90(bdy_IS, bdy_faces, ierr); CHKERRQ(ierr)
          num_faces = size(bdy_faces)
          do iface = 1, num_faces
             f = bdy_faces(iface)
             call DMPlexGetSupport(self%dm, f, cells, ierr); CHKERRQ(ierr)
             call global_section_offset(y_section, cells(2), &
                  y_range_start, y_offset, ierr); CHKERRQ(ierr)
             do i = 1, 2
                call global_section_offset(rock_section, cells(i), &
                     rock_range_start, rock_offsets(i), ierr)
                CHKERRQ(ierr)
             end do
             ! Set primary variables:
             cell_primary => y_array(y_offset : y_offset + np - 1)
             cell_primary = self%bcs(2: np + 1, ibdy)
             ! Copy rock type data from interior cell to boundary ghost cell:
             n = rock_dof - 1
             rock1 => rock_array(rock_offsets(1) : rock_offsets(1) + n)
             rock2 => rock_array(rock_offsets(2) : rock_offsets(2) + n)
             rock2 = rock1
          end do
       end if
       call ISDestroy(bdy_IS, ierr); CHKERRQ(ierr)
    end do

    call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)

  end subroutine mesh_set_boundary_values

!------------------------------------------------------------------------

  subroutine mesh_order_vector(self, v, v_order)
    !! Reorders vector v to correspond to the natural cell order of the mesh
    !! DM, rather than that of the given v_order index set.

    use cell_module, only: cell_type
    use dm_utils_module, only: global_section_offset, &
         global_vec_section, global_vec_range_start

    class(mesh_type), intent(in) :: self
    Vec, intent(in out) :: v
    IS, intent(in) :: v_order
    ! Locals:
    Vec :: vinitial
    VecScatter :: scatter
    PetscInt :: blocksize
    IS :: v_order_block, self_order_block
    PetscInt, pointer :: indices(:)
    PetscErrorCode :: ierr

    call VecDuplicate(v, vinitial, ierr); CHKERRQ(ierr)
    call VecCopy(v, vinitial, ierr); CHKERRQ(ierr)
    call VecGetBlockSize(v, blocksize, ierr); CHKERRQ(ierr)

    ! Create index sets with the appropriate block size:
    call ISGetIndicesF90(v_order, indices, ierr); CHKERRQ(ierr)
    call ISCreateBlock(mpi%comm, blocksize, size(indices), indices, &
         PETSC_COPY_VALUES, v_order_block, ierr); CHKERRQ(ierr)
    call ISRestoreIndicesF90(v_order, indices, ierr); CHKERRQ(ierr)
    call ISGetIndicesF90(self%natural_order, indices, ierr); CHKERRQ(ierr)
    call ISCreateBlock(mpi%comm, blocksize, size(indices), indices, &
         PETSC_COPY_VALUES, self_order_block, ierr); CHKERRQ(ierr)
    call ISRestoreIndicesF90(self%natural_order, indices, ierr); CHKERRQ(ierr)

    call VecScatterCreate(vinitial, v_order_block, v, self_order_block, &
         scatter, ierr); CHKERRQ(ierr)
    call ISDestroy(v_order_block, ierr); CHKERRQ(ierr)
    call ISDestroy(self_order_block, ierr); CHKERRQ(ierr)

    call VecScatterBegin(scatter, vinitial, v, INSERT_VALUES, &
         SCATTER_FORWARD, ierr); CHKERRQ(ierr)
    call VecScatterEnd(scatter, vinitial, v, INSERT_VALUES, &
         SCATTER_FORWARD, ierr); CHKERRQ(ierr)

    call VecScatterDestroy(scatter, ierr); CHKERRQ(ierr)
    call VecDestroy(vinitial, ierr); CHKERRQ(ierr)

  end subroutine mesh_order_vector

!------------------------------------------------------------------------

end module mesh_module
