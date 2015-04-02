module mesh_module
  !! Module for mesh type.

  use kinds_module
  use mpi_module
  use fson

  implicit none

  private

#include <petsc-finclude/petsc.h90>

  integer, parameter, public :: max_mesh_filename_length = 200
  character(len = 16), public :: open_boundary_label_name = "open boundary"

  type, public :: mesh_type
     !! Mesh type.
     private
     character(max_mesh_filename_length), public :: filename
     DM, public :: dm
     Vec, public :: cell_geom, face_geom
     PetscInt, public :: start_cell, end_cell, end_interior_cell
     PetscInt, public :: num_cells, num_interior_cells
   contains
     procedure :: distribute => mesh_distribute
     procedure :: construct_ghost_cells => mesh_construct_ghost_cells
     procedure :: setup_discretization => mesh_setup_discretization
     procedure :: setup_geometry => mesh_setup_geometry
     procedure :: label_boundaries => mesh_label_boundaries
     procedure :: get_bounds => mesh_get_bounds
     procedure, public :: init => mesh_init
     procedure, public :: destroy => mesh_destroy
  end type mesh_type

contains

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

  subroutine mesh_setup_discretization(self, dof)
    !! Sets up finite volume discretization for the mesh.

    class(mesh_type), intent(in out) :: self
    PetscInt, intent(in) :: dof !! Degrees of freedom
    ! Locals:
    PetscInt :: dim
    PetscFV :: fvm
    PetscDS :: ds
    PetscErrorCode :: ierr

    call PetscFVCreate(mpi%comm, fvm, ierr); CHKERRQ(ierr)
    call PetscFVSetNumComponents(fvm, dof, ierr); CHKERRQ(ierr)
    call DMGetDimension(self%dm, dim, ierr); CHKERRQ(ierr)
    call PetscFVSetSpatialDimension(fvm, dim, ierr); CHKERRQ(ierr)

    call DMGetDS(self%dm, ds, ierr); CHKERRQ(ierr)
    call PetscDSAddDiscretization(ds, fvm, ierr); CHKERRQ(ierr)

  end subroutine mesh_setup_discretization

!------------------------------------------------------------------------

  subroutine mesh_setup_geometry(self)
    !! Sets up global vectors containing geometry data (e.g. cell volumes,
    !! cell centroids, face areas, face-to-centroid distances) for the mesh.

    use face_module

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    Vec :: petsc_face_geom, local_face_geom
    PetscInt :: dim
    DM :: dm_cell, dm_face, petsc_dm_face
    PetscSection :: face_section, petsc_face_section, cell_section
    PetscInt :: fstart, fend, f, face_dof, ghost, i
    PetscInt :: face_offset, petsc_face_offset
    PetscInt :: cell_offset(2)
    type(face_type) :: face
    type(petsc_face_type) :: petsc_face
    PetscScalar, pointer :: face_geom_array(:), petsc_face_geom_array(:)
    PetscScalar, pointer :: cell_geom_array(:)
    DMLabel :: ghost_label
    PetscInt, pointer :: cells(:)

    ! First call PETSc geometry routine- we use the cell geometry vector but need to 
    ! create our own face geometry vector, containing additional parameters:
    call DMPlexTSGetGeometryFVM(self%dm, petsc_face_geom, self%cell_geom, &
         PETSC_NULL_REAL, ierr); CHKERRQ(ierr)

    ! Set up face geometry vector:
    call DMGetDimension(self%dm, dim, ierr); CHKERRQ(ierr)
    call DMClone(self%dm, dm_face, ierr); CHKERRQ(ierr)
    call PetscSectionCreate(mpi%comm, face_section, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(self%dm, 1, fStart, fEnd, ierr); CHKERRQ(ierr)
    call PetscSectionSetChart(face_section, fStart, fEnd, ierr); CHKERRQ(ierr)
    face_dof = face%dof()
    do f = fstart, fend - 1
       call PetscSectionSetDof(face_section, f, face_dof, ierr); CHKERRQ(ierr)
    end do
    call PetscSectionSetUp(face_section, ierr); CHKERRQ(ierr)
    call DMSetDefaultSection(dm_face, face_section, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(dm_face, self%face_geom, ierr); CHKERRQ(ierr)

    call VecGetDM(petsc_face_geom, petsc_dm_face, ierr); CHKERRQ(ierr)
    call DMGetDefaultSection(petsc_dm_face, petsc_face_section, ierr); CHKERRQ(ierr)

    call DMGetLocalVector(dm_face, local_face_geom, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(local_face_geom, face_geom_array, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(petsc_face_geom, petsc_face_geom_array, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(self%cell_geom, cell_geom_array, ierr); CHKERRQ(ierr)

    call VecGetDM(self%cell_geom, dm_cell, ierr); CHKERRQ(ierr)
    call DMGetDefaultSection(dm_cell, cell_section, ierr); CHKERRQ(ierr)

    call DMPlexGetLabel(self%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)

    do f = fstart, fend - 1

       call DMLabelGetValue(ghost_label, f, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then

          call PetscSectionGetOffset(face_section, f, face_offset, ierr); CHKERRQ(ierr)
          call PetscSectionGetOffset(petsc_face_section, f, petsc_face_offset, ierr)
          CHKERRQ(ierr)
          face_offset = face_offset + 1; petsc_face_offset = petsc_face_offset + 1

          call DMPlexGetSupport(self%dm, f, cells, ierr); CHKERRQ(ierr)
          do i = 1, 2
             call PetscSectionGetOffset(cell_section, cells(i), cell_offset(i), ierr)
             CHKERRQ(ierr)
          end do
          cell_offset = cell_offset + 1

          call petsc_face%assign(petsc_face_geom_array, petsc_face_offset)
          call face%assign(face_geom_array, face_offset, cell_geom_array, cell_offset)

          face%centroid = petsc_face%centroid
          face%area = norm2(petsc_face%area_normal)
          face%normal = petsc_face%area_normal / face%area
          do i = 1, 2
             face%distance(i) = norm2(face%centroid - face%cell(i)%centroid)
          end do

       end if

    end do

    call VecRestoreArrayF90(self%cell_geom, cell_geom_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(local_face_geom, face_geom_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(petsc_face_geom, petsc_face_geom_array, ierr)
    CHKERRQ(ierr)
    call DMLocalToGlobalBegin(dm_face, local_face_geom, INSERT_VALUES, &
         self%face_geom, ierr); CHKERRQ(ierr)
    call DMLocalToGlobalEnd(dm_face, local_face_geom, INSERT_VALUES, &
         self%face_geom, ierr); CHKERRQ(ierr)
    call DMRestoreLocalVector(dm_face, local_face_geom, ierr); CHKERRQ(ierr)

    call PetscSectionDestroy(face_section, ierr); CHKERRQ(ierr)

  end subroutine mesh_setup_geometry

!------------------------------------------------------------------------

  subroutine mesh_label_boundaries(self)
    !! Labels boundary faces as appropriate for assigning boundary
    !! conditions.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    PetscBool :: has_label

    call DMPlexHasLabel(self%dm, open_boundary_label_name, has_label, &
         ierr); CHKERRQ(ierr)
    if (.not.(has_label)) then
       call DMPlexCreateLabel(self%dm, open_boundary_label_name, &
            ierr); CHKERRQ(ierr)
       ! could read boundary faces from input here if needed- i.e. if labels
       ! not present in mesh file
    end if

  end subroutine mesh_label_boundaries

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

    self%num_cells = self%end_cell - self%start_cell
    self%num_interior_cells = self%end_interior_cell - self%start_cell

  end subroutine mesh_get_bounds

!------------------------------------------------------------------------

  subroutine mesh_init(self, json, dof)
    !! Initializes mesh, reading filename from JSON input file.
    !! If the filename is not present, an error is raised.
    !! Otherwise, the PETSc DM is read in.

    class(mesh_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON file pointer
    PetscInt, intent(in) :: dof !! Degrees of freedom on cells
    ! Locals:
    PetscErrorCode :: ierr
    type(fson_value), pointer :: mesh

    if (mpi%rank == mpi%input_rank) then
       call fson_get(json, "mesh", mesh)
       if (associated(mesh)) then
          call fson_get(mesh, ".", self%filename)
       else
          write (*,'(a,a)') "Stopping: mesh not found in input."
          stop
       end if
    end if

    call MPI_bcast(self%filename, max_mesh_filename_length, MPI_CHARACTER, &
         mpi%input_rank, mpi%comm, ierr)

    ! Read in DM:
    call DMPlexCreateFromFile(mpi%comm, self%filename, PETSC_TRUE, self%dm, ierr)
    CHKERRQ(ierr)

    ! Set up adjacency for finite volume mesh:
    call DMPlexSetAdjacencyUseCone(self%dm, PETSC_TRUE, ierr); CHKERRQ(ierr)
    call DMPlexSetAdjacencyUseClosure(self%dm, PETSC_FALSE, ierr); CHKERRQ(ierr)

    call self%distribute()
    call self%label_boundaries()
    call self%construct_ghost_cells()
    call self%setup_discretization(dof)
    call self%setup_geometry()
    call self%get_bounds()

  end subroutine mesh_init

!------------------------------------------------------------------------

  subroutine mesh_destroy(self)
    !! Destroys mesh object.

    class(mesh_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    
    call VecDestroy(self%face_geom, ierr); CHKERRQ(ierr)
    call DMDestroy(self%dm, ierr); CHKERRQ(ierr)

  end subroutine mesh_destroy

!------------------------------------------------------------------------

end module mesh_module
