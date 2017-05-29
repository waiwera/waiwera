module rock_mechanics_module

#include <petsc/finclude/petsc.h90>

  use petsc
  use kinds_module

  implicit none

  private

  PetscInt, parameter :: max_primary_variable_name_length = 12
  PetscInt, parameter :: max_gauss_order = 2
  PetscInt, parameter :: max_spatial_dimension = 3
  
  type, public :: rock_mechanics_type
     private
     DM :: dm_rock
     PetscInt :: start_cell, end_cell, end_interior_cell
     character(max_primary_variable_name_length), allocatable :: primary_variable_names(:)
     
   contains
     private
     procedure :: mesh_init => rock_mechanics_mesh_init
     procedure, public :: init => rock_mechanics_init
     procedure, public :: run => rock_mechanics_run

  end type rock_mechanics_type

contains
  !------------------------------------------------------------------------
  ! public routines
  !------------------------------------------------------------------------
  subroutine rock_mechanics_init(self,flow_dm,logfile)
    use logfile_module
    !use rock_parameters_module, only : setup_rocktype_labels
    
    class(rock_mechanics_type), intent(in out) :: self
    DM, intent(in) :: flow_dm
    type(logfile_type), intent(in out), optional :: logfile
    
    print *,"in rock mechanics init"
    call self%mesh_init(flow_dm,logfile)
    
  end subroutine rock_mechanics_init
  !------------------------------------------------------------------------
  subroutine rock_mechanics_run(self,fluid,rock,fluid_vector,fluid_range_start,rock_vector,rock_range_start)
    use fluid_module, only: fluid_type
    use rock_module, only: rock_type
    use dm_utils_module, only: global_section_offset, global_vec_section

    class(rock_mechanics_type), intent(in out) :: self
    type(fluid_type), intent (in out) :: fluid
    type(rock_type), intent(in out) :: rock
    Vec, intent(in) :: fluid_vector
    PetscInt, intent(in) :: fluid_range_start
    Vec, intent(in) :: rock_vector
    PetscInt, intent(in) :: rock_range_start

    ! **************  locals **************
    PetscInt :: e, i, ndof
    
    ! stuff for material parameters
    PetscInt :: fluid_offset, rock_offset
    PetscSection :: fluid_section, rock_section
    PetscReal, pointer, contiguous :: fluid_array(:)
    PetscReal, pointer, contiguous :: rock_array(:)
    
    ! stuff for fe quadrature 
    PetscFE :: fem
    PetscQuadrature :: quad  
    PetscInt :: numele
    PetscInt :: quad_dim, num_quad
    PetscReal, target, dimension(max_spatial_dimension*(max_gauss_order**max_spatial_dimension)) :: q_points
    PetscReal, target, dimension(max_gauss_order**max_spatial_dimension) :: q_weights
    PetscReal, pointer :: pq_points(:)
    PetscReal, pointer :: pq_weights(:)

    ! stuff for integration
    Vec :: global_coordinates
    PetscSection :: coord_section
    PetscReal, target, dimension(max_spatial_dimension*(max_gauss_order**max_spatial_dimension)) :: local_coordinates
    PetscReal, pointer :: pcoords(:)
    
    PetscReal, target, dimension(100) :: v0
    PetscReal, target, dimension(100) :: J
    PetscReal, target, dimension(100) :: invJ
    PetscReal, target :: detJ
    PetscReal, pointer :: pv0(:),pJ(:),pinvJ(:)
    PetscReal, pointer :: pdetJ
   
    PetscErrorCode :: ierr

    pq_points => q_points
    pq_weights => q_weights
    pcoords => local_coordinates
    pv0 => v0
    pJ => J
    pinvJ => invJ
    pdetJ => detJ
    
    print *,"in rock mechanics run"

    ! setting up things for temperature/pressure retrieval
    call VecGetArrayReadF90(fluid_vector,fluid_array,ierr); CHKERRQ(ierr)
    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayReadF90(rock_vector,rock_array,ierr); CHKERRQ(ierr)
    call global_vec_section(rock_vector, rock_section)

    ! numerical quadrature setup stuff
    numele = self%end_interior_cell-self%start_cell
    call DMGetDimension(self%dm_rock, ndof, ierr); CHKERRQ(ierr)
    call PetscFECreateDefault(self%dm_rock,ndof,numele,PETSC_FALSE,"fe",1,fem,ierr); CHKERRQ(ierr)
    call PetscFEGetQuadrature(fem,quad,ierr); CHKERRQ(ierr)
    call PetscQuadratureGetData(quad,quad_dim,num_quad,pq_points,pq_weights,ierr)
    call PetscFEGetDefaultTabulation(fem,pv0,pJ,pinvJ,ierr); CHKERRQ(ierr)
    
    CHKERRQ(ierr)
    print *,'quad weights and points:'
    do i=1,num_quad
       print *,pq_points(3*i-2:3*i)
    end do

    ! set up stuff for element coordinates   
    call DMGetCoordinatesLocal(self%dm_rock,global_coordinates,ierr); CHKERRQ(ierr)
    call DMGetCoordinateSection(self%dm_rock,coord_section,ierr); CHKERRQ(ierr)

    ! loop over elements for integration
    do e=self%start_cell,self%end_interior_cell-1     
       ! assign fluid properties to get pressure and temperature
       call global_section_offset(fluid_section, e, fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)
       call global_section_offset(rock_section, e, rock_range_start, rock_offset, ierr); CHKERRQ(ierr)
       call fluid%assign(fluid_array, fluid_offset)
       call rock%assign(rock_array, rock_offset)

       call DMPlexVecGetClosure(self%dm_rock,coord_section,global_coordinates,e,pcoords,ierr)
       CHKERRQ(ierr)
       !print *," pres:  ",fluid%pressure
       !print *," temp:  ",fluid%temperature
       !print *," young: ",rock%youngs_modulus
       !print *," poisson: ",rock%poissons_ratio
       !print *," density: ",rock%density

       !call DMPlexComputeCellGeometryFEM(self%dm_rock,e,fem,pv0,pJ,pinvJ,pdetJ,ierr); CHKERRQ(ierr)
       
       print *,"************************"
       print *,'global coords:'
       do i=1,8
          print *,pcoords(3*i-2:3*i)
       end do
       print *,"************************"
       !print *,pv0(1),pv0(2)
       print *,"************************"    


       !call DMPlexComputeCellGeometryFEM(self%dm_rock,e,fem,pv0,pJ,pinvJ,pdetJ,ierr); CHKERRQ(ierr)
       !print *,ele%pcoords
    end do

    call VecRestoreArrayReadF90(fluid_vector,fluid_array,ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(rock_vector,rock_array,ierr); CHKERRQ(ierr)
    call VecDestroy(global_coordinates,ierr); CHKERRQ(ierr)
    call PetscFEDestroy(fem,ierr); CHKERRQ(ierr)
    nullify(pq_points)
    nullify(pq_weights)
  end subroutine rock_mechanics_run
  
  !------------------------------------------------------------------------
  ! private routines
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  
   subroutine rock_mechanics_mesh_init(self,flow_dm,logfile)
    use logfile_module
    
    class(rock_mechanics_type), intent (in out) :: self
    DM, intent(in) :: flow_dm
    type(logfile_type), intent(in out), optional :: logfile
    
    ! locals
    PetscErrorCode :: ierr
    
    ! Clone flow DM:
    call DMClone(flow_dm, self%dm_rock, ierr); CHKERRQ(ierr)   
    self%primary_variable_names = ["displacement"]
   
    call DMPlexGetHeightStratum(self%dm_rock, 0, self%start_cell, self%end_cell,ierr)
    CHKERRQ(ierr)
    
    call DMPlexGetHybridBounds(self%dm_rock, self%end_interior_cell, &
         PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
         ierr)
    CHKERRQ(ierr)

  end subroutine rock_mechanics_mesh_init
  !------------------------------------------------------------------------
end module rock_mechanics_module
