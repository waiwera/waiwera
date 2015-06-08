module face_module
  !! Defines type for accessing local quantities defined on a mesh face.

  use cell_module

  implicit none
  private

#include <petsc-finclude/petscdef.h>

  type face_type
     !! Type for accessing local face properties.
     private
     PetscReal, pointer, public :: area !! face area
     PetscReal, pointer, public :: distance(:) !! cell centroid distances on either side of the face
     PetscReal, pointer, public :: normal(:) !! normal vector to face
     PetscReal, pointer, public :: centroid(:) !! centroid of face
     type(cell_type), allocatable, public :: cell(:) !! cells on either side of face
   contains
     private
     procedure, public :: init => face_init
     procedure, public :: assign => face_assign
     procedure, public :: dof => face_dof
     procedure, public :: destroy => face_destroy
  end type face_type

  type petsc_face_type
     !! Type for accessing face geometry parameters calculated by
     !! PETSc DMPlexTSGetGeometryFVM().
     private
     PetscReal, pointer, public :: area_normal(:) !! normal vector multiplied by area
     PetscReal, pointer, public :: centroid(:) !! centroid of face
   contains
     private
     procedure, public :: assign => petsc_face_assign
     procedure, public :: destroy => petsc_face_destroy
  end type petsc_face_type

  public :: face_type, petsc_face_type

contains

!------------------------------------------------------------------------

  subroutine face_init(self, num_components, num_phases)
    !! Initialises a face.

    class(face_type), intent(in out) :: self
    PetscInt, intent(in) :: num_components !! Number of fluid components
    PetscInt, intent(in) :: num_phases     !! Number of fluid phases
    ! Locals:
    PetscInt, parameter :: num_cells = 2
    PetscInt :: i

    allocate(self%cell(num_cells))
    do i = 1, num_cells
       call self%cell(i)%init(num_components, num_phases)
    end do

  end subroutine face_init

!------------------------------------------------------------------------

  subroutine face_assign(self, face_geom_data, face_geom_offset, cell_geom_data, &
       cell_geom_offsets, cell_rock_data, cell_rock_offsets, cell_fluid_data, &
       cell_fluid_offsets)
    !! Assigns pointers in a face to elements of the specified data arrays,
    !! starting from the given offsets.

    class(face_type), intent(in out) :: self
    PetscReal, target, intent(in) :: face_geom_data(:)  !! array with face geometry data
    PetscInt, intent(in)  :: face_geom_offset  !! face geometry array offset for this face
    PetscReal, target, intent(in), optional :: cell_geom_data(:)  !! array with cell geometry data
    PetscInt, intent(in), optional  :: cell_geom_offsets(:)  !! cell geometry array offsets for the face cells
    PetscReal, target, intent(in), optional :: cell_rock_data(:)  !! array with cell rock data
    PetscInt, intent(in), optional  :: cell_rock_offsets(:)  !! cell rock array offsets for the face cells
    PetscReal, target, intent(in), optional :: cell_fluid_data(:)  !! array with cell fluid data
    PetscInt, intent(in), optional  :: cell_fluid_offsets(:)  !! cell fluid array offsets for the face cells
    ! Locals:
    PetscInt :: i
    
    self%area => face_geom_data(face_geom_offset)
    self%distance => face_geom_data(face_geom_offset + 1: face_geom_offset + 2)
    self%normal => face_geom_data(face_geom_offset + 3: face_geom_offset + 5)
    self%centroid => face_geom_data(face_geom_offset + 6: face_geom_offset + 8)

    if ((present(cell_geom_data)).and.(present(cell_geom_offsets))) then

       if ((present(cell_fluid_data)).and.(present(cell_fluid_offsets))) then

          if ((present(cell_rock_data)).and.(present(cell_rock_offsets))) then
             ! Assign geometry, rock and fluid:

             do i = 1, 2
                call self%cell(i)%assign(cell_geom_data, cell_geom_offsets(i), &
                     cell_rock_data, cell_rock_offsets(i), &
                     cell_fluid_data, cell_fluid_offsets(i))
             end do

          else
             ! Assign geometry and fluid:

             do i = 1, 2
                call self%cell(i)%assign(cell_geom_data, cell_geom_offsets(i), &
                     fluid_data = cell_fluid_data, fluid_offset = cell_fluid_offsets(i))
             end do

          end if

       else
       
          if ((present(cell_rock_data)).and.(present(cell_rock_offsets))) then
             ! Assign geometry and rock:

             do i = 1, 2
                call self%cell(i)%assign(cell_geom_data, cell_geom_offsets(i), &
                     cell_rock_data, cell_rock_offsets(i))
             end do

          else
             ! Assign geometry:

             do i = 1, 2
                call self%cell(i)%assign(cell_geom_data, cell_geom_offsets(i))
             end do

          end if

       end if

    end if

  end subroutine face_assign

!------------------------------------------------------------------------

  PetscInt function face_dof(self)
    !! Returns number of degrees of freedom in a face object.

    class(face_type), intent(in) :: self
    ! Locals:
    PetscInt, parameter :: fixed_dof = 9

    face_dof = fixed_dof

  end function face_dof

!------------------------------------------------------------------------

  subroutine face_destroy(self)
    !! Destroys a face (nullifies all pointer components).

    class(face_type), intent(in out) :: self

    nullify(self%area)
    nullify(self%distance)
    nullify(self%normal)
    nullify(self%centroid)
    if (allocated(self%cell)) then
       deallocate(self%cell)
    end if

  end subroutine face_destroy

!------------------------------------------------------------------------
! petsc_face_type routines
!------------------------------------------------------------------------

  subroutine petsc_face_assign(self, face_geom_data, face_geom_offset)
    !! Assigns pointers in a petsc_face to elements of the specified data
    !! array, starting from the given offset.

    class(petsc_face_type), intent(in out) :: self
    PetscReal, target, intent(in) :: face_geom_data(:)  !! array with face geometry data
    PetscInt, intent(in)  :: face_geom_offset  !! face geometry array offset for this face

    self%area_normal => face_geom_data(face_geom_offset: face_geom_offset + 2)
    self%centroid => face_geom_data(face_geom_offset + 3: face_geom_offset + 5)

  end subroutine petsc_face_assign

!------------------------------------------------------------------------

  subroutine petsc_face_destroy(self)
    !! Destroys a petsc_face (nullifies all pointer components).

    class(petsc_face_type), intent(in out) :: self

    nullify(self%area_normal)
    nullify(self%centroid)

  end subroutine petsc_face_destroy

!------------------------------------------------------------------------

end module face_module
