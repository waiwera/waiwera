module fluid_module

  implicit none
  private

#include <petsc-finclude/petsc.h90>

  PetscInt, parameter, public :: num_phase_variables = 6
  PetscInt, parameter, public :: num_fluid_variables = 3

  type phase_type
     !! Type for accessing local fluid properties for a particular phase.
     PetscReal, pointer :: density    !! Phase density
     PetscReal, pointer :: viscosity  !! Viscosity
     PetscReal, pointer :: saturation !! Phase saturation
     PetscReal, pointer :: relative_permeability !! Relative permeability
     PetscReal, pointer :: specific_enthalpy !! Specific enthalpy
     PetscReal, pointer :: internal_energy !! Internal energy
     PetscReal, pointer :: mass_fraction(:)
   contains
     private
     procedure, public :: init => phase_init
     procedure, public :: destroy => phase_destroy
     procedure, public :: dof => phase_dof
  end type phase_type

  type fluid_type
     !! Type for accessing local fluid properties.
     private
     PetscReal, pointer, public :: pressure    !! Bulk pressure
     PetscReal, pointer, public :: temperature !! Temperature
     PetscReal, pointer, public :: region      !! Thermodynamic region
     type(phase_type), allocatable, public :: phase(:)
   contains
     private
     procedure, public :: init => fluid_init
     procedure, public :: assign => fluid_assign
     procedure, public :: destroy => fluid_destroy
     procedure, public :: dof => fluid_dof
     procedure, public :: component_density => fluid_component_density
  end type fluid_type

  public :: fluid_type, setup_fluid_vector

contains

!------------------------------------------------------------------------
! Phase procedures
!------------------------------------------------------------------------

  subroutine phase_init(self, num_components)
    !! Initialises a phase object.
    
    class(phase_type), intent(in out) :: self
    PetscInt, intent(in) :: num_components !! Number of fluid components

    allocate(self%mass_fraction(num_components))
    
  end subroutine phase_init

!------------------------------------------------------------------------

  subroutine phase_destroy(self)
    !! Destroys a phase object.

    class(phase_type), intent(in out) :: self

    nullify(self%density)
    nullify(self%viscosity)
    nullify(self%saturation)
    nullify(self%relative_permeability)
    nullify(self%specific_enthalpy)
    nullify(self%internal_energy)
    nullify(self%mass_fraction)

  end subroutine phase_destroy

!------------------------------------------------------------------------

  PetscInt function phase_dof(self)
    !! Returns number of degrees of freedom in the phase object.

    class(phase_type), intent(in) :: self

    phase_dof = num_phase_variables + size(self%mass_fraction)

  end function phase_dof

!------------------------------------------------------------------------
! Fluid procedures
!------------------------------------------------------------------------

  subroutine fluid_init(self, num_components, num_phases)
    !! Initialises a fluid object.
    
    class(fluid_type), intent(in out) :: self
    PetscInt, intent(in) :: num_components !! Number of fluid components
    PetscInt, intent(in) :: num_phases     !! Number of fluid phases
    ! Locals:
    PetscInt :: i

    allocate(self%phase(num_phases))
    do i = 1, num_phases
       call self%phase(i)%init(num_components)
    end do

  end subroutine fluid_init
    
!------------------------------------------------------------------------

  subroutine fluid_assign(self, data, offset)
    !! Assigns pointers in a fluid object to elements in the data array,
    !! starting from the specified offset.

    class(fluid_type), intent(in out) :: self
    PetscReal, target, intent(in) :: data(:)  !! fluid data array
    PetscInt, intent(in) :: offset  !! fluid array offset
    ! Locals:
    PetscInt :: i, p, nc

    self%pressure => data(offset)
    self%temperature => data(offset + 1)
    self%region => data(offset + 2)
    
    i = offset + num_fluid_variables
    do p = 1, size(self%phase)
       self%phase(p)%density => data(i)
       self%phase(p)%viscosity => data(i+1)
       self%phase(p)%saturation => data(i+2)
       self%phase(p)%relative_permeability => data(i+3)
       self%phase(p)%specific_enthalpy => data(i+4)
       self%phase(p)%internal_energy => data(i+5)
       nc = size(self%phase(i)%mass_fraction)
       self%phase(p)%mass_fraction => data(i+6: i+6 + nc-1)
       i = i + self%phase(p)%dof()
    end do

  end subroutine fluid_assign

!------------------------------------------------------------------------

  subroutine fluid_destroy(self)
    !! Destroys a fluid object.
    
    class(fluid_type), intent(in out) :: self
    ! Locals:
    PetscInt :: p

    nullify(self%pressure)
    nullify(self%temperature)
    nullify(self%region)

    do p = 1, size(self%phase)
       call self%phase(p)%destroy()
    end do
    deallocate(self%phase)

  end subroutine fluid_destroy

!------------------------------------------------------------------------

  PetscInt function fluid_dof(self)
    !! Returns number of degrees of freedom in a fluid object.

    class(fluid_type), intent(in) :: self
    ! Locals:
    PetscInt :: p

    fluid_dof = num_fluid_variables
    do p = 1, size(self%phase)
       fluid_dof = fluid_dof + self%phase(p)%dof()
    end do

  end function fluid_dof

!------------------------------------------------------------------------

  PetscReal function fluid_component_density(self, comp) result(d)
    !! Returns fluid density for a given mass component.

    use kinds_module

    class(fluid_type), intent(in) :: self
    PetscInt, intent(in) :: comp !! Index of mass component
    ! Locals:
    PetscInt :: p

    d = 0._dp
    do p = 1, size(self%phase)
       d = d + self%phase(p)%saturation * &
            self%phase(p)%density * self%phase(p)%mass_fraction(comp)
    end do

  end function fluid_component_density

!------------------------------------------------------------------------
! Fluid vector setup routine
!------------------------------------------------------------------------

  subroutine setup_fluid_vector(dm, num_phases, num_components, fluid)
    !! Sets up global vector for fluid properties, with specified
    !! numbers of components and phases.

    use dm_utils_module, only: set_dm_data_layout

    DM, intent(in) :: dm
    PetscInt, intent(in) :: num_phases, num_components
    Vec, intent(out) :: fluid
    ! Locals:
    PetscInt :: num_vars
    PetscInt, allocatable :: num_field_components(:), field_dim(:)
    DM :: dm_fluid
    PetscErrorCode :: ierr

    num_vars = num_fluid_variables + num_phases * &
         (num_phase_variables + num_components)

    allocate(num_field_components(num_vars), field_dim(num_vars))

    ! All fluid variables are scalars defined on cells:
    num_field_components = 1
    field_dim = 3

    call DMClone(dm, dm_fluid, ierr); CHKERRQ(ierr)

    call set_dm_data_layout(dm_fluid, num_field_components, field_dim)

    call DMCreateGlobalVector(dm_fluid, fluid, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(fluid, "fluid", ierr); CHKERRQ(ierr)

    deallocate(num_field_components, field_dim)
    call DMDestroy(dm_fluid, ierr); CHKERRQ(ierr)

  end subroutine setup_fluid_vector

!------------------------------------------------------------------------

end module fluid_module
