module fluid_module

  implicit none
  private

#include <petsc-finclude/petscdef.h>

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
  end type fluid_type

  public :: fluid_type

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
    PetscInt :: i, p, nc, np

    self%pressure => data(offset)
    self%temperature => data(offset + 1)
    self%region => data(offset + 2)
    
    i = offset + num_fluid_variables
    np = size(self%phase)
    do p = 1, np
       self%phase(p)%density => data(i)
       self%phase(p)%viscosity => data(i+1)
       self%phase(p)%saturation => data(i+2)
       self%phase(p)%relative_permeability => data(i+3)
       self%phase(p)%specific_enthalpy => data(i+4)
       nc = size(self%phase(i)%mass_fraction)
       self%phase(p)%mass_fraction => data(i+5: i + 5 + nc-1)
       i = i + self%phase(p)%dof()
    end do

  end subroutine fluid_assign

!------------------------------------------------------------------------

  subroutine fluid_destroy(self)
    !! Destroys a fluid object.
    
    class(fluid_type), intent(in out) :: self
    ! Locals:
    PetscInt :: i

    nullify(self%pressure)
    nullify(self%temperature)
    nullify(self%region)

    do i = 1, size(self%phase)
       call self%phase(i)%destroy()
    end do
    deallocate(self%phase)

  end subroutine fluid_destroy

!------------------------------------------------------------------------

  PetscInt function fluid_dof(self)
    !! Returns number of degrees of freedom in a fluid object.

    class(fluid_type), intent(in) :: self
    ! Locals:
    PetscInt :: i

    fluid_dof = num_fluid_variables
    do i = 1, size(self%phase)
       fluid_dof = fluid_dof + self%phase(i)%dof()
    end do

  end function fluid_dof

!------------------------------------------------------------------------

end module fluid_module
