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

module rock_module
  !! Defines type for accessing local rock properties on cells and faces.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use relative_permeability_module
  use capillary_pressure_module

  implicit none
  private

  type rock_type
     !! Local rock properties.
     private
     PetscReal, pointer, contiguous, public :: permeability(:)   !! Permeability
     PetscReal, pointer, public :: wet_conductivity, dry_conductivity !! Heat conductivities
     PetscReal, pointer, public :: porosity          !! Porosity
     PetscReal, pointer, public :: density           !! Grain density
     PetscReal, pointer, public :: specific_heat     !! Specific heat
     class(relative_permeability_type), pointer, &
          public :: relative_permeability !! Relative permeability functions
     class(capillary_pressure_type), pointer, &
          public :: capillary_pressure !! Capillary pressure function
     PetscInt, public :: dof !! Number of degrees of freedom
   contains
     private
     procedure, public :: init => rock_init
     procedure, public :: assign => rock_assign
     procedure, public :: assign_relative_permeability => &
          rock_assign_relative_permeability
     procedure, public :: assign_capillary_pressure => &
          rock_assign_capillary_pressure
     procedure, public :: destroy => rock_destroy
     procedure, public :: energy => rock_energy
  end type rock_type

  PetscInt, parameter :: num_rock_variables = 6
  PetscInt, parameter :: max_rock_variable_name_length = 32
  character(max_rock_variable_name_length), parameter, public :: &
       rock_variable_names(num_rock_variables) = &
       [character(max_rock_variable_name_length):: &
       "permeability", "wet_conductivity", "dry_conductivity", "porosity", &
       "density", "specific_heat"]
  PetscInt, parameter, public :: &
       rock_variable_num_components(num_rock_variables) = &
       [3, 1, 1, 1, 1, 1]
  PetscInt, parameter, public :: max_rockname_length = 24

  ! Default rock properties:
  PetscReal, parameter, public :: default_permeability_scalar = 1.e-13_dp
  PetscReal, parameter, public :: default_permeability(3) = [ &
       default_permeability_scalar, default_permeability_scalar, &
       default_permeability_scalar]
  PetscReal, parameter, public :: default_porosity = 0.1_dp
  PetscReal, parameter, public :: default_density = 2200.0_dp
  PetscReal, parameter, public :: default_specific_heat = 1000._dp
  PetscReal, parameter, public :: default_heat_conductivity = 2.5_dp

  character(len = 9), public :: rock_type_label_name = "rock_type"

  public :: rock_type, num_rock_variables

contains

!------------------------------------------------------------------------

  subroutine rock_init(self)
    !! Initialises rock object.

    class(rock_type), intent(in out) :: self

    self%dof = sum(rock_variable_num_components)

  end subroutine rock_init

!------------------------------------------------------------------------

  subroutine rock_assign(self, data, offset)
    !! Assigns pointers in a rock object to elements of the specified
    !! data array, starting from the given offset.

    class(rock_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! rock data array
    PetscInt, intent(in) :: offset !! rock array offset

    self%permeability => data(offset: offset + 2)
    self%wet_conductivity => data(offset + 3)
    self%dry_conductivity => data(offset + 4)
    self%porosity => data(offset + 5)
    self%density => data(offset + 6)
    self%specific_heat => data(offset + 7)

  end subroutine rock_assign
    
!------------------------------------------------------------------------

  subroutine rock_assign_relative_permeability(self, relative_permeability)
    !! Assigns relative permeability pointer for a rock object.

    class(rock_type), intent(in out) :: self
    class(relative_permeability_type), intent(in), &
         target :: relative_permeability

    self%relative_permeability => relative_permeability

  end subroutine rock_assign_relative_permeability

!------------------------------------------------------------------------

  subroutine rock_assign_capillary_pressure(self, capillary_pressure)
    !! Assigns capillary pressure pointer for a rock object.

    class(rock_type), intent(in out) :: self
    class(capillary_pressure_type), intent(in), &
         target :: capillary_pressure

    self%capillary_pressure => capillary_pressure

  end subroutine rock_assign_capillary_pressure

!------------------------------------------------------------------------

  PetscReal function rock_energy(self, temperature)
    !! Returns rock energy density at a given temperature.

    class(rock_type), intent(in) :: self
    PetscReal, intent(in) :: temperature !! Temperature (deg C)

    rock_energy = self%density * self%specific_heat * temperature

  end function rock_energy

!------------------------------------------------------------------------

  subroutine rock_destroy(self)
    !! Destroys a rock object (nullifies all pointer components).

    class(rock_type), intent(in out) :: self

    self%permeability => null()
    self%wet_conductivity => null()
    self%dry_conductivity => null()
    self%porosity => null()
    self%density => null()
    self%specific_heat => null()
    self%relative_permeability => null()
    self%capillary_pressure => null()

  end subroutine rock_destroy

!------------------------------------------------------------------------

end module rock_module
