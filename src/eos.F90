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

module eos_module
  !! Equations of state for different fluid mass components. The
  !! behaviour of each combination of fluid components is governed by
  !! an EOS object.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use fson
  use thermodynamics_module
  use interpolation_module
  use hdf5io_module, only: max_field_name_length

  implicit none
  private

  PetscInt, parameter, public :: max_eos_name_length = 8
  PetscInt, parameter, public :: max_eos_description_length = 80
  PetscInt, parameter, public :: max_primary_variable_name_length = 32
  PetscInt, parameter, public :: max_phase_name_length = 13
  PetscInt, parameter, public :: max_component_name_length = 8
  character(max_component_name_length), parameter, public :: &
       energy_component_name  = 'energy'

  type, public, abstract :: eos_type
     !! Abstract type for equation of state (EOS) objects.
     private
     PetscReal, allocatable, public :: primary_scale(:, :)
     character(max_eos_name_length), public :: name !! EOS name
     character(max_eos_description_length), public :: description !! EOS description
     character(max_primary_variable_name_length), allocatable, public :: primary_variable_names(:) !! Names of primary variables
     character(max_phase_name_length), allocatable, public :: phase_names(:) !! Names of fluid phases
     character(max_component_name_length), allocatable, public :: component_names(:) !! Names of mass components
     PetscInt, public :: num_primary_variables !! Number of primary variables
     PetscInt, public :: num_phases !! Number of possible phases
     PetscInt, public :: num_components !! Number of mass components
     PetscReal, allocatable, public :: default_primary(:) !! Default primary variable values
     PetscInt, public :: default_region !! Default thermodynamic region
     character(max_field_name_length), allocatable, public :: required_output_fluid_fields(:) !! Fluid output fields that are required in output (for restarting)
     character(max_field_name_length), allocatable, public :: default_output_fluid_fields(:) !! Default fluid fields written to output
     class(thermodynamics_type), pointer, public :: thermo !! Thermodynamic formulation
     PetscBool, public :: isothermal = PETSC_FALSE !! Whether the EOS is restricted to isothermal fluid conditions
     procedure(eos_scale_procedure), pointer, public :: scale => eos_scale_default
     procedure(eos_unscale_procedure), pointer, public :: unscale => eos_unscale_default
   contains
     private
     procedure(eos_init_procedure), public, deferred :: init
     procedure, public :: destroy => eos_destroy
     procedure(eos_transition_procedure), public, deferred :: transition
     procedure(eos_bulk_properties_procedure), public, deferred :: bulk_properties
     procedure, public :: phase_composition => eos_phase_composition
     procedure(eos_phase_properties_procedure), public, deferred :: phase_properties
     procedure(eos_primary_variables_procedure), public, deferred :: primary_variables
     procedure(eos_check_primary_variables_procedure), public, deferred :: check_primary_variables
     procedure, public :: conductivity => eos_conductivity
     procedure, public :: component_index => eos_component_index
  end type eos_type

  type, public, extends(interpolation_table_type) :: primary_variable_interpolator_type
     !! Interpolator for primary variable arrays, including
     !! thermodynamics object for when interpolator is used as a context
     !! for root finding.
     private
     class(thermodynamics_type), pointer, public :: thermo
   contains
     procedure, public :: destroy =>  primary_variable_interpolator_destroy
  end type primary_variable_interpolator_type

  abstract interface

     subroutine eos_init_procedure(self, json, thermo, logfile)
       !! Initialises the EOS object
       use logfile_module
       import :: eos_type, thermodynamics_type, fson_value
       class(eos_type), intent(in out) :: self
       type(fson_value), pointer, intent(in) :: json
       class(thermodynamics_type), intent(in), target :: thermo
       type(logfile_type), intent(in out), optional :: logfile
     end subroutine eos_init_procedure

     subroutine eos_transition_procedure(self, old_primary, primary, &
          old_fluid, fluid, transition, err)
       !! Check primary variables for a cell and make thermodynamic
       !! region transitions if needed.
       use fluid_module, only: fluid_type
       import :: eos_type
       class(eos_type), intent(in out) :: self
       PetscReal, intent(in) :: old_primary(self%num_primary_variables)
       PetscReal, intent(in out) :: primary(self%num_primary_variables)
       type(fluid_type), intent(in) :: old_fluid
       type(fluid_type), intent(in out) :: fluid
       PetscBool, intent(out) :: transition
       PetscErrorCode, intent(out) :: err
     end subroutine eos_transition_procedure

     subroutine eos_bulk_properties_procedure(self, primary, fluid, err)
       !! Calculate bulk fluid properties from primary variables.
       use fluid_module, only: fluid_type
       import :: eos_type
       class(eos_type), intent(in out) :: self
       PetscReal, intent(in) :: primary(self%num_primary_variables)
       type(fluid_type), intent(in out) :: fluid
       PetscErrorCode, intent(out) :: err
     end subroutine eos_bulk_properties_procedure

     subroutine eos_phase_properties_procedure(self, primary, rock, fluid, err)
       !! Calculate phase fluid properties from primary variables.
       use rock_module, only: rock_type
       use fluid_module, only: fluid_type
       import :: eos_type
       class(eos_type), intent(in out) :: self
       PetscReal, intent(in) :: primary(self%num_primary_variables)
       type(rock_type), intent(in out) :: rock
       type(fluid_type), intent(in out) :: fluid
       PetscErrorCode, intent(out) :: err
     end subroutine eos_phase_properties_procedure

     subroutine eos_primary_variables_procedure(self, fluid, primary)
       !! Determine primary variables from fluid properties.
       use fluid_module, only: fluid_type
       import :: eos_type
       class(eos_type), intent(in) :: self
       type(fluid_type), intent(in) :: fluid
       PetscReal, intent(out) :: primary(self%num_primary_variables)
     end subroutine eos_primary_variables_procedure

     subroutine eos_check_primary_variables_procedure(self, &
          fluid, primary, changed, err)
       !! Check if primary variables are in acceptable bounds, and
       !! return error code accordingly. Also set changed to true if
       !! primary variables have been changed during the check.
       use fluid_module, only: fluid_type
       import :: eos_type
       class(eos_type), intent(in) :: self
       type(fluid_type), intent(in) :: fluid
       PetscReal, intent(in out) :: primary(self%num_primary_variables)
       PetscBool, intent(out) :: changed
       PetscErrorCode, intent(out) :: err
     end subroutine eos_check_primary_variables_procedure

     function eos_scale_procedure(self, primary, region) result(scaled_primary)
       !! Non-dimensionalise primary variables by scaling.
       import :: eos_type
       class(eos_type), intent(in) :: self
       PetscReal, intent(in) :: primary(self%num_primary_variables)
       PetscInt, intent(in) :: region
       PetscReal :: scaled_primary(self%num_primary_variables)
     end function eos_scale_procedure

     function eos_unscale_procedure(self, scaled_primary, region) result(primary)
       !! Re-dimensionalise scaled primary variables.
       import :: eos_type
       class(eos_type), intent(in) :: self
       PetscReal, intent(in) :: scaled_primary(self%num_primary_variables)
       PetscInt, intent(in) :: region
       PetscReal :: primary(self%num_primary_variables)
     end function eos_unscale_procedure

  end interface

contains

!------------------------------------------------------------------------

  function eos_scale_default(self, primary, region) result(scaled_primary)
    !! Default function for non-dimensionalising primary variables by scaling.

    class(eos_type), intent(in) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables)
    PetscInt, intent(in) :: region
    PetscReal :: scaled_primary(self%num_primary_variables)

    scaled_primary = primary / self%primary_scale(:, region)

  end function eos_scale_default

!------------------------------------------------------------------------

  function eos_unscale_default(self, scaled_primary, region) result(primary)
    !! Default function for re-dimensionalising scaled primary variables.

    class(eos_type), intent(in) :: self
    PetscReal, intent(in) :: scaled_primary(self%num_primary_variables)
    PetscInt, intent(in) :: region
    PetscReal :: primary(self%num_primary_variables)

    primary = scaled_primary * self%primary_scale(:, region)

  end function eos_unscale_default

!------------------------------------------------------------------------

  subroutine eos_phase_composition(self, fluid, err)
    !! Determines fluid phase composition from bulk properties and
    !! thermodynamic region.

    use fluid_module, only: fluid_type

    class(eos_type), intent(in out) :: self
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: region, phases

    region = nint(fluid%region)
    phases = self%thermo%phase_composition(region, fluid%pressure, &
         fluid%temperature)
    if (phases > 0) then
       fluid%phase_composition = dble(phases)
       err = 0
    else
       err = 1
    end if

  end subroutine eos_phase_composition

!------------------------------------------------------------------------

  PetscReal function eos_conductivity(self, rock, fluid) result(cond)
    !! Returns effective rock heat conductivity for given fluid properties.
    !! This uses a square-root dependence on liquid saturation.

    use rock_module, only: rock_type
    use fluid_module, only: fluid_type

    class(eos_type), intent(in) :: self
    type(rock_type), intent(in) :: rock !! Rock object
    type(fluid_type), intent(in) :: fluid !! Fluid object
    ! Locals:
    PetscReal :: sl

    sl = fluid%phase(1)%saturation
    cond = rock%dry_conductivity + sqrt(sl) * &
         (rock%wet_conductivity - rock%dry_conductivity)

  end function eos_conductivity

!------------------------------------------------------------------------

  PetscInt function eos_component_index(self, component_name) result (index)
    !! Returns index of specified component name (or -1 if no such
    !! component exists).

    use utils_module, only: str_to_lower, str_array_index

    class(eos_type), intent(in) :: self
    character(len = *), intent(in) :: component_name
    ! Locals:
    character(len = len(component_name)) :: lowercase_name

    lowercase_name = str_to_lower(component_name)

    if ((trim(lowercase_name) == trim(energy_component_name)) .and. &
         (.not. self%isothermal)) then
       index = self%num_primary_variables
    else
       index = str_array_index(lowercase_name, self%component_names)
    end if

  end function eos_component_index

!------------------------------------------------------------------------

  subroutine eos_destroy(self)
    !! Destroys EOS.

    class(eos_type), intent(in out) :: self

    deallocate(self%primary_variable_names)
    deallocate(self%phase_names, self%component_names)
    deallocate(self%default_primary)
    self%thermo => null()

  end subroutine eos_destroy

!------------------------------------------------------------------------
! Primary variable interpolator
!------------------------------------------------------------------------

  subroutine primary_variable_interpolator_destroy(self)
    !! Destroys primary variable interpolator.

    class(primary_variable_interpolator_type), intent(in out) :: self

    call self%interpolation_table_type%destroy()
    self%thermo => null()

  end subroutine primary_variable_interpolator_destroy

!------------------------------------------------------------------------

end module eos_module
