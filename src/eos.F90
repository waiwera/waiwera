module eos_module
  !! Equations of state for different fluid mass components. The
  !! behaviour of each combination of fluid components is governed by
  !! an EOS object.

  use kinds_module
  use fson
  use thermodynamics_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>
#include <petsc/finclude/petscsys.h>

  PetscInt, parameter, public :: max_eos_name_length = 8
  PetscInt, parameter, public :: max_eos_description_length = 80
  PetscInt, parameter, public :: max_primary_variable_name_length = 32
  PetscInt, parameter, public :: max_phase_name_length = 13
  PetscInt, parameter, public :: max_component_name_length = 8

  type, public, abstract :: eos_type
     !! Abstract type for equation of state (EOS) objects.
     private
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
     class(thermodynamics_type), pointer, public :: thermo !! Thermodynamic formulation
     PetscBool, public :: isothermal = PETSC_FALSE !! Whether the EOS is restricted to isothermal fluid conditions
   contains
     private
     procedure(eos_init_procedure), public, deferred :: init
     procedure(eos_destroy_procedure), public, deferred :: destroy
     procedure(eos_transition_procedure), public, deferred :: transition
     procedure(eos_bulk_properties_procedure), public, deferred :: bulk_properties
     procedure(eos_phase_composition_procedure), public, deferred :: phase_composition
     procedure(eos_phase_properties_procedure), public, deferred :: phase_properties
     procedure(eos_primary_variables_procedure), public, deferred :: primary_variables
     procedure(eos_check_primary_variables_procedure), public, deferred :: check_primary_variables
     procedure(eos_conductivity_procedure), public, deferred :: conductivity
  end type eos_type

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

     subroutine eos_destroy_procedure(self)
       !! Destroy the EOS object
       import :: eos_type
       class(eos_type), intent(in out) :: self
     end subroutine eos_destroy_procedure

     subroutine eos_transition_procedure(self, primary, old_fluid, fluid, &
          transition, err)
       !! Check primary variables for a cell and make thermodynamic
       !! region transitions if needed.
       use fluid_module, only: fluid_type
       import :: eos_type
       class(eos_type), intent(in out) :: self
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

     subroutine eos_phase_composition_procedure(self, fluid, err)
       !! Calculate fluid phase composition from bulk properties.
       use fluid_module, only: fluid_type
       import :: eos_type
       class(eos_type), intent(in out) :: self
       type(fluid_type), intent(in out) :: fluid
       PetscErrorCode, intent(out) :: err
     end subroutine eos_phase_composition_procedure

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

     PetscErrorCode function eos_check_primary_variables_procedure(self, fluid, primary)
       !! Check if primary variables are in acceptable bounds, and return
       !! error code accordingly.
       use fluid_module, only: fluid_type
       import :: eos_type
       class(eos_type), intent(in) :: self
       type(fluid_type), intent(in) :: fluid
       PetscReal, intent(in) :: primary(self%num_primary_variables)
     end function eos_check_primary_variables_procedure

     PetscReal function eos_conductivity_procedure(self, rock, fluid)
       !! Calculate effective rock heat conductivity for given fluid properties.
       use rock_module, only: rock_type
       use fluid_module, only: fluid_type
       import :: eos_type
       class(eos_type), intent(in) :: self
       type(rock_type), intent(in) :: rock
       type(fluid_type), intent(in) :: fluid
     end function eos_conductivity_procedure

  end interface

end module eos_module
