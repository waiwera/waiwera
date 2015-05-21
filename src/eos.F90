module eos_module
  !! Abstract type for equation of state (EOS) objects.

  use kinds_module
  use thermodynamics_module

  implicit none
  private

  integer, parameter, public :: max_eos_name_length = 8
  integer, parameter, public :: max_eos_description_length = 80
  integer, parameter, public :: max_primary_variable_name_length = 16

  type, public, abstract :: eos_type
     private
     character(max_eos_name_length), public :: name
     character(max_eos_description_length), public :: description
     character(max_primary_variable_name_length), allocatable, public :: primary_variable_names(:)
     integer, public :: num_primary_variables
     integer, public :: num_phases
     integer, public :: num_components
     class(thermodynamics_type), pointer, public :: thermo
   contains
     private
     procedure(eos_transition_procedure), deferred :: transition
     procedure(eos_init_procedure), public, deferred :: init
     procedure(eos_check_primary_procedure), public, deferred :: check_primary
     procedure(eos_fluid_procedure), public, deferred :: fluid_properties
  end type eos_type

  abstract interface

     subroutine eos_init_procedure(self, thermo)
       !! Initialise EOS object
       import :: eos_type, thermodynamics_type
       class(eos_type), intent(in out) :: self
       class(thermodynamics_type), intent(in), target :: thermo
     end subroutine eos_init_procedure

     subroutine eos_transition_procedure(self, &
          region1, region2, primary)
       !! Perform transitions between thermodynamic regions
       import :: eos_type, dp
       class(eos_type), intent(in) :: self
       integer, intent(in) :: region1, region2
       real(dp), intent(in out), target :: primary(self%num_primary_variables)
     end subroutine eos_transition_procedure

     subroutine eos_check_primary_procedure(self, region, primary)
       !! Check primary variables for current region and make
       !! transition if needed
       import :: eos_type, dp
       class(eos_type), intent(in) :: self
       integer, intent(in out) :: region
       real(dp), intent(in out), target :: primary(self%num_primary_variables)
     end subroutine eos_check_primary_procedure

     subroutine eos_fluid_procedure(self, region, primary, fluid)
       !! Calculate fluid properties from region and primary variables
       use fluid_module, only: fluid_type
       import :: eos_type, dp
       class(eos_type), intent(in out) :: self
       integer, intent(in) :: region
       real(dp), intent(in), target :: primary(self%num_primary_variables)
       type(fluid_type), intent(out) :: fluid
     end subroutine eos_fluid_procedure

  end interface

end module eos_module
