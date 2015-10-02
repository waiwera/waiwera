module eos_module
  !! Equations of state.

  ! All EOSes are in here at present. This is to work around a
  ! gfortran 4.7 bug (free_pi_tree(): Unresolved fixup) which occurs
  ! when they are in their own modules. When stable compilers no
  ! longer have this problem we can give each EOS its own module
  ! again.

  use kinds_module
  use fson
  use thermodynamics_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>

  PetscInt, parameter, public :: max_eos_name_length = 8
  PetscInt, parameter, public :: max_eos_description_length = 80
  PetscInt, parameter, public :: max_primary_variable_name_length = 16
  PetscInt, parameter, public :: max_phase_name_length = 6

  type, public, abstract :: eos_type
     !! Abstract type for equation of state (EOS) objects.
     private
     character(max_eos_name_length), public :: name
     character(max_eos_description_length), public :: description
     character(max_primary_variable_name_length), allocatable, public :: primary_variable_names(:)
     character(max_phase_name_length), allocatable, public :: phase_names(:)
     PetscInt, allocatable, public :: phase_index(:)
     PetscInt, public :: num_primary_variables
     PetscInt, public :: num_phases
     PetscInt, public :: num_components
     PetscBool, public :: isothermal = .false.
     class(thermodynamics_type), pointer, public :: thermo
   contains
     private
     procedure(eos_transition_procedure), deferred :: transition
     procedure(eos_init_procedure), public, deferred :: init
     procedure(eos_destroy_procedure), public, deferred :: destroy
     procedure(eos_check_primary_procedure), public, deferred :: check_primary
     procedure(eos_fluid_procedure), public, deferred :: fluid_properties
  end type eos_type

  type, public, extends(eos_type) :: eos_w_type
     !! Isothermal pure water equation of state type.
     private
     PetscReal, public :: temperature  !! Constant temperature
   contains
     private
     procedure :: transition => eos_w_transition
     procedure, public :: init => eos_w_init
     procedure, public :: destroy => eos_w_destroy
     procedure, public :: check_primary => eos_w_check_primary
     procedure, public :: fluid_properties => eos_w_fluid_properties
  end type eos_w_type

  type, public, extends(eos_w_type) :: eos_we_type
     !! Pure water and energy equation of state type.
     private
   contains
     private
     procedure :: transition => eos_we_transition
     procedure, public :: init => eos_we_init
     procedure, public :: check_primary => eos_we_check_primary
     procedure, public :: fluid_properties => eos_we_fluid_properties
  end type eos_we_type

  abstract interface

     subroutine eos_init_procedure(self, json, thermo)
       !! Initialise EOS object
       import :: eos_type, thermodynamics_type, fson_value
       class(eos_type), intent(in out) :: self
       type(fson_value), pointer, intent(in) :: json
       class(thermodynamics_type), intent(in), target :: thermo
     end subroutine eos_init_procedure

     subroutine eos_destroy_procedure(self)
       !! Destroy EOS object
       import :: eos_type
       class(eos_type), intent(in out) :: self
     end subroutine eos_destroy_procedure

     subroutine eos_transition_procedure(self, &
          region1, region2, primary)
       !! Perform transitions between thermodynamic regions
       import :: eos_type, dp
       class(eos_type), intent(in) :: self
       PetscInt, intent(in) :: region1, region2
       PetscReal, intent(in out), target :: primary(self%num_primary_variables)
     end subroutine eos_transition_procedure

     subroutine eos_check_primary_procedure(self, region, primary)
       !! Check primary variables for current region and make
       !! transition if needed
       import :: eos_type, dp
       class(eos_type), intent(in) :: self
       PetscInt, intent(in out) :: region
       PetscReal, intent(in out), target :: primary(self%num_primary_variables)
     end subroutine eos_check_primary_procedure

     subroutine eos_fluid_procedure(self, region, primary, fluid)
       !! Calculate fluid properties from region and primary variables
       use fluid_module, only: fluid_type
       import :: eos_type, dp
       class(eos_type), intent(in out) :: self
       PetscInt, intent(in) :: region
       PetscReal, intent(in), target :: primary(self%num_primary_variables)
       type(fluid_type), intent(in out) :: fluid
     end subroutine eos_fluid_procedure

  end interface

public :: setup_eos

contains

!------------------------------------------------------------------------
! EOS setup routine
!------------------------------------------------------------------------

  subroutine setup_eos(json, thermo, eos)
    !! Reads equation of state from JSON input file.
    !! If not present, a default value is assigned.

    use utils_module, only : str_to_lower
    use fson_mpi_module, only: fson_get_mpi

    type(fson_value), pointer, intent(in) :: json
    class(thermodynamics_type), intent(in) :: thermo
    class(eos_type), allocatable, intent(in out) :: eos
    ! Locals:
    character(max_eos_name_length), parameter :: &
         default_eos_name = "we"
    character(max_eos_name_length) :: eos_name

    call fson_get_mpi(json, "eos.name", default_eos_name, eos_name)
    eos_name = str_to_lower(eos_name)

    select case (eos_name)
    case ("w")
       allocate(eos_w_type :: eos)
    case ("we")
       allocate(eos_we_type :: eos)
    case default
       allocate(eos_we_type :: eos)
    end select

    call eos%init(json, thermo)

  end subroutine setup_eos

!------------------------------------------------------------------------
! eos_w
!------------------------------------------------------------------------

  subroutine eos_w_init(self, json, thermo)
    !! Initialise isothermal pure water EOS.

    use fson_mpi_module, only: fson_get_mpi

    class(eos_w_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    ! Locals:
    PetscReal, parameter :: default_temperature = 20._dp ! deg C

    self%name = "w"
    self%description = "Isothermal pure water"
    self%primary_variable_names = ["Pressure"]
    self%phase_names = ["Liquid"]
    self%phase_index = [1]

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 1
    self%num_components = 1
    self%isothermal = .true.

    self%thermo => thermo

    call fson_get_mpi(json, "eos.temperature", default_temperature, &
         self%temperature)

  end subroutine eos_w_init

!------------------------------------------------------------------------

  subroutine eos_w_destroy(self)
    !! Destroy isothermal pure water EOS.

    class(eos_w_type), intent(in out) :: self

    deallocate(self%primary_variable_names)
    deallocate(self%phase_names, self%phase_index)

  end subroutine eos_w_destroy

!------------------------------------------------------------------------
  
  subroutine eos_w_transition(self, region1, region2, primary)
    !! Perform transitions between thermodynamic regions for isothermal
    !! pure water

    class(eos_w_type), intent(in) :: self
    PetscInt, intent(in) :: region1, region2
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)

    continue ! no transitions needed

  end subroutine eos_w_transition

!------------------------------------------------------------------------

  subroutine eos_w_check_primary(self, region, primary)
    !! Check primary variables for current region and make
    !! transition if needed for isothermal pure water

    class(eos_w_type), intent(in) :: self
    PetscInt, intent(in out) :: region
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)

    continue ! no checks needed

  end subroutine eos_w_check_primary

!------------------------------------------------------------------------

  subroutine eos_w_fluid_properties(self, region, primary, fluid)
    !! Calculate fluid properties from region and primary variables
    !! for isothermal pure water.

    use fluid_module, only: fluid_type
    class(eos_w_type), intent(in out) :: self
    PetscInt, intent(in) :: region !! Thermodynamic region index
    PetscReal, intent(in), target :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object

    fluid%pressure = primary(1)
    fluid%temperature = self%temperature
    fluid%region = region

    call fluid%phase_properties(self%thermo)

  end subroutine eos_w_fluid_properties

!------------------------------------------------------------------------
! eos_we
!------------------------------------------------------------------------

  subroutine eos_we_init(self, json, thermo)
    !! Initialise pure water and energy EOS.

    use fson_mpi_module, only: fson_get_mpi

    class(eos_we_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object

    self%name = "we"
    self%description = "Pure water and energy"
    self%primary_variable_names = ["Pressure   ", "Temperature"]
    self%phase_names = ["Liquid", "Vapour"]
    self%phase_index = [1, 2]

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 2
    self%num_components = 1

    self%thermo => thermo

  end subroutine eos_we_init

!------------------------------------------------------------------------
  
  subroutine eos_we_transition(self, region1, region2, primary)
    !! Perform transitions between thermodynamic regions.

    class(eos_we_type), intent(in) :: self
    PetscInt, intent(in) :: region1, region2
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)

    continue ! TODO

  end subroutine eos_we_transition

!------------------------------------------------------------------------

  subroutine eos_we_check_primary(self, region, primary)
    !! Check primary variables for current region and make
    !! transition if needed.

    class(eos_we_type), intent(in) :: self
    PetscInt, intent(in out) :: region
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)

    continue ! TODO

  end subroutine eos_we_check_primary

!------------------------------------------------------------------------

  subroutine eos_we_fluid_properties(self, region, primary, fluid)
    !! Calculate fluid properties from region and primary variables.

    use fluid_module, only: fluid_type
    class(eos_we_type), intent(in out) :: self
    PetscInt, intent(in) :: region !! Thermodynamic region index
    PetscReal, intent(in), target :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object

    fluid%pressure = primary(1)
    fluid%temperature = primary(2)
    fluid%region = region

    call fluid%phase_properties(self%thermo)

  end subroutine eos_we_fluid_properties

!------------------------------------------------------------------------

end module eos_module
