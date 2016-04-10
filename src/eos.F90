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
#include <petsc/finclude/petscsys.h>

  PetscInt, parameter, public :: max_eos_name_length = 8
  PetscInt, parameter, public :: max_eos_description_length = 80
  PetscInt, parameter, public :: max_primary_variable_name_length = 32
  PetscInt, parameter, public :: max_phase_name_length = 13
  PetscInt, parameter, public :: max_component_name_length = 8

  type, public, abstract :: eos_type
     !! Abstract type for equation of state (EOS) objects.
     private
     character(max_eos_name_length), public :: name
     character(max_eos_description_length), public :: description
     character(max_primary_variable_name_length), allocatable, public :: primary_variable_names(:)
     character(max_phase_name_length), allocatable, public :: phase_names(:)
     character(max_component_name_length), allocatable, public :: component_names(:)
     PetscInt, public :: num_primary_variables
     PetscInt, public :: num_phases
     PetscInt, public :: num_components
     PetscReal, allocatable, public :: default_primary(:)
     PetscInt, public :: default_region
     class(thermodynamics_type), pointer, public :: thermo
     PetscBool, public :: isothermal = PETSC_FALSE
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

  type, public, extends(eos_type) :: eos_w_type
     !! Isothermal pure water equation of state type.
     private
     PetscReal, public :: temperature  !! Constant temperature
   contains
     private
     procedure, public :: init => eos_w_init
     procedure, public :: destroy => eos_w_destroy
     procedure, public :: transition => eos_w_transition
     procedure, public :: bulk_properties => eos_w_bulk_properties
     procedure, public :: phase_composition => eos_w_phase_composition
     procedure, public :: phase_properties => eos_w_phase_properties
     procedure, public :: primary_variables => eos_w_primary_variables
     procedure, public :: check_primary_variables => eos_w_check_primary_variables
     procedure, public :: conductivity => eos_w_conductivity
  end type eos_w_type

  type, public, extends(eos_w_type) :: eos_we_type
     !! Pure water and energy equation of state type.
     private
   contains
     private
     procedure, public :: init => eos_we_init
     procedure, public :: transition => eos_we_transition
     procedure, public :: transition_to_single_phase => eos_we_transition_to_single_phase
     procedure, public :: transition_to_two_phase => eos_we_transition_to_two_phase
     procedure, public :: bulk_properties => eos_we_bulk_properties
     procedure, public :: phase_properties => eos_we_phase_properties
     procedure, public :: primary_variables => eos_we_primary_variables
     procedure :: phase_saturations => eos_we_phase_saturations
     procedure, public :: check_primary_variables => eos_we_check_primary_variables
     procedure, public :: conductivity => eos_we_conductivity
  end type eos_we_type

  abstract interface

     subroutine eos_init_procedure(self, json, thermo, logfile)
       !! Initialise EOS object
       use logfile_module
       import :: eos_type, thermodynamics_type, fson_value
       class(eos_type), intent(in out) :: self
       type(fson_value), pointer, intent(in) :: json
       class(thermodynamics_type), intent(in), target :: thermo
       type(logfile_type), intent(in out), optional :: logfile
     end subroutine eos_init_procedure

     subroutine eos_destroy_procedure(self)
       !! Destroy EOS object
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
       PetscReal, intent(in out), target :: primary(self%num_primary_variables)
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
       PetscReal, intent(in), target :: primary(self%num_primary_variables)
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
       PetscReal, intent(in), target :: primary(self%num_primary_variables)
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

public :: setup_eos

contains

!------------------------------------------------------------------------
! EOS setup routine
!------------------------------------------------------------------------

  subroutine setup_eos(json, thermo, eos, logfile)
    !! Reads equation of state from JSON input file.
    !! If not present, a default value is assigned.

    use utils_module, only : str_to_lower
    use fson_mpi_module, only: fson_get_mpi
    use logfile_module

    type(fson_value), pointer, intent(in) :: json
    class(thermodynamics_type), intent(in) :: thermo
    class(eos_type), allocatable, intent(in out) :: eos
    type(logfile_type), intent(in out) :: logfile
    ! Locals:
    character(max_eos_name_length), parameter :: &
         default_eos_name = "we"
    character(max_eos_name_length) :: eos_name

    call fson_get_mpi(json, "eos.name", default_eos_name, &
         eos_name, logfile)
    eos_name = str_to_lower(eos_name)

    select case (eos_name)
    case ("w")
       allocate(eos_w_type :: eos)
    case ("we")
       allocate(eos_we_type :: eos)
    case default
       allocate(eos_we_type :: eos)
    end select

    call eos%init(json, thermo, logfile)

  end subroutine setup_eos

!------------------------------------------------------------------------
! eos_w
!------------------------------------------------------------------------

  subroutine eos_w_init(self, json, thermo, logfile)
    !! Initialise isothermal pure water EOS.

    use fson_mpi_module, only: fson_get_mpi
    use logfile_module

    class(eos_w_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscReal, parameter :: default_pressure = 1.0e5_dp
    PetscReal, parameter :: default_temperature = 20._dp ! deg C

    self%name = "w"
    self%description = "Isothermal pure water"
    self%primary_variable_names = ["Pressure"]

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 1
    self%phase_names = ["liquid"]
    self%num_components = 1
    self%component_names = ["water"]
    self%isothermal = PETSC_TRUE

    self%default_primary = [default_pressure]
    self%default_region = 1

    self%thermo => thermo

    call fson_get_mpi(json, "eos.temperature", default_temperature, &
         self%temperature, logfile)

  end subroutine eos_w_init

!------------------------------------------------------------------------

  subroutine eos_w_destroy(self)
    !! Destroy isothermal pure water EOS.

    class(eos_w_type), intent(in out) :: self

    deallocate(self%primary_variable_names)
    deallocate(self%phase_names, self%component_names)
    deallocate(self%default_primary)
    nullify(self%thermo)

  end subroutine eos_w_destroy

!------------------------------------------------------------------------
  
  subroutine eos_w_transition(self, primary, old_fluid, fluid, &
       transition, err)
    !! For eos_w, check primary variables for a cell and make
    !! thermodynamic region transitions if needed

    use fluid_module, only: fluid_type

    class(eos_w_type), intent(in out) :: self
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)
    type(fluid_type), intent(in) :: old_fluid
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err

    ! no transitions needed
    err = 0
    transition = PETSC_FALSE

  end subroutine eos_w_transition

!------------------------------------------------------------------------

  subroutine eos_w_bulk_properties(self, primary, fluid, err)
    !! Calculate fluid bulk properties from region and primary variables
    !! for isothermal pure water.

    use fluid_module, only: fluid_type

    class(eos_w_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err

    err = 0
    fluid%pressure = primary(1)
    fluid%temperature = self%temperature

    fluid%phase(1)%saturation = 1._dp
    call self%phase_composition(fluid, err)

  end subroutine eos_w_bulk_properties

!------------------------------------------------------------------------

  subroutine eos_w_phase_composition(self, fluid, err)
    !! Determines fluid phase composition from bulk properties and
    !! thermodynamic region.

    use fluid_module, only: fluid_type

    class(eos_w_type), intent(in out) :: self
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

  end subroutine eos_w_phase_composition

!------------------------------------------------------------------------

  subroutine eos_w_phase_properties(self, primary, rock, fluid, err)
    !! Calculate fluid phase properties from region and primary variables
    !! for isothermal pure water.
    !! Bulk properties need to be calculated before calling this routine.

    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(eos_w_type), intent(in out) :: self
    PetscReal, intent(in), target :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: region, p
    PetscReal :: properties(2)

    err = 0
    region = nint(fluid%region)

    p = region

    call self%thermo%region(region)%ptr%properties( &
         [fluid%pressure, fluid%temperature], &
         properties, err)

    if (err == 0) then

       fluid%phase(p)%density = properties(1)
       fluid%phase(p)%internal_energy = properties(2)
       fluid%phase(p)%specific_enthalpy = &
            fluid%phase(p)%internal_energy + &
            fluid%pressure / fluid%phase(p)%density

       fluid%phase(p)%relative_permeability = 1._dp
       fluid%phase(p)%mass_fraction(1) = 1._dp

       call self%thermo%region(region)%ptr%viscosity( &
            fluid%temperature, fluid%pressure, &
            fluid%phase(p)%density, fluid%phase(p)%viscosity)

    end if

  end subroutine eos_w_phase_properties

!------------------------------------------------------------------------

  subroutine eos_w_primary_variables(self, fluid, primary)
    !! Determine primary variables from fluid properties.

    use fluid_module, only: fluid_type

    class(eos_w_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(out) :: primary(self%num_primary_variables)

    primary(1) = fluid%pressure

  end subroutine eos_w_primary_variables

!------------------------------------------------------------------------

  PetscErrorCode function eos_w_check_primary_variables(self, fluid, &
       primary) result(err)
    !! Check if primary variables are in acceptable bounds, and return
    !! error code accordingly.

    use fluid_module, only: fluid_type

    class(eos_w_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(in), target :: primary(self%num_primary_variables)
    ! Locals:
    PetscReal, pointer :: p

    p => primary(1)
    if ((p < 0._dp) .or. (p > 100.e6_dp)) then
       err = 1
    else
       err = 0
    end if

  end function eos_w_check_primary_variables

!------------------------------------------------------------------------

  PetscReal function eos_w_conductivity(self, rock, fluid) result(cond)
    !! Returns effective rock heat conductivity for given fluid properties.

    use rock_module, only: rock_type
    use fluid_module, only: fluid_type

    class(eos_w_type), intent(in) :: self
    type(rock_type), intent(in) :: rock !! Rock object
    type(fluid_type), intent(in) :: fluid !! Fluid object

    cond = rock%wet_conductivity

  end function eos_w_conductivity

!------------------------------------------------------------------------
! eos_we
!------------------------------------------------------------------------

  subroutine eos_we_init(self, json, thermo, logfile)
    !! Initialise pure water and energy EOS.

    use logfile_module

    class(eos_we_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscReal, parameter :: default_pressure = 1.0e5_dp
    PetscReal, parameter :: default_temperature = 20._dp ! deg C

    self%name = "we"
    self%description = "Pure water and energy"
    self%primary_variable_names = ["Pressure                     ", &
         "Temperature/Vapour saturation"]

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 2
    self%phase_names = ["liquid", "vapour"]
    self%num_components = 1
    self%component_names = ["water"]

    self%default_primary = [default_pressure, default_temperature]
    self%default_region = 1

    self%thermo => thermo

  end subroutine eos_we_init

!------------------------------------------------------------------------

  subroutine eos_we_transition_to_single_phase(self, old_fluid, &
       new_region, primary, fluid, transition, err)
    !! For eos_we, make transition from two-phase to single-phase with
    !! specified region.

    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in out) :: self
    type(fluid_type), intent(in) :: old_fluid
    PetscInt, intent(in) :: new_region
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal, pointer :: pressure, temperature
    PetscReal :: old_saturation_pressure, factor
    PetscReal, parameter :: small = 1.e-6_dp

    err = 0
    pressure => primary(1)
    temperature => primary(2)

    call self%thermo%saturation%pressure(old_fluid%temperature, &
         old_saturation_pressure, err)

    if (err == 0) then

       if (new_region == 1) then
          factor = 1._dp + small
       else
          factor = 1._dp - small
       end if

       pressure = factor * old_saturation_pressure
       temperature = old_fluid%temperature

       fluid%region = dble(new_region)
       transition = PETSC_TRUE

    end if

  end subroutine eos_we_transition_to_single_phase

!------------------------------------------------------------------------

  subroutine eos_we_transition_to_two_phase(self, saturation_pressure, &
       old_fluid, primary, fluid, transition, err)
    !! For eos_we, make transition from single-phase to two-phase.

    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in out) :: self
    PetscReal, intent(in) :: saturation_pressure
    type(fluid_type), intent(in) :: old_fluid
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal, pointer :: pressure, vapour_saturation
    PetscInt :: old_region
    PetscReal, parameter :: small = 1.e-6_dp

    err = 0
    pressure => primary(1)
    vapour_saturation => primary(2)
    old_region = nint(old_fluid%region)

    pressure = saturation_pressure

    if (old_region == 1) then
       vapour_saturation = small
    else
       vapour_saturation = 1._dp - small
    end if

    fluid%region = dble(4)
    transition = PETSC_TRUE

  end subroutine eos_we_transition_to_two_phase

!------------------------------------------------------------------------

  subroutine eos_we_transition(self, primary, old_fluid, fluid, &
       transition, err)
    !! For eos_we, check primary variables for a cell and make
    !! thermodynamic region transitions if needed.

    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in out) :: self
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)
    type(fluid_type), intent(in) :: old_fluid
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: old_region
    PetscReal, pointer :: pressure, temperature, vapour_saturation
    PetscReal :: saturation_pressure

    err = 0
    transition = PETSC_FALSE
    old_region = nint(old_fluid%region)

    if (old_region == 4) then  ! Two-phase

       vapour_saturation => primary(2)

       if (vapour_saturation < 0._dp) then
          call self%transition_to_single_phase(old_fluid, 1, primary, &
               fluid, transition, err)
       else if (vapour_saturation > 1._dp) then
          call self%transition_to_single_phase(old_fluid, 2, primary, &
               fluid, transition, err)
       end if

    else  ! Single-phase

       pressure => primary(1)
       temperature => primary(2)

       call self%thermo%saturation%pressure(temperature, &
            saturation_pressure, err)

       if (err == 0) then
          if (((old_region == 1) .and. (pressure < saturation_pressure)) .or. &
               ((old_region == 2) .and. (pressure > saturation_pressure))) then
             call self%transition_to_two_phase(saturation_pressure, old_fluid, &
                  primary, fluid, transition, err)
          end if
       end if

    end if

  end subroutine eos_we_transition

!------------------------------------------------------------------------

  subroutine eos_we_bulk_properties(self, primary, fluid, err)
    !! Calculate fluid bulk properties from region and primary variables
    !! for non-isothermal pure water.

    use fluid_module, only: fluid_type
    class(eos_we_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscInt :: region

    err = 0
    fluid%pressure = primary(1)
    region = nint(fluid%region)

    if (region == 4) then
       ! Two-phase
       call self%thermo%saturation%temperature(fluid%pressure, &
            fluid%temperature, err)
    else
       ! Single-phase
       fluid%temperature = primary(2)
    end if

    if (err == 0) then
       call self%phase_composition(fluid, err)
       if (err == 0) then
          call self%phase_saturations(primary, fluid)
       end if
    end if

  end subroutine eos_we_bulk_properties

!------------------------------------------------------------------------

  subroutine eos_we_phase_saturations(self, primary, fluid)
    !! Assigns fluid phase saturations from fluid region and primary variables.

    use fluid_module, only: fluid_type
    class(eos_we_type), intent(in out) :: self
    PetscReal, intent(in), target :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    ! Locals:
    PetscInt :: region

    region = nint(fluid%region)

    select case (region)
    case (1)
       fluid%phase(1)%saturation = 1._dp
       fluid%phase(2)%saturation = 0._dp
    case (2)
       fluid%phase(1)%saturation = 0._dp
       fluid%phase(2)%saturation = 1._dp
    case (4)
       fluid%phase(1)%saturation = 1._dp - primary(2)
       fluid%phase(2)%saturation = primary(2)
    end select

  end subroutine eos_we_phase_saturations

!------------------------------------------------------------------------

  subroutine eos_we_phase_properties(self, primary, rock, fluid, err)
    !! Calculate fluid phase properties from updated fluid region and primary variables.
    !! Bulk properties need to be calculated before calling this routine.

    use rock_module, only: rock_type
    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in out) :: self
    PetscReal, intent(in), target :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscInt :: p, phases, region
    PetscReal :: properties(2), sl, relative_permeability(2)

    err = 0
    region = nint(fluid%region)
    phases = nint(fluid%phase_composition)

    sl = fluid%phase(1)%saturation
    relative_permeability = rock%relative_permeability%values(sl)

    do p = 1, self%num_phases

       if (btest(phases, p - 1)) then

          call self%thermo%region(p)%ptr%properties( &
               [fluid%pressure, fluid%temperature], &
               properties, err)

          if (err == 0) then

             fluid%phase(p)%density = properties(1)
             fluid%phase(p)%internal_energy = properties(2)
             fluid%phase(p)%specific_enthalpy = &
                  fluid%phase(p)%internal_energy + &
                  fluid%pressure / fluid%phase(p)%density

             fluid%phase(p)%mass_fraction(1) = 1._dp
             fluid%phase(p)%relative_permeability = relative_permeability(p)

             call self%thermo%region(p)%ptr%viscosity( &
                  fluid%temperature, fluid%pressure, &
                  fluid%phase(p)%density, fluid%phase(p)%viscosity)

          else
             exit
          end if

       else

          fluid%phase(p)%density = 0._dp
          fluid%phase(p)%internal_energy = 0._dp
          fluid%phase(p)%specific_enthalpy = 0._dp
          fluid%phase(p)%relative_permeability = 0._dp
          fluid%phase(p)%viscosity = 0._dp
          fluid%phase(p)%mass_fraction(1) = 0._dp

       end if

    end do

  end subroutine eos_we_phase_properties

!------------------------------------------------------------------------

  subroutine eos_we_primary_variables(self, fluid, primary)
    !! Determine primary variables from fluid properties.

    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(out) :: primary(self%num_primary_variables)
    ! Locals:
    PetscInt :: region

    primary(1) = fluid%pressure

    region = nint(fluid%region)
    if (region == 4) then
       primary(2) = fluid%phase(2)%saturation
    else
       primary(2) = fluid%temperature
    end if

  end subroutine eos_we_primary_variables

!------------------------------------------------------------------------

  PetscErrorCode function eos_we_check_primary_variables(self, fluid, &
       primary) result(err)
    !! Check if primary variables are in acceptable bounds, and return error
    !! code accordingly.

    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(in), target :: primary(self%num_primary_variables)
    ! Locals:
    PetscInt :: region
    PetscReal, pointer :: p, t, vapour_saturation

    err = 0
    p => primary(1)
    if ((p < 0._dp) .or. (p > 100.e6_dp)) then
       err = 1
    else
       region = nint(fluid%region)
       if (region == 4) then
          vapour_saturation => primary(2)
          if ((vapour_saturation < -1._dp) .or. &
               (vapour_saturation > 2._dp)) then
             err = 1
          end if
       else
          t => primary(2)
          if ((t < 0._dp) .or. (t > 800._dp)) then
             err = 1
          end if
       end if
    end if

  end function eos_we_check_primary_variables

!------------------------------------------------------------------------

  PetscReal function eos_we_conductivity(self, rock, fluid) result(cond)
    !! Returns effective rock heat conductivity for given fluid properties.

    use rock_module, only: rock_type
    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in) :: self
    type(rock_type), intent(in) :: rock !! Rock object
    type(fluid_type), intent(in) :: fluid !! Fluid object
    ! Locals:
    PetscReal :: sl

    sl = fluid%phase(1)%saturation
    cond = rock%dry_conductivity + sqrt(sl) * &
         (rock%wet_conductivity - rock%dry_conductivity)

  end function eos_we_conductivity

!------------------------------------------------------------------------

end module eos_module
