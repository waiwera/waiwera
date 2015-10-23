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
     PetscBool, public :: isothermal = .false.
     class(thermodynamics_type), pointer, public :: thermo
   contains
     private
     procedure(eos_init_procedure), public, deferred :: init
     procedure(eos_destroy_procedure), public, deferred :: destroy
     procedure(eos_transition_procedure), public, deferred :: transition
     procedure(eos_bulk_properties_procedure), public, deferred :: bulk_properties
     procedure(eos_phase_composition_procedure), public, deferred :: phase_composition
     procedure(eos_phase_properties_procedure), public, deferred :: phase_properties
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
     procedure :: phase_saturations => eos_we_phase_saturations
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

     subroutine eos_transition_procedure(self, primary, fluid)
       !! Check primary variables for a cell and make thermodynamic
       !! region transitions if needed.
       use fluid_module, only: fluid_type
       import :: eos_type, dp
       class(eos_type), intent(in out) :: self
       PetscReal, intent(in out), target :: primary(self%num_primary_variables)
       type(fluid_type), intent(in out) :: fluid
     end subroutine eos_transition_procedure

     subroutine eos_bulk_properties_procedure(self, primary, fluid)
       !! Calculate bulk fluid properties from primary variables.
       use fluid_module, only: fluid_type
       import :: eos_type, dp
       class(eos_type), intent(in out) :: self
       PetscReal, intent(in) :: primary(self%num_primary_variables)
       type(fluid_type), intent(in out) :: fluid
     end subroutine eos_bulk_properties_procedure

     subroutine eos_phase_composition_procedure(self, fluid)
       !! Calculate fluid phase composition from bulk properties.
       use fluid_module, only: fluid_type
       import :: eos_type, dp
       class(eos_type), intent(in out) :: self
       type(fluid_type), intent(in out) :: fluid
     end subroutine eos_phase_composition_procedure

     subroutine eos_phase_properties_procedure(self, primary, rock, fluid)
       !! Calculate phase fluid properties from primary variables.
       use rock_module, only: rock_type
       use fluid_module, only: fluid_type
       import :: eos_type, dp
       class(eos_type), intent(in out) :: self
       PetscReal, intent(in), target :: primary(self%num_primary_variables)
       type(rock_type), intent(in out) :: rock
       type(fluid_type), intent(in out) :: fluid
     end subroutine eos_phase_properties_procedure

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

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 1
    self%phase_names = ["liquid"]
    self%num_components = 1
    self%component_names = ["water"]
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
    deallocate(self%phase_names, self%component_names)
    nullify(self%thermo)

  end subroutine eos_w_destroy

!------------------------------------------------------------------------
  
  subroutine eos_w_transition(self, primary, fluid)
    !! For eos_w, check primary variables for a cell and make
    !! thermodynamic region transitions if needed

    use fluid_module, only: fluid_type

    class(eos_w_type), intent(in out) :: self
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid

    continue ! no transitions needed

  end subroutine eos_w_transition

!------------------------------------------------------------------------

  subroutine eos_w_bulk_properties(self, primary, fluid)
    !! Calculate fluid bulk properties from region and primary variables
    !! for isothermal pure water.

    use fluid_module, only: fluid_type

    class(eos_w_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object

    fluid%pressure = primary(1)
    fluid%temperature = self%temperature

  end subroutine eos_w_bulk_properties

!------------------------------------------------------------------------

  subroutine eos_w_phase_composition(self, fluid)
    !! Determines fluid phase composition from bulk properties and
    !! thermodynamic region.

    use fluid_module, only: fluid_type

    class(eos_w_type), intent(in out) :: self
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    ! Locals:
    PetscInt :: region, phases

    region = nint(fluid%region)
    phases = self%thermo%phase_composition(region, fluid%pressure, &
         fluid%temperature)
    fluid%phase_composition = dble(phases)

  end subroutine eos_w_phase_composition

!------------------------------------------------------------------------

  subroutine eos_w_phase_properties(self, primary, rock, fluid)
    !! Calculate fluid phase properties from region and primary variables
    !! for isothermal pure water.

    use fluid_module, only: fluid_type
    use rock_module, only: rock_type

    class(eos_w_type), intent(in out) :: self
    PetscReal, intent(in), target :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    ! Locals:
    PetscInt :: region, p, ierr
    PetscReal :: properties(2)

    region = nint(fluid%region)

    p = region

    call self%thermo%region(region)%ptr%properties( &
         [fluid%pressure, fluid%temperature], &
         properties, ierr)

    fluid%phase(p)%density = properties(1)
    fluid%phase(p)%internal_energy = properties(2)
    fluid%phase(p)%specific_enthalpy = &
         fluid%phase(p)%internal_energy + &
         fluid%pressure / fluid%phase(p)%density

    fluid%phase(p)%saturation = 1._dp
    fluid%phase(p)%relative_permeability = 1._dp
    fluid%phase(p)%mass_fraction(1) = 1._dp

    call self%thermo%region(region)%ptr%viscosity( &
         fluid%temperature, fluid%pressure, &
         fluid%phase(p)%density, fluid%phase(p)%viscosity)

  end subroutine eos_w_phase_properties

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

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 2
    self%phase_names = ["liquid", "vapour"]
    self%num_components = 1
    self%component_names = ["water"]

    self%thermo => thermo

  end subroutine eos_we_init

!------------------------------------------------------------------------

  subroutine eos_we_transition_to_single_phase(self, primary, fluid, &
       new_region)
    !! For eos_we, make transition from two-phase to single-phase with
    !! specified region.

    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in out) :: self
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscInt, intent(in) :: new_region
    ! Locals:
    PetscReal, pointer :: pressure, temperature
    PetscReal :: old_saturation_pressure, factor
    PetscInt :: ierr
    PetscReal, parameter :: small = 1.e-6_dp

    pressure => primary(1)
    temperature => primary(2)

    call self%thermo%saturation%pressure(fluid%temperature, &
         old_saturation_pressure, ierr)

    if (ierr >= 0) then

       if (new_region == 1) then
          factor = 1._dp + small
       else
          factor = 1._dp - small
       end if

       pressure = factor * old_saturation_pressure
       temperature = fluid%temperature

       fluid%region = dble(new_region)
       call self%phase_composition(fluid)

    else
       fluid%phase_composition = -1
    end if

  end subroutine eos_we_transition_to_single_phase

!------------------------------------------------------------------------

  subroutine eos_we_transition_to_two_phase(self, primary, fluid, &
       saturation_pressure)
    !! For eos_we, make transition from single-phase to two-phase.

    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in out) :: self
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscReal, intent(in) :: saturation_pressure
    ! Locals:
    PetscReal, pointer :: pressure, vapour_saturation
    PetscInt :: old_region
    PetscReal, parameter :: small = 1.e-6_dp

    pressure => primary(1)
    vapour_saturation => primary(2)
    old_region = nint(fluid%region)

    pressure = saturation_pressure

    if (old_region == 1) then
       vapour_saturation = small
    else
       vapour_saturation = 1._dp - small
    end if

    fluid%region = dble(4)
    call self%phase_composition(fluid)

  end subroutine eos_we_transition_to_two_phase

!------------------------------------------------------------------------

  subroutine eos_we_transition(self, primary, fluid)
    !! For eos_we, check primary variables for a cell and make
    !! thermodynamic region transitions if needed.

    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in out) :: self
    PetscReal, intent(in out), target :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    ! Locals:
    PetscInt :: region, ierr
    PetscReal, pointer :: pressure, temperature, vapour_saturation
    PetscReal :: saturation_pressure

    region = nint(fluid%region)

    if (region == 4) then  ! Two-phase

       vapour_saturation => primary(2)

       if (vapour_saturation < 0._dp) then
          call self%transition_to_single_phase(primary, fluid, 1)
       else if (vapour_saturation > 1._dp) then
          call self%transition_to_single_phase(primary, fluid, 2)
       end if

    else  ! Single-phase

       pressure => primary(1)
       temperature => primary(2)

       call self%thermo%saturation%pressure(temperature, &
            saturation_pressure, ierr)

       if (ierr >= 0) then

          if (((region == 1) .and. (pressure < saturation_pressure)) .or. &
               ((region == 2) .and. (pressure > saturation_pressure))) then
             call self%transition_to_two_phase(primary, fluid, &
                  saturation_pressure)
          end if

       else
          fluid%phase_composition = -1
       end if

    end if

  end subroutine eos_we_transition

!------------------------------------------------------------------------

  subroutine eos_we_bulk_properties(self, primary, fluid)
    !! Calculate fluid bulk properties from region and primary variables
    !! for non-isothermal pure water.

    use fluid_module, only: fluid_type
    class(eos_we_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    ! Locals:
    PetscInt :: region, ierr

    fluid%pressure = primary(1)
    region = nint(fluid%region)

    if (region == 4) then
       ! Two-phase
       call self%thermo%saturation%temperature(fluid%pressure, &
            fluid%temperature, ierr)
    else
       ! Single-phase
       fluid%temperature = primary(2)
    end if

  end subroutine eos_we_bulk_properties

!------------------------------------------------------------------------

  subroutine eos_we_phase_saturations(self, region, primary, saturation)
    !! Calculates fluid phase saturations from region and primary variables.

    use fluid_module, only: fluid_type
    class(eos_we_type), intent(in out) :: self
    PetscInt, intent(in) :: region
    PetscReal, intent(in), target :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    PetscReal, intent(out), target :: saturation(self%num_primary_variables) !! Phase saturation

    select case (region)
    case (1)
       saturation = [1._dp, 0._dp]
    case (2)
       saturation = [0._dp, 1._dp]
    case (4)
       saturation = [1._dp - primary(2), primary(2)]
    end select

  end subroutine eos_we_phase_saturations

!------------------------------------------------------------------------

  subroutine eos_we_phase_properties(self, primary, rock, fluid)
    !! Calculate fluid phase properties from updated fluid region and primary variables.

    use rock_module, only: rock_type
    use fluid_module, only: fluid_type

    class(eos_we_type), intent(in out) :: self
    PetscReal, intent(in), target :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    ! Locals:
    PetscInt :: p, ierr, phases, region
    PetscReal :: properties(2), saturation(2), relative_permeability(2)

    region = nint(fluid%region)
    phases = nint(fluid%phase_composition)

    if (phases >= 0) then

       call self%phase_saturations(region, primary, saturation)
       relative_permeability = rock%relative_permeability%values(saturation(1))

       do p = 1, self%num_phases

          fluid%phase(p)%saturation = saturation(p)

          if (btest(phases, p - 1)) then

             call self%thermo%region(p)%ptr%properties( &
                  [fluid%pressure, fluid%temperature], &
                  properties, ierr)

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
             fluid%phase(p)%density = 0._dp
             fluid%phase(p)%internal_energy = 0._dp
             fluid%phase(p)%specific_enthalpy = 0._dp
             fluid%phase(p)%relative_permeability = 0._dp
             fluid%phase(p)%viscosity = 0._dp
             fluid%phase(p)%mass_fraction(1) = 0._dp
          end if

       end do

    else
       fluid%pressure = qnan_dp
    end if

  end subroutine eos_we_phase_properties

!------------------------------------------------------------------------

end module eos_module
