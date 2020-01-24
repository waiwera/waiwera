module eos_wge_module
  !! Equation of state for non-isothermal water and non-condensible gas.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use eos_module
  use ncg_thermodynamics_module
  use root_finder_module

  implicit none
  private

  type, public, extends(eos_type) :: eos_wge_type
     !! Pure water, non-condensible gas and energy equation of state type.
     private
     class(ncg_thermodynamics_type), allocatable, public :: gas
     type(root_finder_type) :: saturation_line_finder
     type(primary_variable_interpolator_type), pointer :: &
          primary_variable_interpolator
   contains
     private
     procedure, public :: init => eos_wge_init
     procedure, public :: destroy => eos_wge_destroy
     procedure, public :: transition => eos_wge_transition
     procedure, public :: transition_to_single_phase => eos_wge_transition_to_single_phase
     procedure, public :: transition_to_two_phase => eos_wge_transition_to_two_phase
     procedure, public :: bulk_properties => eos_wge_bulk_properties
     procedure, public :: phase_properties => eos_wge_phase_properties
     procedure, public :: primary_variables => eos_wge_primary_variables
     procedure, public :: phase_saturations => eos_wge_phase_saturations
     procedure, public :: check_primary_variables => eos_wge_check_primary_variables
  end type eos_wge_type

contains

!------------------------------------------------------------------------

  subroutine eos_wge_init(self, json, thermo, logfile)
    !! Initialise pure water, non-condensible gas and energy EOS.

    use fson
    use fson_mpi_module, only: fson_get_mpi, fson_type_mpi
    use fson_value_m, only: TYPE_STRING, TYPE_REAL, TYPE_NULL
    use logfile_module
    use thermodynamics_module

    class(eos_wge_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    procedure(root_finder_function), pointer :: f
    class(*), pointer :: pinterp
    PetscReal, allocatable :: data(:, :)
    PetscReal :: pressure_scale, temperature_scale, partial_pressure_scale
    PetscInt :: scale_type
    PetscErrorCode :: err
    PetscReal, parameter :: default_pressure = 1.0e5_dp
    PetscReal, parameter :: default_temperature = 20._dp ! deg C
    PetscReal, parameter :: default_gas_partial_pressure = 0._dp
    PetscReal, parameter :: default_pressure_scale = 1.e6_dp !! Default scale factor for non-dimensionalising pressure
    PetscReal, parameter :: default_temperature_scale = 1.e2_dp !! Default scale factor for non-dimensionalising temperature
    PetscReal, parameter :: default_partial_pressure_scale = 1.e6_dp !! Default scale factor for non-dimensionalising partial pressure

    self%name = "wge"
    self%description = "Water, non-condensible gas and energy"
    self%primary_variable_names = [ &
         "pressure                     ", &
         "temperature/vapour_saturation", &
         "gas partial pressure         "]

    self%num_primary_variables = size(self%primary_variable_names)
    self%num_phases = 2
    self%phase_names = ["liquid", "vapour"]
    self%num_components = 2
    self%component_names = ["water", "gas  "]

    self%default_primary = [default_pressure, default_temperature, &
         default_gas_partial_pressure]
    self%default_region = 1
    self%required_output_fluid_fields = [ &
         "pressure            ", "temperature         ", &
         "region              ", "gas_partial_pressure", &
         "vapour_saturation   "]
    self%default_output_fluid_fields = [ &
         "pressure              ", "temperature           ", &
         "region                ", "gas_partial_pressure  ", &
         "vapour_saturation     "]

    call fson_get_mpi(json, "eos.primary.scale.pressure", default_pressure_scale, &
         pressure_scale, logfile)
    call fson_get_mpi(json, "eos.primary.scale.temperature", default_temperature_scale, &
         temperature_scale, logfile)

    scale_type = fson_type_mpi(json, "eos.primary.scale.partial_pressure")
    select case (scale_type)
    case (TYPE_STRING, TYPE_NULL)
       self%scale => eos_wge_scale_adaptive
       self%unscale => eos_wge_unscale_adaptive
       partial_pressure_scale = 0._dp
    case (TYPE_REAL)
       call fson_get_mpi(json, "eos.primary.scale.partial_pressure", &
            default_partial_pressure_scale, partial_pressure_scale, logfile)
    end select
    allocate(self%primary_scale(3, 4))
    self%primary_scale = reshape([ &
          pressure_scale, temperature_scale, partial_pressure_scale, &
          pressure_scale, temperature_scale, partial_pressure_scale, &
          0._dp, 0._dp, 0._dp, &
          pressure_scale, 1._dp, partial_pressure_scale], [3, 4])

    self%thermo => thermo

    ! Set up saturation line finder:
    allocate(primary_variable_interpolator_type :: &
         self%primary_variable_interpolator)
    allocate(data(2, 1 + self%num_primary_variables))
    data = 0._dp
    data(:, 1) = [0._dp, 1._dp]
    call self%primary_variable_interpolator%init(data, err)
    deallocate(data)
    self%primary_variable_interpolator%thermo => self%thermo
    f => eos_wge_saturation_difference
    pinterp => self%primary_variable_interpolator
    call self%saturation_line_finder%init(f, context = pinterp)

  end subroutine eos_wge_init

!------------------------------------------------------------------------

  subroutine eos_wge_destroy(self)
    !! Destroy pure water, non-condensible gas and energy EOS.

    class(eos_wge_type), intent(in out) :: self

    deallocate(self%primary_variable_names)
    deallocate(self%phase_names, self%component_names)
    deallocate(self%default_primary)
    deallocate(self%primary_scale)
    self%thermo => null()

    call self%saturation_line_finder%destroy()
    call self%primary_variable_interpolator%destroy()
    deallocate(self%primary_variable_interpolator)
    call self%gas%destroy()

  end subroutine eos_wge_destroy

!------------------------------------------------------------------------

  subroutine eos_wge_transition_to_single_phase(self, old_primary, old_fluid, &
       new_region, primary, fluid, transition, err)
    !! For eos_wge, make transition from two-phase to single-phase with
    !! specified region.

    use fluid_module, only: fluid_type

    class(eos_wge_type), intent(in out) :: self
    type(fluid_type), intent(in) :: old_fluid
    PetscInt, intent(in) :: new_region
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: old_saturation_pressure, pressure_factor
    PetscReal :: saturation_bound, xi, interpolated_water_pressure
    PetscReal :: interpolated_primary(self%num_primary_variables)
    PetscReal, parameter :: small = 1.e-6_dp

    err = 0
    transition = PETSC_FALSE

    if (new_region == 1) then
       saturation_bound = 0._dp
       pressure_factor = 1._dp + small
    else
       saturation_bound = 1._dp
       pressure_factor = 1._dp - small
    end if

    associate (pressure => primary(1), temperature => primary(2), &
         partial_pressure => primary(3))

      partial_pressure = max(0._dp, min(partial_pressure, pressure))
      self%primary_variable_interpolator%val(:, 1) = old_primary
      self%primary_variable_interpolator%val(:, 2) = primary
      call self%primary_variable_interpolator%find_component_at_index(&
           saturation_bound, 2, xi, err)

      if (err == 0) then

         interpolated_primary = self%primary_variable_interpolator%interpolate(xi)
         associate(interpolated_pressure => interpolated_primary(1), &
              interpolated_partial_pressure => interpolated_primary(3))
           interpolated_water_pressure = interpolated_pressure - &
                interpolated_partial_pressure
           pressure = pressure_factor * interpolated_water_pressure + &
                interpolated_partial_pressure
           partial_pressure = interpolated_partial_pressure
         end associate

         call self%thermo%saturation%temperature(interpolated_water_pressure, &
              temperature, err)
         if (err == 0) then
            fluid%region = dble(new_region)
            transition = PETSC_TRUE
         end if

      else

         call self%thermo%saturation%pressure(old_fluid%temperature, &
              old_saturation_pressure, err)
         if (err == 0) then
            pressure = pressure_factor * old_saturation_pressure + &
                 partial_pressure
            temperature = old_fluid%temperature
            fluid%region = dble(new_region)
            transition = PETSC_TRUE
         end if

      end if

    end associate

  end subroutine eos_wge_transition_to_single_phase

!------------------------------------------------------------------------

  subroutine eos_wge_transition_to_two_phase(self, saturation_pressure, &
       old_primary, old_fluid, primary, fluid, transition, err)
    !! For eos_wge, make transition from single-phase to two-phase.

    use fluid_module, only: fluid_type

    class(eos_wge_type), intent(in out) :: self
    PetscReal, intent(in) :: saturation_pressure
    type(fluid_type), intent(in) :: old_fluid
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: old_region
    PetscReal :: interpolated_primary(self%num_primary_variables)
    PetscReal :: xi
    PetscReal, parameter :: small = 1.e-6_dp

    err = 0
    associate (pressure => primary(1), vapour_saturation => primary(2), &
         partial_pressure => primary(3))

      partial_pressure = max(0._dp, min(partial_pressure, pressure))
      self%primary_variable_interpolator%val(:, 1) = old_primary
      self%primary_variable_interpolator%val(:, 2) = primary
      call self%saturation_line_finder%find()

      if (self%saturation_line_finder%err == 0) then
         xi = self%saturation_line_finder%root
         interpolated_primary = self%primary_variable_interpolator%interpolate(xi)
         associate(interpolated_pressure => interpolated_primary(1), &
              interpolated_partial_pressure => interpolated_primary(3))
           pressure = interpolated_pressure
           partial_pressure = interpolated_partial_pressure
         end associate
      else
         pressure = saturation_pressure + partial_pressure
      end if

      old_region = nint(old_fluid%region)
      if (old_region == 1) then
         vapour_saturation = small
      else
         vapour_saturation = 1._dp - small
      end if

      fluid%region = dble(4)
      transition = PETSC_TRUE

    end associate

  end subroutine eos_wge_transition_to_two_phase

!------------------------------------------------------------------------

  subroutine eos_wge_transition(self, old_primary, primary, &
       old_fluid, fluid, transition, err)
    !! Check primary variables for eos_wge and make thermodynamic
    !! region transitions if needed.

    use fluid_module, only: fluid_type
    
    class(eos_wge_type), intent(in out) :: self
    PetscReal, intent(in) :: old_primary(self%num_primary_variables)
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in) :: old_fluid
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: old_region
    PetscReal :: saturation_pressure
    PetscReal :: water_pressure

    err = 0
    transition = PETSC_FALSE
    old_region = nint(old_fluid%region)

    if (old_region == 4) then  ! Two-phase
       associate (vapour_saturation => primary(2))

         if (vapour_saturation < 0._dp) then
            call self%transition_to_single_phase(old_primary, old_fluid, &
                 1, primary, fluid, transition, err)
         else if (vapour_saturation > 1._dp) then
            call self%transition_to_single_phase(old_primary, old_fluid, &
                 2, primary, fluid, transition, err)
         end if

     end associate
    else  ! Single-phase
       associate (pressure => primary(1), temperature => primary(2), &
            partial_pressure => primary(3))

         call self%thermo%saturation%pressure(temperature, &
              saturation_pressure, err)

         if (err == 0) then
            water_pressure = pressure - partial_pressure
            if (((old_region == 1) .and. (water_pressure < saturation_pressure)) .or. &
                 ((old_region == 2) .and. (water_pressure > saturation_pressure))) then
               call self%transition_to_two_phase(saturation_pressure, &
                    old_primary, old_fluid, primary, fluid, transition, err)
            end if
         end if

       end associate
    end if

  end subroutine eos_wge_transition

!------------------------------------------------------------------------

  subroutine eos_wge_bulk_properties(self, primary, fluid, err)
    !! Calculate fluid bulk properties from region and primary variables
    !! for non-isothermal water and non-condensible gas.

    use fluid_module, only: fluid_type

    class(eos_wge_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscInt :: region

    err = 0
    fluid%pressure = primary(1)
    region = nint(fluid%region)

    associate(partial_pressure => primary(3))
      fluid%partial_pressure(1) = fluid%pressure - partial_pressure
      fluid%partial_pressure(2) = partial_pressure
    end associate

    if (region == 4) then
       ! Two-phase
       call self%thermo%saturation%temperature( &
            fluid%partial_pressure(1), fluid%temperature, err)
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

  end subroutine eos_wge_bulk_properties

!------------------------------------------------------------------------

  subroutine eos_wge_phase_saturations(self, primary, fluid)
    !! Assigns fluid phase saturations from fluid region and primary variables.

    use fluid_module, only: fluid_type
    class(eos_wge_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
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

  end subroutine eos_wge_phase_saturations

!------------------------------------------------------------------------

  subroutine eos_wge_phase_properties(self, primary, rock, fluid, err)
    !! Calculate fluid phase properties from updated fluid region and primary variables.
    !! Bulk properties need to be calculated before calling this routine.

    use rock_module, only: rock_type
    use fluid_module, only: fluid_type

    class(eos_wge_type), intent(in out) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables) !! Primary thermodynamic variables
    type(rock_type), intent(in out) :: rock !! Rock object
    type(fluid_type), intent(in out) :: fluid !! Fluid object
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscInt :: p, phases
    PetscReal :: sl, xg
    PetscReal :: henrys_constant
    PetscReal :: water_properties(2), water_viscosity, water_enthalpy
    PetscReal :: relative_permeability(2), capillary_pressure(2)
    PetscReal :: water_pressure(2), energy_solution(2)
    PetscReal :: gas_properties(2), effective_gas_properties(2)

    err = 0
    phases = nint(fluid%phase_composition)

    sl = fluid%phase(1)%saturation
    relative_permeability = rock%relative_permeability%values(sl)
    capillary_pressure = 0._dp
    henrys_constant = 0._dp
    energy_solution = 0._dp

    call self%gas%properties(fluid%partial_pressure(2), fluid%temperature, &
         gas_properties, err)

    if (err == 0) then

       ! effective water pressure in each phase:
       water_pressure = [fluid%pressure, fluid%partial_pressure(1)]

       if (btest(phases, 0)) then
          capillary_pressure(1) = rock%capillary_pressure%value(sl, &
               fluid%temperature)
          call self%gas%henrys_constant(fluid%temperature, henrys_constant, err)
          if (err == 0) then
             call self%gas%energy_solution(fluid%temperature, henrys_constant, &
                  energy_solution(1), err)
          end if
       end if

       if (err == 0) then

          do p = 1, self%num_phases
             associate(phase => fluid%phase(p), region => self%thermo%region(p)%ptr)

               if (btest(phases, p - 1)) then

                  call region%properties([water_pressure(p), fluid%temperature], &
                       water_properties, err)

                  if (err == 0) then

                     call self%gas%effective_properties(gas_properties, p, &
                          effective_gas_properties)

                     associate(water_density => water_properties(1), &
                          water_internal_energy => water_properties(2), &
                          gas_density => effective_gas_properties(1), &
                          gas_enthalpy => effective_gas_properties(2))

                       call self%gas%mass_fraction(fluid%partial_pressure(2), &
                            fluid%temperature, p, gas_density, water_density, &
                            henrys_constant, xg, err)

                       if (err == 0) then

                          call region%viscosity(fluid%temperature, fluid%pressure, &
                               water_density, water_viscosity)
                          call self%gas%mixture_viscosity(water_viscosity, &
                               fluid%temperature, fluid%partial_pressure(2), xg, p, &
                               phase%viscosity, err)

                          if (err == 0) then
                             phase%density = water_density + gas_density
                             phase%mass_fraction = [1._dp - xg, xg]
                             phase%relative_permeability = relative_permeability(p)
                             phase%capillary_pressure =  capillary_pressure(p)
                             water_enthalpy = water_internal_energy &
                                  + water_pressure(p) / water_density
                             phase%specific_enthalpy = water_enthalpy * (1._dp - xg) &
                                  + (gas_enthalpy + energy_solution(p)) * xg
                             phase%internal_energy = phase%specific_enthalpy &
                                  - fluid%pressure / phase%density
                          else
                             exit
                          end if
                       else
                          exit
                       end if

                     end associate
                  else
                     exit
                  end if

               else
                  phase%density = 0._dp
                  phase%internal_energy = 0._dp
                  phase%specific_enthalpy = 0._dp
                  phase%relative_permeability = 0._dp
                  phase%capillary_pressure = 0._dp
                  phase%viscosity = 0._dp
                  phase%mass_fraction = 0._dp
               end if

             end associate
          end do
       end if

    end if

  end subroutine eos_wge_phase_properties

!------------------------------------------------------------------------

  subroutine eos_wge_primary_variables(self, fluid, primary)
    !! Determine primary variables from fluid properties.

    use fluid_module, only: fluid_type

    class(eos_wge_type), intent(in) :: self
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

    primary(3) = fluid%partial_pressure(2)

  end subroutine eos_wge_primary_variables

!------------------------------------------------------------------------

  subroutine eos_wge_check_primary_variables(self, fluid, &
       primary, changed, err)
    !! Check if primary variables are in acceptable bounds, and return error
    !! code accordingly. If gas partial pressure drops below zero, it
    !! is reset to zero, but an error is not raised.

    use fluid_module, only: fluid_type

    class(eos_wge_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    PetscBool, intent(out) :: changed
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: region
    PetscReal :: p
    PetscReal, parameter :: small = 1.e-6_dp

    changed = PETSC_FALSE
    err = 0

    associate (total_pressure => primary(1), partial_pressure => primary(3))

      if (total_pressure > 0._dp) then

         associate(max_partial_pressure => (1._dp - small) * total_pressure)
           if (partial_pressure > max_partial_pressure) then
              partial_pressure = max_partial_pressure
              changed = PETSC_TRUE
           else if (partial_pressure < 0._dp) then
              partial_pressure = 0._dp
              changed = PETSC_TRUE
           end if
         end associate

         p = total_pressure - partial_pressure
         if (p > 100.e6_dp) then
            err = 1
         else
            region = nint(fluid%region)
            if (region == 4) then
               associate (vapour_saturation => primary(2))
                 if ((vapour_saturation < -1._dp) .or. &
                      (vapour_saturation > 2._dp)) then
                    err = 1
                 end if
               end associate
            else
               associate (t => primary(2))
               if ((t < 0._dp) .or. (t > 800._dp)) then
                  err = 1
               end if
             end associate
          end if
       end if

    else
       err = 1
    end if

  end associate

  end subroutine eos_wge_check_primary_variables

!------------------------------------------------------------------------

  function eos_wge_scale_adaptive(self, primary, region) result(scaled_primary)
    !! Non-dimensionalise eos_wge primary variables by scaling. The
    !! first two variables (pressure and temperature or saturation)
    !! are scaled by fixed constants. The third variable, NCG partial
    !! pressure, is scaled adaptively by total pressure in the cell.

    class(eos_type), intent(in) :: self
    PetscReal, intent(in) :: primary(self%num_primary_variables)
    PetscInt, intent(in) :: region
    PetscReal :: scaled_primary(self%num_primary_variables)

    scaled_primary(1:2) = primary(1:2) / self%primary_scale(1:2, region)
    associate(scaled_partial_pressure => scaled_primary(3), &
         pressure => primary(1), partial_pressure => primary(3))
      scaled_partial_pressure = partial_pressure / pressure
    end associate

  end function eos_wge_scale_adaptive

!------------------------------------------------------------------------

  function eos_wge_unscale_adaptive(self, scaled_primary, region) result(primary)
    !! Re-dimensionalise eos_wge scaled primary variables.

    class(eos_type), intent(in) :: self
    PetscReal, intent(in) :: scaled_primary(self%num_primary_variables)
    PetscInt, intent(in) :: region
    PetscReal :: primary(self%num_primary_variables)

    primary(1:2) = scaled_primary(1:2) * self%primary_scale(1:2, region)
    associate(scaled_partial_pressure => scaled_primary(3), &
         pressure => primary(1), partial_pressure => primary(3))
      partial_pressure = scaled_partial_pressure * pressure
    end associate

  end function eos_wge_unscale_adaptive

!------------------------------------------------------------------------

  PetscReal function eos_wge_saturation_difference(x, context) result(dp)
    !! Returns difference between saturation pressure and water
    !! pressure at normalised point 0 <= x <= 1 along line between
    !! start and end primary variables.

    PetscReal, intent(in) :: x
    class(*), pointer, intent(in out) :: context
    ! Locals:
    PetscReal, allocatable :: var(:)
    PetscReal :: Ps
    PetscInt :: err

    select type (context)
    type is (primary_variable_interpolator_type)
       allocate(var(context%dim))
       var = context%interpolate_at_index(x)
       associate(P => var(1), T => var(2), Pg => var(3))
         call context%thermo%saturation%pressure(T, Ps, err)
         dp = P - Pg - Ps
       end associate
       deallocate(var)
    end select

  end function eos_wge_saturation_difference

!------------------------------------------------------------------------

end module eos_wge_module
