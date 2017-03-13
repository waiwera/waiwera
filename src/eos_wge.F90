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
    use fson_mpi_module, only: fson_get_mpi
    use logfile_module
    use thermodynamics_module

    class(eos_wge_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json !! JSON input object
    class(thermodynamics_type), intent(in), target :: thermo !! Thermodynamics object
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    procedure(root_finder_function), pointer :: f
    class(*), pointer :: pinterp
    PetscReal, parameter :: default_pressure = 1.0e5_dp
    PetscReal, parameter :: default_temperature = 20._dp ! deg C
    PetscReal, parameter :: default_gas_partial_pressure = 0._dp

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

    self%thermo => thermo

    ! Set up saturation line finder:
    allocate(primary_variable_interpolator_type :: &
         self%primary_variable_interpolator)
    call self%primary_variable_interpolator%init(self%num_primary_variables)
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
    self%thermo => null()

    call self%saturation_line_finder%destroy()
    call self%primary_variable_interpolator%destroy()
    deallocate(self%primary_variable_interpolator)

  end subroutine eos_we_destroy

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
    PetscReal :: old_saturation_pressure, factor
    PetscReal :: saturation_bound, xi
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

    call self%primary_variable_interpolator%assign(old_primary, primary)
    call self%primary_variable_interpolator%find(2, saturation_bound, xi, err)

    associate (pressure => primary(1), temperature => primary(2), &
         partial_pressure => primary(3), &
         interpolated_pressure => interpolated_primary(1), &
         interpolated_partial_pressure => interpolated_primary(3))

      if (err == 0) then

         interpolated_primary = self%primary_variable_interpolator%interpolate(xi)
         pressure = pressure_factor * interpolated_pressure
         partial_pressure = max(interpolated_partial_pressure, 0._dp)
         call self%thermo%saturation%temperature(interpolated_pressure - &
              interpolated_partial_pressure, temperature, err)
         if (err == 0) then
            fluid%region = dble(new_region)
            transition = PETSC_TRUE
         end if

      else

         call self%thermo%saturation%pressure(old_fluid%temperature, &
              old_saturation_pressure, err)
         if (err == 0) then
            pressure = pressure_factor * (old_saturation_pressure + &
                 partial_pressure)
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
         partial_pressure => primary(3), &
         interpolated_pressure => interpolated_primary(1), &
         interpolated_partial_pressure => interpolated_primary(3))

      call self%primary_variable_interpolator%assign(old_primary, primary)
      call self%saturation_line_finder%find()

      if (self%saturation_line_finder%err == 0) then
         xi = self%saturation_line_finder%root
         interpolated_primary = self%primary_variable_interpolator%interpolate(xi)
         pressure = interpolated_pressure
         partial_pressure = interpolated_partial_pressure
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
    PetscReal :: pressure_water

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
            pressure_water = pressure - partial_pressure
            if (((old_region == 1) .and. (pressure_water < saturation_pressure)) .or. &
                 ((old_region == 2) .and. (pressure_water > saturation_pressure))) then
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
    PetscReal :: pressure_water

    err = 0
    fluid%pressure = primary(1)
    region = nint(fluid%region)

    if (region == 4) then
       ! Two-phase
       associate(partial_pressure => primary(3))
         pressure_water = fluid%pressure - partial_pressure
         call self%thermo%saturation%temperature(pressure_water, &
              fluid%temperature, err)
       end associate
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
    PetscReal :: h2o_properties(2), sl, relative_permeability(2)
    PetscReal :: pressure_water, gas_properties(2), xg
    PetscReal :: h_solution

    err = 0
    phases = nint(fluid%phase_composition)

    sl = fluid%phase(1)%saturation
    relative_permeability = rock%relative_permeability%values(sl)

    associate(partial_pressure => primary(3))

      do p = 1, self%num_phases
         associate(phase => fluid%phase(p), region => self%thermo%region(p)%ptr)

           if (btest(phases, p - 1)) then

              if (p == 1) then
                 pressure_water = fluid%pressure
              else
                 pressure_water = fluid%pressure - partial_pressure
              end if

              call region%properties([pressure_water, fluid%temperature], &
                   h2o_properties, err)

              if (err == 0) then

                 associate(h2o_density => h2o_properties(1), &
                      h2o_internal_energy => h2o_properties(2))

                   call self%gas%properties(partial_pressure, fluid%temperature, &
                        p, h2o_density, gas_properties, xg, err)

                   if (err == 0) then

                      phase%mass_fraction(1) = 1._dp - xg
                      phase%mass_fraction(2) = xg

                      phase%relative_permeability = relative_permeability(p)

                      associate(gas_density => gas_properties(1), &
                           gas_enthalpy => gas_properties(2))

                        if (p == 2) then
                           phase%density = h2o_density + gas_density
                           h_solution = 0._dp
                           call self%gas%vapour_mixture_viscosity(fluid%pressure, &
                                fluid%temperature, partial_pressure, region, xg, &
                                h2o_density, phase%viscosity, err)
                           if (err > 0) then
                              exit
                           end if
                        else
                           phase%density = h2o_density
                           call self%gas%energy_solution(fluid%temperature, h_solution, err)
                           if (err > 0) then
                              exit
                           end if
                           call region%viscosity(fluid%temperature, fluid%pressure, &
                                phase%density, phase%viscosity)
                        end if

                        phase%specific_enthalpy = (h2o_internal_energy &
                             + pressure_water / h2o_density) * (1._dp - xg) &
                             + (gas_enthalpy + h_solution) * xg
                        phase%internal_energy = phase%specific_enthalpy &
                             - fluid%pressure / phase%density

                      end associate
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
              phase%viscosity = 0._dp
              phase%mass_fraction = 0._dp
           end if

         end associate
      end do

    end associate

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
    PetscReal :: xg, xmole, hc, partial_pressure
    PetscErrorCode :: err

    primary(1) = fluid%pressure

    region = nint(fluid%region)
    if (region == 4) then
       primary(2) = fluid%phase(2)%saturation
    else
       primary(2) = fluid%temperature
    end if

    xg = fluid%component_mass_fraction(2)
    xmole = self%gas%mole_fraction(xg)
    call self%gas%henrys_constant(fluid%temperature, hc, err)
    partial_pressure = xmole / hc
    primary(3) = partial_pressure

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

    changed = PETSC_FALSE
    err = 0

    associate (total_pressure => primary(1), partial_pressure => primary(3))
      p = total_pressure - partial_pressure
      if ((p < 0._dp) .or. (p > 100.e6_dp)) then
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
         if (err == 0) then
            if (partial_pressure > total_pressure) then
               err = 1
            else if (partial_pressure < 0._dp) then
               partial_pressure = 0._dp
               changed = PETSC_TRUE
            end if
         end if
      end if
    end associate

  end subroutine eos_wge_check_primary_variables

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
       var = context%interpolate(x)
       associate(P => var(1), T => var(2), Pg => var(3))
         call context%thermo%saturation%pressure(T, Ps, err)
         dp = P - Pg - Ps
       end associate
    end select
    deallocate(var)

  end function eos_wge_saturation_difference

!------------------------------------------------------------------------

end module eos_wge_module
