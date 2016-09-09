module eos_wge_module
  !! Equation of state for non-isothermal water and non-condensible gas.

  use kinds_module
  use eos_module
  use ncg_thermodynamics_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>
#include <petsc/finclude/petscsys.h>

  type, public, extends(eos_type) :: eos_wge_type
     !! Pure water, non-condensible gas and energy equation of state type.
     private
     class(ncg_thermodynamics_type), allocatable, public :: gas
   contains
     private
     procedure, public :: init => eos_wge_init
     procedure, public :: transition => eos_wge_transition
     procedure, public :: transition_to_single_phase => eos_wge_transition_to_single_phase
     procedure, public :: transition_to_two_phase => eos_wge_transition_to_two_phase
     procedure, public :: bulk_properties => eos_wge_bulk_properties
     procedure, public :: phase_properties => eos_wge_phase_properties
     procedure, public :: primary_variables => eos_wge_primary_variables
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

  end subroutine eos_wge_init

!------------------------------------------------------------------------

  subroutine eos_wge_transition_to_single_phase(self, old_fluid, &
       new_region, primary, fluid, transition, err)
    !! For eos_wge, make transition from two-phase to single-phase with
    !! specified region.

    use fluid_module, only: fluid_type

    class(eos_wge_type), intent(in out) :: self
    type(fluid_type), intent(in) :: old_fluid
    PetscInt, intent(in) :: new_region
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: old_saturation_pressure, factor
    PetscReal, parameter :: small = 1.e-6_dp

    err = 0

    associate (pressure => primary(1), temperature => primary(2), &
         partial_pressure => primary(3))

      call self%thermo%saturation%pressure(old_fluid%temperature, &
           old_saturation_pressure, err)

      if (err == 0) then

         if (new_region == 1) then
            factor = 1._dp + small
         else
            factor = 1._dp - small
         end if

         pressure = factor * (old_saturation_pressure + partial_pressure)
         temperature = old_fluid%temperature

         fluid%region = dble(new_region)
         transition = PETSC_TRUE

      end if

    end associate

  end subroutine eos_wge_transition_to_single_phase

!------------------------------------------------------------------------

  subroutine eos_wge_transition_to_two_phase(self, saturation_pressure, &
       old_fluid, primary, fluid, transition, err)
    !! For eos_wge, make transition from single-phase to two-phase.

    use fluid_module, only: fluid_type

    class(eos_wge_type), intent(in out) :: self
    PetscReal, intent(in) :: saturation_pressure
    type(fluid_type), intent(in) :: old_fluid
    PetscReal, intent(in out) :: primary(self%num_primary_variables)
    type(fluid_type), intent(in out) :: fluid
    PetscBool, intent(out) :: transition
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: old_region
    PetscReal, parameter :: small = 1.e-6_dp

    err = 0
    associate (pressure => primary(1), vapour_saturation => primary(2), &
         partial_pressure => primary(3))

      old_region = nint(old_fluid%region)

      pressure = saturation_pressure + partial_pressure

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

  subroutine eos_wge_transition(self, primary, old_fluid, fluid, &
          transition, err)
    !! Check primary variables for eos_wge and make thermodynamic
    !! region transitions if needed.

    use fluid_module, only: fluid_type
    
    class(eos_wge_type), intent(in out) :: self
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
            call self%transition_to_single_phase(old_fluid, 1, primary, &
                 fluid, transition, err)
         else if (vapour_saturation > 1._dp) then
            call self%transition_to_single_phase(old_fluid, 2, primary, &
                 fluid, transition, err)
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
               call self%transition_to_two_phase(saturation_pressure, old_fluid, &
                    primary, fluid, transition, err)
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
    PetscReal :: pressure, gas_properties(2), xg
    PetscReal :: gas_viscosity, h_solution
    
    err = 0
    phases = nint(fluid%phase_composition)

    sl = fluid%phase(1)%saturation
    relative_permeability = rock%relative_permeability%values(sl)

    associate(partial_pressure => primary(3))

      do p = 1, self%num_phases
         associate(phase => fluid%phase(p), region => self%thermo%region(p)%ptr)

           if (btest(phases, p - 1)) then

              if (p == 1) then
                 ! for liquid use total pressure and ignore ncg for density and internal_energy
                 pressure = fluid%pressure
              else
                 pressure = fluid%pressure - partial_pressure
              end if
            
              call region%properties([pressure, fluid%temperature], &
                   h2o_properties, err)

              if (err == 0) then

                 associate(h2o_density => h2o_properties(1), &
                      h2o_internal_energy => h2o_properties(2))

                   call self%gas%properties(partial_pressure, fluid%temperature, &
                        p, h2o_density, gas_properties, xg, err)

                   if (err == 0) then

                      associate(gas_density => gas_properties(1), &
                           gas_internal_energy => gas_properties(2))

                        phase%density = h2o_density + gas_density

                        if (p == 2) then
                           call self%gas%energy_solution(fluid%temperature, h_solution, err)
                           if (err > 0) then
                              exit
                           end if
                        else
                           h_solution = 0._dp
                        end if

                        phase%internal_energy = h2o_internal_energy
                        phase%specific_enthalpy = (h2o_internal_energy &
                             + fluid%pressure / h2o_density) * (1._dp - xg) &
                             + (gas_internal_energy + h_solution) * xg

                        phase%mass_fraction(1) = 1._dp - xg
                        phase%mass_fraction(2) = xg
                        phase%relative_permeability = relative_permeability(p)

                        call region%viscosity(fluid%temperature, fluid%pressure, &
                             phase%density, phase%viscosity)
                        if (p == 2) then
                           call self%gas%viscosity(partial_pressure, fluid%temperature, &
                                region, xg, phase%density, gas_viscosity, err)
                           if (err == 0) then
                              phase%viscosity = phase%viscosity * (1._dp - xg) &
                                   + gas_viscosity * xg
                           else
                              exit
                           end if
                        end if

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

    primary(1) = fluid%pressure

    region = nint(fluid%region)
    if (region == 4) then
       primary(2) = fluid%phase(2)%saturation
    else
       primary(2) = fluid%temperature
    end if

    primary(3) = 0._dp ! TODO: need to calculate

  end subroutine eos_wge_primary_variables

!------------------------------------------------------------------------

  PetscErrorCode function eos_wge_check_primary_variables(self, fluid, &
       primary) result(err)
    !! Check if primary variables are in acceptable bounds, and return error
    !! code accordingly.

    use fluid_module, only: fluid_type

    class(eos_wge_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscReal, intent(in) :: primary(self%num_primary_variables)
    ! Locals:
    PetscInt :: region
    PetscReal :: p
    PetscReal, parameter :: co2_partial_pressure_tol = 1.e-6_dp

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
            if ((partial_pressure < -co2_partial_pressure_tol) .or. &
                 (partial_pressure > total_pressure)) then
               err = 1
            end if
         end if
      end if
    end associate

  end function eos_wge_check_primary_variables

!------------------------------------------------------------------------

end module eos_wge_module
