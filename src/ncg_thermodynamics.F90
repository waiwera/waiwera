module ncg_thermodynamics_module
  !! Module for abstract non-condensible gas thermodynamics type, from
  !! which specific NCG thermodynamics types can be derived.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use thermodynamics_module

  implicit none
  private

  PetscInt, parameter, public :: max_ncg_name_length = 8

  type, public, abstract :: ncg_thermodynamics_type
     !! Non-condensible gas thermodynamics type
     private
     character(max_ncg_name_length), public :: name !! NCG name
     PetscReal, public :: molecular_weight !! NCG molecular weight
     PetscReal, public :: deviation_factor = 1._dp !! Gas deviation (compressibility) factor, to account for non-ideal gas behaviour
   contains
     private
     procedure(ncg_init_procedure), public, deferred :: init
     procedure(ncg_properties_procedure), public, deferred :: properties
     procedure(ncg_henrys_constant_procedure), public, deferred :: henrys_constant
     procedure(ncg_henrys_constant_salt_procedure), public, deferred :: henrys_constant_salt
     procedure(ncg_henrys_derivative_procedure), public, deferred :: henrys_derivative
     procedure(ncg_henrys_derivative_salt_procedure), public, deferred :: henrys_derivative_salt
     procedure(ncg_viscosity_procedure), public, deferred :: viscosity
     procedure(ncg_mixture_viscosity_procedure), public, deferred :: mixture_viscosity
     procedure, public :: energy_solution_internal => ncg_energy_solution_internal
     procedure, public :: energy_solution => ncg_energy_solution
     procedure, public :: energy_solution_salt => ncg_energy_solution_salt
     procedure, public :: partial_pressure => ncg_partial_pressure
     procedure, public :: mass_fraction => ncg_mass_fraction
     procedure, public :: mole_to_mass_fraction => ncg_thermodynamics_mole_to_mass_fraction
     procedure, public :: mass_to_mole_fraction => ncg_thermodynamics_mass_to_mole_fraction
     procedure, public :: effective_properties => ncg_effective_properties
     procedure, public :: destroy => ncg_thermodynamics_destroy
   end type ncg_thermodynamics_type

  abstract interface

     subroutine ncg_init_procedure(self)
       !! Initialise NCG thermodynamics object.
       import :: ncg_thermodynamics_type
       class(ncg_thermodynamics_type), intent(in out) :: self
     end subroutine ncg_init_procedure

     subroutine ncg_properties_procedure(self, partial_pressure, &
          temperature, props, err)
       !! Calculate NCG fluid properties (density and enthalpy).
       import :: ncg_thermodynamics_type
       class(ncg_thermodynamics_type), intent(in) :: self
       PetscReal, intent(in) :: partial_pressure
       PetscReal, intent(in) :: temperature
       PetscReal, intent(out) :: props(:)
       PetscErrorCode, intent(out) :: err
     end subroutine ncg_properties_procedure

     subroutine ncg_henrys_constant_procedure(self, temperature, &
          henrys_constant, err)
       !! Calculate NCG Henry's constant, for calculating dissolution
       !! of gas into water.
       import :: ncg_thermodynamics_type
       class(ncg_thermodynamics_type), intent(in) :: self
       PetscReal, intent(in) :: temperature
       PetscReal, intent(out) :: henrys_constant
       PetscErrorCode, intent(out) :: err
     end subroutine ncg_henrys_constant_procedure

     subroutine ncg_henrys_constant_salt_procedure(self, temperature, &
          salt_mass_fraction, henrys_constant_0, henrys_constant, err)
       !! Calculate NCG Henry's constant, for calculating dissolution
       !! of gas into brine with given salt mass fraction.
       import :: ncg_thermodynamics_type
       class(ncg_thermodynamics_type), intent(in) :: self
       PetscReal, intent(in) :: temperature
       PetscReal, intent(in) :: salt_mass_fraction
       PetscReal, intent(out) :: henrys_constant_0 !! For zero salt
       PetscReal, intent(out) :: henrys_constant
       PetscErrorCode, intent(out) :: err
     end subroutine ncg_henrys_constant_salt_procedure

     subroutine ncg_henrys_derivative_procedure(self, temperature, &
          henrys_constant, henrys_derivative, err)
       !! Calculate derivative of the natural logarithm of Henry's
       !! constant with respect to temperature (used for computing
       !! energy of solution).
       import :: ncg_thermodynamics_type
       class(ncg_thermodynamics_type), intent(in) :: self
       PetscReal, intent(in) :: temperature
       PetscReal, intent(in) :: henrys_constant
       PetscReal, intent(out) :: henrys_derivative
       PetscErrorCode, intent(out) :: err
     end subroutine ncg_henrys_derivative_procedure

     subroutine ncg_henrys_derivative_salt_procedure(self, temperature, &
          salt_mass_fraction, henrys_constant_0, henrys_derivative, err)
       !! Calculate derivative of the natural logarithm of Henry's
       !! constant with respect to temperature (used for computing
       !! energy of solution), in brine with given salt mass fraction.
       import :: ncg_thermodynamics_type
       class(ncg_thermodynamics_type), intent(in) :: self
       PetscReal, intent(in) :: temperature
       PetscReal, intent(in) :: salt_mass_fraction
       PetscReal, intent(in) :: henrys_constant_0
       PetscReal, intent(out) :: henrys_derivative
       PetscErrorCode, intent(out) :: err
     end subroutine ncg_henrys_derivative_salt_procedure

     subroutine ncg_viscosity_procedure(self, partial_pressure, &
          temperature, viscosity, err)
       !! Calculate NCG viscosity.
       import :: ncg_thermodynamics_type
       import :: region_type
       class(ncg_thermodynamics_type), intent(in out) :: self
       PetscReal, intent(in)  :: partial_pressure
       PetscReal, intent(in)  :: temperature
       PetscReal, intent(out) :: viscosity
       PetscErrorCode, intent(out)  :: err
     end subroutine ncg_viscosity_procedure

     subroutine ncg_mixture_viscosity_procedure(self, water_viscosity, &
          temperature, partial_pressure, xg, phase, viscosity, err)
       !! Calculate water-NCG mixture viscosity for given phase.
       import :: ncg_thermodynamics_type
       import :: region_type
       class(ncg_thermodynamics_type), intent(in out) :: self
       PetscReal, intent(in)  :: water_viscosity
       PetscReal, intent(in)  :: temperature
       PetscReal, intent(in)  :: partial_pressure
       PetscReal, intent(in)  :: xg
       PetscInt, intent(in)   :: phase
       PetscReal, intent(out) :: viscosity
       PetscErrorCode, intent(out)  :: err
     end subroutine ncg_mixture_viscosity_procedure

  end interface

contains

!------------------------------------------------------------------------

  PetscReal function ncg_thermodynamics_mole_to_mass_fraction(self, xmole) &
       result(xg)
    !! Calculates NCG mass fraction from mole fraction.

    class(ncg_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: xmole !! NCG mole fraction
    ! Locals:
    PetscReal :: w

    w = xmole * self%molecular_weight
    xg = w / (w + (1._dp - xmole) * water_molecular_weight)

  end function ncg_thermodynamics_mole_to_mass_fraction

!------------------------------------------------------------------------

  PetscReal function ncg_thermodynamics_mass_to_mole_fraction(self, xg) &
       result(xmole)
    !! Calculates NCG mole fraction from mass fraction.

    class(ncg_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: xg !! NCG mass fraction
    ! Locals:
    PetscReal :: w

    w = xg / self%molecular_weight
    xmole = w / (w + (1._dp - xg) / water_molecular_weight)

  end function ncg_thermodynamics_mass_to_mole_fraction

!------------------------------------------------------------------------

  PetscReal function ncg_energy_solution_internal(self, temperature, &
       henrys_derivative) result(energy_solution)
    !! Calculates NCG energy of solution from the given temperature
    !! and Henry's derivative (from Himmelblau, 1959).

    use thermodynamics_module, only: tc_k, gas_constant

    class(ncg_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature
    PetscReal, intent(in) :: henrys_derivative

    associate(tk => temperature + tc_k)
      energy_solution = -1.e3_dp * gas_constant * tk * tk * &
           henrys_derivative / self%molecular_weight
    end associate

  end function ncg_energy_solution_internal

!------------------------------------------------------------------------

  subroutine ncg_energy_solution(self, temperature, henrys_constant, &
       energy_solution, err)
    !! Calculates NCG energy of solution from the given temperature
    !! and Henry's constant.

    use thermodynamics_module, only: tc_k, gas_constant

    class(ncg_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature
    PetscReal, intent(in) :: henrys_constant
    PetscReal, intent(out) :: energy_solution
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: henrys_derivative

    err = 0
    call self%henrys_derivative(temperature, henrys_constant, &
         henrys_derivative, err)
    if (err == 0) then
       energy_solution = self%energy_solution_internal(temperature, &
            henrys_derivative)
    end if

  end subroutine ncg_energy_solution

!------------------------------------------------------------------------

  subroutine ncg_energy_solution_salt(self, temperature, salt_mass_fraction, &
       henrys_constant_0, energy_solution, err)
    !! Calculates NCG energy of solution for brine from the given
    !! temperature, salt mass fraction and Henry's constant (for
    !! zero-salt case).

    use thermodynamics_module, only: tc_k, gas_constant

    class(ncg_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature
    PetscReal, intent(in) :: salt_mass_fraction
    PetscReal, intent(in) :: henrys_constant_0
    PetscReal, intent(out) :: energy_solution
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: henrys_derivative

    err = 0
    call self%henrys_derivative_salt(temperature, salt_mass_fraction, &
         henrys_constant_0, henrys_derivative, err)
    if (err == 0) then
       energy_solution = self%energy_solution_internal(temperature, &
            henrys_derivative)
    end if

  end subroutine ncg_energy_solution_salt

!------------------------------------------------------------------------

  PetscReal function ncg_partial_pressure(self, &
       temperature, total_density, xg) result(partial_pressure)
    !! Calculate NCG partial pressure from mass fraction xg.

    use thermodynamics_module, only: tc_k, gas_constant

    class(ncg_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in)  :: temperature !! Temperature
    PetscReal, intent(in)  :: total_density !! Total mixture density
    PetscReal, intent(in)  :: xg !! NCG mass fraction

    associate(tk => temperature + tc_k, gas_density => total_density * xg)
      partial_pressure = gas_density / self%molecular_weight * &
           (1.e3_dp * gas_constant * self%deviation_factor * tk)
    end associate

  end function ncg_partial_pressure

!------------------------------------------------------------------------

  subroutine ncg_mass_fraction(self, partial_pressure, &
       temperature, phase, gas_density, water_density, &
       henrys_constant, xg, err)
    !! Calculate NCG mass fraction from partial pressure.

    class(ncg_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: partial_pressure !! NCG partial pressure
    PetscReal, intent(in) :: temperature !! Temperature
    PetscInt, intent(in) :: phase !! Phase index
    PetscReal, intent(in) :: gas_density !! NCG density in this phase
    PetscReal, intent(in) :: water_density !! Water density in this phase
    PetscReal, intent(in) :: henrys_constant !! Henry's constant
    PetscReal, intent(out) :: xg !! NCG mass fraction
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscReal :: xmole
    PetscReal :: total_density
    PetscReal, parameter :: small = 1.e-30_dp

    err = 0
    if (phase == 1) then
       xmole =  partial_pressure / henrys_constant
       xg = self%mole_to_mass_fraction(xmole)
    else
       total_density = gas_density + water_density
       if (total_density < small) then
          xg = 0._dp
       else
          xg = gas_density / total_density
       end if
    end if

  end subroutine ncg_mass_fraction

!------------------------------------------------------------------------

  subroutine ncg_effective_properties(self, props, phase, &
       effective_props)
    !! Returns effective NCG properties for specified phase.
    !! NCG density is treated as effectively zero in the liquid phase.

    class(ncg_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: props(:) !! NCG properties (density, enthalpy)
    PetscInt, intent(in) :: phase !! Phase index
    PetscReal, intent(out) :: effective_props(:) !! Effective properties

    effective_props = props

    if (phase == 1) then
       associate(ncg_density => effective_props(1))
         ncg_density = 0._dp
       end associate
    end if

  end subroutine ncg_effective_properties

!------------------------------------------------------------------------

  subroutine ncg_thermodynamics_destroy(self)
    !! Destroys NCG thermodynamics. Dummy routine to be overridden by
    !! derived types.

    class(ncg_thermodynamics_type), intent(in out) :: self

    continue

  end subroutine ncg_thermodynamics_destroy

!------------------------------------------------------------------------

end module ncg_thermodynamics_module
