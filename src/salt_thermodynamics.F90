module salt_thermodynamics_module
  !! Module for salt thermodynamics routines.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use thermodynamics_module
  use utils_module, only: polynomial

  implicit none
  private
  PetscReal, parameter :: salt_molecular_weight = 58.448_dp ! g/mol
  PetscReal, parameter :: halite_solubility_data(3) = &
       [2.6218e-1_dp, 7.2e-2_dp, 1.06_dp]
  PetscReal, parameter :: brine_sat_pressure_a_data(4) = &
       [0._dp, 5.93582e-1_dp, -5.19386_dp, 1.23156_dp]
  PetscReal, parameter :: brine_sat_pressure_b_data(6) = &
       [0._dp, 1.15420_dp, 1.41254_dp, -1.92476_dp, &
       -1.70717_dp, 1.05390_dp]

  public :: halite_solubility, brine_saturation_pressure

contains

!------------------------------------------------------------------------

  subroutine halite_solubility(temperature, solubility, err)
    !! Equilibrium solubility of salt in water as a function of
    !! temperature. From Chou, I.M. (1987).

    PetscReal, intent(in) :: temperature
    PetscReal, intent(out) :: solubility
    PetscErrorCode, intent(out) :: err

    if ((0._dp <= temperature) .and. (temperature <= 382._dp)) then
       solubility = polynomial(halite_solubility_data, temperature / 1.e3_dp)
       err = 0
    else
       solubility = 0._dp
       err = 1
    end if

  end subroutine halite_solubility

!------------------------------------------------------------------------

  subroutine brine_saturation_pressure(temperature, salt_mass_fraction, &
       thermo, saturation_pressure, err)
    !! Saturation pressure of brine at given temperature and salt mass
    !! fraction, using the specified pure water thermodynamics. Based
    !! on Haas (1976).

    PetscReal, intent(in) :: temperature !! Temperature (\(^\circ C\))
    PetscReal, intent(in) :: salt_mass_fraction !! Salt mass fraction
    class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics
    PetscReal, intent(out) :: saturation_pressure !! Saturation pressure
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscReal :: smol
    PetscReal :: a, b, tk, t_effective

    smol = 1.e3_dp * salt_mass_fraction / (salt_molecular_weight * &
         (1._dp - salt_mass_fraction))
    a = 1._dp + 1.e-5_dp * polynomial(brine_sat_pressure_a_data, smol)
    b = 1.e-5_dp * polynomial(brine_sat_pressure_b_data, 0.1_dp * smol)

    tk = temperature + tc_k
    t_effective = exp(log(tk)/(a + b * tk)) - tc_k

    call thermo%saturation%pressure(t_effective, saturation_pressure, err)

  end subroutine brine_saturation_pressure

!------------------------------------------------------------------------

end module salt_thermodynamics_module
