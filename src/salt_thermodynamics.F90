module salt_thermodynamics_module
  !! Module for salt thermodynamics routines.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use thermodynamics_module
  use utils_module, only: polynomial, newton1d

  implicit none
  private
  PetscReal, parameter :: salt_molecular_weight = 58.448_dp ! g/mol
  PetscReal, parameter :: halite_solubility_data(3) = &
       [2.6218e-1_dp, 7.2e-2_dp, 1.06_dp]
  PetscReal, parameter :: halite_solubility_two_phase_data(5) = &
       [0.2876823_dp, 0.30122157_dp, -0.39877656_dp, 0.31352381_dp, -0.09062578_dp]
  PetscReal, parameter :: brine_sat_pressure_a_data(4) = &
       [0._dp, 5.93582e-1_dp, -5.19386_dp, 1.23156_dp]
  PetscReal, parameter :: brine_sat_pressure_b_data(6) = &
       [0._dp, 1.15420_dp, 1.41254_dp, -1.92476_dp, &
       -1.70717_dp, 1.05390_dp]

  public :: halite_solubility, halite_solubility_two_phase
  public :: brine_saturation_pressure, brine_saturation_temperature

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

  subroutine halite_solubility_two_phase(pressure, thermo, solubility, err)
    !! Returns solubility of salt in liquid water on the brine
    !! saturation line, given the total pressure.

    PetscReal, intent(in) :: pressure
    class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics
    PetscReal, intent(out) :: solubility
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: xs
    PetscInt, parameter :: maxit = 100
    PetscReal, parameter :: tol = 1.e-10_dp
    PetscReal, parameter :: inc = 1.e-8_dp

    ! initial estimate from polynomial fit:
    xs = polynomial(halite_solubility_two_phase_data, pressure / 1.e7_dp)

    call newton1d(f, xs, inc, tol, maxit, err)
    solubility = xs

  contains

    PetscReal function f(x, err)
      PetscReal, intent(in) :: x
      PetscErrorCode, intent(out) :: err
      ! Locals:
      PetscReal :: temperature, solubility

      call brine_saturation_temperature(pressure, x, thermo, &
           temperature, err)
      if (err == 0) then
         call halite_solubility(temperature, solubility, err)
         f = x - solubility
      else
         f = -1._dp
      end if

    end function f

  end subroutine halite_solubility_two_phase

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

  subroutine brine_saturation_temperature(pressure, salt_mass_fraction, &
       thermo, saturation_temperature, err)
    !! Saturation temperature of brine at given pressure and salt mass
    !! fraction, using the specified water thermodynamics.

    PetscReal, intent(in) :: pressure !! Pressure (Pa)
    PetscReal, intent(in) :: salt_mass_fraction !! Salt mass fraction
    class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics
    PetscReal, intent(out) :: saturation_temperature !! Saturation temperature
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscReal :: t
    PetscInt, parameter :: maxit = 100
    PetscReal, parameter :: tol = 1.e-10_dp
    PetscReal, parameter :: inc = 1.e-8_dp

    ! Initial estimate from pure water thermodynamics:
    call thermo%saturation%temperature(pressure, t, err)
    if (err == 0) then
       call newton1d(f, t, inc, tol * pressure, maxit, err)
       saturation_temperature = t
    end if

  contains

    PetscReal function f(x, err)
      PetscReal, intent(in) :: x
      PetscErrorCode, intent(out) :: err
      ! Locals:
      PetscReal :: saturation_pressure

      call brine_saturation_pressure(x, salt_mass_fraction, thermo, &
          saturation_pressure, err)
      f = pressure - saturation_pressure

    end function f

  end subroutine brine_saturation_temperature

!------------------------------------------------------------------------

end module salt_thermodynamics_module
