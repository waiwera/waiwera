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
  PetscReal, parameter :: halite_enthalpy_data(3) = &
       [-0.120453_dp, 12.0453_dp, 1.95e-3_dp]
  PetscReal, parameter :: brine_sat_pressure_a_data(4) = &
       [0._dp, 5.93582e-1_dp, -5.19386_dp, 1.23156_dp]
  PetscReal, parameter :: brine_sat_pressure_b_data(6) = &
       [0._dp, 1.15420_dp, 1.41254_dp, -1.92476_dp, &
       -1.70717_dp, 1.05390_dp]
  PetscReal, parameter :: brine_viscosity_data(4) = &
       [1._dp, 0.0816_dp, 0.0122_dp, 1.28e-4_dp]
  PetscReal, parameter :: brine_critical_data(4) = &
       [-92.6824818148_dp, 0.43077335_dp, &
       -6.2561155e-4_dp, 3.6441625e-7_dp]
  PetscReal, parameter :: brine_critical_initial_data(4) = &
       [374.15000416_dp, 762.09488272_dp, 3391.76173286_dp, -6629.23422183_dp]

  public :: halite_solubility, halite_solubility_two_phase
  public :: halite_density, halite_enthalpy
  public :: brine_saturation_pressure, brine_saturation_temperature
  public :: brine_viscosity, brine_critical_temperature

contains

!------------------------------------------------------------------------
! Halite
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
    PetscReal, parameter :: ftol = 1.e-10_dp, xtol = 1.e-10_dp
    PetscReal, parameter :: inc = 1.e-8_dp

    ! initial estimate from polynomial fit:
    xs = polynomial(halite_solubility_two_phase_data, pressure / 1.e7_dp)

    call newton1d(f, xs, ftol, xtol, maxit, inc, err)
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

  subroutine halite_density(pressure, temperature, density, err)
    !! Returns density of halite as a function of pressure and
    !! temperature. From Silvester and Pitzer (1976).

    PetscReal, intent(in):: pressure !! Pressure
    PetscReal, intent(in):: temperature !! Temperature
    PetscReal, intent(out) :: density !! Halite density
    PetscErrorCode, intent(out) :: err !! Error code

    err = 0
    density = 2165._dp * exp(4.e-11_dp * pressure - 1.2e-4_dp * temperature)

  end subroutine halite_density

!------------------------------------------------------------------------

  subroutine halite_enthalpy(temperature, enthalpy, err)
    !! Returns enthalpy of halite as a function of temperature. From
    !! Silvester and Pitzer (1976).

    PetscReal, intent(in):: temperature !! Temperature
    PetscReal, intent(out) :: enthalpy !! Halite enthalpy
    PetscErrorCode, intent(out) :: err !! Error code

    err = 0
    enthalpy = 4184._dp *  polynomial(halite_enthalpy_data, temperature) / &
         salt_molecular_weight

  end subroutine halite_enthalpy

!------------------------------------------------------------------------
! Brine
!------------------------------------------------------------------------

  PetscReal function salt_mole_fraction(salt_mass_fraction)
    !! Returns salt mole fraction for given salt mass fraction.

    PetscReal, intent(in) :: salt_mass_fraction

    salt_mole_fraction = 1.e3_dp * salt_mass_fraction / (salt_molecular_weight * &
         (1._dp - salt_mass_fraction))

  end function salt_mole_fraction

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

    smol = salt_mole_fraction(salt_mass_fraction)
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
    PetscReal, parameter :: ftol = 1.e-10_dp, xtol = 1.e-10_dp
    PetscReal, parameter :: inc = 1.e-8_dp

    ! Initial estimate from pure water thermodynamics:
    call thermo%saturation%temperature(pressure, t, err)
    if (err == 0) then
       call newton1d(f, t, ftol * pressure, xtol, maxit, inc, err)
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


    PetscReal, intent(in) :: salt_mass_fraction
    PetscReal, intent(out) :: critical_temperature
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: t
    PetscReal :: poly(size(brine_critical_data))
    PetscInt, parameter :: maxit = 30
    PetscReal, parameter :: ftol = 1.e-10_dp, xtol = 1.e-10_dp



!------------------------------------------------------------------------

  subroutine brine_viscosity(pressure, temperature, density, &
       salt_mass_fraction, thermo, viscosity, err)
    !! Viscosity of brine as a function of pressure, temperature,
    !! density and salt mass fraction. Pure water viscosity is
    !! calculated using the specified thermodynamics. From Phillips,
    !! Igbene, Fair, Ozbek and Tavana (1981).

    PetscReal, intent(in) :: pressure !! Pressure
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(in) :: density !! Density
    PetscReal, intent(in) :: salt_mass_fraction !! Salt mass fraction
    class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics
    PetscReal, intent(out) :: viscosity !! Brine viscosity
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: smol, factor
    PetscReal :: water_viscosity

    err = 0
    smol = salt_mole_fraction(salt_mass_fraction)
    factor = polynomial(brine_viscosity_data, smol) + &
         6.29e-4_dp * temperature * (1._dp - exp(-0.7_dp * smol))
    call thermo%water%viscosity(temperature, pressure, &
         density, water_viscosity)
    viscosity = factor * water_viscosity

  end subroutine brine_viscosity

!------------------------------------------------------------------------

  subroutine brine_critical_temperature(salt_mass_fraction, &
       critical_temperature, err)
    !! Returns critical temperature of brine for given salt mass
    !! fraction.

    PetscReal, intent(in) :: salt_mass_fraction
    PetscReal, intent(out) :: critical_temperature
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: t
    PetscReal :: poly(size(brine_critical_data))
    PetscInt, parameter :: maxit = 30
    PetscReal, parameter :: ftol = 1.e-10_dp, xtol = 1.e-10_dp

    poly = brine_critical_data
    poly(1) = poly(1) - salt_mass_fraction * 100._dp
    t = polynomial(brine_critical_initial_data, salt_mass_fraction)
    call newton1d(poly, t, ftol, xtol, maxit, err)
    critical_temperature = max(t, 374.15_dp)

  end subroutine brine_critical_temperature

!------------------------------------------------------------------------

end module salt_thermodynamics_module
