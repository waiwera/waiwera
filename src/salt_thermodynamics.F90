module salt_thermodynamics_module
  !! Module for salt thermodynamics routines.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use thermodynamics_module
  use utils_module, only: polynomial, newton1d

  implicit none
  private
  PetscReal, parameter :: salt_molecular_weight = 58.443_dp ! g/mol
  PetscReal, parameter :: halite_density_data(3) = &
       [2.1704e3_dp, -2.4599e-1_dp, -9.5797e-5_dp]
  PetscReal, parameter :: halite_enthalpy_data(4) = &
       [-5.615174e5_dp, 8.766380e2_dp, 6.413881e-2_dp, 8.810112e-5_dp]
  PetscReal, parameter :: halite_solubility_data(3) = &
       [2.6218e-1_dp, 7.2e-2_dp, 1.06_dp]
  PetscReal, parameter :: halite_solubility_two_phase_data(5) = &
       [0.2876823_dp, 0.30122157_dp, -0.39877656_dp, 0.31352381_dp, -0.09062578_dp]
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
  PetscReal, parameter :: brine_enthalpy_data(4) = &
       [0._dp, -25.9293_dp, 0.16792_dp, -8.3624e-4_dp]
  PetscReal, parameter :: brine_enthalpy_data2(3,4) = &
       reshape([ &
       -9.6336e3_dp, -4.0800e3_dp, 2.8649e2_dp, &
       1.6658e2_dp, 6.8577e1_dp, -4.6856_dp, &
       -0.90963_dp, -0.36524_dp, 2.49667e-2_dp, &
       1.7965e-3_dp, 7.1924e-4_dp, -4.900e-5_dp], [3,4])

  public :: halite_solubility, halite_solubility_two_phase
  public :: halite_properties
  public :: brine_saturation_pressure, brine_saturation_temperature
  public :: brine_viscosity, brine_critical_temperature
  public :: brine_properties

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

  subroutine halite_properties(pressure, temperature, props, err)
    !! Returns properties (density and internal energy) of halite as a
    !! function of pressure and temperature. From Driesner (2007).

    PetscReal, intent(in):: pressure !! Pressure
    PetscReal, intent(in):: temperature !! Temperature
    PetscReal, intent(out):: props(:) !! Properties (density and internal energy)
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscReal, parameter :: l3 = 5.727e-3_dp, l4 = 2.715e-3_dp, l5 = 733.4_dp
    PetscReal :: density0, l, enthalpy_1bar, enthalpy

    err = 0
    associate(pbar => pressure / 1.e5_dp, &
         density => props(1), internal_energy => props(2))

      density0 = polynomial(halite_density_data, temperature)
      l = l3 + l4 * exp(temperature / l5)
      density = density0 + l * pbar

      enthalpy_1bar = polynomial(halite_enthalpy_data, temperature)
      enthalpy = enthalpy_1bar + 44.14_dp * (pbar - 1._dp)
      internal_energy = enthalpy - pressure / density

    end associate

  end subroutine halite_properties

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

  subroutine brine_properties(pressure, temperature, salt_mass_fraction, &
       thermo, props, water_density, err)
    !! Returns properties (density and internal energy) of brine at
    !! the given pressure, temperature and salt mass fraction, using
    !! the specified pure water thermodynamics. Also returns density
    !! of pure water.

    PetscReal, intent(in) :: pressure !! Pressure
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(in) :: salt_mass_fraction !! Salt mass fraction
    class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics
    PetscReal, intent(out):: props(:) !! Properties (density and internal energy)
    PetscReal, intent(out):: water_density !! Pure water density
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscReal :: Pws, Ps, smol, v0, sk, fi, sat_brine_density
    PetscReal :: hws, delp, Ps1, hws1, factor, tc, hbs
    PetscReal :: tau, xmol, c, hsalt, dhsalt, dppsi
    PetscReal :: sat_water_props(2), sat_water_props1(2)
    PetscReal :: b(4), brine_enthalpy
    PetscInt :: i

    err = 0

    associate(brine_density => props(1), brine_internal_energy => props(2))

      call thermo%saturation%pressure(temperature, Pws, err)
      if (err == 0) then
         call thermo%water%properties([Pws, temperature], sat_water_props, err)
         if (err == 0) then
            associate(dws => sat_water_props(1), uws => sat_water_props(2))

              hws = uws + Pws / dws

              ! From Haas (1976):
              smol = salt_mole_fraction(salt_mass_fraction)
              v0 = 1.e3_dp / dws
              sk = (-13.644_dp + 13.97_dp * v0) * (3.1975_dp / (3.1975_dp - v0))**2
              fi = -167.219_dp + 448.55_dp * v0 - 261.07_dp * v0 * v0 + sk*sqrt(smol)
              sat_brine_density = (1.e3_dp + smol * salt_molecular_weight) / &
                   (1000._dp * v0 + smol * fi) * 1.e3_dp

              call brine_saturation_pressure(temperature, salt_mass_fraction, &
                   thermo, Ps, err)
              if (err == 0) then

                 ! Compressibility from Andersen et al. (1992)
                 call brine_critical_temperature(salt_mass_fraction, tc, err)
                 if (err == 0) then

                    tau = 1._dp - (temperature + 273.15_dp) / (tc + 273.15_dp)
                    xmol = smol / (1.e3_dp / water_molecular_weight + smol)
                    c = -1.14e-6_dp / (tau**1.25_dp - 5.6_dp * xmol**1.5_dp + 0.005_dp)
                    dppsi = 14.50377e-5_dp * (pressure - Ps)
                    brine_density = sat_brine_density / (1._dp + c * dppsi)

                    ! Enthalpy of vapour saturated brine, from Michaelides (1982):
                    hsalt = 4.184_dp  * polynomial(brine_enthalpy_data, temperature) / &
                         salt_molecular_weight
                    do i = 1, 4
                       b(i) = polynomial(brine_enthalpy_data2(:, i), smol)
                    end do
                    dhsalt = polynomial(b, temperature) * (4.184_dp / &
                         (1.e3_dp + salt_molecular_weight *smol))
                    hbs = (1._dp - salt_mass_fraction) * hws + &
                         1.e3_dp * (salt_mass_fraction * hsalt + smol * dhsalt)

                    ! Compressed brine enthalpy:
                    delp = max(pressure - Pws, 1.e3_dp)
                    Ps1 = Pws + delp
                    call thermo%water%properties([Ps1, temperature], sat_water_props1, err)
                    if (err == 0) then
                       water_density = sat_water_props1(1)
                       associate(uws1 => sat_water_props1(2))
                         hws1 = uws1 + Ps1 / water_density
                         factor = (hws1 - hws) / hws / delp
                         brine_enthalpy = hbs * (1._dp + factor * (pressure - Ps))
                         brine_internal_energy = brine_enthalpy - pressure / brine_density
                       end associate
                    end if

                 end if

              end if
            end associate
         end if
      end if
    end associate

  end subroutine brine_properties

!------------------------------------------------------------------------

  subroutine brine_viscosity(temperature, pressure, water_density, &
       salt_mass_fraction, thermo, viscosity, err)
    !! Viscosity of brine as a function of temperature, pressure,
    !! water density and salt mass fraction. Pure water viscosity is
    !! calculated using the specified thermodynamics. From Phillips,
    !! Igbene, Fair, Ozbek and Tavana (1981).

    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(in) :: pressure !! Pressure
    PetscReal, intent(in) :: water_density !! Density of pure water
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
         water_density, water_viscosity)
    viscosity = factor * water_viscosity

  end subroutine brine_viscosity

!------------------------------------------------------------------------

  subroutine brine_critical_temperature(salt_mass_fraction, &
       critical_temperature, err)
    !! Returns critical temperature of brine for given salt mass
    !! fraction. From Sourirayan and Kennedy (1962).

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
