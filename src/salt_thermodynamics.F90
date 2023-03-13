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
  PetscReal, parameter :: halite_solubility_data(7) = &
       [0.2627980_dp, 3.130833e-2_dp, 2.136495_dp, &
       -9.371763_dp, 3.083588e1_dp, -3.959050e1_dp, &
       1.711302e1_dp]
  PetscReal, parameter :: halite_solubility_two_phase_data(5) = &
       [0.2876823_dp, 0.30122157_dp, -0.39877656_dp, 0.31352381_dp, -0.09062578_dp]
  PetscReal, parameter :: brine_sat_pressure_a_data(4) = &
       [0._dp, 5.93582e-1_dp, -5.19386_dp, 1.23156_dp]
  PetscReal, parameter :: brine_sat_pressure_b_data(6) = &
       [0._dp, 1.15420_dp, 1.41254_dp, -1.92476_dp, &
       -1.70717_dp, 1.05390_dp]
  PetscReal, parameter :: brine_viscosity_data(4) = &
       [1._dp, 0.0816_dp, 0.0122_dp, 1.28e-4_dp]

  public :: halite_solubility, halite_solubility_two_phase
  public :: halite_properties
  public :: brine_saturation_pressure, brine_saturation_temperature
  public :: brine_properties, brine_viscosity
  public :: salt_mole_fraction

contains

!------------------------------------------------------------------------
! Halite
!------------------------------------------------------------------------

  subroutine halite_solubility(temperature, solubility, err)
    !! Equilibrium solubility of salt in water as a function of
    !! temperature. From the regression given by Battistelli (2012),
    !! based on Driesner (2007).

    PetscReal, intent(in) :: temperature
    PetscReal, intent(out) :: solubility
    PetscErrorCode, intent(out) :: err

    if (0._dp <= temperature) then
       solubility = polynomial(halite_solubility_data, temperature * 1.e-3_dp)
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
       thermo, props, err)
    !! Returns properties (density and internal energy) of brine at
    !! the given pressure, temperature and salt mass fraction, using
    !! the specified pure water thermodynamics. From Driesner (2007).

    PetscReal, intent(in) :: pressure !! Pressure
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(in) :: salt_mass_fraction !! Salt mass fraction
    class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics
    PetscReal, intent(out):: props(:) !! Properties (density and internal energy)
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscReal :: brine_molecular_weight
    PetscReal :: f, xmol, xmol1, xmol12, tstar_v, ts
    PetscReal :: n1, n2, n10, n11, n12, n20, n21, n22, n23
    PetscReal :: n1x1, n2x1
    PetscReal :: props_star(2)
    PetscReal :: q1, q2, q10, q11, q12, q20, q21, q22, q23
    PetscReal :: q1x1, q2x1
    PetscReal :: tstar_h, hb
    PetscBool :: extrapolate

    err = 0

    associate(brine_density => props(1), brine_internal_energy => props(2), &
         pbar => pressure / 1.e5_dp)

      f = 1._dp / (salt_mass_fraction + (1._dp - salt_mass_fraction) * &
           salt_molecular_weight / water_molecular_weight)
      xmol = salt_mass_fraction * f
      xmol1 = 1._dp - xmol
      xmol12 = xmol1 * xmol1
      brine_molecular_weight = salt_molecular_weight * f

      ! density:
      n11 = -54.2958_dp - 45.7623_dp * exp(-9.44785e-4_dp * pbar)
      n21 = -2.6142_dp - 0.000239092_dp * pbar
      n22 = polynomial([0.0356828_dp, 4.37235e-3_dp, 2.0566e-3_dp], pbar / 1.e3_dp)

      n1x1 = polynomial([330.47_dp + 0.942876_dp * sqrt(pbar), 8.17193_dp, &
           -2.47556e-4_dp, 3.45052e-4_dp], pbar / 1.e2_dp)
      n2x1 = polynomial([-0.0370751_dp + 0.00237723_dp * sqrt(pbar), 5.42049e-1_dp, &
           5.84709e-1_dp, -5.99373e-1_dp], pbar / 1.e4_dp)

      n10 = n1x1
      n20 = 1._dp - n21 * sqrt(n22)
      n12 = -n11 - n10
      n23 = n2x1 - n20 - n21 * sqrt(1._dp + n22)

      n1 = n10 + n11 * xmol1 + n12 * xmol12
      n2 = n20 + n21 * sqrt(xmol + n22) + n23 * xmol

      tstar_v = n1 + n2 * temperature + deviation(pbar, temperature, xmol)
      if (pressure <= pcritical) then
         call thermo%saturation%temperature(pressure, ts, err)
         if (err == 0) then
            extrapolate = (tstar_v > ts)
         end if
      else
         extrapolate = PETSC_FALSE
      end if
      if (err == 0) then
         if (extrapolate) then
            brine_density = extrapolation(pressure, pbar, ts, tstar_v, &
                 brine_molecular_weight, err)
         else
            call thermo%water%properties([pressure, tstar_v], props_star, err)
            if (err == 0) then
               brine_density = props_star(1) * &
                    brine_molecular_weight / water_molecular_weight
            end if
         end if
      end if

      if (err == 0) then
         ! internal energy:

         q11 = -32.1724_dp + 0.0621255_dp * pbar
         q21 = polynomial([-1.69513_dp, -4.52781_dp, -6.04279_dp], pbar / 1.e4_dp)
         q22 = 0.0612567_dp + 1.88082e-5_dp * pbar

         q1x1 = polynomial([47.9048_dp, -9.36994_dp, 6.51059_dp], pbar / 1.e3_dp)
         q2x1 = polynomial([0.241022_dp, 3.45087e-1_dp, -4.28356e-1_dp], pbar / 1.e4_dp)

         q10 = q1x1
         q20 = 1._dp - q21 * sqrt(q22)
         q12 = -q11 - q10
         q23 = q2x1 - q20 - q21 * sqrt(1._dp + q22)

         q1 = q10 + q11 * xmol1 + q12 * xmol12
         q2 = q20 + q21 * sqrt(xmol + q22) + q23 * xmol

         tstar_h = q1 + q2 * temperature
         call thermo%water%properties([pressure, tstar_h], props_star, err)
         if (err == 0) then
            associate(dw => props_star(1), uw => props_star(2))
              hb = uw + pressure / dw
              brine_internal_energy = hb - pressure / brine_density
            end associate
         end if

      end if

    end associate

  contains

    PetscReal function deviation(pbar, temperature, xmol)
      !! Temperature deviation from eq. 14

      PetscReal, intent(in) :: pbar, temperature, xmol
      ! Locals:
      PetscReal :: pp
      PetscReal :: n30, n31, n300, n301, n302, n310, n311, n312

      pp = pbar + 472.051_dp
      n300 = 7.60664e6 / (pp * pp)
      n301 = -50.0_dp - 86.1446_dp * exp(-6.21128e-4_dp * pbar)
      n302 = 294.318_dp * exp(-5.66735e-3_dp * pbar)
      n310 = -0.0732761_dp * exp(-2.3772e-3_dp * pbar) - 5.2948e-5_dp * pbar
      n311 = -47.2747_dp + 24.3653_dp * exp(-1.25533e-3_dp * pbar)
      n312 = -0.278529_dp - 0.00081381_dp * pbar
      n30 = n300 * (exp(n301 * xmol) - 1._dp) + n302 * xmol
      n31 = n310 * exp(n311 * xmol) + n312 * xmol
      deviation = n30 * exp(n31 * temperature)

    end function deviation

!........................................................................

    PetscReal function extrapolation(pressure, pbar, ts, tstar_v, &
         brine_molecular_weight, err)
      !! Brine density from eq. 17 extrapolation

      PetscReal, intent(in) :: pressure, pbar, ts, tstar_v
      PetscReal, intent(in) :: brine_molecular_weight
      PetscErrorCode, intent(in out) :: err
      ! Locals:
      PetscReal :: vws, vws1, dvdt, logp, vb, ts2
      PetscReal :: o0, o1, o2
      PetscReal :: props_s(2), props_s1(2)
      PetscReal, parameter :: dt = 0.2_dp

      call thermo%water%properties([pressure, ts], props_s, err)
      if (err == 0) then
         associate(dws => props_s(1))
           vws = 1.e3_dp * water_molecular_weight / dws
         end associate
         call thermo%water%properties([pressure, ts - dt], props_s1, err)
         if (err == 0) then
            associate(dws1 => props_s1(1))
              vws1 = 1.e3_dp * water_molecular_weight / dws1
            end associate
            dvdt = (vws - vws1) / dt
            logp = log(pbar)
            o2 = polynomial([2.0125e-7_dp + 3.29977e-9_dp * exp(-4.31279_dp * logp), &
                 -1.17748e-7_dp, 7.58009e-8_dp], logp)
            ts2 = ts * ts
            o1 = dvdt - 3._dp * o2 * ts2
            o0 = vws - ts * (o1 + o2 * ts2)
            vb = polynomial([o0, o1, 0._dp, o2], tstar_v)
            extrapolation = 1.e3_dp * brine_molecular_weight / vb
         end if
      end if

    end function extrapolation

  end subroutine brine_properties

!------------------------------------------------------------------------

  subroutine brine_viscosity(temperature, pressure, &
       salt_mass_fraction, thermo, viscosity, err)
    !! Viscosity of brine as a function of temperature, pressure
    !! and salt mass fraction. Pure water viscosity is
    !! calculated using the specified thermodynamics. From Phillips,
    !! Igbene, Fair, Ozbek and Tavana (1981).

    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(in) :: pressure !! Pressure
    PetscReal, intent(in) :: salt_mass_fraction !! Salt mass fraction
    class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics
    PetscReal, intent(out) :: viscosity !! Brine viscosity
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: smol, factor
    PetscReal :: props(2), water_viscosity

    err = 0
    smol = salt_mole_fraction(salt_mass_fraction)
    factor = polynomial(brine_viscosity_data, smol) + &
         6.29e-4_dp * temperature * (1._dp - exp(-0.7_dp * smol))
    call thermo%water%properties([pressure, temperature], props, err)
    if (err == 0) then
       associate(water_density => props(1))
         call thermo%water%viscosity(temperature, pressure, &
              water_density, water_viscosity)
       end associate
       viscosity = factor * water_viscosity
    end if

  end subroutine brine_viscosity

!------------------------------------------------------------------------

end module salt_thermodynamics_module
