module ncg_air_thermodynamics_module
  !! Module for thermodynamics of air non-condensible gas.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use ncg_thermodynamics_module

  implicit none
  private

  PetscReal, parameter, public :: air_molecular_weight = 28.96_dp ! g/mol
  PetscReal, parameter :: enthalpy_data(4) = [&
       1.20740_dp, 9.24502_dp, 0.115984_dp, -5.63568e-4_dp]
  PetscReal, parameter :: henry_weight(2) = [0.79_dp, 0.21_dp]
  PetscReal, parameter :: henry_p0(2) = [1.01325e5_dp, 1.e5_dp]
  PetscReal, parameter :: henry_data(2, 7) = reshape([&
       0.513726_dp, 0.26234_dp, &
       1.58603_dp, 0.610628_dp, &
       -5.9378e-1_dp, 7.00732e-1_dp, &
       -6.98282e-1_dp, -0.139299e1_dp, &
       5.10330e-1_dp, 7.13850e-1_dp, &
       -1.21388e-1_dp, -1.54216e-1_dp, &
       1.00041e-2_dp, 1.23190e-2_dp], &
       [2, 7])
  PetscReal, parameter :: tscale = 100._dp

  type, public, extends(ncg_thermodynamics_type) :: ncg_air_thermodynamics_type
     !! Type for air NCG thermodynamics.
     private
     PetscReal :: fair = 97.0_dp
     PetscReal :: fwat = 363.0_dp
     PetscReal :: cair = 3.617_dp
     PetscReal :: cwat = 2.655_dp
     PetscReal :: fmix, cmix
     PetscReal :: enthalpy_shift
     PetscReal :: henry_derivative_data(2, 6)
   contains
     private
     procedure, public :: init => ncg_air_init
     procedure, public :: properties => ncg_air_properties
     procedure, public :: henrys_constant => ncg_air_henrys_constant
     procedure, public :: henrys_constant_salt => ncg_air_henrys_constant_salt
     procedure, public :: henrys_derivative => ncg_air_henrys_derivative
     procedure, public :: henrys_derivative_salt => ncg_air_henrys_derivative_salt
     procedure, public :: viscosity => ncg_air_viscosity
     procedure, public :: mixture_viscosity => ncg_air_mixture_viscosity
  end type ncg_air_thermodynamics_type

contains

!------------------------------------------------------------------------

  subroutine ncg_air_init(self)
    !! Initialises air NCG thermodynamics object.

    use thermodynamics_module, only: ttriple, tc_k
    use utils_module, only: polynomial, polynomial_derivative

    class(ncg_air_thermodynamics_type), intent(in out) :: self

    self%name = "Air"
    self%molecular_weight = air_molecular_weight

    self%fmix = sqrt(self%fair * self%fwat)
    self%cmix = 0.5_dp * (self%cair + self%cwat)

    ! Enthalpy shift is calculated so that enthalpy at triple point of
    ! water is zero:
    associate(tk => ttriple + tc_k)
      self%enthalpy_shift = polynomial(enthalpy_data, tk / tscale)
    end associate

    self%henry_derivative_data = polynomial_derivative(henry_data)

  end subroutine ncg_air_init

!------------------------------------------------------------------------  

  subroutine ncg_air_properties(self, partial_pressure, temperature, &
       props, err)
    !! Calculates air NCG density and enthalpy. Density is calculated
    !! from the real gas law. Enthalpy is from Irvine and Liley
    !! (1984), "Steam and gas tables with computer equations".

    use thermodynamics_module, only: tc_k, gas_constant
    use utils_module, only: polynomial

    class(ncg_air_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: partial_pressure !! Air partial pressure
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(out):: props(:) !! Properties (density and enthalpy)
    PetscErrorCode, intent(out) :: err !! Error code

    err = 0
    associate(tk => temperature + tc_k, air_density => props(1), &
         air_enthalpy => props(2))
      air_density = partial_pressure * self%molecular_weight / &
           (1.e3_dp * gas_constant * self%deviation_factor * tk)
      air_enthalpy = 1.e4_dp * (polynomial(enthalpy_data, tk / tscale) - &
           self%enthalpy_shift)
    end associate

  end subroutine ncg_air_properties

!------------------------------------------------------------------------

  subroutine ncg_air_henrys_constant(self, temperature, &
       henrys_constant, err)
    !! Henry's constant for air NCG. The formulation is based on
    !! D'Amore and Truesdell (1988), Cramer (1982) and Cygan (1991).

    use utils_module, only: polynomial

    class(ncg_air_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(out) :: henrys_constant !! Henry's constant
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscReal :: hinv(2)

    err = 0
    hinv = polynomial(henry_data, temperature / tscale)
    henrys_constant = 1.e5_dp * sum(henry_weight * henry_p0 * hinv)

  end subroutine ncg_air_henrys_constant

!------------------------------------------------------------------------

  subroutine ncg_air_henrys_constant_salt(self, temperature, &
       salt_mass_fraction, henrys_constant_0, henrys_constant, err)
    !! Henry's constant for air NCG in brine.

    use utils_module, only: polynomial
    use salt_thermodynamics_module, only: salt_mole_fraction

    class(ncg_air_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(in) :: salt_mass_fraction !! Salt mass fraction
    PetscReal, intent(out) :: henrys_constant_0 !! Henry's constant for zero salt
    PetscReal, intent(out) :: henrys_constant !! Henry's constant
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    ! PetscReal :: m, kb

    err = 0
    call self%henrys_constant(temperature, henrys_constant_0, err)
    henrys_constant = henrys_constant_0
    ! if (err == 0) then
    !    m = salt_mole_fraction(salt_mass_fraction)
    !    kb = polynomial(henry_salt_data, temperature / tscale)
    !    henrys_constant = henrys_constant * 10._dp ** (m * kb)
    ! end if

  end subroutine ncg_air_henrys_constant_salt

!------------------------------------------------------------------------

  subroutine ncg_air_henrys_derivative(self, temperature, &
       henrys_constant, henrys_derivative, err)
    !! Returns derivative of natural logarithm of Henry's constant
    !! with respect to temperature.

    use utils_module, only: polynomial

    class(ncg_air_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(in) :: henrys_constant !! Henry's constant
    PetscReal, intent(out) :: henrys_derivative !! Henry's derivative
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscReal :: dhinv(2)

    err = 0
    dhinv = polynomial(self%henry_derivative_data, temperature / tscale)
    henrys_derivative = 1.e5_dp * sum(henry_weight * henry_p0 * dhinv) &
         / (henrys_constant * tscale)

  end subroutine ncg_air_henrys_derivative

!------------------------------------------------------------------------

  subroutine ncg_air_henrys_derivative_salt(self, temperature, &
       salt_mass_fraction, henrys_constant_0, henrys_derivative, err)
    !! Returns derivative of natural logarithm of Henry's constant
    !! with respect to temperature, for brine with given salt mass
    !! fraction.

    use utils_module, only: polynomial
    use salt_thermodynamics_module, only: salt_mole_fraction

    class(ncg_air_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(in) :: salt_mass_fraction !! Salt mass fraction
    PetscReal, intent(in) :: henrys_constant_0 !! Henry's constant for zero salt
    PetscReal, intent(out) :: henrys_derivative !! Henry's derivative
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    ! PetscReal :: m, dkb_dt

    err = 0
    call self%henrys_derivative(temperature, henrys_constant_0, &
         henrys_derivative, err)
    ! if (err == 0) then
    !    m = salt_mole_fraction(salt_mass_fraction)
    !    dkb_dt = polynomial(self%henry_salt_derivative_data, &
    !         temperature / tscale) / tscale
    !    henrys_derivative = henrys_derivative + log(10._dp) * m * dkb_dt
    ! end if

  end subroutine ncg_air_henrys_derivative_salt

!------------------------------------------------------------------------
  
  subroutine ncg_air_viscosity(self, partial_pressure, temperature, &
       viscosity, err)
    !! Viscosity for air, given partial pressure and temperature.
    !! This is not used.

    class(ncg_air_thermodynamics_type), intent(in out) :: self
    PetscReal, intent(in) :: partial_pressure !! Air partial pressure
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(out):: viscosity !! Air viscosity
    PetscInt, intent(out)  :: err !! Error code

    err = 0
    viscosity = 0._dp

  end subroutine ncg_air_viscosity

!------------------------------------------------------------------------

  subroutine ncg_air_mixture_viscosity(self, water_viscosity, &
       temperature, partial_pressure, xg, phase, viscosity, err)
    !! Calculates viscosity for water-air mixture in the given phase.
    !!
    !! Uses a modified version of a formulation based on kinetic
    !! gas theory, as given by J.O. Hirschfelder, C.F. Curtiss, and
    !! R.B. Bird, Molecular Theory of Gases and Liquids, John Wiley
    !! & Sons, 1954, pp. 528-530.
    !!
    !! The modification made to the Hirschfelder et al. expressions is
    !! that for vapour viscosity accurate (empirical) values are used,
    !! rather than the first order expression of kinetic theory.
    !!
    !! The formulation matches experimental data on viscosities of
    !! vapour-air mixtures in the temperature range from 100 to 150
    !! deg. C, for all compositions, to better than 4%.
    
    use thermodynamics_module, only: water_molecular_weight, tc_k

    class(ncg_air_thermodynamics_type), intent(in out) :: self
    PetscReal, intent(in) :: water_viscosity !! Viscosity of water
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(in) :: partial_pressure !! Air partial pressure (not used)
    PetscReal, intent(in) :: xg !! Air mass fraction
    PetscInt, intent(in)  :: phase !! Phase index
    PetscReal, intent(out):: viscosity !! Mixture viscosity
    PetscInt, intent(out)  :: err !! Error code
    ! Locals:
    PetscReal :: x1, x2, ard
    PetscReal :: e, g, h
    PetscReal :: ome1, ome3, rm3
    PetscReal :: trd1, trd3, vis1, vis2, vis3, z1, z2, z3

    err = 0

    if (phase == 1) then

       viscosity = water_viscosity

    else

       associate(rm1 => self%molecular_weight, rm2 => water_molecular_weight)

         x1 = self%mass_to_mole_fraction(xg)
         x2 = 1._dp - x1

         associate (tk => temperature + tc_k)
           trd1 = tk / self%fair
           trd3 = tk / self%fmix
         end associate
         ome1 = (1.188_dp - 0.051_dp * trd1) / trd1
         ome3 = (1.48_dp  - 0.412_dp * log(trd3)) / trd3
         ard = 1.095_dp / trd3
         rm3 = 2._dp * rm1 * rm2 / (rm1 + rm2)
         vis1 = covis(trd1, self%cair, ome1, rm1, self%fair)
         vis2 = 10._dp * water_viscosity
         vis3 = covis(trd3, self%cmix, ome3, rm3, self%fmix)
         z1 = x1 * x1 / vis1 + 2._dp * x2 * x1 / vis3 + x2 * x2 / vis2
         g = x1 * x1 * rm1 / rm2
         h = x2 * x2 * rm2 / rm1
         e = (2._dp * x1 * x2 * rm1 * rm2 / (rm3 * rm3)) * &
              vis3 / (vis1 * vis2)
         z2 = 0.6_dp * ard * (g / vis1 + e + h / vis2)
         z3 = 0.6_dp * ard * (g + e * (vis1 + vis2) - 2._dp * x1 * x2 + h)
         viscosity = 0.1_dp * (1._dp + z3) / (z1 + z2)

       end associate

    end if

  contains

    PetscReal function covis(trd, c, ome, rm, f)
      ! Coefficient of viscosity to the first approximation.
      PetscReal, intent(in)  :: trd, c, ome, rm, f
      covis = 266.93e-7_dp * sqrt(rm * trd * f) / (c * c * ome * trd)
    end function covis

  end subroutine ncg_air_mixture_viscosity

!------------------------------------------------------------------------

end module ncg_air_thermodynamics_module
