module ncg_co2_thermodynamics_module
  !! Module for thermodynamics of CO2 non-condensible gas.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use ncg_thermodynamics_module
  use interpolation_module, only: interpolation_table_type

  implicit none
  private

  PetscReal, parameter, public :: co2_molecular_weight = 44.01_dp ! g/mol
  PetscReal, parameter :: henry_data(6) = [&
       0.783666_dp, 1.96025_dp, 8.20574_dp, &
       -7.40674_dp, 2.18380_dp, -0.220999_dp]
  PetscReal, parameter :: viscosity_data(5, 6) = reshape([ &
       0._dp, 10._dp, 15._dp, 20._dp, 30._dp, &
       1.3578_dp, 3.9189_dp, 9.6607_dp, 13.1566_dp, 14.7968_dp, &
       4.9227e-3_dp, -35.984e-3_dp, -135.479e-3_dp, -179.352e-3_dp, -160.731e-3_dp, &
       -2.9661e-6_dp, 0.25825e-3_dp, 0.90087e-3_dp, 1.12474e-3_dp, 0.850257e-3_dp, &
       2.8529e-9_dp, -7.1178e-7_dp, -2.4727e-6_dp, -2.98864e-6_dp, -1.99076e-6_dp, &
       -2.1829e-12_dp, 6.9578e-10_dp, 2.4156e-9_dp, 2.85911e-9_dp, 1.73423e-9_dp], &
       [5, 6])
  PetscReal, parameter :: tscale = 100._dp

  type, public, extends(ncg_thermodynamics_type) :: ncg_co2_thermodynamics_type
     !! Type for CO2 NCG thermodynamics.
     private
     type(interpolation_table_type) :: viscosity_table
     PetscReal :: henry_derivative_data(5)
   contains
     private
     procedure, public :: init => ncg_co2_init
     procedure, public :: destroy => ncg_co2_destroy
     procedure, public :: properties => ncg_co2_properties
     procedure, public :: henrys_constant => ncg_co2_henrys_constant
     procedure, public :: henrys_derivative => ncg_co2_henrys_derivative
     procedure, public :: viscosity => ncg_co2_viscosity
     procedure, public :: mixture_viscosity => ncg_co2_mixture_viscosity
  end type ncg_co2_thermodynamics_type

contains

!------------------------------------------------------------------------

  subroutine ncg_co2_init(self)
    !! Initialises CO2 NCG thermodynamics object.

    use utils_module, only: polynomial_derivative
    
    class(ncg_co2_thermodynamics_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: err

    self%name = "CO2"
    self%molecular_weight = co2_molecular_weight
    call self%viscosity_table%init(viscosity_data, err)
    self%henry_derivative_data = 10._dp * polynomial_derivative(henry_data)

  end subroutine ncg_co2_init

!------------------------------------------------------------------------

  subroutine ncg_co2_destroy(self)
    !! Destroys CO2 thermodynamics.

    class(ncg_co2_thermodynamics_type), intent(in out) :: self

    call self%viscosity_table%destroy()

  end subroutine ncg_co2_destroy

!------------------------------------------------------------------------  

  subroutine ncg_co2_properties(self, partial_pressure, temperature, &
       props, err)
    !! Calculates CO2 NCG density and enthalpy.

    use thermodynamics_module, only: tc_k

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: partial_pressure !! CO2 partial pressure
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(out):: props(:) !! Properties (density and enthalpy)
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscReal :: pp, tc, hci, vc

    err = 0

    associate(tk => temperature + tc_k, co2_density => props(1), &
         co2_enthalpy => props(2))

      pp = partial_pressure * 1.0e-6_dp
      tc = (0.01_dp * tk) ** 3.3333333333_dp
      hci = 1.667_dp + 0.001542_dp * tk - 0.7948_dp * log10(tk) - 41.35_dp / tk
      co2_enthalpy = 1.e6_dp * (hci - 0.3571_dp * pp * &
           (1._dp + 0.07576_dp * pp) / tc)
      vc = 0.00018882_dp * tk - pp * (0.0824_dp + 0.01249_dp * pp) / tc
      co2_density = pp / vc

    end associate

  end subroutine ncg_co2_properties

!------------------------------------------------------------------------

  subroutine ncg_co2_henrys_constant(self, temperature, henrys_constant, err)
    !! Henry's constant for CO2 NCG. The formulation is from
    !! Batistelli et al. (1997).

    use utils_module, only: polynomial

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(out) :: henrys_constant !! Henry's constant
    PetscErrorCode, intent(out) :: err !! Error code

    henrys_constant = 1.e8_dp * polynomial(henry_data, &
         temperature / tscale)
    err = 0

  end subroutine ncg_co2_henrys_constant

!------------------------------------------------------------------------

  subroutine ncg_co2_henrys_derivative(self, temperature, &
       henrys_constant, henrys_derivative, err)
    !! Returns derivative of natural logarithm of Henry's constant
    !! with respect to temperature.

    use utils_module, only: polynomial

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(in) :: henrys_constant !! Henry's constant
    PetscReal, intent(out) :: henrys_derivative !! Henry's derivative
    PetscErrorCode, intent(out) :: err !! Error code

    henrys_derivative = 1.e7_dp * &
         polynomial(self%henry_derivative_data, temperature / tscale) / &
         (henrys_constant * tscale)
    err = 0

  end subroutine ncg_co2_henrys_derivative

!------------------------------------------------------------------------
  
  subroutine ncg_co2_viscosity(self, partial_pressure, temperature, &
       viscosity, err)
    !! Calculates viscosity for gas phase given partial pressure and
    !! temperature. Formulation from Pritchett et al. (1982).

    use utils_module, only: polynomial

    class(ncg_co2_thermodynamics_type), intent(in out) :: self
    PetscReal, intent(in) :: partial_pressure !! CO2 partial pressure
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(out):: viscosity !! CO2 viscosity
    PetscInt, intent(out)  :: err !! Error code
    ! Locals:
    PetscReal :: coefs(5)
    PetscReal, parameter :: pscale = 1.e6_dp

    if (partial_pressure <= 300.e5_dp) then
       coefs = self%viscosity_table%interpolate(partial_pressure / pscale)
       viscosity = 1.e-5_dp * polynomial(coefs, temperature)
       err = 0
    else
       err = 1
    end if

  end subroutine ncg_co2_viscosity

!------------------------------------------------------------------------

  subroutine ncg_co2_mixture_viscosity(self, water_viscosity, &
       temperature, partial_pressure, xg, phase, viscosity, err)
    !! Calculates viscosity for water-CO2 mixture in the given phase.

    class(ncg_co2_thermodynamics_type), intent(in out) :: self
    PetscReal, intent(in) :: water_viscosity !! Viscosity of water
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(in) :: partial_pressure !! Partial pressure of CO2
    PetscReal, intent(in) :: xg !! CO2 mass fraction
    PetscInt, intent(in)  :: phase !! Phase index
    PetscReal, intent(out):: viscosity !! Mixture viscosity
    PetscInt, intent(out)  :: err !! Error code
    ! Locals:
    PetscReal :: gas_viscosity

    err = 0

    if (phase == 1) then
       viscosity = water_viscosity
    else
       call self%viscosity(partial_pressure, temperature, &
            gas_viscosity, err)
       if (err == 0) then
          viscosity = water_viscosity * (1._dp - xg) &
               + gas_viscosity * xg
       end if
    end if

  end subroutine ncg_co2_mixture_viscosity

!------------------------------------------------------------------------

end module ncg_co2_thermodynamics_module
