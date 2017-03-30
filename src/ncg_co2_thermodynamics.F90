module ncg_co2_thermodynamics_module
  !! Module for thermodynamics of CO2 non-condensible gas.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use ncg_thermodynamics_module

  implicit none
  private

  PetscReal, parameter, public :: co2_molecular_weight = 44.01_dp ! g/mol

  type, public, extends(ncg_thermodynamics_type) :: ncg_co2_thermodynamics_type
     !! Type for CO2 NCG thermodynamics.
     private
     PetscReal :: viscosity_A(5) = [1357.8_dp, 4.9227_dp, &
          -2.9661e-3_dp, 2.8529e-6_dp, -2.1829e-9_dp]
     PetscReal :: viscosity_B(5) = [3918.9_dp, -35.984_dp, &
          2.5825e-1_dp, -7.1178e-4_dp, 6.9578e-7_dp]
     PetscReal :: viscosity_BA(5)
   contains
     private
     procedure, public :: init => ncg_co2_init
     procedure, public :: properties => ncg_co2_properties
     procedure, public :: effective_properties => ncg_co2_effective_properties
     procedure, public :: henrys_constant => ncg_co2_henrys_constant
     procedure, public :: energy_solution => ncg_co2_energy_solution
     procedure, public :: viscosity => ncg_co2_viscosity
     procedure, public :: mixture_viscosity => ncg_co2_mixture_viscosity
  end type ncg_co2_thermodynamics_type

contains

!------------------------------------------------------------------------

  subroutine ncg_co2_init(self)
    !! Initialises CO2 NCG thermodynamics object.

    class(ncg_co2_thermodynamics_type), intent(in out) :: self

    self%name = "CO2"
    self%molecular_weight = co2_molecular_weight
    self%viscosity_BA = self%viscosity_B - self%viscosity_A

  end subroutine ncg_co2_init

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

  subroutine ncg_co2_effective_properties(self, props, phase, &
       effective_props)
    !! Returns effective CO2 NCG properties for specified phase.
    !! CO2 density is treated as effectively zero in the liquid phase.

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: props(:) !! CO2 NCG properties (density, enthalpy)
    PetscInt, intent(in) :: phase !! Phase index
    PetscReal, intent(out) :: effective_props(:) !! Effective NCG properties

    effective_props = props

    if (phase == 1) then
       associate(co2_density => effective_props(1))
         co2_density = 0._dp
       end associate
    else
    end if

  end subroutine ncg_co2_effective_properties

!------------------------------------------------------------------------

  subroutine ncg_co2_henrys_constant(self, temperature, henrys_constant, err)
    !! Henry's constant for CO2 NCG.

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(out) :: henrys_constant !! Henry's constant
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscReal :: T2, T3, RKH

    associate(T => temperature)
      if (T <= 300._dp) then
         T2 = T * T
         T3 = T2 * T
         RKH = 783.666e5_dp + 19.6025e5_dp * T + 0.820574e5_dp * T2 &
              -T3 * (7.40674e2_dp - 2.18380_dp * T + 2.20999e-3_dp * T2)
         henrys_constant = 1.0_dp / RKH
         err = 0
      else
         err = 1
      end if
    end associate

  end subroutine ncg_co2_henrys_constant

!------------------------------------------------------------------------

  subroutine ncg_co2_energy_solution(self, temperature, energy_solution, err)
    !! Calculates enthalpy of CO2 dissolution in liquid.

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(out):: energy_solution !! Energy of solution
    PetscInt, intent(out) :: err     !! error code
    ! Locals:
    PetscReal :: T, T2, T3, T4

    err = 0

    T = 0.01_dp * temperature
    T2 = T * T
    T3 = T * T2
    T4 = T * T3
    energy_solution = -0.549491e6_dp + 0.456571e6_dp * T &
         - 0.070404e6_dp * T2 - 0.031035e6_dp * T3 + 0.014121e6_dp * T4

  end subroutine ncg_co2_energy_solution

!------------------------------------------------------------------------
  
  subroutine ncg_co2_viscosity(self, partial_pressure, temperature, &
       viscosity, err)
    !! Calculates viscosity for gas phase given partial pressure and
    !! temperature.

    use thermodynamics_module, only: region_type

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: partial_pressure !! CO2 partial pressure
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(out):: viscosity !! CO2 viscosity
    PetscInt, intent(out)  :: err !! Error code
    ! Locals:
    PetscReal :: T2, T3, T4, C(5)

    associate(pscale => partial_pressure * 1.0e-7_dp, T => temperature)
      if (pscale <= 1._dp) then
         C = self%viscosity_A + pscale * self%viscosity_BA
         T2 = T * T
         T3 = T * T2
         T4 = T2 * T2
         viscosity = 1.0e-8_dp * (C(1) + C(2) * T + &
              C(3) * T2 + C(4) * T3 + C(5) * T4)
         err = 0
      else
         err = 1
      end if
    end associate

  end subroutine ncg_co2_viscosity

!------------------------------------------------------------------------

  subroutine ncg_co2_mixture_viscosity(self, temperature, pressure, &
       partial_pressure, region, xg, total_density, &
       phase, viscosity, err)
    !! Calculates viscosity for water-CO2 mixture in the given phase.

    use thermodynamics_module, only: region_type

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(in) :: pressure !! Pressure
    PetscReal, intent(in) :: partial_pressure !! Partial pressure of CO2
    class(region_type), pointer :: region !! Thermodynamic region
    PetscReal, intent(in) :: xg !! CO2 mass fraction
    PetscReal, intent(in) :: total_density !! Total density
    PetscInt, intent(in)  :: phase !! Phase index
    PetscReal, intent(out):: viscosity !! Mixture viscosity
    PetscInt, intent(out)  :: err !! Error code
    ! Locals:
    PetscReal :: gas_viscosity, water_viscosity

    err = 0

    call region%viscosity(temperature, pressure, total_density, &
         water_viscosity)

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
