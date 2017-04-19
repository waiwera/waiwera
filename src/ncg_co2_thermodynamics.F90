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
  PetscReal, parameter :: viscosity_data(5, 6) = reshape([ &
       0._dp, 100.e5_dp, 150.e5_dp, 200.e5_dp, 300.e5_dp, &
       1357.8_dp, 3918.9_dp, 9660.7_dp, 1.31566e4_dp, 1.47968e4_dp, &
       4.9227_dp, -35.984_dp, -135.479_dp, -179.352_dp, -160.731_dp, &
       -2.9661e-3_dp, 0.25825_dp, 0.90087_dp, 1.12474_dp, 0.850257_dp, &
       2.8529e-6_dp, -7.1178e-4_dp, -2.4727e-3_dp, -2.98864e-3_dp, -1.99076e-3_dp, &
       -2.1829e-9_dp, 6.9578e-7_dp, 2.4156e-6_dp, 2.85911e-6_dp, 1.73423e-6_dp], &
       [5, 6])

  type, public, extends(ncg_thermodynamics_type) :: ncg_co2_thermodynamics_type
     !! Type for CO2 NCG thermodynamics.
     private
     type(interpolation_table_type) :: viscosity_table
   contains
     private
     procedure, public :: init => ncg_co2_init
     procedure, public :: destroy => ncg_co2_destroy
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
    call self%viscosity_table%init(viscosity_data)

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
    !! temperature. Formulation from Pritchett et al. (1982).

    use thermodynamics_module, only: region_type

    class(ncg_co2_thermodynamics_type), intent(in out) :: self
    PetscReal, intent(in) :: partial_pressure !! CO2 partial pressure
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(out):: viscosity !! CO2 viscosity
    PetscInt, intent(out)  :: err !! Error code
    ! Locals:
    PetscReal :: a(5)

    if (partial_pressure <= 300.e5_dp) then
       a = self%viscosity_table%interpolate(partial_pressure)
       associate(t => temperature)
         viscosity = a(1) + t * (a(2) + t * (a(3) + t * (a(4) + t * a(5))))
       end associate
       viscosity = 1.e-8_dp * viscosity
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
