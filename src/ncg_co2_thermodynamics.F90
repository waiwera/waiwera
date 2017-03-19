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
     procedure, public :: henrys_constant => ncg_co2_henrys_constant
     procedure, public :: energy_solution => ncg_co2_energy_solution
     procedure, public :: viscosity => ncg_co2_viscosity
     procedure, public :: vapour_mixture_viscosity => ncg_co2_vapour_mixture_viscosity
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
       phase, water_density, props, xg, err)
    !! Calculates density and internal energy of mixture of a fluid and and a ncg
    !! pressure (Pa) and temperature (deg C).

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: partial_pressure !! CO2 partial pressure
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(in) :: water_density !! Water density
    PetscInt, intent(in)  :: phase !! Fluid phase
    PetscReal, intent(out):: props(:) !! Properties (density and internal energy)
    PetscReal, intent(out):: xg !! Mass fraction of the ncg in this phase
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscReal :: total_density, hc, xmole
    PetscReal, parameter :: small = 1.e-30_dp

    err = 0

    associate(co2_density => props(1), co2_enthalpy => props(2))

      call co2_rho_h(partial_pressure, temperature, co2_density, &
           co2_enthalpy, err)

      if (err == 0) then
         if (phase == 1) then
            ! liquid
            co2_density = 0._dp    ! not used for mixture density
            call self%henrys_constant(temperature, hc, err)
            if (err == 0) then
               xmole = hc * partial_pressure
               xg = self%mass_fraction(xmole)
            end if
         else
            ! vapour
            total_density = co2_density + water_density
            if (total_density < small) then
               xg = 0._dp
            else
               xg = co2_density / total_density
            end if
         end if
      end if

    end associate

  contains

    subroutine co2_rho_h(PP, T, DC, HC, err)

      use thermodynamics_module, only: tc_k

      PetscReal, intent(in)  :: T, PP
      PetscReal, intent(out) :: DC, HC
      PetscInt, intent(out)  :: err
      ! Locals:
      PetscReal :: PPb, TA, TB, TC
      PetscReal :: HCI, VC1, VC2

      err = 0

      PPb = PP * 1.0e-6_dp
      TA = T + tc_k
      TB = 0.01_dp * TA
      TC = TB ** 3.3333333333_dp
      HCI = 1.667_dp + 0.001542_dp * TA - 0.7948_dp * log10(TA) - 41.35_dp / TA
      HC = HCI - 0.3571_dp * PPb * (1._dp + 0.07576_dp * PPb) / TC
      HC = HC * 1.e6_dp

      VC1 = 0.00018882_dp * TA
      VC2= - PPb * (0.0824_dp + 0.01249_dp * PPb) / TC
      DC = PPb / (VC1 + VC2)

    end subroutine co2_rho_h

  end subroutine ncg_co2_properties

!------------------------------------------------------------------------

  subroutine ncg_co2_henrys_constant(self, temperature, henrys_constant, err)
    !! Henry's constant for CO2 NCG.

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature
    PetscReal, intent(out) :: henrys_constant
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: T2, T3, RKH

    associate(T => temperature)
      if (T <= 300._dp) then
         T2 = T * T
         T3 = T2 * T
         RKH = (783.666_dp + 19.6025_dp * T + 0.820574_dp * T2) * 1.0e5_dp &
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
    !! Calculates enthalpy of CO2 dissolution in liquid

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature
    PetscReal, intent(out):: energy_solution
    PetscInt, intent(out) :: err     !! error code
    ! Locals:
    PetscReal :: T, T2, T3, T4

    err = 0

    T = 0.01_dp * temperature
    T2 = T * T
    T3 = T * T2
    T4 = T * T3
    energy_solution = 1.e6_dp * (-0.549491_dp + 0.456571_dp * T &
         - 0.070404_dp * T2 - 0.031035_dp * T3 + 0.014121_dp * T4)

  end subroutine ncg_co2_energy_solution

!------------------------------------------------------------------------
  
  subroutine ncg_co2_viscosity(self, partial_pressure, temperature, &
       region, xg, density_g, viscosity, err)
    !! Calculates viscosity for gas phase given partial pressure and temperature
    ! Currently using TOUGH2 formulation

    use thermodynamics_module, only: region_type

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: partial_pressure 
    PetscReal, intent(in) :: temperature
    class(region_type), pointer :: region
    PetscReal, intent(in) :: xg
    PetscReal, intent(in) :: density_g
    PetscReal, intent(out):: viscosity
    PetscInt, intent(out)  :: err
    ! Locals:
    PetscReal :: pbar, pscale, T2, T3, T4
    PetscReal :: C(5)

    pbar = partial_pressure * 1.0e-5_dp
    pscale = 0.01_dp * pbar
    if (pscale <= 1._dp) then
       C = self%viscosity_A + pscale * self%viscosity_BA
       associate(T => temperature)
         T2 = T * T
         T3 = T * T2
         T4 = T2 * T2
         viscosity = 1.0e-8_dp * (C(1) + C(2) * T + &
              C(3) * T2 + C(4) * T3 + C(5) * T4)
       end associate
       err = 0
    else
       err = 1
    end if

  end subroutine ncg_co2_viscosity

!------------------------------------------------------------------------

  subroutine ncg_co2_vapour_mixture_viscosity(self, pressure, &
       temperature, partial_pressure, region, xg, total_density, &
       viscosity, err)
    !! Calculates viscosity for gas phase mixture given partial
    !! pressure and temperature.

    use thermodynamics_module, only: region_type

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: pressure
    PetscReal, intent(in) :: temperature
    PetscReal, intent(in) :: partial_pressure
    class(region_type), pointer :: region
    PetscReal, intent(in) :: xg
    PetscReal, intent(in) :: total_density
    PetscReal, intent(out):: viscosity
    PetscInt, intent(out)  :: err
    ! Locals:
    PetscReal :: gas_viscosity, water_viscosity

    call region%viscosity(temperature, pressure, total_density, &
         water_viscosity)

    call self%viscosity(partial_pressure, temperature, &
         region, xg, xg * total_density, gas_viscosity, err)

    if (err == 0) then
       viscosity = water_viscosity * (1._dp - xg) &
            + gas_viscosity * xg
    end if

  end subroutine ncg_co2_vapour_mixture_viscosity

!------------------------------------------------------------------------

end module ncg_co2_thermodynamics_module
