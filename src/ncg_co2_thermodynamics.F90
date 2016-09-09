module ncg_co2_thermodynamics_module
  !! Module for thermodynamics of CO2 non-condensible gas.

  use kinds_module
  use ncg_thermodynamics_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>
#include <petsc/finclude/petscsys.h>

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
     procedure, public :: henrys_law => ncg_co2_henrys_law
     procedure, public :: energy_solution => ncg_co2_energy_solution
     procedure, public :: viscosity => ncg_co2_viscosity
  end type ncg_co2_thermodynamics_type

contains

!------------------------------------------------------------------------

  subroutine ncg_co2_init(self)
    !! Initialises CO2 NCG thermodynamics object.

    class(ncg_co2_thermodynamics_type), intent(in out) :: self

    self%name = "CO2"

    self%viscosity_BA = self%viscosity_B - self%viscosity_A

  end subroutine ncg_co2_init

!------------------------------------------------------------------------  

  subroutine ncg_co2_properties(self, partial_pressure, temperature, &
       phase, density_water, props, xg, err)
    !! Calculates density and internal energy of mixture of a fluid and and a ncg
    !! pressure (Pa) and temperature (deg C).

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: partial_pressure
    PetscReal, intent(in) :: temperature
    PetscReal, intent(in) :: density_water 
    PetscInt, intent(in)  :: phase !! Fluid phase
    PetscReal, intent(out):: props(:) !! Properties (density and internal energy)
    PetscReal, intent(out):: xg !! Mass fraction of the ncg in this phase
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscReal :: total_density
    PetscReal, parameter :: small = 1.e-30_dp

    err = 0

    associate(density_co2 => props(1), internal_energy_co2 => props(2))

      call co2_rho_h(partial_pressure, temperature, density_co2, &
           internal_energy_co2, err)

      if (err == 0) then
         if (phase == 1) then
            ! liquid
            density_co2 = 0._dp    ! not used for mixture density
            call self%henrys_law(partial_pressure, temperature, xg, err)
         else
            ! vapour
            total_density = density_co2 + density_water
            if (total_density < small) then
               xg = 0._dp
            else
               xg = density_co2 / total_density
            end if
         end if
      end if

    end associate

  contains

    subroutine co2_rho_h(PP, T, DC, HC, err)

      use thermodynamics_module, only: tc_k

      implicit none

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

  subroutine ncg_co2_henrys_law(self, partial_pressure, temperature, xg, err)
    !! Henry's Law for CO2 NCG.

    use thermodynamics_module, only: h2o_molecular_weight

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: partial_pressure
    PetscReal, intent(in) :: temperature
    PetscReal, intent(out) :: xg
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: rkh, xmole
    PetscReal :: T2, T3, T4, T5

    err = 0

    T2 = temperature * temperature
    T3 = T2 * temperature
    T4 = T2 * T2
    T5 = T2 * T3

    rkh = (783.666_dp + 19.6025_dp * temperature + 0.820574_dp * T2 &
         - 7.40674e-3_dp * T3 + 2.18380e-5_dp * T4 - 2.20999e-8_dp * T5) &
         * 1.e5_dp
    xmole = partial_pressure / rkh

    xg = xmole * co2_molecular_weight / &
         (xmole * co2_molecular_weight + (1._dp - xmole) * h2o_molecular_weight)

  end subroutine ncg_co2_henrys_law

!------------------------------------------------------------------------

  subroutine ncg_co2_energy_solution(self, temperature, h_solution, err)
    !! Calculates enthalpy of CO2 dissolution in liquid

    class(ncg_co2_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature
    PetscReal, intent(out):: h_solution
    PetscInt, intent(out) :: err     !! error code
    ! Locals:
    PetscReal :: T, T2, T3, T4

    err = 0

    T = 0.01_dp * temperature
    T2 = T * T
    T3 = T * T2
    T4 = T * T3
    h_solution = 1.e6_dp * (-0.073696_dp - 0.56405_dp * T + &
         0.70363_dp * T2 - 0.27882_dp * T3 + 0.042579_dp * T4)

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
       T2 = temperature * temperature
       T3 = temperature * T2
       T4 = T2 * T2
       viscosity = 1.0e-8_dp * (C(1) + C(2) * temperature + &
            C(3) * T2 + C(4) * T3 + C(5) * T4)
       err = 0
    else
       err = 1
    end if

  end subroutine ncg_co2_viscosity

!------------------------------------------------------------------------

end module ncg_co2_thermodynamics_module
