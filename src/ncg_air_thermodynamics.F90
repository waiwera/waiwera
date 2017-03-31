module ncg_air_thermodynamics_module
  !! Module for thermodynamics of air non-condensible gas.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use ncg_thermodynamics_module

  implicit none
  private

  PetscReal, parameter, public :: air_molecular_weight = 28.96_dp ! g/mol
  PetscReal, parameter, public :: air_specific_heat = 733.0_dp ! J/kg/K

  type, public, extends(ncg_thermodynamics_type) :: ncg_air_thermodynamics_type
     !! Type for air NCG thermodynamics.
     private
     PetscReal :: specific_heat
     PetscReal :: fair = 97.0_dp
     PetscReal :: fwat = 363.0_dp
     PetscReal :: cair = 3.617_dp
     PetscReal :: cwat = 2.655_dp
   contains
     private
     procedure, public :: init => ncg_air_init
     procedure, public :: properties => ncg_air_properties
     procedure, public :: effective_properties => ncg_air_effective_properties
     procedure, public :: henrys_constant => ncg_air_henrys_constant
     procedure, public :: energy_solution => ncg_air_energy_solution
     procedure, public :: viscosity => ncg_air_viscosity
     procedure, public :: mixture_viscosity => ncg_air_mixture_viscosity
  end type ncg_air_thermodynamics_type

contains

!------------------------------------------------------------------------

  subroutine ncg_air_init(self)
    !! Initialises air NCG thermodynamics object.

    class(ncg_air_thermodynamics_type), intent(in out) :: self

    self%name = "Air"
    self%molecular_weight = air_molecular_weight
    self%specific_heat = air_specific_heat

  end subroutine ncg_air_init

!------------------------------------------------------------------------  

  subroutine ncg_air_properties(self, partial_pressure, temperature, &
       props, err)
    !! Calculates air NCG density and enthalpy.

    use thermodynamics_module, only: tc_k, gas_constant

    class(ncg_air_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: partial_pressure !! Air partial pressure
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(out):: props(:) !! Properties (density and enthalpy)
    PetscErrorCode, intent(out) :: err !! Error code

    err = 0

    associate(tk => temperature + tc_k, air_density => props(1), &
         air_enthalpy => props(2))

      air_density = partial_pressure * self%molecular_weight / &
           (1.e3 * gas_constant * self%deviation_factor * tk)
      if (air_density > 0._dp) then
         air_enthalpy = self%specific_heat * temperature + &
              partial_pressure / air_density
      else
         air_enthalpy = 0._dp
      end if

    end associate

  end subroutine ncg_air_properties

!------------------------------------------------------------------------

  subroutine ncg_air_effective_properties(self, props, phase, &
       effective_props)
    !! Returns effective air NCG properties for specified phase.
    !! Air density and enthalpy are treated as effectively zero in the
    !! liquid phase.

    class(ncg_air_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: props(:) !! Air NCG properties (density, enthalpy)
    PetscInt, intent(in) :: phase !! Phase index
    PetscReal, intent(out) :: effective_props(:) !! Effective NCG properties

    effective_props = props

    if (phase == 1) then
       associate(air_density => effective_props(1), &
            air_enthalpy => effective_props(2))
         air_density = 0._dp
         air_enthalpy = 0._dp
       end associate
    end if

  end subroutine ncg_air_effective_properties

!------------------------------------------------------------------------

  subroutine ncg_air_henrys_constant(self, temperature, &
       henrys_constant, err)
    !! Henry's constant for air NCG.

    class(ncg_air_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(out) :: henrys_constant !! Henry's constant
    PetscErrorCode, intent(out) :: err !! Error code

    err = 0
    henrys_constant = 1.e-10_dp

  end subroutine ncg_air_henrys_constant

!------------------------------------------------------------------------

  subroutine ncg_air_energy_solution(self, temperature, energy_solution, err)
    !! Enthalpy of air dissolution in liquid.

    class(ncg_air_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(out):: energy_solution !! Energy of solution
    PetscInt, intent(out) :: err     !! Error code

    err = 0
    energy_solution = 0._dp

  end subroutine ncg_air_energy_solution

!------------------------------------------------------------------------
  
  subroutine ncg_air_viscosity(self, partial_pressure, temperature, &
       viscosity, err)
    !! Viscosity for air, given partial pressure and temperature.
    !! This is not used.

    use thermodynamics_module, only: region_type

    class(ncg_air_thermodynamics_type), intent(in) :: self
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
    
    use thermodynamics_module

    class(ncg_air_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: water_viscosity !! Viscosity of water
    PetscReal, intent(in) :: temperature !! Temperature
    PetscReal, intent(in) :: partial_pressure !! Air partial pressure (not used)
    PetscReal, intent(in) :: xg !! Air mass fraction
    PetscInt, intent(in)  :: phase !! Phase index
    PetscReal, intent(out):: viscosity !! Mixture viscosity
    PetscInt, intent(out)  :: err !! Error code
    ! Locals:
    PetscReal :: x1, x2, ard, cmix
    PetscReal :: e, fmix, g, h
    PetscReal :: ome1, ome3, rm1, rm2, rm3
    PetscReal :: trd1, trd3, vis1, vis2, vis3, z1, z2, z3

    err = 0

    if (phase == 1) then

       viscosity = water_viscosity

    else

       rm1 = self%molecular_weight
       rm2 = water_molecular_weight

       fmix = sqrt(self%fair * self%fwat)
       cmix = 0.5_dp * (self%cair + self%cwat)

       x1 = self%mass_to_mole_fraction(xg)
       x2 = 1._dp - x1

       associate (tk => temperature + tc_k)
         trd1 = tk / self%fair
         trd3 = tk / fmix
       end associate
       ome1 = (1.188_dp - 0.051_dp * trd1) / trd1
       ome3 = (1.48_dp  - 0.412_dp * log(trd3)) / trd3
       ard = 1.095_dp / trd3
       rm3 = 2._dp * rm1 * rm2 / (rm1 + rm2)
       vis1 = covis(trd1, self%cair, ome1, rm1, self%fair)
       vis2 = 10._dp * water_viscosity
       vis3 = covis(trd3, cmix, ome3, rm3, fmix)
       z1 = x1 * x1 / vis1 + 2._dp * x2 * x1 / vis3 + x2 * x2 / vis2
       g = x1 * x1 * rm1 / rm2
       h = x2 * x2 * rm2 / rm1
       e = (2._dp * x1 * x2 * rm1 * rm2 / (rm3 * rm3)) * &
            vis3 / (vis1 * vis2)
       z2 = 0.6_dp * ard * (g / vis1 + e + h / vis2)
       z3 = 0.6_dp * ard * (g + e * (vis1 + vis2) - 2._dp * x1 * x2 + h)
       viscosity = 0.1_dp * (1._dp + z3) / (z1 + z2)

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
