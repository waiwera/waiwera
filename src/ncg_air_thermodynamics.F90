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
     procedure, public :: henrys_constant => ncg_air_henrys_constant
     procedure, public :: energy_solution => ncg_air_energy_solution
     procedure, public :: viscosity => ncg_air_viscosity
     procedure, public :: vapour_mixture_viscosity => ncg_air_vapour_mixture_viscosity
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
       phase, density_water, props, xg, err)
    !! Calculates density and internal energy of mixture of a fluid and and a ncg
    !! pressure (Pa) and temperature (deg C).

    class(ncg_air_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: partial_pressure
    PetscReal, intent(in) :: temperature
    PetscReal, intent(in) :: density_water 
    PetscInt, intent(in)  :: phase !! Fluid phase
    PetscReal, intent(out):: props(:) !! Properties (density and internal energy)
    PetscReal, intent(out):: xg !! Mass fraction of the ncg in this phase
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscReal :: total_density, hc, xmole
    PetscReal, parameter :: small = 1.e-30_dp

    err = 0

    associate(density_air => props(1), enthalpy_air => props(2))

      call air_rho_h(partial_pressure, temperature, density_air, &
           enthalpy_air, err)

      if (err == 0) then
         if (phase == 1) then
            ! liquid
            density_air = 0._dp    ! not used for mixture density
            enthalpy_air = 0._dp
            call self%henrys_constant(temperature, hc, err)
            if (err == 0) then
               xmole = hc * partial_pressure
               xg = self%mass_fraction(xmole)
            end if
         else
            ! vapour
            total_density = density_air + density_water
            if (total_density < small) then
               xg = 0._dp
            else
               xg = density_air / total_density
            end if
         end if
      end if

    end associate

  contains

    subroutine air_rho_h(partial_pressure, temperature, density_air, &
           enthalpy_air, err)

      use thermodynamics_module, only: tc_k, gas_constant

      PetscReal, intent(in)  :: partial_pressure, temperature
      PetscReal, intent(out) :: density_air, enthalpy_air
      PetscInt, intent(out)  :: err

      err = 0

      associate(TK => temperature + tc_k)
        density_air = partial_pressure * self%molecular_weight / &
             (1.e3 * gas_constant * self%deviation_factor * TK)
      end associate

      if (density_air > 0._dp) then 
         enthalpy_air = self%specific_heat * temperature + &
              partial_pressure / density_air
      else
         enthalpy_air = 0._dp
      end if

    end subroutine air_rho_h

  end subroutine ncg_air_properties

!------------------------------------------------------------------------

  subroutine ncg_air_henrys_constant(self, temperature, hc, err)
    !! Henry's constant for air NCG.

    class(ncg_air_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature
    PetscReal, intent(out) :: hc
    PetscErrorCode, intent(out) :: err

    err = 0
    hc = 1.e-10_dp

  end subroutine ncg_air_henrys_constant

!------------------------------------------------------------------------

  subroutine ncg_air_energy_solution(self, temperature, h_solution, err)
    !! Enthalpy of air dissolution in liquid.

    class(ncg_air_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: temperature
    PetscReal, intent(out):: h_solution
    PetscInt, intent(out) :: err     !! error code

    err = 0
    h_solution = 0._dp

  end subroutine ncg_air_energy_solution

!------------------------------------------------------------------------
  
  subroutine ncg_air_viscosity(self, partial_pressure, temperature, &
       region, xg, density_g, viscosity, err)
    !! Viscosity for air, given partial pressure and temperature.

    use thermodynamics_module, only: region_type

    class(ncg_air_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: partial_pressure 
    PetscReal, intent(in) :: temperature
    class(region_type), pointer :: region
    PetscReal, intent(in) :: xg
    PetscReal, intent(in) :: density_g
    PetscReal, intent(out):: viscosity
    PetscInt, intent(out)  :: err

    err = 0
    viscosity = 0._dp

  end subroutine ncg_air_viscosity

!------------------------------------------------------------------------

  subroutine ncg_air_vapour_mixture_viscosity(self, pressure, &
       temperature, partial_pressure, region, xg, density, viscosity, err)
    !! Calculates viscosity for the gas phase mixture, given partial
    !! pressure and temperature.

    ! This routine computes the viscosity of vapor-air mixtures.
    ! It uses a modified version of a formulation based on kinetic
    ! gas theory, as given by J.O. Hirschfelder, C.F. Curtiss, and
    ! R.B. Bird, Molecular Theory of Gases and Liquids, John Wiley
    ! & Sons, 1954, pp. 528-530.
    !
    ! The modification made to the Hirschfelder et al. expressions is
    ! that for vapor viscosity accurate (empirical) values are used,
    ! rather than the first order expression of kinetic theory.
    !
    ! The formulation matches experimental data on viscosities of
    ! vapor-air mixtures in the temperature range from 100 to 150
    ! deg. C, for all compositions, to better than 4%.
    
    use thermodynamics_module

    class(ncg_air_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: pressure 
    PetscReal, intent(in) :: temperature
    PetscReal, intent(in) :: partial_pressure 
    class(region_type), pointer :: region
    PetscReal, intent(in) :: xg
    PetscReal, intent(in) :: density
    PetscReal, intent(out):: viscosity
    PetscInt, intent(out)  :: err
    ! Locals:
    PetscReal :: vs, x1, x2, ard, cmix
    PetscReal :: e, fmix, g, h
    PetscReal :: ome1, ome3, rm1, rm2, rm3
    PetscReal :: trd1, trd3, vis1, vis2, vis3, z1, z2, z3

    err = 0

    rm1 = self%molecular_weight
    rm2 = water_molecular_weight

    fmix = sqrt(self%fair * self%fwat)
    cmix = 0.5_dp * (self%cair + self%cwat)

    x1 = self%mole_fraction(xg)
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

    call region%viscosity(temperature, pressure, density, vs)

    vis2 = 10._dp * vs
    vis3 = covis(trd3, cmix, ome3, rm3, fmix)
    z1 = x1 * x1 / vis1 + 2._dp * x2 * x1 / vis3 + x2 * x2 / vis2
    g = x1 * x1 * rm1 / rm2
    h = x2 * x2 * rm2 / rm1
    e = (2._dp * x1 * x2 * rm1 * rm2 / (rm3 * rm3)) * &
         vis3 / (vis1 * vis2)
    z2 = 0.6_dp * ard * (g / vis1 + e + h / vis2)
    z3 = 0.6_dp * ard * (g + e * (vis1 + vis2) - 2._dp * x1 * x2 + h)
    viscosity = 0.1_dp * (1._dp + z3) / (z1 + z2)

  contains

    PetscReal function covis(trd, c, ome, rm, f)
      ! Coefficient of viscosity to the first approximation.
      PetscReal, intent(in)  :: trd, c, ome, rm, f
      covis = 266.93e-7_dp * sqrt(rm * trd * f) / (c * c * ome * trd)
    end function covis

  end subroutine ncg_air_vapour_mixture_viscosity

!------------------------------------------------------------------------

end module ncg_air_thermodynamics_module
