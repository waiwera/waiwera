module ncg_thermodynamics_module
  !! Module for abstract non-condensible gas thermodynamics type, from
  !! which specific NCG thermodynamics types can be derived.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use thermodynamics_module

  implicit none
  private

  PetscInt, parameter, public :: max_ncg_name_length = 8

  type, public, abstract :: ncg_thermodynamics_type
     !! Non-condensible gas thermodynamics type
     private
     character(max_ncg_name_length), public :: name !! NCG name
     PetscReal, public :: molecular_weight !! NCG molecular weight
   contains
     private
     procedure(ncg_init_procedure), public, deferred :: init
     procedure(ncg_properties_procedure), public, deferred :: properties
     procedure(ncg_henrys_constant_procedure), public, deferred :: henrys_constant
     procedure(ncg_energy_solution_procedure), public, deferred :: energy_solution
     procedure(ncg_viscosity_procedure), public, deferred :: viscosity
     procedure, public :: mass_fraction => ncg_thermodynamics_mass_fraction
     procedure, public :: mole_fraction => ncg_thermodynamics_mole_fraction
   end type ncg_thermodynamics_type

  abstract interface

     subroutine ncg_init_procedure(self)
       !! Initialise NCG thermodynamics object.
       import :: ncg_thermodynamics_type
       class(ncg_thermodynamics_type), intent(in out) :: self
     end subroutine ncg_init_procedure

     subroutine ncg_properties_procedure(self, partial_pressure, &
          temperature, phase, density_water, props, xg, err)
       !! Calculate NCG fluid properties (density and internal energy)
       !! and mass fraction.
       import :: ncg_thermodynamics_type
       class(ncg_thermodynamics_type), intent(in) :: self
       PetscReal, intent(in) :: partial_pressure
       PetscReal, intent(in) :: temperature
       PetscInt, intent(in)  :: phase
       PetscReal, intent(in) :: density_water
       PetscReal, intent(out) :: props(:)
       PetscReal, intent(out) :: xg
       PetscErrorCode, intent(out) :: err
     end subroutine ncg_properties_procedure

     subroutine ncg_henrys_constant_procedure(self, temperature, hc, err)
       !! Calculate NCG Henry's constant, for calculating dissolution
       !! of gas into water.
       import :: ncg_thermodynamics_type
       class(ncg_thermodynamics_type), intent(in) :: self
       PetscReal, intent(in) :: temperature
       PetscReal, intent(out) :: hc
       PetscErrorCode, intent(out) :: err
     end subroutine ncg_henrys_constant_procedure

     subroutine ncg_energy_solution_procedure(self, temperature, &
          h_solution, err)
       !! Calculate NCG energy of solution.
       import :: ncg_thermodynamics_type
       class(ncg_thermodynamics_type), intent(in) :: self
       PetscReal, intent(in) :: temperature
       PetscReal, intent(out) :: h_solution
       PetscErrorCode, intent(out) :: err
     end subroutine ncg_energy_solution_procedure

     subroutine ncg_viscosity_procedure(self, partial_pressure, &
          temperature, region, xg, density_g, viscosity, err)
       !! Calculate NCG viscosity.
       import :: ncg_thermodynamics_type
       import :: region_type
       class(ncg_thermodynamics_type), intent(in) :: self 
       PetscReal, intent(in)  :: partial_pressure
       PetscReal, intent(in)  :: temperature
       class(region_type), pointer :: region
       PetscReal, intent(in)  :: xg
       PetscReal, intent(in)  :: density_g
       PetscReal, intent(out) :: viscosity
       PetscErrorCode, intent(out)  :: err
     end subroutine ncg_viscosity_procedure

  end interface

contains

!------------------------------------------------------------------------

  PetscReal function ncg_thermodynamics_mass_fraction(self, xmole) &
       result(xmass)
    !! Calculates NCG mass fraction from mole fraction.
    class(ncg_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: xmole !! NCG mole fraction
    ! Locals:
    PetscReal :: w

    w = xmole * self%molecular_weight
    xmass = w / (w + (1._dp - xmole) * h2o_molecular_weight)

  end function ncg_thermodynamics_mass_fraction

!------------------------------------------------------------------------

  PetscReal function ncg_thermodynamics_mole_fraction(self, xmass) &
       result(xmole)
    !! Calculates NCG mole fraction from mass fraction.
    class(ncg_thermodynamics_type), intent(in) :: self
    PetscReal, intent(in) :: xmass !! NCG mass fraction
    ! Locals:
    PetscReal :: w

    w = xmass / self%molecular_weight
    xmole = w / (w + (1._dp - xmass) / h2o_molecular_weight)

  end function ncg_thermodynamics_mole_fraction

!------------------------------------------------------------------------

end module ncg_thermodynamics_module
