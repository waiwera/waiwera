module ncg_thermodynamics_module
  !! Module for abstract non-condensible gas thermodynamics type, from
  !! which specific NCG thermodynamics types can be derived.

  use thermodynamics_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>

  PetscInt, parameter, public :: max_ncg_name_length = 8

  type, public, abstract :: ncg_thermodynamics_type
     !! Non-condensible gas thermodynamics type
     private
     character(max_ncg_name_length), public :: name !! NCG name
   contains
     private
     procedure(ncg_init_procedure), public, deferred :: init
     procedure(ncg_properties_procedure), public, deferred :: properties
     procedure(ncg_henrys_law_procedure), public, deferred :: henrys_law
     procedure(ncg_energy_solution_procedure), public, deferred :: energy_solution
     procedure(ncg_viscosity_procedure), public, deferred :: viscosity
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

     subroutine ncg_henrys_law_procedure(self, partial_pressure, &
          temperature, xg, err)
       !! Calculate NCG mass fraction from Henry's Law.
       import :: ncg_thermodynamics_type
       class(ncg_thermodynamics_type), intent(in) :: self
       PetscReal, intent(in) :: partial_pressure
       PetscReal, intent(in) :: temperature
       PetscReal, intent(out) :: xg
       PetscErrorCode, intent(out) :: err
     end subroutine ncg_henrys_law_procedure

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

end module ncg_thermodynamics_module
