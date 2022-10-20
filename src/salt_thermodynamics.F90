module salt_thermodynamics_module
  !! Module for salt thermodynamics routines.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use thermodynamics_module
  use utils_module, only: polynomial

  implicit none
  private
  PetscReal :: halite_solubility_data(3) = [2.6218e-1_dp, 7.2e-2_dp, 1.06_dp]

  public :: halite_solubility

contains

!------------------------------------------------------------------------

  subroutine halite_solubility(temperature, solubility, err)
    !! Equilibrium solubility of salt in water as a function of
    !! temperature. From Chou, I.M. (1987).

    PetscReal, intent(in) :: temperature
    PetscReal, intent(out) :: solubility
    PetscErrorCode, intent(out) :: err

    if ((0._dp <= temperature) .and. (temperature <= 382._dp)) then
       solubility = polynomial(halite_solubility_data, temperature / 1.e3_dp)
       err = 0
    else
       solubility = 0._dp
       err = 1
    end if

  end subroutine halite_solubility

!------------------------------------------------------------------------

end module salt_thermodynamics_module
