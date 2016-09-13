module ncg_co2_thermodynamics_test

  ! Tests for CO2 NCG thermodynamics module

  use kinds_module
  use mpi_module
  use ncg_co2_thermodynamics_module
  use fruit

  implicit none
  private

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscdef.h>

  public :: test_ncg_co2_henrys_constant

contains

!------------------------------------------------------------------------

  subroutine test_ncg_co2_henrys_constant

    ! CO2 Henry's constant tests.
    ! Expected values are from AUTOUGH2.

    type(ncg_co2_thermodynamics_type) :: gas
    PetscReal :: temperature, expected, hc
    PetscErrorCode :: err
    character(33) :: s = "CO2 Henry's constant, temperature"
    PetscReal, parameter :: tol = 1.e-15_dp

    if (mpi%rank == mpi%output_rank) then

       temperature = 20._dp
       expected = 0.690552871945e-08_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(expected, hc, tol, s // " 20 deg C")

       temperature = 100._dp
       expected = 0.181629386327e-08_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(expected, hc, tol, s // " 100 deg C")

       temperature = 240._dp
       expected = 0.191626750106e-08_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(expected, hc, tol, s // " 240 deg C")

       temperature = 300._dp
       expected = 0.268879436880e-08_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(expected, hc, tol, s // " 300 deg C")

    end if

  end subroutine test_ncg_co2_henrys_constant

!------------------------------------------------------------------------

end module ncg_co2_thermodynamics_test
