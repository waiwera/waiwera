module ncg_co2_thermodynamics_test

  ! Tests for CO2 NCG thermodynamics module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use ncg_co2_thermodynamics_module
  use fruit

  implicit none
  private

  public :: test_ncg_co2_henrys_constant, test_ncg_co2_energy_solution

contains

!------------------------------------------------------------------------

  subroutine test_ncg_co2_henrys_constant

    ! CO2 Henry's constant tests.
    ! Expected values are from AUTOUGH2.

    type(ncg_co2_thermodynamics_type) :: gas
    PetscReal :: temperature, expected, hc
    PetscErrorCode :: err
    character(33) :: s = "CO2 Henry's constant, temperature"
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: tol = 1.e-15_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       temperature = 20._dp
       expected = 0.690552871945e-08_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(0, err, s // " 20 deg C error")
       call assert_equals(expected, hc, tol, s // " 20 deg C")

       temperature = 100._dp
       expected = 0.181629386327e-08_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(0, err, s // " 100 deg C error")
       call assert_equals(expected, hc, tol, s // " 100 deg C")

       temperature = 240._dp
       expected = 0.191626750106e-08_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(0, err, s // " 240 deg C error")
       call assert_equals(expected, hc, tol, s // " 240 deg C")

       temperature = 300._dp
       expected = 0.268879436880e-08_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(0, err, s // " 300 deg C error")
       call assert_equals(expected, hc, tol, s // " 300 deg C")

    end if

  end subroutine test_ncg_co2_henrys_constant

!------------------------------------------------------------------------

  subroutine test_ncg_co2_energy_solution

    type(ncg_co2_thermodynamics_type) :: gas
    PetscReal :: temperature, expected, hs
    PetscErrorCode :: err
    character(35) :: s = "CO2 energy of solution, temperature"
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: tol = 1.e-4_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       temperature = 20._dp
       expected = -461218.6464_dp
       call gas%energy_solution(temperature, hs, err)
       call assert_equals(0, err, s // " 20 deg C error")
       call assert_equals(expected, hs, tol, s // " 20 deg C")

       temperature = 100._dp
       expected = -180238.0_dp
       call gas%energy_solution(temperature, hs, err)
       call assert_equals(0, err, s // " 100 deg C error")
       call assert_equals(expected, hs, tol, s // " 100 deg C")

       temperature = 240._dp
       expected = 180225.4096_dp
       call gas%energy_solution(temperature, hs, err)
       call assert_equals(0, err, s // " 240 deg C error")
       call assert_equals(expected, hs, tol, s // " 240 deg C")

       temperature = 300._dp
       expected = 492442.0_dp
       call gas%energy_solution(temperature, hs, err)
       call assert_equals(0, err, s // " 300 deg C error")
       call assert_equals(expected, hs, tol, s // " 300 deg C")

    end if

  end subroutine test_ncg_co2_energy_solution

!------------------------------------------------------------------------

end module ncg_co2_thermodynamics_test
