module ncg_co2_thermodynamics_test

  ! Tests for CO2 NCG thermodynamics module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use ncg_co2_thermodynamics_module
  use fruit

  implicit none
  private

  public :: test_ncg_co2_henrys_constant, test_ncg_co2_energy_solution, &
       test_ncg_co2_viscosity

contains

!------------------------------------------------------------------------

  subroutine test_ncg_co2_henrys_constant

    ! CO2 Henry's constant.
    ! Expected values are from AUTOUGH2.

    type(ncg_co2_thermodynamics_type) :: gas
    PetscReal :: temperature, expected, hc
    PetscErrorCode :: err
    character(33) :: s = "CO2 Henry's constant, temperature"
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: tol = 1.e-15_dp

    call gas%init()

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
    ! CO2 energy of solution

    type(ncg_co2_thermodynamics_type) :: gas
    PetscReal :: temperature, expected, hs
    PetscErrorCode :: err
    character(35) :: s = "CO2 energy of solution, temperature"
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: tol = 1.e-4_dp

    call gas%init()

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

  subroutine test_ncg_co2_viscosity
    ! CO2 viscosity

    type(ncg_co2_thermodynamics_type) :: gas
    PetscMPIInt :: rank
    PetscInt :: ierr, err
    PetscInt :: ip, it
    PetscReal :: visc
    character(40) :: s
    PetscReal, parameter :: tol = 1.e-15_dp
    PetscReal, parameter :: pc(4) = [0._dp, 1.e5_dp, 40.e5_dp, 90.e5_dp]
    PetscReal, parameter :: t(4) = [20._dp, 100._dp, 240._dp, 300._dp]
    PetscReal, parameter :: expected_visc(4,4) = reshape( &
         [1.4550900339359998e-05_dp, &
         1.4735085044446403e-05_dp, &
         2.1918288542816003e-05_dp, &
         3.1127523797136011e-05_dp, &
         1.8230436100000001e-05_dp, &
         1.8274211538999999e-05_dp, &
         1.9981453660000002e-05_dp, &
         2.2170225610000008e-05_dp, &
         2.4005967912959995e-05_dp, &
         2.4028580667110404e-05_dp, &
         2.4910478078976012e-05_dp, &
         2.6041115786496016e-05_dp, &
         2.6270078099999996e-05_dp, &
         2.6285773119000002e-05_dp, &
         2.689787886e-05_dp, &
         2.7682629810000054e-05_dp], &
         [4, 4])

    call gas%init()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       do ip = 1, 4
          do it = 1, 4
             call gas%viscosity(pc(ip), t(it), visc, err)
             write(s, '(a, e8.3, a, f6.2)') 'p = ', &
                  pc(ip), ', t = ', t(it)
             call assert_equals(0, err, s // ' error')
             call assert_equals(expected_visc(ip, it), visc, tol, s)
          end do
       end do

    end if

  end subroutine test_ncg_co2_viscosity

!------------------------------------------------------------------------

end module ncg_co2_thermodynamics_test
