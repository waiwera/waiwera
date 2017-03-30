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
       test_ncg_co2_viscosity, test_ncg_co2_properties

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
       call assert_equals(0, err, trim(s) // " 20 deg C error")
       call assert_equals(expected, hc, tol, trim(s) // " 20 deg C")

       temperature = 100._dp
       expected = 0.181629386327e-08_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(0, err, trim(s) // " 100 deg C error")
       call assert_equals(expected, hc, tol, trim(s) // " 100 deg C")

       temperature = 240._dp
       expected = 0.191626750106e-08_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(0, err, trim(s) // " 240 deg C error")
       call assert_equals(expected, hc, tol, trim(s) // " 240 deg C")

       temperature = 300._dp
       expected = 0.268879436880e-08_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(0, err, trim(s) // " 300 deg C error")
       call assert_equals(expected, hc, tol, trim(s) // " 300 deg C")

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
       call assert_equals(0, err, trim(s) // " 20 deg C error")
       call assert_equals(expected, hs, tol, trim(s) // " 20 deg C")

       temperature = 100._dp
       expected = -180238.0_dp
       call gas%energy_solution(temperature, hs, err)
       call assert_equals(0, err, trim(s) // " 100 deg C error")
       call assert_equals(expected, hs, tol, trim(s) // " 100 deg C")

       temperature = 240._dp
       expected = 180225.4096_dp
       call gas%energy_solution(temperature, hs, err)
       call assert_equals(0, err, trim(s) // " 240 deg C error")
       call assert_equals(expected, hs, tol, trim(s) // " 240 deg C")

       temperature = 300._dp
       expected = 492442.0_dp
       call gas%energy_solution(temperature, hs, err)
       call assert_equals(0, err, trim(s) // " 300 deg C error")
       call assert_equals(expected, hs, tol, trim(s) // " 300 deg C")

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
             call assert_equals(0, err, trim(s) // ' error')
             call assert_equals(expected_visc(ip, it), visc, tol, s)
          end do
       end do

    end if

  end subroutine test_ncg_co2_viscosity

!------------------------------------------------------------------------

  subroutine test_ncg_co2_properties
    ! CO2 properties

    type(ncg_co2_thermodynamics_type) :: gas
    PetscMPIInt :: rank
    PetscInt :: ierr, err, i
    PetscReal :: props(2)
    character(40) :: s
    PetscReal, parameter :: htol = 1.e-3_dp, dtol = 1.e-7_dp
    PetscInt, parameter :: num_cases = 14
    PetscReal, parameter :: data(4, num_cases) = reshape([ &
         0.0_dp, 20.0_dp, 17140.18077231938_dp, 0.0_dp, &
         100000.0_dp, 20.0_dp, 16142.247883091828_dp, 1.8142044368713437_dp, &
         0.0_dp, 100.0_dp, 87450.99131436742_dp, 0.0_dp, &
         100000.0_dp, 100.0_dp, 87004.524163092_dp, 1.4213754811567743_dp, &
         4000000.0_dp, 100.0_dp, 64355.3813832885_dp, 62.608990505735434_dp, &
         9000000.0_dp, 100.0_dp, 20379.357776952613_dp, 184.7959892299282_dp, &
         0.0_dp, 240.0_dp, 223594.37705727902_dp, 0.0_dp, &
         100000.0_dp, 240.0_dp, 223439.99865083068_dp, 1.0324489144812645_dp, &
         4000000.0_dp, 240.0_dp, 215608.4290498441_dp, 42.27375154306431_dp, &
         9000000.0_dp, 240.0_dp, 200402.49860929986_dp, 100.70459422220841_dp, &
         0.0_dp, 300.0_dp, 286380.4950504236_dp, 0.0_dp, &
         100000.0_dp, 300.0_dp, 286273.71092985675_dp, 0.9242369906584087_dp, &
         4000000.0_dp, 300.0_dp, 280856.58497462136_dp, 37.5055455044134_dp, &
         9000000.0_dp, 300.0_dp, 270338.58607276645_dp, 87.3627658128452_dp &
         ], [4, num_cases])

    call gas%init()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       do i = 1, num_cases
          associate(pc => data(1, i), t => data(2, i), &
               expected_h => data(3, i), expected_d => data(4, i), &
               d => props(1), h => props(2))
            call gas%properties(pc, t, props, err)
            write(s, '(a, e8.3, a, f6.2)') 'pc = ', pc, ', t = ', t
            call assert_equals(0, err, trim(s) // ' error')
            call assert_equals(expected_h, h, htol, trim(s) // ' enthalpy')
            call assert_equals(expected_d, d, dtol, trim(s) // ' density')
          end associate
       end do

    end if

  end subroutine test_ncg_co2_properties

!------------------------------------------------------------------------

end module ncg_co2_thermodynamics_test
