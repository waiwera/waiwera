module ncg_co2_thermodynamics_test

  ! Tests for CO2 NCG thermodynamics module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use ncg_co2_thermodynamics_module
  use zofu

  implicit none
  private

  public :: setup, teardown
  public :: test_ncg_co2_henrys_constant, test_ncg_co2_energy_solution, &
       test_ncg_co2_viscosity, test_ncg_co2_properties

contains

!------------------------------------------------------------------------

  subroutine setup()

    use profiling_module, only: init_profiling

    ! Locals:
    PetscErrorCode :: ierr

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call init_profiling()

  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown()

    PetscErrorCode :: ierr

    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

  subroutine test_ncg_co2_henrys_constant(test)

    ! CO2 Henry constant.
    ! Expected values are from AUTOUGH2.

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(ncg_co2_thermodynamics_type) :: gas
    PetscReal :: temperature, expected, hc, chc(1)
    PetscErrorCode :: err
    character(31) :: s = "CO2 Henry constant, temperature"
    PetscMPIInt :: rank
    PetscInt :: ierr

    call gas%init()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       temperature = 20._dp
       expected = 1.44811504032e+08_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call test%assert(0, err, trim(s) // " 20 deg C error")
       call test%assert(expected, hc, trim(s) // " 20 deg C")

       call gas%henrys_constant_salt(temperature, 0._dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 20 deg C zero salt error")
       call test%assert(expected, hc, trim(s) // " 20 deg C zero salt")

       call gas%henrys_constant_salt(temperature, 0.1_dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 20 deg C 0.1 salt error")
       call test%assert(1.599700034044_dp * expected, hc, trim(s) // " 20 deg C 0.1 salt")

       call gas%henrys_constant_salt(temperature, 0.25_dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 20 deg C 0.25 salt error")
       call test%assert(4.093696693331_dp * expected, hc, trim(s) // " 20 deg C 0.25 salt")

       temperature = 100._dp
       expected = 5.50571700000e+08_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call test%assert(0, err, trim(s) // " 100 deg C error")
       call test%assert(expected, hc, trim(s) // " 100 deg C")

       call gas%henrys_constant_salt(temperature, 0._dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 100 deg C zero salt error")
       call test%assert(expected, hc, trim(s) // " 100 deg C zero salt")

       call gas%henrys_constant_salt(temperature, 0.1_dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 100 deg C 0.1 salt error")
       call test%assert(1.470515461734_dp * expected, hc, trim(s) // " 100 deg C 0.1 salt")

       call gas%henrys_constant_salt(temperature, 0.3_dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 100 deg C 0.3 salt error")
       call test%assert(4.425416524041_dp * expected, hc, trim(s) // " 100 deg C 0.3 salt")

       temperature = 240._dp
       expected = 5.21847810624e+08_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call test%assert(0, err, trim(s) // " 240 deg C error")
       call test%assert(expected, hc, trim(s) // " 240 deg C")

       temperature = 300._dp
       expected = 3.71913900000e+08_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call test%assert(0, err, trim(s) // " 300 deg C error")
       call test%assert(expected, hc, trim(s) // " 300 deg C")

       call gas%henrys_constant_salt(temperature, 0._dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 300 deg C zero salt error")
       call test%assert(expected, hc, trim(s) // " 300 deg C zero salt")

       call gas%henrys_constant_salt(temperature, 0.1_dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 300 deg C 0.1 salt error")
       call test%assert(1.981753363144_dp * expected, hc, trim(s) // " 300 deg C 0.1 salt")

       call gas%henrys_constant_salt(temperature, 0.3_dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 300 deg C 0.3 salt error")
       call test%assert(13.988229784674_dp * expected, hc, trim(s) // " 300 deg C 0.3 salt")

       temperature = 350._dp
       expected = 2.23454746875e+08_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call test%assert(0, err, trim(s) // " 350 deg C error")
       call test%assert(expected, hc, trim(s) // " 350 deg C")

    end if

  end subroutine test_ncg_co2_henrys_constant

!------------------------------------------------------------------------

  subroutine test_ncg_co2_energy_solution(test)
    ! CO2 energy of solution

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(ncg_co2_thermodynamics_type) :: gas
    PetscReal :: temperature, expected, hs, hc, chc(1)
    PetscErrorCode :: err
    character(35) :: s = "CO2 energy of solution, temperature"
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal :: xs

    call gas%init()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       temperature = 20._dp
       expected = -495750.87299689_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call gas%energy_solution(temperature, chc, hs, err)
       call test%assert(0, err, trim(s) // " 20 deg C error")
       call test%assert(expected, hs, trim(s) // " 20 deg C")

       xs = 0._dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 20 deg C 0 salt error")
       call test%assert(expected, hs, trim(s) // " 20 deg C 0 salt")

       xs = 0.1_dp
       expected = -457912.15278475045_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 20 deg C 0.1 salt error")
       call test%assert(expected, hs, trim(s) // " 20 deg C 0.1 salt")

       xs = 0.2_dp
       expected = -410613.7525195794_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 20 deg C 0.2 salt error")
       call test%assert(expected, hs, trim(s) // " 20 deg C 0.2 salt")

       xs = 0.3_dp
       expected = -349801.52360721654_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 20 deg C 0.3 salt error")
       call test%assert(expected, hs, trim(s) // " 20 deg C 0.3 salt")

       temperature = 100._dp
       expected = -180685.98723494_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call gas%energy_solution(temperature, chc, hs, err)
       call test%assert(0, err, trim(s) // " 100 deg C error")
       call test%assert(expected, hs, trim(s) // " 100 deg C")

       xs = 0.1_dp
       expected = -180881.77711991896_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 100 deg C 0.1 salt error")
       call test%assert(expected, hs, trim(s) // " 100 deg C 0.1 salt")

       xs = 0.2_dp
       expected = -181126.5144761395_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 100 deg C 0.2 salt error")
       call test%assert(expected, hs, trim(s) // " 100 deg C 0.2 salt")

       xs = 0.3_dp
       expected = -181441.17679128022_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 100 deg C 0.3 salt error")
       call test%assert(expected, hs, trim(s) // " 100 deg C 0.3 salt")

       temperature = 240._dp
       expected = 242741.64505202_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call gas%energy_solution(temperature, chc, hs, err)
       call test%assert(0, err, trim(s) // " 240 deg C error")
       call test%assert(expected, hs, trim(s) // " 240 deg C")

       temperature = 300._dp
       expected = 407409.27618764_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call gas%energy_solution(temperature, chc, hs, err)
       call test%assert(0, err, trim(s) // " 300 deg C error")
       call test%assert(expected, hs, trim(s) // " 300 deg C")

       xs = 0.1_dp
       expected = 241434.91186903339_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 300 deg C 0.1 salt error")
       call test%assert(expected, hs, trim(s) // " 300 deg C 0.1 salt")

       xs = 0.2_dp
       expected = 33966.95647077112_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 300 deg C 0.2 salt error")
       call test%assert(expected, hs, trim(s) // " 300 deg C 0.2 salt")

       xs = 0.3_dp
       expected = -232777.5576127092_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 300 deg C 0.3 salt error")
       call test%assert(expected, hs, trim(s) // " 300 deg C 0.3 salt")

    end if

  end subroutine test_ncg_co2_energy_solution

!------------------------------------------------------------------------

  subroutine test_ncg_co2_viscosity(test)
    ! CO2 viscosity

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(ncg_co2_thermodynamics_type) :: gas
    PetscMPIInt :: rank
    PetscInt :: ierr, err
    PetscInt :: ip, it
    PetscReal :: Pi, Ti, visc
    character(40) :: s
    PetscReal, parameter :: tol = 1.e-8_dp
    PetscReal, parameter :: pc(6) = [0.1e6_dp, 1.e6_dp, 5.e6_dp, &
         10.e6_dp, 20.e6_dp, 30.e6_dp]
    PetscReal, parameter :: t(5) = [20._dp, 100._dp, 200._dp, 300._dp, 350._dp]
    PetscReal, parameter :: expected_visc(6, 5) = reshape([ &
         1.47350850e-5_dp, 1.63927474e-5_dp, 2.37601356e-5_dp, &
         3.29693708e-5_dp, 9.99600434e-5_dp, 1.19066342e-4_dp, &
         1.82742115e-5_dp, 1.86681905e-5_dp, 2.04192081e-5_dp, &
         2.26079800e-5_dp, 3.76607100e-5_dp, 5.40893300e-5_dp, &
         2.24530737e-5_dp, 2.26583470e-5_dp, 2.35706728e-5_dp, &
         2.47110800e-5_dp, 2.94125600e-5_dp, 3.50956800e-5_dp, &
         2.62857731e-5_dp, 2.64270283e-5_dp, 2.70548291e-5_dp, &
         2.78395800e-5_dp, 3.04311100e-5_dp, 3.39737300e-5_dp, &
         2.80772517e-5_dp, 2.81462344e-5_dp, 2.84528241e-5_dp, &
         2.88360612e-5_dp, 2.93062944e-5_dp, 3.36788644e-5_dp], &
         [6, 5])

    call gas%init()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       do ip = 1, 6
          Pi = pc(ip)
          do it = 1, 5
             Ti = t(it)
             call gas%viscosity(Pi, Ti, visc, err)
             write(s, '(a, e8.3, a, f6.2)') 'p = ', &
                  Pi, ', t = ', Ti
             call test%assert(0, err, trim(s) // ' error')
             call test%assert(expected_visc(ip, it), visc, s, tol)
          end do
       end do

    end if

  end subroutine test_ncg_co2_viscosity

!------------------------------------------------------------------------

  subroutine test_ncg_co2_properties(test)
    ! CO2 properties

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(ncg_co2_thermodynamics_type) :: gas
    PetscMPIInt :: rank
    PetscInt :: ierr, err, i
    PetscReal :: props(2)
    character(40) :: s
    PetscReal, parameter :: tol = 1.e-9_dp
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
            call test%assert(0, err, trim(s) // ' error')
            call test%assert(expected_h, h, trim(s) // ' enthalpy', tol)
            call test%assert(expected_d, d, trim(s) // ' density', tol)
          end associate
       end do

    end if

  end subroutine test_ncg_co2_properties

!------------------------------------------------------------------------

end module ncg_co2_thermodynamics_test
