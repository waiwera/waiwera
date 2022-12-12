module ncg_air_thermodynamics_test

  ! Tests for air NCG thermodynamics module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use ncg_air_thermodynamics_module
  use zofu

  implicit none
  private

  public :: setup, teardown
  public :: test_ncg_air_enthalpy, test_ncg_air_henry, &
       test_ncg_air_energy_solution, test_ncg_air_mixture_viscosity

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

  subroutine test_ncg_air_enthalpy(test)
    ! Air enthalpy

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(ncg_air_thermodynamics_type) :: gas
    PetscReal :: temperature, expected, props(2)
    PetscErrorCode :: err
    character(33) :: s = "Air enthalpy, temperature"
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: pa = 1.e5_dp

    call gas%init()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       temperature = 20._dp
       expected = 19766.68740112_dp
       call gas%properties(pa, temperature, props, err)
       call test%assert(0, err, trim(s) // " 20 deg C error")
       call test%assert(expected, props(2), trim(s) // " 20 deg C")

       temperature = 100._dp
       expected = 99758.4176379_dp
       call gas%properties(pa, temperature, props, err)
       call test%assert(0, err, trim(s) // " 20 deg C error")
       call test%assert(expected, props(2), trim(s) // " 20 deg C")

       temperature = 240._dp
       expected = 243111.52506676_dp
       call gas%properties(pa, temperature, props, err)
       call test%assert(0, err, trim(s) // " 20 deg C error")
       call test%assert(expected, props(2), trim(s) // " 20 deg C")

       temperature = 300._dp
       expected = 305841.67918959_dp
       call gas%properties(pa, temperature, props, err)
       call test%assert(0, err, trim(s) // " 20 deg C error")
       call test%assert(expected, props(2), trim(s) // " 20 deg C")

       temperature = 350._dp
       expected = 358701.73310281_dp
       call gas%properties(pa, temperature, props, err)
       call test%assert(0, err, trim(s) // " 20 deg C error")
       call test%assert(expected, props(2), trim(s) // " 20 deg C")

    end if

  end subroutine test_ncg_air_enthalpy

!------------------------------------------------------------------------

  subroutine test_ncg_air_henry(test)
    ! Air Henry constant

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(ncg_air_thermodynamics_type) :: gas
    PetscReal :: temperature, expected, hc, chc(2)
    PetscErrorCode :: err
    character(33) :: s = "Air Henry constant, temperature"
    PetscMPIInt :: rank
    PetscInt :: ierr

    call gas%init()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       temperature = 20._dp
       expected = 7.26786761e+09_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call test%assert(0, err, trim(s) // " 20 deg C error")
       call test%assert(expected, hc, trim(s) // " 20 deg C")

       call gas%henrys_constant_salt(temperature, 0._dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 20 deg C zero salt error")
       call test%assert(expected, hc, trim(s) // " 20 deg C zero salt")

       call gas%henrys_constant_salt(temperature, 0.1_dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 20 deg C 0.1 salt error")
       call test%assert(0.13689413e11_dp, hc, trim(s) // " 20 deg C 0.1 salt")

       call gas%henrys_constant_salt(temperature, 0.25_dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 20 deg C 0.25 salt error")
       call test%assert(0.48571791e11_dp, hc, trim(s) // " 20 deg C 0.25 salt")

       temperature = 100._dp
       expected = 1.12393541e+10_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call test%assert(0, err, trim(s) // " 100 deg C error")
       call test%assert(expected, hc, trim(s) // " 100 deg C")

       call gas%henrys_constant_salt(temperature, 0._dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 100 deg C zero salt error")
       call test%assert(expected, hc, trim(s) // " 100 deg C zero salt")

       call gas%henrys_constant_salt(temperature, 0.1_dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 100 deg C 0.1 salt error")
       call test%assert(0.19031197e11_dp, hc, trim(s) // " 100 deg C 0.1 salt")

       call gas%henrys_constant_salt(temperature, 0.3_dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 100 deg C 0.3 salt error")
       call test%assert(0.86620445e11_dp, hc, trim(s) // " 100 deg C 0.3 salt")

       temperature = 240._dp
       expected = 3.95721621e+09_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call test%assert(0, err, trim(s) // " 240 deg C error")
       call test%assert(expected, hc, trim(s) // " 240 deg C")

       temperature = 300._dp
       expected = 1.90169063e+09_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call test%assert(0, err, trim(s) // " 300 deg C error")
       call test%assert(expected, hc, trim(s) // " 300 deg C")

       call gas%henrys_constant_salt(temperature, 0._dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 300 deg C zero salt error")
       call test%assert(expected, hc, trim(s) // " 300 deg C zero salt")

       call gas%henrys_constant_salt(temperature, 0.1_dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 300 deg C 0.1 salt error")
       call test%assert(0.91717193e10_dp, hc, trim(s) // " 300 deg C 0.1 salt")

       call gas%henrys_constant_salt(temperature, 0.3_dp, hc, chc, err)
       call test%assert(0, err, trim(s) // " 300 deg C 0.3 salt error")
       call test%assert(0.99492925e12_dp, hc, trim(s) // " 300 deg C 0.3 salt")

       temperature = 350._dp
       expected = 6.10714550e+08_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call test%assert(0, err, trim(s) // " 300 deg C error")
       call test%assert(expected, hc, trim(s) // " 300 deg C")

    end if

  end subroutine test_ncg_air_henry

!------------------------------------------------------------------------

  subroutine test_ncg_air_energy_solution(test)
    ! Air energy of solution

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(ncg_air_thermodynamics_type) :: gas
    PetscReal :: temperature, expected, hs, hc, chc(2), xs
    PetscErrorCode :: err
    character(35) :: s = "Air energy of solution, temperature"
    PetscMPIInt :: rank
    PetscInt :: ierr

    call gas%init()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       temperature = 20._dp
       expected = -.40693012e6_dp
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
       expected = -.25858067e6_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 20 deg C 0.1 salt error")
       call test%assert(expected, hs, trim(s) // " 20 deg C 0.1 salt")

       xs = 0.2_dp
       expected = -.73143866e5_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 20 deg C 0.2 salt error")
       call test%assert(expected, hs, trim(s) // " 20 deg C 0.2 salt")

       xs = 0.3_dp
       expected = 0.16527489e6_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 20 deg C 0.3 salt error")
       call test%assert(expected, hs, trim(s) // " 20 deg C 0.3 salt")

       temperature = 100._dp
       expected = 53825.89_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call gas%energy_solution(temperature, chc, hs, err)
       call test%assert(0, err, trim(s) // " 100 deg C error")
       call test%assert(expected, hs, trim(s) // " 100 deg C")

       xs = 0.1_dp
       expected = -.23307100e5_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 100 deg C 0.1 salt error")
       call test%assert(expected, hs, trim(s) // " 100 deg C 0.1 salt")

       xs = 0.2_dp
       expected = -.11972333e6_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 100 deg C 0.2 salt error")
       call test%assert(expected, hs, trim(s) // " 100 deg C 0.2 salt")

       xs = 0.3_dp
       expected = -.24368706E+06_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 100 deg C 0.3 salt error")
       call test%assert(expected, hs, trim(s) // " 100 deg C 0.3 salt")

       temperature = 240._dp
       expected = 850547.65_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call gas%energy_solution(temperature, chc, hs, err)
       call test%assert(0, err, trim(s) // " 240 deg C error")
       call test%assert(expected, hs, trim(s) // " 240 deg C")

       temperature = 300._dp
       expected = 1378067.69_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call gas%energy_solution(temperature, chc, hs, err)
       call test%assert(0, err, trim(s) // " 300 deg C error")
       call test%assert(expected, hs, trim(s) // " 300 deg C")

       xs = 0.1_dp
       expected = 0.73588495e6_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 300 deg C 0.1 salt error")
       call test%assert(expected, hs, trim(s) // " 300 deg C 0.1 salt")

       xs = 0.2_dp
       expected = -.66843484e5_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 300 deg C 0.2 salt error")
       call test%assert(expected, hs, trim(s) // " 300 deg C 0.2 salt")

       xs = 0.3_dp
       expected = -.10989229e7_dp
       call gas%henrys_constant_salt(temperature, xs, hc, chc, err)
       call gas%energy_solution_salt(temperature, xs, chc, hs, err)
       call test%assert(0, err, trim(s) // " 300 deg C 0.3 salt error")
       call test%assert(expected, hs, trim(s) // " 300 deg C 0.3 salt")

       temperature = 350._dp
       expected = 4498309.62_dp
       call gas%henrys_constant(temperature, hc, chc, err)
       call gas%energy_solution(temperature, chc, hs, err)
       call test%assert(0, err, trim(s) // " 350 deg C error")
       call test%assert(expected, hs, trim(s) // " 350 deg C")

    end if

  end subroutine test_ncg_air_energy_solution

!------------------------------------------------------------------------

    subroutine test_ncg_air_mixture_viscosity(test)
    ! Water-air mixture viscosity

    use IFC67_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(ncg_air_thermodynamics_type) :: gas
    PetscMPIInt :: rank
    PetscInt :: ierr, err
    PetscInt :: i
    PetscReal :: visc
    character(60) :: s
    type(IFC67_type) :: thermo
    PetscReal, parameter :: Pair = 0.e5_dp ! not used
    PetscInt, parameter :: phase = 2 ! vapour phase
    PetscInt, parameter :: num_cases = 3
    PetscReal, parameter :: t(num_cases) = [240._dp, 120._dp, 20._dp]
    PetscReal, parameter :: xg(num_cases) = [0.1_dp, 0.5_dp, 0.8_dp]
    PetscReal, parameter :: water_viscosity(num_cases) = [ &
         0.171595480e-4_dp, 0.128139659e-4_dp, 8.73278989112e-6_dp]
    PetscReal, parameter :: expected_visc(num_cases) = [ &
         1.81800535828e-05_dp, 0.178797007e-4_dp, 1.67537163543e-5_dp]
    PetscReal, parameter :: tol = 1.e-5_dp

    call thermo%init()
    call gas%init()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       do i = 1, num_cases
          call gas%mixture_viscosity(water_viscosity(i), t(i), Pair, &
               xg(i), phase, visc, err)
          write(s, '(a, i2)') 'case ', i
          call test%assert(0, err, trim(s) // ' error')
          call test%assert(expected_visc(i), visc, s, tol)
       end do

    end if

    call thermo%destroy()

  end subroutine test_ncg_air_mixture_viscosity

!------------------------------------------------------------------------

end module ncg_air_thermodynamics_test
