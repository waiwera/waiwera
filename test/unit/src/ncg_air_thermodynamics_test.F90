module ncg_air_thermodynamics_test

  ! Tests for air NCG thermodynamics module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use ncg_air_thermodynamics_module
  use fruit

  implicit none
  private

  public :: test_ncg_air_enthalpy, test_ncg_air_henry, &
       test_ncg_air_energy_solution, test_ncg_air_mixture_viscosity

contains

!------------------------------------------------------------------------

  subroutine test_ncg_air_enthalpy
    ! Air enthalpy

    type(ncg_air_thermodynamics_type) :: gas
    PetscReal :: temperature, expected, props(2)
    PetscErrorCode :: err
    character(33) :: s = "Air enthalpy, temperature"
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: pa = 1.e5_dp
    PetscReal, parameter :: tol = 1.e-6_dp

    call gas%init()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       temperature = 20._dp
       expected = 19766.68296464_dp
       call gas%properties(pa, temperature, props, err)
       call assert_equals(0, err, trim(s) // " 20 deg C error")
       call assert_equals(expected, props(2), tol, trim(s) // " 20 deg C")

       temperature = 100._dp
       expected = 99758.41320142_dp
       call gas%properties(pa, temperature, props, err)
       call assert_equals(0, err, trim(s) // " 20 deg C error")
       call assert_equals(expected, props(2), tol, trim(s) // " 20 deg C")

       temperature = 240._dp
       expected = 243111.52063028_dp
       call gas%properties(pa, temperature, props, err)
       call assert_equals(0, err, trim(s) // " 20 deg C error")
       call assert_equals(expected, props(2), tol, trim(s) // " 20 deg C")

       temperature = 300._dp
       expected = 305841.67475311_dp
       call gas%properties(pa, temperature, props, err)
       call assert_equals(0, err, trim(s) // " 20 deg C error")
       call assert_equals(expected, props(2), tol, trim(s) // " 20 deg C")

       temperature = 350._dp
       expected = 358701.72866633_dp
       call gas%properties(pa, temperature, props, err)
       call assert_equals(0, err, trim(s) // " 20 deg C error")
       call assert_equals(expected, props(2), tol, trim(s) // " 20 deg C")

    end if

  end subroutine test_ncg_air_enthalpy

!------------------------------------------------------------------------

  subroutine test_ncg_air_henry
    ! Air Henry's constant

    type(ncg_air_thermodynamics_type) :: gas
    PetscReal :: temperature, expected, hc
    PetscErrorCode :: err
    character(33) :: s = "Air Henry's constant, temperature"
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: tol = 1.e-15_dp

    call gas%init()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       temperature = 20._dp
       expected = 1.37591939e-10_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(0, err, trim(s) // " 20 deg C error")
       call assert_equals(expected, hc, tol, trim(s) // " 20 deg C")

       temperature = 100._dp
       expected = 8.89730840e-11_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(0, err, trim(s) // " 100 deg C error")
       call assert_equals(expected, hc, tol, trim(s) // " 100 deg C")

       temperature = 240._dp
       expected = 2.52702897e-10_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(0, err, trim(s) // " 240 deg C error")
       call assert_equals(expected, hc, tol, trim(s) // " 240 deg C")

       temperature = 300._dp
       expected = 5.25847886e-10_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(0, err, trim(s) // " 300 deg C error")
       call assert_equals(expected, hc, tol, trim(s) // " 300 deg C")

       temperature = 350._dp
       expected = 1.63742619e-09_dp
       call gas%henrys_constant(temperature, hc, err)
       call assert_equals(0, err, trim(s) // " 300 deg C error")
       call assert_equals(expected, hc, tol, trim(s) // " 300 deg C")

    end if

  end subroutine test_ncg_air_henry

!------------------------------------------------------------------------

  subroutine test_ncg_air_energy_solution
    ! Air energy of solution

    type(ncg_air_thermodynamics_type) :: gas
    PetscReal :: temperature, expected, hs, hc
    PetscErrorCode :: err
    character(35) :: s = "Air energy of solution, temperature"
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: tol = 1.e-4_dp

    call gas%init()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       temperature = 20._dp
       expected = -400992.95145511_dp
       call gas%henrys_constant(temperature, hc, err)
       call gas%energy_solution(temperature, hc, hs, err)
       call assert_equals(0, err, trim(s) // " 20 deg C error")
       call assert_equals(expected, hs, tol, trim(s) // " 20 deg C")

       temperature = 100._dp
       expected = 58149.1888649_dp
       call gas%henrys_constant(temperature, hc, err)
       call gas%energy_solution(temperature, hc, hs, err)
       call assert_equals(0, err, trim(s) // " 100 deg C error")
       call assert_equals(expected, hs, tol, trim(s) // " 100 deg C")

       temperature = 240._dp
       expected = 845139.95259237_dp
       call gas%henrys_constant(temperature, hc, err)
       call gas%energy_solution(temperature, hc, hs, err)
       call assert_equals(0, err, trim(s) // " 240 deg C error")
       call assert_equals(expected, hs, tol, trim(s) // " 240 deg C")

       temperature = 300._dp
       expected = 1358812.70708159_dp
       call gas%henrys_constant(temperature, hc, err)
       call gas%energy_solution(temperature, hc, hs, err)
       call assert_equals(0, err, trim(s) // " 300 deg C error")
       call assert_equals(expected, hs, tol, trim(s) // " 300 deg C")

       temperature = 350._dp
       expected = 4390404.09402286_dp
       call gas%henrys_constant(temperature, hc, err)
       call gas%energy_solution(temperature, hc, hs, err)
       call assert_equals(0, err, trim(s) // " 350 deg C error")
       call assert_equals(expected, hs, tol, trim(s) // " 350 deg C")

    end if

  end subroutine test_ncg_air_energy_solution

!------------------------------------------------------------------------

    subroutine test_ncg_air_mixture_viscosity
    ! Water-air mixture viscosity

    use IFC67_module

    type(ncg_air_thermodynamics_type) :: gas
    PetscMPIInt :: rank
    PetscInt :: ierr, err
    PetscInt :: i
    PetscReal :: visc
    character(60) :: s
    type(IFC67_type) :: thermo
    PetscReal, parameter :: Pair = 0.e5_dp ! not used
    PetscInt, parameter :: phase = 2 ! vapour phase
    PetscReal, parameter :: tol = 1.e-9_dp
    PetscInt, parameter :: num_cases = 3
    PetscReal, parameter :: t(num_cases) = [240._dp, 120._dp, 20._dp]
    PetscReal, parameter :: xg(num_cases) = [0.1_dp, 0.5_dp, 0.8_dp]
    PetscReal, parameter :: water_viscosity(num_cases) = [ &
         0.171595480e-4_dp, 0.128139659e-4_dp, 8.73278989112e-6_dp]
    PetscReal, parameter :: expected_visc(num_cases) = [ &
         1.81800535828e-05_dp, 0.178797007e-4_dp, 1.67537163543e-5_dp]

    call thermo%init()
    call gas%init()

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       do i = 1, num_cases
          call gas%mixture_viscosity(water_viscosity(i), t(i), Pair, &
               xg(i), phase, visc, err)
          write(s, '(a, i2)') 'case ', i
          call assert_equals(0, err, trim(s) // ' error')
          call assert_equals(expected_visc(i), visc, tol, s)
       end do

    end if

    call thermo%destroy()

  end subroutine test_ncg_air_mixture_viscosity

!------------------------------------------------------------------------

end module ncg_air_thermodynamics_test
