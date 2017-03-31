module ncg_air_thermodynamics_test

  ! Tests for air NCG thermodynamics module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use ncg_air_thermodynamics_module
  use fruit

  implicit none
  private

  public :: test_ncg_air_mixture_viscosity

contains

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
