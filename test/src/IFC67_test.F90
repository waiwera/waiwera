module IFC67_test

  ! Tests for IAPWS thermodynamics module

  use kinds_module
  use IFC67_module
  use thermodynamics_module, only: tc_k
  use fruit

  implicit none
  private

#include <petsc-finclude/petscsys.h>

  real(dp), parameter :: density_tol = 1.e-5_dp, energy_tol = 1.e-2_dp
  real(dp), parameter :: pressure_tol = 1.e-1_dp, temperature_tol = 1.e-6_dp
  real(dp), parameter :: viscosity_tol = 1.e-12_dp

  public :: test_IFC67_region1, test_IFC67_region2, test_IFC67_saturation

  contains

!------------------------------------------------------------------------

    subroutine test_IFC67_region1

      ! IFC-67 region 1 tests

      integer, parameter :: n = 3 
      real(dp) :: params(n,2) = reshape([ &
           3.e6_dp, 80.e6_dp, 3.e6_dp, &
           300._dp, 300._dp,  500._dp], [n,2])
      real(dp), parameter :: rho(n) = [997.95721560998174_dp, 1029.7256888266911_dp, &
           831.84196191567298_dp]
      real(dp), parameter ::  u(n) = [112247.43313085975_dp, 106310.47344628950_dp, &
           971985.91117384087_dp]
      integer :: i, err
      real(dp) :: param(2), props(2)
      PetscMPIInt :: rank
      PetscErrorCode :: ierr

      call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr); CHKERRQ(ierr)
      if (rank == 0) then
         params(:,2) = params(:,2) - tc_k  ! convert temperatures to Celcius
         do i = 1, n
            param = params(i,:)
            call IFC67%water%properties(param, props, err)
            call assert_equals(rho(i), props(1), density_tol, 'density')
            call assert_equals(u(i), props(2), energy_tol, 'energy')
            call assert_equals(0, err, 'error')
         end do
      end if

    end subroutine test_IFC67_region1

!------------------------------------------------------------------------

    subroutine test_IFC67_region2

      ! IFC-67 region 2 tests

      integer, parameter :: n = 3 
      real(dp) :: params(n,2) = reshape([ &
           0.0035e6_dp, 0.0035e6_dp, 30.e6_dp, &
           300._dp, 700._dp,  700._dp], [n,2])
      real(dp), parameter :: rho(n) = [2.5316826343790743e-2_dp, 1.0834441421293962e-2_dp, &
           183.90041953968711_dp]
      real(dp), parameter ::  u(n) = [2412405.0932077002_dp, 3012229.4965919587_dp, &
           2474981.3799304822_dp]
      integer :: i, err
      real(dp) :: param(2), props(2)
      real(dp) :: rho_ref(n), u_ref(n)
      PetscMPIInt :: rank
      PetscErrorCode :: ierr

      call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr); CHKERRQ(ierr)
      if (rank == 0) then
         params(:,2) = params(:,2) - tc_k  ! convert temperatures to Celcius
         do i = 1, n
            param = params(i,:)
            call IFC67%steam%properties(param, props, err)
            call assert_equals(rho(i), props(1), density_tol, 'density')
            call assert_equals(u(i), props(2), energy_tol, 'energy')
            call assert_equals(0, err, 'error')
         end do
      end if

    end subroutine test_IFC67_region2

!------------------------------------------------------------------------

    subroutine test_IFC67_saturation

      ! IFC67 saturation curve tests

      integer, parameter :: n = 3
      real(dp), parameter ::  t(n) = [300._dp, 500._dp, 600._dp] - tc_k
      real(dp), parameter :: p(n) = [0.35323426e4_dp, 0.263961572e7_dp, &
           0.123493902e8_dp]
      real(dp) :: ps, ts
      integer :: i, err
      PetscMPIInt :: rank
      PetscErrorCode :: ierr

      call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr); CHKERRQ(ierr)
      if (rank == 0) then
         do i = 1, n
            call IFC67%saturation%pressure(t(i), ps, err)
            call assert_equals(p(i), ps, pressure_tol, 'pressure')
            call assert_equals(0, err, 'pressure error')
            call IFC67%saturation%temperature(ps, ts, err)
            call assert_equals(t(i), ts, temperature_tol, 'temperature')
            call assert_equals(0, err, 'temperature error')
         end do
      end if
      
    end subroutine test_IFC67_saturation

!------------------------------------------------------------------------

end module IFC67_test
