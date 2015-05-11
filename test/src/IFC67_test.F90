module IFC67_test

  ! Tests for IAPWS thermodynamics module

  use kinds_module
  use mpi_module
  use IFC67_module
  use thermodynamics_module, only: tc_k
  use fruit

  implicit none
  private

#include <petsc-finclude/petscsys.h>

  real(dp), parameter :: density_tol = 1.e-5_dp, energy_tol = 1.e-2_dp
  real(dp), parameter :: pressure_tol = 1.e-1_dp, temperature_tol = 1.e-6_dp
  real(dp), parameter :: viscosity_tol = 1.e-9_dp

  public :: setup_IFC67, teardown_IFC67
  public :: test_IFC67_region1, test_IFC67_region2
  public :: test_IFC67_saturation, test_IFC67_viscosity

  contains

!------------------------------------------------------------------------

    subroutine setup_IFC67

      ! Sets up IFC67 tests

      call IFC67%init()

    end subroutine setup_IFC67

!------------------------------------------------------------------------

    subroutine teardown_IFC67

      ! Tears down IFC67 tests

      call IFC67%destroy()

    end subroutine teardown_IFC67

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

      if (mpi%rank == mpi%output_rank) then
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

      if (mpi%rank == mpi%output_rank) then
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

      ! IFC-67 saturation curve tests

      integer, parameter :: n = 3
      real(dp), parameter ::  t(n) = [300._dp, 500._dp, 600._dp] - tc_k
      real(dp), parameter :: p(n) = [0.35323426e4_dp, 0.263961572e7_dp, &
           0.123493902e8_dp]
      real(dp) :: ps, ts
      integer :: i, err

      if (mpi%rank == mpi%output_rank) then
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

    subroutine test_IFC67_viscosity

      ! IFC-67 viscosity tests

      integer, parameter :: n1 = 2, n2 = 2
      ! region 1 tests:
      real(dp), parameter :: t1(n1) = [298.15_dp, 373.15_dp] - tc_k
      real(dp), parameter :: p1(n1) = [1977563.58349_dp, 99834578.2816_dp]
      real(dp), parameter :: visc1(n1) = [8.903129e-04_dp, 2.988268e-04_dp]
      ! region 2 tests:
      real(dp), parameter :: t2(n2) = [873.15_dp, 873.15_dp] - tc_k
      real(dp), parameter :: d2(n2) = [1._dp, 100._dp]
      real(dp), parameter :: visc2(n2) = [3.249537e-05_dp, 3.667671e-05_dp]
      real(dp) :: v, p = 0.0_dp, d = 0.0_dp
      integer :: i

      if (mpi%rank == mpi%output_rank) then
         do i = 1, n1
            call IFC67%water%viscosity(t1(i), p1(i), d, v)
            call assert_equals(visc1(i), v, viscosity_tol)
         end do
         do i = 1, n2
            call IFC67%steam%viscosity(t2(i), p, d2(i), v)
            call assert_equals(visc2(i), v, viscosity_tol)
         end do
      end if

    end subroutine test_IFC67_viscosity

!------------------------------------------------------------------------

end module IFC67_test
