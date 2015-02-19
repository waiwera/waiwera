module IAPWS_test

  ! Tests for IAPWS thermodynamics module

  use kinds_module
  use mpi_module
  use IAPWS_module
  use thermodynamics_module, only: tc_k
  use fruit

  implicit none
  private

#include <petsc-finclude/petscsys.h>

  real(dp), parameter :: density_tol = 1.e-5_dp, energy_tol = 1.e-2_dp
  real(dp), parameter :: pressure_tol = 1.e-1_dp, temperature_tol = 1.e-6_dp
  real(dp), parameter :: viscosity_tol = 1.e-12_dp

  public :: test_IAPWS_region1, test_IAPWS_region2, test_IAPWS_region3, &
       test_IAPWS_saturation, test_IAPWS_viscosity, test_IAPWS_boundary23

  contains

!------------------------------------------------------------------------

    subroutine test_IAPWS_region1

      ! IAPWS-97 region 1 tests

      integer, parameter :: n = 3 
      real(dp) :: params(n,2) = reshape([ &
           3.e6_dp, 80.e6_dp, 3.e6_dp, &
           300._dp, 300._dp,  500._dp], [n,2])
      real(dp), parameter :: nu(n) = [0.100215168e-2_dp, 0.971180894e-3_dp, 0.120241800e-2_dp]
      real(dp), parameter ::  u(n) = [0.112324818e6_dp,  0.106448356e6_dp,  0.971934985e6_dp]
      real(dp), parameter :: rho(n) = 1._dp / nu
      integer :: i, err
      real(dp) :: param(2), props(2)

      if (mpi%rank == mpi%output_rank) then
         params(:,2) = params(:,2) - tc_k  ! convert temperatures to Celcius
         do i = 1, n
            param = params(i,:)
            call IAPWS%water%properties(param, props, err)
            call assert_equals(rho(i), props(1), density_tol, 'density')
            call assert_equals(u(i), props(2), energy_tol, 'energy')
            call assert_equals(0, err, 'error')
         end do
      end if

    end subroutine test_IAPWS_region1

!------------------------------------------------------------------------

    subroutine test_IAPWS_region2

      ! IAPWS-97 region 2 tests

      integer, parameter :: n = 3 
      real(dp) :: params(n,2) = reshape([ &
           0.0035e6_dp, 0.0035e6_dp, 30.e6_dp, &
           300._dp, 700._dp,  700._dp], [n,2])
      real(dp), parameter :: nu(n) = [0.394913866e2_dp, 0.923015898e2_dp, 0.542946619e-2_dp]
      real(dp), parameter ::  u(n) = [0.241169160e7_dp, 0.301262819e7_dp, 0.246861076e7_dp]
      real(dp), parameter :: rho(n) = 1._dp / nu
      integer :: i, err
      real(dp) :: param(2), props(2)

      if (mpi%rank == mpi%output_rank) then
         params(:,2) = params(:,2) - tc_k  ! convert temperatures to Celcius
         do i = 1, n
            param = params(i,:)
            call IAPWS%steam%properties(param, props, err)
            call assert_equals(rho(i), props(1), density_tol, 'density')
            call assert_equals(u(i), props(2), energy_tol, 'energy')
            call assert_equals(0, err, 'error')
         end do
      end if

    end subroutine test_IAPWS_region2

!------------------------------------------------------------------------

    subroutine test_IAPWS_region3

      ! IAPWS-97 region 3 tests

      integer, parameter :: n = 3 
      real(dp) :: params(n,2) = reshape([ &
           500._dp, 200._dp, 500._dp, &
           650._dp, 650._dp, 750._dp], [n,2])
      real(dp), parameter :: p(n) = [0.255837018e8_dp, 0.222930643e8_dp, 0.783095639e8_dp]
      real(dp), parameter ::  u(n) = [0.181226279e7_dp, 0.226365868e7_dp, 0.210206932e7_dp]
      integer :: i, err
      real(dp) :: param(2), props(2)

      if (mpi%rank == mpi%output_rank) then
         params(:,2) = params(:,2) - tc_k  ! convert temperatures to Celcius
         do i = 1, n
            param = params(i,:)
            call IAPWS%supercritical%properties(param, props, err)
            call assert_equals(p(i), props(1), pressure_tol, 'pressure')
            call assert_equals(u(i), props(2), energy_tol, 'energy')
            call assert_equals(0, err, 'error')
         end do
      end if

    end subroutine test_IAPWS_region3

!------------------------------------------------------------------------

    subroutine test_IAPWS_saturation

      ! IAPWS-97 saturation curve tests

      integer, parameter :: n = 3
      real(dp), parameter ::  t(n) = [300._dp, 500._dp, 600._dp] - tc_k
      real(dp), parameter :: p(n) = [0.353658941e4_dp, 0.263889776e7_dp, &
           0.123443146e8_dp]
      real(dp) :: ps, ts
      integer :: i, err

      if (mpi%rank == mpi%output_rank) then
         do i = 1, n
            call IAPWS%saturation%pressure(t(i), ps, err)
            call assert_equals(p(i), ps, pressure_tol, 'pressure')
            call assert_equals(0, err, 'pressure error')
            call IAPWS%saturation%temperature(ps, ts, err)
            call assert_equals(t(i), ts, temperature_tol, 'temperature')
            call assert_equals(0, err, 'temperature error')
         end do
      end if

    end subroutine test_IAPWS_saturation

!------------------------------------------------------------------------

    subroutine test_IAPWS_viscosity

      ! IAPWS viscosity tests

      integer, parameter :: n = 11
      real(dp), parameter :: t(n) = [298.15_dp, 298.15_dp, 373.15_dp, 433.15_dp, 433.15_dp, &
           873.15_dp, 873.15_dp, 873.15_dp, 1173.15_dp, 1173.15_dp, 1173.15_dp] - tc_k
      real(dp), parameter :: d(n) = [998._dp, 1200._dp, 1000._dp, 1._dp, &
           1000._dp, 1._dp, 100._dp, 600._dp, 1._dp, 100._dp, 400._dp]
      real(dp), parameter :: visc(n) = [889.735100_dp, 1437.649467_dp, 307.883622_dp, &
           14.538324_dp, 217.685358_dp, 32.619287_dp, 35.802262_dp, 77.430195_dp, &
           44.217245_dp, 47.640433_dp, 64.154608_dp] * 1.e-6_dp
      real(dp) :: v
      integer :: i

      if (mpi%rank == mpi%output_rank) then
         do i = 1, n
            call IAPWS%viscosity(d(i), t(i), v)
            call assert_equals(visc(i), v, viscosity_tol)
         end do
      end if

    end subroutine test_IAPWS_viscosity

!------------------------------------------------------------------------

    subroutine test_IAPWS_boundary23

      ! Boundary 23 test
      
      real(dp), parameter :: t0 = 0.62315e3_dp - tc_k
      real(dp), parameter :: p0 = 0.165291643e8_dp
      real(dp) :: p, t

      if (mpi%rank == mpi%output_rank) then
         call IAPWS%boundary23%pressure(t0, p)
         call assert_equals(p0, p, pressure_tol, "Pressure")
         call IAPWS%boundary23%temperature(p0, t)
         call assert_equals(t0, t, temperature_tol, "Temperature")
      end if

    end subroutine test_IAPWS_boundary23

!------------------------------------------------------------------------

end module IAPWS_test
