module IAPWS_test

  ! Tests for IAPWS thermodynamics module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use IAPWS_module
  use thermodynamics_module, only: tc_k
  use zofu

  implicit none
  private

  type(IAPWS_type) :: IAPWS

  public :: setup, teardown, setup_test
  public :: test_IAPWS_region1, test_IAPWS_region2, test_IAPWS_region3, &
       test_IAPWS_saturation, test_IAPWS_viscosity, test_IAPWS_boundary23, &
       test_IAPWS_phase_composition, test_IAPWS_region3_subbdy, &
       test_IAPWS_region3_dpdd, test_IAPWS_region3_density, &
       test_IAPWS_region3_saturation_density, &
       test_region3_widom, test_region3_pi_liquidlike, &
       test_IAPWS_region1_pressure, test_IAPWS_region2_pressure

  contains

!------------------------------------------------------------------------

  subroutine setup()

    use profiling_module, only: init_profiling
    use fson
    use fson_mpi_module

    ! Locals:
    type(fson_value), pointer :: json
    character(120) :: json_str = &
         '{"thermodynamics": {"widom": {"delta": {"growth": 25, "min": 0.1}}}}'
    PetscErrorCode :: ierr

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call init_profiling()
    json => fson_parse_mpi(str = json_str)
    call IAPWS%init(json)
    call fson_destroy_mpi(json)

  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown()

    PetscErrorCode :: ierr

    call IAPWS%destroy()
    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

  subroutine setup_test(test)

    class(unit_test_type), intent(in out) :: test

    test%tolerance = 1.e-7

  end subroutine setup_test

!------------------------------------------------------------------------

  subroutine test_IAPWS_region1(test)

    ! IAPWS-97 region 1 tests

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: n = 3, nerr = 2
    PetscReal :: params(n,2) = reshape([ &
         3.e6_dp, 80.e6_dp, 3.e6_dp, &
         300._dp, 300._dp,  500._dp], [n,2])
    PetscReal, parameter :: nu(n) = [0.100215168e-2_dp, 0.971180894e-3_dp, 0.120241800e-2_dp]
    PetscReal, parameter ::  u(n) = [0.112324818e6_dp,  0.106448356e6_dp,  0.971934985e6_dp]
    PetscReal, parameter :: rho(n) = 1._dp / nu
    PetscInt :: i, err
    PetscReal :: param(2), props(2)
    PetscReal :: err_params(nerr,2) = reshape([ &
         20.e6_dp, 101.e6_dp, &
         360._dp,    60._dp], [nerr,2])
    PetscInt :: ierr
    PetscMPIInt :: rank

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       params(:,2) = params(:,2) - tc_k  ! convert temperatures to Celcius
       do i = 1, n
          param = params(i,:)
          call IAPWS%water%properties(param, props, err)
          call test%assert(rho(i), props(1), 'density')
          call test%assert(u(i), props(2), 'energy')
          call test%assert(0, err, 'no error')
       end do
       do i = 1, nerr
          param = err_params(i,:)
          call IAPWS%water%properties(param, props, err)
          call test%assert(1, err, 'error')
       end do
    end if

  end subroutine test_IAPWS_region1

!------------------------------------------------------------------------

  subroutine test_IAPWS_region2(test)

    ! IAPWS-97 region 2 tests

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: n = 3, nerr = 2
    PetscReal :: params(n,2) = reshape([ &
         0.0035e6_dp, 0.0035e6_dp, 30.e6_dp, &
         300._dp, 700._dp,  700._dp], [n,2])
    PetscReal, parameter :: nu(n) = [0.394913866e2_dp, 0.923015898e2_dp, 0.542946619e-2_dp]
    PetscReal, parameter ::  u(n) = [0.241169160e7_dp, 0.301262819e7_dp, 0.246861076e7_dp]
    PetscReal, parameter :: rho(n) = 1._dp / nu
    PetscInt :: i, err
    PetscReal :: param(2), props(2)
    PetscReal :: err_params(nerr,2) = reshape([ &
         20.e6_dp, 101.e6_dp, &
         801._dp,    60._dp], [nerr,2])
    PetscInt :: ierr
    PetscMPIInt :: rank

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       params(:,2) = params(:,2) - tc_k  ! convert temperatures to Celcius
       do i = 1, n
          param = params(i,:)
          call IAPWS%steam%properties(param, props, err)
          call test%assert(rho(i), props(1), 'density')
          call test%assert(u(i), props(2), 'energy')
          call test%assert(0, err, 'error')
       end do
       do i = 1, nerr
          param = err_params(i,:)
          call IAPWS%steam%properties(param, props, err)
          call test%assert(1, err, 'error')
       end do
    end if

  end subroutine test_IAPWS_region2

!------------------------------------------------------------------------

  subroutine test_IAPWS_region3(test)

    ! IAPWS-97 region 3 tests

    class(unit_test_type), intent(in out) :: test
    PetscInt, parameter :: n = 3, nerr = 4
    PetscReal :: params(n,2) = reshape([ &
         500._dp, 200._dp, 500._dp, &
         650._dp, 650._dp, 750._dp], [n,2])
    PetscReal, parameter :: p(n) = [0.255837018e8_dp, 0.222930643e8_dp, 0.783095639e8_dp]
    PetscReal, parameter ::  u(n) = [0.181226279e7_dp, 0.226365868e7_dp, 0.210206932e7_dp]
    PetscInt :: i, err
    PetscReal :: param(2), props(2)
    PetscReal :: err_params(nerr,2) = reshape([ &
         800._dp, 500._dp, 100._dp, 300._dp, &
         400._dp, 340._dp, 380._dp, 600._dp], [nerr, 2])
    PetscInt :: ierr
    PetscMPIInt :: rank

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       params(:,2) = params(:,2) - tc_k  ! convert temperatures to Celcius
       do i = 1, n
          param = params(i,:)
          call IAPWS%supercritical%properties(param, props, err)
          call test%assert(p(i), props(1), 'pressure')
          call test%assert(u(i), props(2), 'energy')
          call test%assert(0, err, 'error')
       end do
       do i = 1, nerr
          param = err_params(i,:)
          call IAPWS%supercritical%properties(param, props, err)
          call test%assert(1, err, 'error')
       end do
    end if

  end subroutine test_IAPWS_region3

!------------------------------------------------------------------------

  subroutine test_IAPWS_saturation(test)

    ! IAPWS-97 saturation curve tests
    
    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: n = 3, nerr = 1
    PetscReal, parameter ::  t(n) = [300._dp, 500._dp, 600._dp] - tc_k
    PetscReal, parameter :: p(n) = [0.353658941e4_dp, 0.263889776e7_dp, &
         0.123443146e8_dp]
    PetscReal :: ps, ts
    PetscInt :: i, err
    PetscReal :: terr(nerr) = [380._dp], perr(nerr) = [30.e6_dp]
    PetscInt :: ierr
    PetscMPIInt :: rank

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       do i = 1, n
          call IAPWS%saturation%pressure(t(i), ps, err)
          call test%assert(p(i), ps, 'pressure')
          call test%assert(0, err, 'pressure no error')
          call IAPWS%saturation%temperature(ps, ts, err)
          call test%assert(t(i), ts, 'temperature')
          call test%assert(0, err, 'temperature no error')
       end do
       do i = 1, nerr
          call IAPWS%saturation%pressure(terr(i), ps, err)
          call test%assert(1, err, 'pressure error')
          call IAPWS%saturation%temperature(perr(i), ts, err)
          call test%assert(1, err, 'temperature error')
       end do
    end if

  end subroutine test_IAPWS_saturation

!------------------------------------------------------------------------

  subroutine test_IAPWS_viscosity(test)

    ! IAPWS viscosity tests

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: n = 11
    PetscReal, parameter :: t(n) = [298.15_dp, 298.15_dp, 373.15_dp, 433.15_dp, 433.15_dp, &
         873.15_dp, 873.15_dp, 873.15_dp, 1173.15_dp, 1173.15_dp, 1173.15_dp] - tc_k
    PetscReal, parameter :: d(n) = [998._dp, 1200._dp, 1000._dp, 1._dp, &
         1000._dp, 1._dp, 100._dp, 600._dp, 1._dp, 100._dp, 400._dp]
    PetscReal, parameter :: visc(n) = [889.735100_dp, 1437.649467_dp, 307.883622_dp, &
         14.538324_dp, 217.685358_dp, 32.619287_dp, 35.802262_dp, 77.430195_dp, &
         44.217245_dp, 47.640433_dp, 64.154608_dp] * 1.e-6_dp
    PetscInt, parameter :: reg(n) = [1, 1, 1, 2, 1, 2, 2, 3, 2, 2, 2]
    PetscReal, parameter :: p = 1.e5 ! dummy pressure
    PetscReal :: v
    PetscInt :: i
    PetscInt :: ierr
    PetscMPIInt :: rank

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       do i = 1, n
          call IAPWS%region(reg(i))%ptr%viscosity(t(i), p, d(i), v)
          call test%assert(visc(i), v)
       end do
    end if

  end subroutine test_IAPWS_viscosity

!------------------------------------------------------------------------

  subroutine test_IAPWS_boundary23(test)

    ! Boundary 23 test
      
    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscReal, parameter :: t0 = 0.62315e3_dp - tc_k
    PetscReal, parameter :: p0 = 0.165291643e8_dp
    PetscReal :: p, t
    PetscInt :: ierr
    PetscMPIInt :: rank

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       call IAPWS%boundary23%pressure(t0, p)
       call test%assert(p0, p, "Pressure")
       call IAPWS%boundary23%temperature(p0, t)
       call test%assert(t0, t, "Temperature")
    end if

  end subroutine test_IAPWS_boundary23

!------------------------------------------------------------------------

  subroutine test_IAPWS_phase_composition(test)

    ! Phase composition tests

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt :: phases, expected_phases
    PetscInt :: ierr
    PetscMPIInt :: rank

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       phases = IAPWS%phase_composition(1, 1.e5_dp, 20._dp)
       expected_phases = int(b'001')
       call test%assert(expected_phases, phases, &
            "Region 1 liquid")

       phases = IAPWS%phase_composition(3, 200.e5_dp, 360._dp)
       expected_phases = int(b'001')
       call test%assert(expected_phases, phases, &
            "Region 3 liquid")

       phases = IAPWS%phase_composition(2, 1.e5_dp, 110._dp)
       expected_phases = int(b'010')
       call test%assert(expected_phases, phases, &
            "Region 2 steam below critical temperature")

       phases = IAPWS%phase_composition(2, 175.e5_dp, 360._dp)
       expected_phases = int(b'010')
       call test%assert(expected_phases, phases, &
            "Region 2 steam above critical temperature")

       phases = IAPWS%phase_composition(2, 150.e5_dp, 700._dp)
       expected_phases = int(b'010')
       call test%assert(expected_phases, phases, &
            "High temperature region 2 steam")

       phases = IAPWS%phase_composition(3, 180.e5_dp, 360._dp)
       expected_phases = int(b'010')
       call test%assert(expected_phases, phases, &
            "Region 3 steam below critical temperature")

       phases = IAPWS%phase_composition(3, 210.e5_dp, 380._dp)
       expected_phases = int(b'010')
       call test%assert(expected_phases, phases, &
            "Region 3 steam above critical temperature")

       phases = IAPWS%phase_composition(4, 33.466518715101621e5_dp, 240._dp)
       expected_phases = int(b'011')
       call test%assert(expected_phases, phases, &
            "Two-phase at 240 deg C")

       phases = IAPWS%phase_composition(3, 500.e5_dp, 390._dp)
       expected_phases = int(b'100')
       call test%assert(expected_phases, phases, &
            "Region 3 supercritical, 390 deg C")

       phases = IAPWS%phase_composition(3, 560.e5_dp, 500._dp)
       expected_phases = int(b'100')
       call test%assert(expected_phases, phases, &
            "Region 3 supercritical, 500 deg C")

       phases = IAPWS%phase_composition(2, 300.e5_dp, 700._dp)
       expected_phases = int(b'100')
       call test%assert(expected_phases, phases, &
            "Region 2 supercritical")

    end if

  end subroutine test_IAPWS_phase_composition

!------------------------------------------------------------------------

  subroutine test_IAPWS_region3_subbdy(test)
    ! Region 3 subregion boundary tests
    ! From tables 3 and 11 of IAPWS (2014)

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscReal, parameter :: MPa = 1.e6_dp
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal :: tb

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       select type (region3 => IAPWS%supercritical)
       type is (IAPWS_region3_type)

          call logpoly_test(40._dp, 6.930341408e2_dp, region3, &
               region3%subregion_bdy_n_ab, region3%subregion_bdy_ninv_ab, 'ab')
          call poly_test(25._dp, 6.493659208e2_dp, region3, &
               region3%subregion_bdy_n_cd, 'cd')
          tb = region3%subregion_boundary_3ef(40._dp * MPa)
          call test%assert(7.139593992e2_dp, tb, 'IAPWS region 3 subregion bdy ef')
          call poly_test(23._dp, 6.498873759e2_dp, region3, &
               region3%subregion_bdy_n_gh, 'gh')
          call poly_test(23._dp, 6.515778091e2_dp, region3, &
               region3%subregion_bdy_n_ij, 'ij')
          call poly_test(23._dp, 6.558338344e2_dp, region3, &
               region3%subregion_bdy_n_jk, 'jk')
          call poly_test(22.8_dp, 6.496054133e2_dp, region3, &
               region3%subregion_bdy_n_mn, 'mn')
          call logpoly_test(22.8_dp, 6.500106943e2_dp, region3, &
               region3%subregion_bdy_n_op, region3%subregion_bdy_ninv_op, 'op')
          call poly_test(22._dp, 6.456355027e2_dp, region3, &
               region3%subregion_bdy_n_qu, 'qu')
          call poly_test(22._dp, 6.482622754e2_dp, region3, &
               region3%subregion_bdy_n_rx, 'rx')
          call poly_test(22.3_dp, 6.477996121e2_dp, region3, &
               region3%subregion_bdy_n_uv, 'uv')
          call logpoly_test(22.3_dp, 6.482049480e2_dp, region3, &
               region3%subregion_bdy_n_wx, region3%subregion_bdy_ninv_wx, 'wx')

       end select
    end if

  contains

    subroutine logpoly_test(p, tk, region3, poly, invpoly, name)

      PetscReal, intent(in) :: p, tk
      type(IAPWS_region3_type), intent(in) :: region3
      PetscReal, intent(in) :: poly(:), invpoly(:)
      character(*), intent(in) :: name
      ! Locals:
      PetscReal :: tb

      tb = region3%subregion_boundary_logpoly(poly, invpoly, p * MPa)
      call test%assert(tk, tb, 'IAPWS region 3 subregion bdy ' // name)

    end subroutine logpoly_test

    subroutine poly_test(p, tk, region3, poly, name)

      PetscReal, intent(in) :: p, tk
      type(IAPWS_region3_type), intent(in) :: region3
      PetscReal, intent(in) :: poly(:)
      character(*), intent(in) :: name
      ! Locals:
      PetscReal :: tb

      tb = region3%subregion_boundary_poly(poly, p * MPa)
      call test%assert(tk, tb, 'IAPWS region 3 subregion bdy ' // name)

    end subroutine poly_test

  end subroutine test_IAPWS_region3_subbdy

!------------------------------------------------------------------------

  subroutine test_IAPWS_region3_dpdd(test)
    ! Region 3 dp/dd tests

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       select type (region3 => IAPWS%supercritical)
       type is (IAPWS_region3_type)
          call dpdd_test(region3, [600._dp, 400._dp], 'case 1')
          call dpdd_test(region3, [700._dp, 360._dp], 'case 2')
          call dpdd_test(region3, [400._dp, 500._dp], 'case 3')
          call dpdd_test(region3, [200._dp, 380._dp], 'case 4')
       end select
    end if

  contains

    subroutine dpdd_test(region3, param, name)
      ! Test dp/dd against finite difference approximation

      type(IAPWS_region3_type), intent(in out) :: region3
      PetscReal, intent(in) :: param(2)
      character(*), intent(in) :: name
      ! Locals:
      PetscReal :: dp_drho, dp_drho2, del_rho
      PetscReal :: rho2, props(2), p, p2
      PetscErrorCode :: err
      PetscReal, parameter :: drho = 1e-8_dp

      dp_drho = region3%dpdd(param)

      associate(rho => param(1), t => param(2))
        call region3%properties(param, props, err)
        p = props(1)
        del_rho = drho * rho
        rho2 = rho + del_rho
        call region3%properties([rho2, t], props, err)
        p2 = props(1)
        dp_drho2 = (p2 - p) / del_rho
      end associate

      call test%assert(dp_drho2, dp_drho, &
           'IAPWS region 3 dp/dd ' // name, tol = 1.e-5_dp)

    end subroutine dpdd_test

  end subroutine test_IAPWS_region3_dpdd

!------------------------------------------------------------------------

  subroutine test_IAPWS_region3_density(test)
    ! Region 3 density tests
    ! From tables 5 and 13 of IAPWS (2014)

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: n = 52
    PetscReal :: data(n, 3) = transpose(reshape([ &
         50._dp, 630._dp, 1.470853100e-3_dp, &
         80._dp, 670._dp, 1.503831359e-3_dp, &
         50._dp, 710._dp, 2.204728587e-3_dp, &
         80._dp, 750._dp, 1.973692940e-3_dp, &
         20._dp, 630._dp, 1.761696406e-3_dp, &
         30._dp, 650._dp, 1.819560617e-3_dp, &
         26._dp, 656._dp, 2.245587720e-3_dp, &
         30._dp, 670._dp, 2.506897702e-3_dp, &
         26._dp, 661._dp, 2.970225962e-3_dp, &
         30._dp, 675._dp, 3.004627086e-3_dp, &
         26._dp, 671._dp, 5.019029401e-3_dp, &
         30._dp, 690._dp, 4.656470142e-3_dp, &
         23.6_dp, 649._dp, 2.163198378e-3_dp, &
         24._dp, 650._dp, 2.166044161e-3_dp, &
         23.6_dp, 652._dp, 2.651081407e-3_dp, &
         24._dp, 654._dp, 2.967802335e-3_dp, &
         23.6_dp, 653._dp, 3.273916816e-3_dp, &
         24._dp, 655._dp, 3.550329864e-3_dp, &
         23.5_dp, 655._dp, 4.545001142e-3_dp, &
         24._dp, 660._dp, 5.100267704e-3_dp, &
         23._dp, 660._dp, 6.109525997e-3_dp, &
         24._dp, 670._dp, 6.427325645e-3_dp, &
         22.6_dp, 646._dp, 2.117860851e-3_dp, &
         23._dp, 646._dp, 2.062374674e-3_dp, &
         22.6_dp, 648.6_dp, 2.533063780e-3_dp, &
         22.8_dp, 649.3_dp, 2.572971781e-3_dp, &
         22.6_dp, 649.0_dp, 2.923432711e-3_dp, &
         22.8_dp, 649.7_dp, 2.913311494e-3_dp, &
         22.6_dp, 649.1_dp, 3.131208996e-3_dp, &
         22.8_dp, 649.9_dp, 3.221160278e-3_dp, &
         22.6_dp, 649.4_dp, 3.715596186e-3_dp, &
         22.8_dp, 650.2_dp, 3.664754790e-3_dp, &
         21.1_dp, 640._dp, 1.970999272e-3_dp, &
         21.8_dp, 643._dp, 2.043919161e-3_dp, &
         21.1_dp, 644._dp, 5.251009921e-3_dp, &
         21.8_dp, 648._dp, 5.256844741e-3_dp, &
         19.1_dp, 635._dp, 1.932829079e-3_dp, &
         20._dp, 638._dp, 1.985387227e-3_dp, &
         17._dp, 626._dp, 8.483262001e-3_dp, &
         20._dp, 640._dp, 6.227528101e-3_dp, &
         21.5_dp, 644.6_dp, 2.268366647e-3_dp, &
         22._dp, 646.1_dp,  2.296350553e-3_dp, &
         22.5_dp, 648.6_dp,  2.832373260e-3_dp, &
         22.3_dp, 647.9_dp,  2.811424405e-3_dp, &
         22.15_dp, 647.5_dp,  3.694032281e-3_dp, &
         22.3_dp, 648.1_dp,  3.622226305e-3_dp, &
         22.11_dp, 648._dp,  4.528072649e-3_dp, &
         22.3_dp, 649._dp,  4.556905799e-3_dp, &
         22._dp, 646.84_dp,  2.698354719e-3_dp, &
         22.064_dp, 647.05_dp,  2.717655648e-3_dp, &
         22._dp, 646.89_dp,  3.798732962e-3_dp, &
         22.064_dp, 647.15_dp,  3.701940010e-3_dp], &
         [3, n]))
    PetscInt, parameter :: subregion(n) = [ &
         1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, &
         8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, &
         14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, &
         21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26]
    PetscReal, parameter :: MPa = 1.e6_dp
    PetscInt :: i, err, sr
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal :: param(2), density, p, props(2)
    character(2) :: istr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       data(:,1) = data(:,1) * MPa
       data(:,2) = data(:,2) - tc_k
       do i = 1, n
          write(istr, '(i2)') i
          param = data(i, 1:2)
          select type (region3 => IAPWS%region(3)%ptr)
          type is (IAPWS_region3_type)

             sr = region3%subregion_index(param)
             call test%assert(subregion(i), sr, 'subregion' // istr)
             call region3%density(param, density, err, polish = PETSC_FALSE)
             call test%assert(0, err, 'error' // istr)

             if (err == 0) then
                call test%assert(1._dp / data(i,3), density, 'density' // istr)
                p = data(i, 1)
                param(1) = density
                call region3%properties(param, props, err)
                associate(pback => props(1))
                  call test%assert(p, pback, 'p back' // istr, tol = 3.e-4_dp)
                end associate
             end if

             param = data(i, 1:2)
             call region3%density(param, density, err, polish = PETSC_TRUE)
             call test%assert(0, err, 'polished error' // istr)
             if (err == 0) then
                p = data(i, 1)
                param(1) = density
                call region3%properties(param, props, err)
                associate(pback => props(1))
                  call test%assert(p, pback, 'p back polished' // istr)
                end associate
             end if

          end select
       end do
    end if

  end subroutine test_IAPWS_region3_density

!------------------------------------------------------------------------

  subroutine test_IAPWS_region3_saturation_density(test)
    ! Region 3 saturation density tests

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: n = 10
    PetscReal, parameter :: pressure(n) = &
         [16.52916425260448e6_dp, 19.00881189173929e6_dp, &
         20.5e6_dp, 21.e6_dp, 21.9e6_dp, 21.99e6_dp, &
         22.0639911177352e6_dp, 22.063994473664202e6_dp, &
         22.063991311720084e6_dp, 22.064e6_dp]
    PetscInt, parameter :: liquid_index(n) = &
         [3, 3, 19, 19, 21, 25, 25, 25, 25, 25]
    PetscInt, parameter :: vapour_index(n) = &
         [20, 20, 20, 18, 24, 26, 26, 26, 26, 26]
    PetscInt :: i, ierr, sr
    character(2) :: istr
    PetscMPIInt :: rank
    PetscReal :: Ts, density, props(2)
    PetscErrorCode :: err

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       do i = 1, n
          write(istr, '(i2)') i
          call IAPWS%saturation%temperature(pressure(i), Ts, err)
          call density_case(i, PETSC_TRUE)
          call density_case(i, PETSC_FALSE)
       end do

    end if

  contains

    subroutine density_case(i, liquid)

      PetscInt, intent(in) :: i
      PetscBool, intent(in) :: liquid
      ! Locals:
      character(6) :: phase_str
      PetscInt :: expected_sr

      if (liquid) then
         phase_str = 'liquid'
         expected_sr = liquid_index(i)
      else
         phase_str = 'vapour'
         expected_sr = vapour_index(i)
      end if

      select type (region3 => IAPWS%region(3)%ptr)
      type is (IAPWS_region3_type)

         sr = region3%saturation_subregion_index(pressure(i), liquid)
         call test%assert(expected_sr, sr, ' ' // phase_str // ' sr ' // istr)
         call region3%saturation_density([pressure(i), Ts], &
              liquid, density, err, polish = PETSC_TRUE)
         call test%assert(0, err, ' ' // phase_str // ' err ' // istr)
         if (err == 0) then
            call region3%properties([density, Ts], props, err)
            associate(P2 => props(1))
              call test%assert(pressure(i), P2, ' ' // phase_str // &
                   ' pressure ' // istr)
            end associate
         end if
      end select

    end subroutine density_case

  end subroutine test_IAPWS_region3_saturation_density

!------------------------------------------------------------------------

  subroutine test_region3_widom(test)
    ! Region 3 Widom line and delta tests

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call widom_case(IAPWS%critical%pressure, &
            [373.9308294690435_dp, 373.96117053095657_dp], 0, 'case 1')
       call widom_case(40.e6_dp, &
            [423.18886692525234_dp, 443.541905594127_dp], 0, 'case 2')
       call widom_case(75.e6_dp, &
            [466.14304658909657_dp, 526.1534456640626_dp], 0, 'case 3')
       call widom_case(100.e6_dp, &
            [480.712211436166_dp, 569.049296515483_dp], 0, 'case 4')

    end if

  contains

    PetscReal function widom_p(t)
      ! Widom line from Banuti et al. (2017)
      PetscReal, intent(in) :: t

       select type (region3 => IAPWS%region(3)%ptr)
       type is (IAPWS_region3_type)
          widom_p = IAPWS%critical%pressure * exp(region3%widom_slope * &
               ((t + tc_k) / IAPWS%critical%temperature_k - 1._dp))
       end select

     end function widom_p

     subroutine widom_case(p, expected_delta, expected_err, name)

       PetscReal, intent(in) :: p, expected_delta(2)
       PetscErrorCode, intent(in) :: expected_err
       character(*), intent(in) :: name
       ! Locals:
       PetscReal :: t, delta(2)
       PetscErrorCode :: err

       select type (region3 => IAPWS%region(3)%ptr)
       type is (IAPWS_region3_type)
          call region3%widom(p, t, err)
          call test%assert(expected_err, err, name // ' widom error')
          if (err == 0) then
             call test%assert(p, widom_p(t), name // ' widom temperature')
             call region3%widom_delta(p, delta, err)
             call test%assert(expected_err, err, name // ' widom delta error')
             if (err == 0) then
                call test%assert(expected_delta, delta, name // ' widom delta')
             end if
          end if
       end select

     end subroutine widom_case

  end subroutine test_region3_widom

!------------------------------------------------------------------------

  subroutine test_region3_pi_liquidlike(test)
    ! Region 3 pi_liquidlike tests

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call pi_liquidlike_case(IAPWS%critical%pressure, &
            IAPWS%critical%temperature, 0._dp, 0.5_dp, 3, 0, 'case 1')
       call pi_liquidlike_case(60.e6_dp, IAPWS%critical%temperature, &
            0._dp, 1.0_dp, 1, 0, 'case 2')
       call pi_liquidlike_case(IAPWS%critical%pressure, 700._dp, 0._dp, &
            0.0_dp, 2, 0, 'case 3')
       call pi_liquidlike_case(64.e6_dp, 475._dp, 0._dp, &
            0.6646566520355781_dp, 3, 0, 'case 4')
       call pi_liquidlike_case(64.e6_dp, 465._dp, 0._dp, &
            0.9161795159368163_dp, 3, 0, 'case 5')
       call pi_liquidlike_case(64.e6_dp, 495._dp, 0._dp, &
            0.095494398752873355_dp, 3, 0, 'case 6')
       call pi_liquidlike_case(40.e6_dp, 377._dp, 0._dp, &
            1.0_dp, 1, 0, 'case 7')
       call pi_liquidlike_case(30.e6_dp, 500._dp, 0._dp, &
            0.0_dp, 2, 0, 'case 8')
       call pi_liquidlike_case(40.45712306840381e6_dp, 360._dp, 650._dp, &
            1.0_dp, 0, 0, 'case 9')
       call pi_liquidlike_case(20.16407678415291_dp, 360._dp, 550._dp, &
            1.0_dp, 0, 0, 'case 10')
       call pi_liquidlike_case(21.51502247903555_dp, 380._dp, 150._dp, &
            0.0_dp, 0, 0, 'case 11')
       call pi_liquidlike_case(18.843641755076366_dp, 365._dp, 130._dp, &
            0.0_dp, 0, 0, 'case 12')
       call pi_liquidlike_case(21.976407368553966e6_dp, 373.6172901183901_dp, &
            370.0243306702589_dp, 1.0_dp, 0, 0, 'case 13')

    end if

  contains

    subroutine pi_liquidlike_case(pressure, temperature, density, &
         expected_pi, expected_phases, expected_err, name)

      PetscReal, intent(in) :: pressure, temperature, density
      PetscReal, intent(in) :: expected_pi
      PetscInt, intent(in) :: expected_phases
      PetscErrorCode, intent(in) :: expected_err
      character(*), intent(in) :: name
      ! Locals:
      PetscReal :: pi_liq
      PetscInt :: phases
      PetscErrorCode :: err

      select type (region3 => IAPWS%region(3)%ptr)
      type is (IAPWS_region3_type)
         call region3%pi_liquidlike(pressure, temperature, density, &
              pi_liq, phases, err)
         call test%assert(expected_err, err, name // ' pi_liq error')
         if (err == 0) then
            call test%assert(expected_pi, pi_liq, name // ' pi_liq')
            call test%assert(expected_phases, phases, name // ' phases')
         end if
      end select

    end subroutine pi_liquidlike_case

  end subroutine test_region3_pi_liquidlike

!------------------------------------------------------------------------

  subroutine test_IAPWS_region1_pressure(test)
    ! Region 1 pressure tests

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call pressure_case([990._dp, 100._dp], 60.e6_dp, 0, 'case 1')
       call pressure_case([800._dp, 300._dp], 60.e6_dp, 0, 'case 2')
       call pressure_case([700._dp, 320._dp], 20.e6_dp, 0, 'case 3')
       call pressure_case([700._dp, 320._dp], 80.e6_dp, 0, 'case 4')
       call pressure_case([600._dp, 340._dp], 20.e6_dp, 0, 'case 5')

    end if

  contains

    subroutine pressure_case(param, P0, expected_err, name)

      PetscReal, intent(in) :: param(2), P0
      PetscErrorCode, intent(in) :: expected_err
      character(*), intent(in) :: name
      ! Locals:
      PetscReal :: pressure, props(2)
      PetscErrorCode :: err

      select type (region1 => IAPWS%region(1)%ptr)
      type is (IAPWS_region1_type)
         pressure = P0
         call region1%pressure(param, pressure, err)
         call test%assert(expected_err, err, name // ' error')
         if (err == 0) then
            associate(t => param(2), expected_density => param(1), &
                 d => props(1))
              call region1%properties([pressure, t], props, err)
              call test%assert(expected_err, err, name // ' error')
              if (err == 0) then
                 call test%assert(expected_density, d, name // ' density')
              end if
            end associate
         end if
      end select

    end subroutine pressure_case

  end subroutine test_IAPWS_region1_pressure

!------------------------------------------------------------------------

  subroutine test_IAPWS_region2_pressure(test)
    ! Region 2 pressure tests

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call pressure_case([0.01_dp, 50._dp], 1.e3_dp, 0, 'case 1')
       call pressure_case([10._dp, 300._dp], 10.e5_dp, 0, 'case 2')
       call pressure_case([100._dp, 400._dp], 25.e6_dp, 0, 'case 3')
       call pressure_case([250._dp, 500._dp], 70.e6_dp, 0, 'case 4')
       call pressure_case([200._dp, 700._dp], 50.e6_dp, 0, 'case 5')
       call pressure_case([10._dp, 750._dp], 1.e6_dp, 0, 'case 6')

    end if

  contains

    subroutine pressure_case(param, P0, expected_err, name)

      PetscReal, intent(in) :: param(2), P0
      PetscErrorCode, intent(in) :: expected_err
      character(*), intent(in) :: name
      ! Locals:
      PetscReal :: pressure, props(2)
      PetscErrorCode :: err

      select type (region2 => IAPWS%region(2)%ptr)
      type is (IAPWS_region2_type)
         pressure = P0
         call region2%pressure(param, pressure, err)
         call test%assert(expected_err, err, name // ' error')
         if (err == 0) then
            associate(t => param(2), expected_density => param(1), &
                 d => props(1))
              call region2%properties([pressure, t], props, err)
              call test%assert(expected_err, err, name // ' error')
              if (err == 0) then
                 call test%assert(expected_density, d, name // ' density')
              end if
            end associate
         end if
      end select

    end subroutine pressure_case

  end subroutine test_IAPWS_region2_pressure

!------------------------------------------------------------------------

end module IAPWS_test
