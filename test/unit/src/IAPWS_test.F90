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
       test_IAPWS_phase_composition

  contains

!------------------------------------------------------------------------

  subroutine setup()

    use profiling_module, only: init_profiling

    ! Locals:
    PetscErrorCode :: ierr

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call init_profiling()
    call IAPWS%init()

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
    PetscInt, parameter :: n = 3, nerr = 1
    PetscReal :: params(n,2) = reshape([ &
         500._dp, 200._dp, 500._dp, &
         650._dp, 650._dp, 750._dp], [n,2])
    PetscReal, parameter :: p(n) = [0.255837018e8_dp, 0.222930643e8_dp, 0.783095639e8_dp]
    PetscReal, parameter ::  u(n) = [0.181226279e7_dp, 0.226365868e7_dp, 0.210206932e7_dp]
    PetscInt :: i, err
    PetscReal :: param(2), props(2)
    PetscReal :: err_params(nerr,2) = reshape([ &
         800._dp, &
         400._dp], [nerr,2])
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

end module IAPWS_test
