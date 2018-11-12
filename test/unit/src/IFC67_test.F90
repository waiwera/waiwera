module IFC67_test

  ! Tests for IAPWS thermodynamics module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use IFC67_module
  use thermodynamics_module, only: tc_k
  use zofu

  implicit none
  private

  type(IFC67_type) :: IFC67

  public :: setup, teardown, setup_test
  public :: test_IFC67_region1, test_IFC67_region2
  public :: test_IFC67_saturation, test_IFC67_viscosity
  public :: test_IFC67_phase_composition

contains

!------------------------------------------------------------------------

  subroutine setup()

    use profiling_module, only: init_profiling

    ! Locals:
    PetscErrorCode :: ierr

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call init_profiling()
    call IFC67%init()

  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown()

    PetscErrorCode :: ierr

    call IFC67%destroy()
    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

  subroutine setup_test(test)

    class(unit_test_type), intent(in out) :: test

    test%tolerance = 1.e-7

  end subroutine setup_test

!------------------------------------------------------------------------

  subroutine test_IFC67_region1(test)

    ! IFC-67 region 1 tests

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: n = 3, nerr = 2
    PetscReal :: params(n,2) = reshape([ &
         3.e6_dp, 80.e6_dp, 3.e6_dp, &
         300._dp, 300._dp,  500._dp], [n,2])
    PetscReal, parameter :: rho(n) = [997.95721560998174_dp, 1029.7256888266911_dp, &
         831.84196191567298_dp]
    PetscReal, parameter ::  u(n) = [112247.43313085975_dp, 106310.47344628950_dp, &
         971985.91117384087_dp]
    PetscInt :: i, err
    PetscReal :: param(2), props(2)
    PetscReal :: err_params(nerr,2) = reshape([ &
         20.e6_dp, 101.e6_dp, &
         360._dp,    60._dp], [nerr,2])
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       params(:,2) = params(:,2) - tc_k  ! convert temperatures to Celcius
       do i = 1, n
          param = params(i,:)
          call IFC67%water%properties(param, props, err)
          call test%assert(rho(i), props(1), 'density')
          call test%assert(u(i), props(2), 'energy')
          call test%assert(0, err, 'error')
       end do
       do i = 1, nerr
          param = err_params(i,:)
          call IFC67%water%properties(param, props, err)
          call test%assert(1, err, 'error')
       end do
    end if

  end subroutine test_IFC67_region1

!------------------------------------------------------------------------

  subroutine test_IFC67_region2(test)

    ! IFC-67 region 2 tests

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: n = 3, nerr = 2
    PetscReal :: params(n,2) = reshape([ &
         0.0035e6_dp, 0.0035e6_dp, 30.e6_dp, &
         300._dp, 700._dp,  700._dp], [n,2])
    PetscReal, parameter :: rho(n) = [2.5316826343790743e-2_dp, 1.0834441421293962e-2_dp, &
         183.90041953968711_dp]
    PetscReal, parameter ::  u(n) = [2412405.0932077002_dp, 3012229.4965919587_dp, &
         2474981.3799304822_dp]
    PetscInt :: i, err
    PetscReal :: param(2), props(2)
    PetscReal :: err_params(nerr,2) = reshape([ &
         20.e6_dp, 101.e6_dp, &
         801._dp,    60._dp], [nerr,2])
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       params(:,2) = params(:,2) - tc_k  ! convert temperatures to Celcius
       do i = 1, n
          param = params(i,:)
          call IFC67%steam%properties(param, props, err)
          call test%assert(rho(i), props(1), 'density')
          call test%assert(u(i), props(2), 'energy')
          call test%assert(0, err, 'error')
       end do
       do i = 1, nerr
          param = err_params(i,:)
          call IFC67%steam%properties(param, props, err)
          call test%assert(1, err, 'error')
       end do
    end if

  end subroutine test_IFC67_region2

!------------------------------------------------------------------------

  subroutine test_IFC67_saturation(test)

    ! IFC-67 saturation curve tests

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: n = 3, nerr = 1
    PetscReal, parameter ::  t(n) = [300._dp, 500._dp, 600._dp] - tc_k
    PetscReal, parameter :: p(n) = [0.35323426e4_dp, 0.263961572e7_dp, &
         0.123493902e8_dp]
    PetscReal :: ps, ts, ps1, ts1
    PetscInt :: i, err
    PetscReal :: terr(nerr) = [380._dp], perr(nerr) = [30.e6_dp]
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       do i = 1, n

          call IFC67%saturation%pressure(t(i), ps, err)
          call test%assert(p(i), ps, 'pressure')
          call test%assert(0, err, 'pressure error')
          call IFC67%saturation%temperature(ps, ts, err)
          call test%assert(t(i), ts, 'temperature')
          call test%assert(0, err, 'temperature error')

          ! Test region 1 saturation object:
          call IFC67%region(1)%ptr%saturation%pressure(t(i), ps1, err)
          call test%assert(p(i), ps1, 'region 1 pressure')
          call test%assert(0, err, 'region 1 pressure error')
          call IFC67%region(1)%ptr%saturation%temperature(ps1, ts1, err)
          call test%assert(t(i), ts1, 'region 1 temperature')
          call test%assert(0, err, 'region 1 temperature error')

       end do
       do i = 1, nerr
          call IFC67%saturation%pressure(terr(i), ps, err)
          call test%assert(1, err, 'error')
          call IFC67%saturation%temperature(perr(i), ts, err)
          call test%assert(1, err, 'temperature error')
       end do

    end if

  end subroutine test_IFC67_saturation

!------------------------------------------------------------------------

  subroutine test_IFC67_viscosity(test)

    ! IFC-67 viscosity tests
    
    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: n1 = 2, n2 = 2
    ! region 1 tests:
    PetscReal, parameter :: t1(n1) = [298.15_dp, 373.15_dp] - tc_k
    PetscReal, parameter :: p1(n1) = [1977563.58349_dp, 99834578.2816_dp]
    PetscReal, parameter :: visc1(n1) = [8.903129e-04_dp, 2.988268e-04_dp]
    ! region 2 tests:
    PetscReal, parameter :: t2(n2) = [873.15_dp, 873.15_dp] - tc_k
    PetscReal, parameter :: d2(n2) = [1._dp, 100._dp]
    PetscReal, parameter :: visc2(n2) = [3.249537e-05_dp, 3.667671e-05_dp]
    PetscReal :: v, p = 0.0_dp, d = 0.0_dp
    PetscInt :: i
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       do i = 1, n1
          call IFC67%water%viscosity(t1(i), p1(i), d, v)
          call test%assert(visc1(i), v)
       end do
       do i = 1, n2
          call IFC67%steam%viscosity(t2(i), p, d2(i), v)
          call test%assert(visc2(i), v)
       end do
    end if

  end subroutine test_IFC67_viscosity

!------------------------------------------------------------------------

  subroutine test_IFC67_phase_composition(test)

    ! IFC-67 phase composition tests

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt :: phases, expected_phases
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       phases = IFC67%phase_composition(1, 1.e5_dp, 20._dp)
       expected_phases = b'01'
       call test%assert(expected_phases, phases, &
            "Region 1 liquid")

       phases = IFC67%phase_composition(2, 1.e5_dp, 110._dp)
       expected_phases = b'10'
       call test%assert(expected_phases, phases, &
            "Region 2 steam")

       phases = IFC67%phase_composition(4, 33.466518715101621e5_dp, 240._dp)
       expected_phases = b'11'
       call test%assert(expected_phases, phases, &
            "Two-phase at 240 deg C")

    end if

  end subroutine test_IFC67_phase_composition

!------------------------------------------------------------------------

end module IFC67_test
