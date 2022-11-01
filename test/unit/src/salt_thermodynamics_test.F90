module salt_thermodynamics_test

  ! Tests for salt thermodynamics module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use salt_thermodynamics_module
  use IFC67_module
  use zofu

  implicit none
  private

  public :: setup, teardown
  public :: test_halite_solubility, test_brine_saturation_pressure
  public :: test_brine_viscosity

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

  subroutine test_halite_solubility(test)
    ! Halite solubility

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then
       call solubility_case("-1 deg C", -1._dp, 0._dp, 1)
       call solubility_case("20 deg C", 20._dp, 0.264044_dp, 0)
       call solubility_case("100 deg C", 100._dp, 0.27998_dp, 0)
       call solubility_case("200 deg C", 200._dp, 0.31898_dp, 0)
       call solubility_case("300 deg C", 300._dp, 0.37918_dp, 0)
       call solubility_case("350 deg C", 350._dp, 0.41723_dp, 0)
       call solubility_case("400 deg C", 400._dp, 0._dp, 1)
    end if

  contains

    subroutine solubility_case(name, temperature, &
         expected_val, expected_err)

      character(*), intent(in) :: name
      PetscReal, intent(in) :: temperature
      PetscReal, intent(in) :: expected_val
      PetscErrorCode, intent(in) :: expected_err
      ! Locals:
      PetscReal :: s
      PetscErrorCode :: err

      call halite_solubility(temperature, s, err)
      if (expected_err == 0) then
         call test%assert(expected_val, s, trim(name) // " value")
      end if
      call test%assert(expected_err, err, trim(name) // " error")

    end subroutine solubility_case

  end subroutine test_halite_solubility

!------------------------------------------------------------------------

  subroutine test_brine_saturation_pressure(test)
    ! Brine saturation pressure

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(IFC67_type) :: thermo
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       call thermo%init()

       call zero_salt_case("20", 20._dp)
       call zero_salt_case("100", 100._dp)
       call zero_salt_case("200", 200._dp)
       call zero_salt_case("300", 300._dp)
       call zero_salt_case("350", 350._dp)

       call sat_case("20, 0",   20._dp, 0._dp, 2.33656155e3_dp, 0)
       call sat_case("100, 0", 100._dp, 0._dp, 1.01325262e5_dp, 0)
       call sat_case("200, 0", 200._dp, 0._dp, 1.55488024e6_dp, 0)
       call sat_case("300, 0", 300._dp, 0._dp,  8.59269200e6_dp, 0)
       call sat_case("350, 0", 350._dp, 0._dp, 1.65351241e7_dp, 0)

       call sat_case("20, 0.1",    20._dp, 0.1_dp, 2.18333836e3_dp, 0)
       call sat_case("100, 0.1",  100._dp, 0.1_dp, 9.46576988e4_dp, 0)
       call sat_case("200, 0.1",  200._dp, 0.1_dp, 1.45230253e6_dp, 0)
       call sat_case("300, 0.1",  300._dp, 0.1_dp, 8.00790120e6_dp, 0)
       call sat_case("350, 0.1",  350._dp, 0.1_dp, 1.53587484e7_dp, 0)

       call sat_case("20, 0.2",    20._dp, 0.2_dp, 1.96358923e3_dp, 0)
       call sat_case("100, 0.2",  100._dp, 0.2_dp, 8.55903276e4_dp, 0)
       call sat_case("200, 0.2",  200._dp, 0.2_dp, 1.31845596e6_dp, 0)
       call sat_case("300, 0.2",  300._dp, 0.2_dp, 7.26818306e6_dp, 0)
       call sat_case("350, 0.2",  350._dp, 0.2_dp, 1.38995352e7_dp, 0)

       call sat_case("20, 0.3",    20._dp, 0.3_dp, 1.57331654e3_dp, 0)
       call sat_case("100, 0.3",  100._dp, 0.3_dp, 7.21180799e4_dp, 0)
       call sat_case("200, 0.3",  200._dp, 0.3_dp, 1.15342823e6_dp, 0)
       call sat_case("300, 0.3",  300._dp, 0.3_dp, 6.49165260e6_dp, 0)
       call sat_case("350, 0.3",  350._dp, 0.3_dp, 1.24828644e7_dp, 0)

       call thermo%destroy()

    end if

  contains

    subroutine sat_case(name, temperature, salt_mass_fraction, &
         expected_val, expected_err)

      character(*), intent(in) :: name
      PetscReal, intent(in) :: temperature
      PetscReal, intent(in) :: salt_mass_fraction
      PetscReal, intent(in) :: expected_val
      PetscErrorCode, intent(in) :: expected_err
      ! Locals:
      PetscReal :: Ps, Ts
      PetscErrorCode :: err

      call brine_saturation_pressure(temperature, salt_mass_fraction, &
           thermo, Ps, err)
      call test%assert(expected_err, err, trim(name) // " error")
      if (expected_err == 0) then
         call test%assert(expected_val, Ps, trim(name) // " value")
         if (err == 0) then
            call brine_saturation_temperature(Ps, salt_mass_fraction, &
                 thermo, Ts, err)
            call test%assert(temperature, Ts, trim(name) // " inverse value")
         end if
      end if

    end subroutine sat_case

    subroutine zero_salt_case(name, temperature)

      character(*), intent(in) :: name
      PetscReal, intent(in) :: temperature
      ! Locals:
      PetscReal :: Ps, Ps0
      PetscReal, parameter :: salt_mass_fraction = 0._dp
      PetscErrorCode :: err

      call thermo%saturation%pressure(temperature, Ps0, err)
      call brine_saturation_pressure(temperature, salt_mass_fraction, &
           thermo, Ps, err)
      call test%assert(Ps0, Ps, trim(name) // " zero salt value")

    end subroutine zero_salt_case

  end subroutine test_brine_saturation_pressure

!------------------------------------------------------------------------

  subroutine test_brine_viscosity(test)
    ! Brine viscosity

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(IFC67_type) :: thermo
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       call thermo%init()

       call visc_case("20, 0",   20._dp, 0._dp, 1.00174876e-03_dp, 0)
       call visc_case("100, 0", 100._dp, 0._dp, 2.78979538e-04_dp, 0)
       call visc_case("200, 0", 200._dp, 0._dp, 1.33827807e-04_dp, 0)
       call visc_case("300, 0", 300._dp, 0._dp, 9.01208869e-05_dp, 0)
       call visc_case("350, 0", 350._dp, 0._dp, 7.86360089e-05_dp, 0)

       call visc_case("20, 0.1",    20._dp, 0.1_dp, 1.21146267e-03_dp, 0)
       call visc_case("100, 0.1",  100._dp, 0.1_dp, 3.47711409e-04_dp, 0)
       call visc_case("200, 0.1",  200._dp, 0.1_dp, 1.72991869e-04_dp, 0)
       call visc_case("300, 0.1",  300._dp, 0.1_dp, 1.20664792e-04_dp, 0)
       call visc_case("350, 0.1",  350._dp, 0.1_dp, 1.07106934e-04_dp, 0)

       call visc_case("20, 0.2",    20._dp, 0.2_dp, 1.59698593e-03_dp, 0)
       call visc_case("100, 0.2",  100._dp, 0.2_dp, 4.58083838e-04_dp, 0)
       call visc_case("200, 0.2",  200._dp, 0.2_dp, 2.27741188e-04_dp, 0)
       call visc_case("300, 0.2",  300._dp, 0.2_dp, 1.58747744e-04_dp, 0)
       call visc_case("350, 0.2",  350._dp, 0.2_dp, 1.40866404e-04_dp, 0)

       call visc_case("20, 0.3",    20._dp, 0.3_dp, 2.32129913e-03_dp, 0)
       call visc_case("100, 0.3",  100._dp, 0.3_dp, 6.60419870e-04_dp, 0)
       call visc_case("200, 0.3",  200._dp, 0.3_dp, 3.25174645e-04_dp, 0)
       call visc_case("300, 0.3",  300._dp, 0.3_dp, 2.24610781e-04_dp, 0)
       call visc_case("350, 0.3",  350._dp, 0.3_dp, 1.98445212e-04_dp, 0)

       call thermo%destroy()

    end if

  contains

    subroutine visc_case(name, temperature, salt_mass_fraction, &
         expected_val, expected_err)

      character(*), intent(in) :: name
      PetscReal, intent(in) :: temperature
      PetscReal, intent(in) :: salt_mass_fraction
      PetscReal, intent(in) :: expected_val
      PetscErrorCode, intent(in) :: expected_err
      ! Locals:
      PetscReal :: Ps, d, visc
      PetscErrorCode :: err

      d = 0._dp
      call thermo%saturation%pressure(temperature, Ps, err)
      call brine_viscosity(Ps, temperature, d, salt_mass_fraction, &
           thermo, visc, err)
      call test%assert(expected_err, err, trim(name) // " error")
      if (expected_err == 0) then
         call test%assert(expected_val, visc, trim(name) // " value")
      end if

    end subroutine visc_case

  end subroutine test_brine_viscosity

!------------------------------------------------------------------------

  end module salt_thermodynamics_test
