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
  public :: test_halite_solubility
  public :: test_halite_density, test_halite_enthalpy
  public :: test_brine_saturation_pressure
  public :: test_brine_viscosity, test_brine_critical_temperature
  public :: test_brine_properties

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

  subroutine test_halite_density(test)
    ! Halite density

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr
    type(IFC67_type) :: thermo

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then
       call thermo%init()
       call density_case("20 deg C",   20._dp, 2.15981043e+03_dp, 0)
       call density_case("100 deg C", 100._dp, 2.13918393e+03_dp, 0)
       call density_case("200 deg C", 200._dp, 2.11379003e+03_dp, 0)
       call density_case("300 deg C", 300._dp, 2.08916417e+03_dp, 0)
       call density_case("350 deg C", 350._dp, 2.07732657e+03_dp, 0)
       call thermo%destroy()
    end if

  contains

    subroutine density_case(name, temperature, &
         expected_val, expected_err)

      character(*), intent(in) :: name
      PetscReal, intent(in) :: temperature
      PetscReal, intent(in) :: expected_val
      PetscErrorCode, intent(in) :: expected_err
      ! Locals:
      PetscReal :: Ps, density
      PetscErrorCode :: err

      call thermo%saturation%pressure(temperature, Ps, err)
      call halite_density(Ps, temperature, density, err)
      if (expected_err == 0) then
         call test%assert(expected_val, density, trim(name) // " value")
      end if
      call test%assert(expected_err, err, trim(name) // " error")

    end subroutine density_case

  end subroutine test_halite_density

!------------------------------------------------------------------------

  subroutine test_halite_enthalpy(test)
    ! Halite enthalpy

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr
    type(IFC67_type) :: thermo

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then
       call thermo%init()
       call enthalpy_case("20 deg C",   20._dp, 1.72924694e+04_dp, 0)
       call enthalpy_case("100 deg C", 100._dp, 8.76135632e+04_dp, 0)
       call enthalpy_case("200 deg C", 200._dp, 1.78027564e+05_dp, 0)
       call enthalpy_case("300 deg C", 300._dp, 2.71233380e+05_dp, 0)
       call enthalpy_case("350 deg C", 350._dp, 3.18883218e+05_dp, 0)
       call thermo%destroy()
    end if

  contains

    subroutine enthalpy_case(name, temperature, &
         expected_val, expected_err)

      character(*), intent(in) :: name
      PetscReal, intent(in) :: temperature
      PetscReal, intent(in) :: expected_val
      PetscErrorCode, intent(in) :: expected_err
      ! Locals:
      PetscReal :: h
      PetscErrorCode :: err

      call halite_enthalpy(temperature, h, err)
      if (expected_err == 0) then
         call test%assert(expected_val, h, trim(name) // " value")
      end if
      call test%assert(expected_err, err, trim(name) // " error")

    end subroutine enthalpy_case

  end subroutine test_halite_enthalpy

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
      call brine_viscosity(temperature, Ps, d, salt_mass_fraction, &
           thermo, visc, err)
      call test%assert(expected_err, err, trim(name) // " error")
      if (expected_err == 0) then
         call test%assert(expected_val, visc, trim(name) // " value")
      end if

    end subroutine visc_case

  end subroutine test_brine_viscosity

!------------------------------------------------------------------------

  subroutine test_brine_critical_temperature(test)
    ! Brine critical temperature

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then
       call critical_temperature_case("0",      0._dp, 374.15_dp, 0)
       call critical_temperature_case("0.05", 0.05_dp, 421.191707134445_dp, 0)
       call critical_temperature_case("0.1",   0.1_dp, 476.973186669071_dp, 0)
       call critical_temperature_case("0.15", 0.15_dp, 541.4997586262932_dp, 0)
       call critical_temperature_case("0.2",   0.2_dp, 609.8012731271594_dp, 0)
       call critical_temperature_case("0.25", 0.25_dp, 673.5698534299003_dp, 0)
       call critical_temperature_case("0.3",   0.3_dp, 728.4174135605035_dp, 0)
    end if
  contains

    subroutine critical_temperature_case(name, salt_mass_fraction, &
         expected_val, expected_err)

      character(*), intent(in) :: name
      PetscReal, intent(in) :: salt_mass_fraction
      PetscReal, intent(in) :: expected_val
      PetscErrorCode, intent(in) :: expected_err
      ! Locals:
      PetscReal :: tc
      PetscErrorCode :: err

      call brine_critical_temperature(salt_mass_fraction, tc, err)
      call test%assert(expected_err, err, trim(name) // " error")
      if (expected_err == 0) then
         call test%assert(expected_val, tc, trim(name) // " value")
      end if

    end subroutine critical_temperature_case

  end subroutine test_brine_critical_temperature

!------------------------------------------------------------------------

  subroutine test_brine_properties(test)
    ! Brine properties

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr, i, j, k
    type(IFC67_type) :: thermo
    character(3) :: name
    PetscReal, parameter ::  p(3) = [1.e5_dp, 10.e5_dp, 100.e5_dp]
    PetscReal, parameter :: t(3) = [100._dp, 200._dp, 300._dp]
    PetscReal, parameter :: xs(4) = [0._dp, 0.1_dp, 0.2_dp, 0.3_dp]
    PetscReal, parameter :: expected_density(4, 3, 3) = reshape([ &
         0.95812163E+03_dp, 0.10282636E+04_dp, 0.11035250E+04_dp, 0.11845521E+04_dp, &
         0.86363535E+03_dp, 0.94310446E+03_dp, 0.10229377E+04_dp, 0.11013288E+04_dp, &
         0.69852800E+03_dp, 0.81231985E+03_dp, 0.92003650E+03_dp, 0.10128337E+04_dp, &
         0.95853309E+03_dp, 0.10286517E+04_dp, 0.11039345E+04_dp, 0.11850704E+04_dp, &
         0.86428163E+03_dp, 0.94364599E+03_dp, 0.10234795E+04_dp, 0.11020369E+04_dp, &
         0.69995391E+03_dp, 0.81320323E+03_dp, 0.92084705E+03_dp, 0.10140055E+04_dp, &
         0.96266716E+03_dp, 0.10325492E+04_dp, 0.11080466E+04_dp, 0.11902792E+04_dp, &
         0.87079804E+03_dp, 0.94909575E+03_dp, 0.10289289E+04_dp, 0.11091680E+04_dp, &
         0.71453985E+03_dp, 0.82214387E+03_dp, 0.92903176E+03_dp, 0.10258747E+04_dp], &
         [4, 3, 3])
    PetscReal, parameter :: expected_enthalpy(4, 3, 3) = reshape([ &
         0.41906369E+06_dp, 0.36102859E+06_dp, 0.30161341E+06_dp, 0.24685456E+06_dp, &
         0.85178076E+06_dp, 0.75114204E+06_dp, 0.66777960E+06_dp, 0.59343873E+06_dp, &
         0.13556369E+07_dp, 0.11964168E+07_dp, 0.11000645E+07_dp, 0.10235516E+07_dp, &
         0.41973864E+06_dp, 0.36161006E+06_dp, 0.30209918E+06_dp, 0.24725213E+06_dp, &
         0.85214586E+06_dp, 0.75146399E+06_dp, 0.66806580E+06_dp, 0.59369305E+06_dp, &
         0.13545154E+07_dp, 0.11954264E+07_dp, 0.10991533E+07_dp, 0.10227032E+07_dp, &
         0.42650354E+06_dp, 0.36743804E+06_dp, 0.30696797E+06_dp, 0.25123689E+06_dp, &
         0.85592476E+06_dp, 0.75479761E+06_dp, 0.67103086E+06_dp, 0.59632955E+06_dp, &
         0.13433636E+07_dp, 0.11856024E+07_dp, 0.10901413E+07_dp, 0.10143385E+07_dp], &
         [4, 3, 3])

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       call thermo%init()

       do i = 1, 3
          do j = 1, 3
             do k = 1, 4
                write(name, '(i1, i1, i1)') i, j, k
                call props_case(name, p(i), t(j), xs(k), &
                     [expected_density(k, j, i), expected_enthalpy(k, j, i)], 0)
             end do
          end do
       end do

       call thermo%destroy()

    end if
  contains

    subroutine props_case(name, p, t, salt_mass_fraction, &
         expected_props, expected_err)

      character(*), intent(in) :: name
      PetscReal, intent(in) :: p, t
      PetscReal, intent(in) :: salt_mass_fraction
      PetscReal, intent(in) :: expected_props(2)
      PetscErrorCode, intent(in) :: expected_err
      ! Locals:
      PetscReal :: props(2)
      PetscErrorCode :: err

      call brine_properties(p, t, salt_mass_fraction, thermo, props, err)
      call test%assert(expected_err, err, trim(name) // " error")
      if (expected_err == 0) then
         call test%assert(expected_props, props, trim(name) // " value")
      end if

    end subroutine props_case

  end subroutine test_brine_properties

!------------------------------------------------------------------------

end module salt_thermodynamics_test
