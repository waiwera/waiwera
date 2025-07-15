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
  public :: test_halite_solubility_two_phase
  public :: test_halite_properties
  public :: test_brine_saturation_pressure
  public :: test_brine_saturation_temperature
  public :: test_brine_viscosity
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
       call solubility_case("-1 deg C",   -1._dp, 0._dp, 1)
       call solubility_case("20 deg C",   20._dp, 0.26420860_dp, 0)
       call solubility_case("100 deg C", 100._dp, 0.28062682_dp, 0)
       call solubility_case("200 deg C", 200._dp, 0.31730904_dp, 0)
       call solubility_case("300 deg C", 300._dp, 0.37747855_dp, 0)
       call solubility_case("400 deg C", 400._dp, 0.47145444_dp, 0)
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

  subroutine test_halite_solubility_two_phase(test)
    ! Two-phase halite solubility

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    type(IFC67_type) :: thermo
    PetscInt :: i
    PetscInt, parameter :: n = 20
    character(16) :: case
    PetscReal :: p
    PetscReal, parameter :: p0 = 1.e5_dp, p1 = 13.5e6_dp
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then
       call thermo%init()
       call solubility_case("-1 bar",  -1.e5_dp, 1)
       call solubility_case("23 MPa",  23.e6_dp, 1)
       do i = 1, n
          p = p0 + (i - 1._dp) / (n - 1._dp) * (p1 - p0)
          write(case, '(a4, f6.2, a4)') 'P = ', p / 1.e6_dp, ' MPa'
          call solubility_case(case, p, 0)
       end do
       call thermo%destroy()
    end if

  contains

    subroutine solubility_case(name, pressure, expected_err)

      character(*), intent(in) :: name
      PetscReal, intent(in) :: pressure
      PetscErrorCode, intent(in) :: expected_err
      ! Locals:
      PetscReal :: s, temperature, expected_solubility
      PetscErrorCode :: err, err2

      err2 = 0
      call halite_solubility_two_phase(pressure, thermo, s, err)
      if (expected_err == 0) then
         call brine_saturation_temperature(pressure, s, thermo, &
              temperature, err)
         if (err == 0) then
            call halite_solubility(temperature, expected_solubility, err2)
            if (err2 == 0) then
               call test%assert(expected_solubility, s, trim(name) // " value")
            end if
         end if
      end if
      call test%assert(expected_err, err, trim(name) // " error")
      call test%assert(0, err2, trim(name) // " error2")

    end subroutine solubility_case

  end subroutine test_halite_solubility_two_phase

!------------------------------------------------------------------------

  subroutine test_halite_properties(test)
    ! Halite properties

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    PetscInt :: ierr
    type(IFC67_type) :: thermo

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then
       call thermo%init()
       call properties_case("20 deg C",   20._dp, [2.16544208e3_dp, -5.44002467e5_dp], 0)
       call properties_case("100 deg C", 100._dp, [2.14485199e3_dp, -4.73170767e5_dp], 0)
       call properties_case("200 deg C", 200._dp, [2.11751462e3_dp, -3.83011549e5_dp], 0)
       call properties_case("300 deg C", 300._dp, [2.08882457e3_dp, -2.90739752e5_dp], 0)
       call properties_case("350 deg C", 350._dp, [2.07423883e3_dp, -2.43776955e5_dp], 0)
       call thermo%destroy()
    end if

  contains

    subroutine properties_case(name, temperature, &
         expected_props, expected_err)

      character(*), intent(in) :: name
      PetscReal, intent(in) :: temperature
      PetscReal, intent(in) :: expected_props(2)
      PetscErrorCode, intent(in) :: expected_err
      ! Locals:
      PetscReal :: Ps, props(2)
      PetscErrorCode :: err

      call thermo%saturation%pressure(temperature, Ps, err)
      call halite_properties(Ps, temperature, props, err)
      if (expected_err == 0) then
         call test%assert(expected_props, props, trim(name) // " value")
      end if
      call test%assert(expected_err, err, trim(name) // " error")

    end subroutine properties_case

  end subroutine test_halite_properties

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
       call sat_case("300, 0", 300._dp, 0._dp, 8.59269200e6_dp, 0)
       call sat_case("350, 0", 350._dp, 0._dp, 1.65351241e7_dp, 0)

       call sat_case("20, 0.1",    20._dp, 0.1_dp, 2.18332495e3_dp, 0)
       call sat_case("100, 0.1",  100._dp, 0.1_dp, 9.46571127e4_dp, 0)
       call sat_case("200, 0.1",  200._dp, 0.1_dp, 1.45229349e6_dp, 0)
       call sat_case("300, 0.1",  300._dp, 0.1_dp, 8.00784960e6_dp, 0)
       call sat_case("350, 0.1",  350._dp, 0.1_dp, 1.53586450e7_dp, 0)

       call sat_case("20, 0.2",    20._dp, 0.2_dp, 1.96355069e3_dp, 0)
       call sat_case("100, 0.2",  100._dp, 0.2_dp, 8.55888487e4_dp, 0)
       call sat_case("200, 0.2",  200._dp, 0.2_dp, 1.31843548e6_dp, 0)
       call sat_case("300, 0.2",  300._dp, 0.2_dp, 7.26807523e6_dp, 0)
       call sat_case("350, 0.2",  350._dp, 0.2_dp, 1.38993274e7_dp, 0)

       call sat_case("20, 0.3",    20._dp, 0.3_dp, 1.57322067e3_dp, 0)
       call sat_case("100, 0.3",  100._dp, 0.3_dp, 7.21150101e4_dp, 0)
       call sat_case("200, 0.3",  200._dp, 0.3_dp, 1.15339496e6_dp, 0)
       call sat_case("300, 0.3",  300._dp, 0.3_dp, 6.49151820e6_dp, 0)
       call sat_case("350, 0.3",  350._dp, 0.3_dp, 1.24826404e7_dp, 0)

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

  subroutine test_brine_saturation_temperature(test)
    ! Brine saturation temperature

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(IFC67_type) :: thermo
    PetscMPIInt :: rank
    PetscInt, parameter :: m = 50, n = 10
    character(24) :: case
    PetscReal :: p, x
    PetscInt :: i, j
    PetscReal, parameter :: p0 = 1.e5_dp, p1 = 16.5e6_dp
    PetscReal, parameter :: x0 = 0._dp, x1 = 0.3_dp
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       call thermo%init()

       do i = 1, m
          p = p0 + (i - 1._dp) / (m - 1._dp) * (p1 - p0)
          do j = 1, n
             x = x0 + (j - 1._dp) / (n - 1._dp) * (x1 - x0)
             write(case, '(f6.2, a4, f6.2)') p / 1.e6_dp, ' MPa, ', x
             call sat_case(case, p, x, 0)
          end do
       end do

       call thermo%destroy()

    end if

  contains

    subroutine sat_case(name, pressure, salt_mass_fraction, &
         expected_err)

      character(*), intent(in) :: name
      PetscReal, intent(in) :: pressure
      PetscReal, intent(in) :: salt_mass_fraction
      PetscErrorCode, intent(in) :: expected_err
      ! Locals:
      PetscReal :: Ps, Ts
      PetscErrorCode :: err, err2

      err2 = 0

      call brine_saturation_temperature(pressure, salt_mass_fraction, &
           thermo, Ts, err)
      call test%assert(expected_err, err, trim(name) // " error")
      if (expected_err == 0) then
         call brine_saturation_pressure(Ts, salt_mass_fraction, thermo, Ps, err2)
         if (err2 == 0) then
            call test%assert(Ps, pressure, trim(name) // " value")
         end if
         call test%assert(0, err2, trim(name) // " error2")
      end if

    end subroutine sat_case

  end subroutine test_brine_saturation_temperature

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

       call visc_case("20, 0.1",   20._dp, 0.1_dp, 1.21148413e-3_dp, 0)
       call visc_case("100, 0.1", 100._dp, 0.1_dp, 3.47717807e-4_dp, 0)
       call visc_case("200, 0.1", 200._dp, 0.1_dp, 1.72995191e-4_dp, 0)
       call visc_case("300, 0.1", 300._dp, 0.1_dp, 1.20667200e-4_dp, 0)
       call visc_case("350, 0.1", 350._dp, 0.1_dp, 1.07109110e-4_dp, 0)

       call visc_case("20, 0.2",   20._dp, 0.2_dp, 1.59705684e-3_dp, 0)
       call visc_case("100, 0.2", 100._dp, 0.2_dp, 4.58103766e-4_dp, 0)
       call visc_case("200, 0.2", 200._dp, 0.2_dp, 2.27750855e-4_dp, 0)
       call visc_case("300, 0.2", 300._dp, 0.2_dp, 1.58754327e-4_dp, 0)
       call visc_case("350, 0.2", 350._dp, 0.2_dp, 1.40872179e-4_dp, 0)

       call visc_case("20, 0.3",   20._dp, 0.3_dp, 2.32147586e-3_dp, 0)
       call visc_case("100, 0.3", 100._dp, 0.3_dp, 6.60469123e-4_dp, 0)
       call visc_case("200, 0.3", 200._dp, 0.3_dp, 3.25198294e-4_dp, 0)
       call visc_case("300, 0.3", 300._dp, 0.3_dp, 2.24626721e-4_dp, 0)
       call visc_case("350, 0.3", 350._dp, 0.3_dp, 1.98459127e-4_dp, 0)

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
      PetscReal :: Ps, visc
      PetscErrorCode :: err

      call thermo%saturation%pressure(temperature, Ps, err)
      call brine_viscosity(temperature, Ps, salt_mass_fraction, &
           thermo, visc, err)
      call test%assert(expected_err, err, trim(name) // " error")
      if (expected_err == 0) then
         call test%assert(expected_val, visc, trim(name) // " value")
      end if

    end subroutine visc_case

  end subroutine test_brine_viscosity

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
    PetscReal, parameter :: t(4) = [10._dp, 100._dp, 200._dp, 300._dp]
    PetscReal, parameter :: xs(4) = [0._dp, 0.1_dp, 0.2_dp, 0.25_dp]
    PetscReal, parameter :: expected_density(size(xs), size(t), size(p)) = reshape([ &
         0.99979351E+03_dp, 0.10739598E+04_dp, 0.11529153E+04_dp, 0.11933788E+04_dp, &
         0.95812176E+03_dp, 0.10272547E+04_dp, 0.11027768E+04_dp, 0.11433353E+04_dp, &
         0.85626589E+03_dp, 0.93720194E+03_dp, 0.10191556E+04_dp, 0.10620991E+04_dp, &
         0.70015006E+03_dp, 0.80298857E+03_dp, 0.90027464E+03_dp, 0.94946347E+03_dp, &
         0.10002189E+04_dp, 0.10743439E+04_dp, 0.11531331E+04_dp, 0.11935841E+04_dp, &
         0.95855591E+03_dp, 0.10275680E+04_dp, 0.11029480E+04_dp, 0.11434872E+04_dp, &
         0.86275186E+03_dp, 0.94343137E+03_dp, 0.10250794E+04_dp, 0.10680578E+04_dp, &
         0.68154614E+03_dp, 0.79377210E+03_dp, 0.89628366E+03_dp, 0.94706615E+03_dp, &
         0.10044352E+04_dp, 0.10781109E+04_dp, 0.11561697E+04_dp, 0.11969145E+04_dp, &
         0.96282473E+03_dp, 0.10314023E+04_dp, 0.11060135E+04_dp, 0.11461556E+04_dp, &
         0.87107596E+03_dp, 0.94917659E+03_dp, 0.10293831E+04_dp, 0.10716744E+04_dp, &
         0.71538348E+03_dp, 0.82526593E+03_dp, 0.92194681E+03_dp, 0.97012926E+03_dp], &
         [size(xs), size(t), size(p)])

    PetscReal, parameter :: expected_enthalpy(size(xs), size(t), size(p)) = reshape([ &
         0.42090543E+05_dp, 0.47403251E+05_dp, 0.54240665E+05_dp, 0.58259798E+05_dp, &
         0.41906369E+06_dp, 0.39316472E+06_dp, 0.37156721E+06_dp, 0.36193592E+06_dp, &
         0.85178481E+06_dp, 0.78747999E+06_dp, 0.73169015E+06_dp, 0.70586576E+06_dp, &
         0.13586546E+07_dp, 0.12268205E+07_dp, 0.11217912E+07_dp, 0.10746141E+07_dp, &
         0.42969298E+05_dp, 0.48328407E+05_dp, 0.55213434E+05_dp, 0.59256807E+05_dp, &
         0.41973864E+06_dp, 0.39389202E+06_dp, 0.37234207E+06_dp, 0.36273464E+06_dp, &
         0.85214645E+06_dp, 0.78794005E+06_dp, 0.73222111E+06_dp, 0.70642830E+06_dp, &
         0.13568585E+07_dp, 0.12262123E+07_dp, 0.11216359E+07_dp, 0.10745992E+07_dp, &
         0.51707550E+05_dp, 0.57543140E+05_dp, 0.64919394E+05_dp, 0.69213759E+05_dp, &
         0.42650354E+06_dp, 0.40118549E+06_dp, 0.38012053E+06_dp, 0.37075751E+06_dp, &
         0.85592476E+06_dp, 0.79266142E+06_dp, 0.73763322E+06_dp, 0.71215154E+06_dp, &
         0.13433636E+07_dp, 0.12215744E+07_dp, 0.11207975E+07_dp, 0.10749918E+07_dp], &
         [size(xs), size(t), size(p)])

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       call thermo%init()

       do i = 1, size(p)
          do j = 1, size(t)
             do k = 1, size(xs)
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
      PetscReal :: props(2), h
      PetscErrorCode :: err

      call brine_properties(p, t, salt_mass_fraction, thermo, props, err)
      call test%assert(expected_err, err, trim(name) // " error")
      if (expected_err == 0) then
         associate(d => props(1), u => props(2))
           h = u + p / d
           call test%assert(expected_props, [d, h], trim(name) // " value")
         end associate
      end if

    end subroutine props_case

  end subroutine test_brine_properties

!------------------------------------------------------------------------

end module salt_thermodynamics_test
