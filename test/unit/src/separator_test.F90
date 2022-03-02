module separator_test

  ! Tests for separator module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use separator_module

  implicit none
  private

  character(len = 512) :: data_path

  public :: setup, teardown
  public :: test_separator, test_separator_2stage

contains

!------------------------------------------------------------------------

  subroutine setup()

    use profiling_module, only: init_profiling

    ! Locals:
    PetscErrorCode :: ierr
    PetscInt :: ios

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call init_profiling()

    call get_environment_variable('WAIWERA_TEST_DATA_PATH', &
         data_path, status = ios)
    if (ios /= 0) data_path = ''

  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown()

    PetscErrorCode :: ierr

    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

  subroutine test_separator(test)
    ! Separator

    use IAPWS_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    type(IAPWS_type) :: thermo
    PetscReal, pointer, contiguous :: data(:)
    type(separator_type) :: sep
    PetscReal :: water_rate, steam_rate
    PetscReal :: water_enthalpy, steam_enthalpy
    PetscReal, parameter :: separator_pressure(1) = [10.e5_dp]
    PetscInt, parameter :: offset = 1
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       allocate(data(8))

       call thermo%init()
       call sep%assign(data, offset)
       call sep%init(separator_pressure, thermo)

       call sep%separate(-10._dp, 500.e3_dp, water_rate, water_enthalpy, &
            steam_rate, steam_enthalpy)
       call separator_test_case(test, sep, water_rate, water_enthalpy, &
            steam_rate, steam_enthalpy, "water", 0._dp, -10._dp, 500.e3_dp, 0._dp, 0._dp)

       call sep%separate(-10._dp, 3000.e3_dp, water_rate, water_enthalpy, &
            steam_rate, steam_enthalpy)
       call separator_test_case(test, sep, water_rate, water_enthalpy, &
            steam_rate, steam_enthalpy, "steam", 1._dp, 0._dp, 0._dp, -10._dp, 3000.e3_dp)

       call sep%separate(-10._dp, 1200.e3_dp, water_rate, water_enthalpy, &
            steam_rate, steam_enthalpy)
       call separator_test_case(test, sep, water_rate, water_enthalpy, &
            steam_rate, steam_enthalpy, "two-phase", 0.21709153586628488_dp, &
            -7.829084641337152_dp, 762682.8443354106_dp, &
            -2.1709153586628487_dp, 2777119.5376846623_dp)

       deallocate(data)
       call thermo%destroy()
       call sep%destroy()

    end if

  end subroutine test_separator

!------------------------------------------------------------------------

  subroutine test_separator_2stage(test)
    ! 2-stage separator

    use IAPWS_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscMPIInt :: rank
    type(IAPWS_type) :: thermo
    PetscReal, pointer, contiguous :: data(:)
    type(separator_type) :: sep
    PetscReal :: water_rate, steam_rate
    PetscReal :: water_enthalpy, steam_enthalpy
    PetscReal, parameter :: separator_pressure(2) = [1.45e6_dp, 0.55e6_dp]
    PetscInt, parameter :: offset = 1
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       allocate(data(8))

       call thermo%init()
       call sep%assign(data, offset)
       call sep%init(separator_pressure, thermo)

       call sep%separate(-10._dp, 500.e3_dp, water_rate, water_enthalpy, &
            steam_rate, steam_enthalpy)
       call separator_test_case(test, sep, water_rate, water_enthalpy, &
            steam_rate, steam_enthalpy, "water", 0._dp, -10._dp, 500.e3_dp, 0._dp, 0._dp)

       call sep%separate(-10._dp, 3000.e3_dp, water_rate, water_enthalpy, &
            steam_rate, steam_enthalpy)
       call separator_test_case(test, sep, water_rate, water_enthalpy, &
            steam_rate, steam_enthalpy, "steam", 1._dp, 0._dp, 0._dp, -10._dp, 3000.e3_dp)

       call sep%separate(-10._dp, 1200.e3_dp, water_rate, water_enthalpy, &
            steam_rate, steam_enthalpy)
       call separator_test_case(test, sep, water_rate, water_enthalpy, &
            steam_rate, steam_enthalpy, "two-phase", 0.256210105124_dp, &
            -7.437898948764703_dp, 655876.6515067405_dp, &
            -2.5621010512352966_dp, 2779615.4799612807_dp)       
       deallocate(data)
       call thermo%destroy()
       call sep%destroy()

    end if

  end subroutine test_separator_2stage

!------------------------------------------------------------------------

  subroutine separator_test_case(test, separator, water_rate, water_enthalpy, &
       steam_rate, steam_enthalpy, name, &
       expected_steam_fraction, expected_water_rate, expected_water_enthalpy, &
       expected_steam_rate, expected_steam_enthalpy)

    class(unit_test_type), intent(in out) :: test
    type(separator_type) :: separator
    PetscReal, intent(in) :: water_rate, water_enthalpy
    PetscReal, intent(in) :: steam_rate, steam_enthalpy
    character(*), intent(in) :: name
    PetscReal, intent(in) :: expected_steam_fraction
    PetscReal, intent(in) :: expected_water_rate, expected_water_enthalpy
    PetscReal, intent(in) :: expected_steam_rate, expected_steam_enthalpy

    call test%assert(expected_steam_fraction, separator%steam_fraction, &
         trim(name) // " steam fraction")
    call test%assert(expected_water_rate, water_rate, trim(name) // " water rate")
    call test%assert(expected_water_enthalpy, water_enthalpy, trim(name) // " water enthalpy")
    call test%assert(expected_steam_rate, steam_rate, trim(name) // " steam rate")
    call test%assert(expected_steam_enthalpy, steam_enthalpy, trim(name) // " steam enthalpy")

  end subroutine separator_test_case

!------------------------------------------------------------------------

end module separator_test
