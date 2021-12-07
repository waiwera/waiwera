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
  public :: test_separator

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
    PetscReal, parameter :: separator_pressure = 10.e5_dp
    PetscInt, parameter :: offset = 1
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       allocate(data(8))

       call thermo%init()
       call sep%assign(data, offset)
       call sep%init(thermo, separator_pressure)

       call sep%separate(rate = -10._dp, enthalpy = 500.e3_dp)
       call separator_test_case("water", 0._dp, -10._dp, 500.e3_dp, 0._dp, 0._dp)

       call sep%separate(rate = -10._dp, enthalpy = 3000.e3_dp)
       call separator_test_case("steam", 1._dp, 0._dp, 0._dp, -10._dp, 3000.e3_dp)

       call sep%separate(rate = -10._dp, enthalpy = 1200.e3_dp)
       call separator_test_case("two-phase", 0.21709153586628488_dp, &
            -7.829084641337152_dp, 762682.8443354106_dp, &
            -2.1709153586628487_dp, 2777119.5376846623_dp)

       deallocate(data)
       call thermo%destroy()

    end if

  contains

    subroutine separator_test_case(name, steam_fraction, water_rate, water_enthalpy, &
         steam_rate, steam_enthalpy)

      character(*), intent(in) :: name
      PetscReal, intent(in) :: steam_fraction
      PetscReal, intent(in) :: water_rate, water_enthalpy
      PetscReal, intent(in) :: steam_rate, steam_enthalpy

      call test%assert(steam_fraction, sep%steam_fraction, trim(name) // " steam fraction")
      call test%assert(water_rate, sep%water_rate, trim(name) // " water rate")
      call test%assert(water_enthalpy, sep%water_enthalpy, trim(name) // " water enthalpy")
      call test%assert(steam_rate, sep%steam_rate, trim(name) // " steam rate")
      call test%assert(steam_enthalpy, sep%steam_enthalpy, trim(name) // " steam enthalpy")

    end subroutine separator_test_case

  end subroutine test_separator

!------------------------------------------------------------------------

end module separator_test
