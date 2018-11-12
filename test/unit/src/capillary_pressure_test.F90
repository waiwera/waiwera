module capillary_pressure_test

  ! Tests for capillary pressure module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use zofu
  use capillary_pressure_module
  use fson

  implicit none
  private

  public :: setup, teardown
  public :: test_capillary_pressure_zero, &
       test_capillary_pressure_linear, &
       test_capillary_pressure_van_genuchten, &
       test_capillary_pressure_table

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

  subroutine test_capillary_pressure_zero(test)

    ! Zero capillary pressure function

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(capillary_pressure_zero_type) :: cp
    type(fson_value), pointer :: json
    character(100), parameter :: json_str = '{"type": "zero"}'
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscErrorCode :: err
    PetscReal, parameter :: t = 20._dp ! dummy temperature

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse(str = json_str)
    call cp%init(json, err = err)
    call fson_destroy(json)

    if (rank == 0) then

       call test%assert(0, err, "error")
       call test%assert("zero", cp%name, "Name")
       call test%assert(0._dp, cp%value(0._dp, t), "0")
       call test%assert(0._dp, cp%value(0.6_dp, t), "0.6")
       call test%assert(0._dp, cp%value(0.9_dp, t), "0.9")

    end if

    call cp%destroy()

  end subroutine test_capillary_pressure_zero

!------------------------------------------------------------------------

  subroutine test_capillary_pressure_linear(test)

    ! Linear capillary pressure function

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(capillary_pressure_linear_type) :: cp
    type(fson_value), pointer :: json
    character(100), parameter :: json_str = &
         '{"type": "linear", "saturation_limits": [0.1, 0.8], "pressure": 0.2e5}'
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscErrorCode :: err
    PetscReal, parameter :: t = 20._dp ! dummy temperature

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse(str = json_str)
    call cp%init(json, err = err)
    call fson_destroy(json)

    if (rank == 0) then

       call test%assert(0, err, "error")
       call test%assert("linear", cp%name, "Name")
       call test%assert(-0.2e5_dp, cp%value(0._dp, t), "0")
       call test%assert(-0.0571428571428e5_dp, cp%value(0.6_dp, t), "0.6")
       call test%assert(0._dp, cp%value(0.9_dp, t), "0.9")

    end if

    call cp%destroy()

  end subroutine test_capillary_pressure_linear

!------------------------------------------------------------------------

  subroutine test_capillary_pressure_van_genuchten(test)

    ! Van Genuchten capillary pressure function

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(capillary_pressure_van_genuchten_type) :: cp
    type(fson_value), pointer :: json
    character(100) :: json_str
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscErrorCode :: err
    PetscReal, parameter :: t = 20._dp ! dummy temperature

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json_str = '{"type": "van Genuchten", "P0": 0.2e5, "lambda": 0.5, ' // &
         '"slr": 0.1, "sls": 0.8}'
    json => fson_parse(str = json_str)
    call cp%init(json, err = err)
    call fson_destroy(json)

    if (rank == 0) then
       call test%assert(0, err, "error")
       call test%assert("van Genuchten", cp%name, "Name")
       call test%assert(-0.2e5_dp, cp%value(0._dp, t), "0")
       call test%assert(-6.99714227381e5_dp, cp%value(0.12_dp, t), "0.12")
       call test%assert(-2.79284800875e5_dp, cp%value(0.15_dp, t), "0.15")
       call test%assert(-0.911652955412e5_dp, cp%value(0.25_dp, t), "0.25")
       call test%assert(-0.195959179423e5_dp, cp%value(0.6_dp, t), "0.6")
       call test%assert(0._dp, cp%value(0.9_dp, t), "0.9")
    end if

    call cp%destroy()

    json_str = '{"type": "van Genuchten", "P0": 0.2e5, "lambda": 0.5, ' // &
         '"slr": 0.1, "sls": 0.8, "Pmax": 6.e5}'
    json => fson_parse(str = json_str)
    call cp%init(json, err = err)
    call fson_destroy(json)

    if (rank == 0) then
       call test%assert(0, err, "error")
       call test%assert(-6.e5_dp, cp%value(0.12_dp, t), "0.12 Pmax")
       call test%assert(-0.195959179423e5_dp, cp%value(0.6_dp, t), "0.6 Pmax")
    end if

    call cp%destroy()

  end subroutine test_capillary_pressure_van_genuchten

!------------------------------------------------------------------------

  subroutine test_capillary_pressure_table(test)

    ! Table capillary pressure function

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(capillary_pressure_table_type) :: cp
    type(fson_value), pointer :: json
    character(100) :: json_str
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscErrorCode :: err
    PetscReal, parameter :: t = 20._dp ! dummy temperature

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json_str = '{"type": "table", "pressure": ' // &
         '[[0, -5.e5], [0.4, -1.e5], [0.7, 0]]}'
    json => fson_parse(str = json_str)
    call cp%init(json, err = err)
    call fson_destroy(json)

    if (rank == 0) then
       call test%assert(0, err, "error")
       call test%assert("table", cp%name, "Name")
       call test%assert(-5.e5_dp, cp%value(0._dp, t), "0")
       call test%assert(-2.e5_dp, cp%value(0.3_dp, t), "0.3")
       call test%assert(0.e5_dp, cp%value(0.7_dp, t), "0.7")
       call test%assert(0.e5_dp, cp%value(0.9_dp, t), "0.9")
    end if

    call cp%destroy() 

  end subroutine test_capillary_pressure_table

!------------------------------------------------------------------------

end module capillary_pressure_test
