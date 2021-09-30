module relative_permeability_test

  ! Tests for relative permeability module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use zofu
  use relative_permeability_module
  use fson

  implicit none
  private

  public :: setup, teardown
  public :: test_relative_permeability_linear, &
       test_relative_permeability_pickens, &
       test_relative_permeability_corey, &
       test_relative_permeability_grant, &
       test_relative_permeability_van_genuchten, &
       test_relative_permeability_table, &
       test_relative_permeability_fully_mobile

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

  subroutine relative_permeability_case(test, sl, relp, expected)

    ! Test one case

    class(unit_test_type), intent(in out) :: test
    PetscReal, intent(in) :: sl
    class(relative_permeability_type), intent(in out) :: relp
    PetscReal, intent(in) :: expected(2)
    ! Locals:
    PetscReal :: rp(2)
    character(12) :: msg

    rp = relp%values(sl)
    write(msg, '(a, f4.2)') ' sl = ', sl
    call test%assert(expected(1), rp(1), "Liquid" // trim(msg))
    call test%assert(expected(2), rp(2), "Vapour" // trim(msg))

  end subroutine relative_permeability_case

!------------------------------------------------------------------------

  subroutine test_relative_permeability_linear(test)

    ! Linear relative permeability functions

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(relative_permeability_linear_type) :: linear
    type(fson_value), pointer :: json
    character(100), parameter :: json_str = &
         '{"type": "linear", "liquid": [0.1, 0.8], "vapour": [0.3, 0.75]}'
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse(str = json_str)
    call linear%init(json)
    call fson_destroy(json)

    if (rank == 0) then

       call test%assert("Linear", linear%name, "Name")

       call relative_permeability_case(test, 0.01_dp, linear, [0._dp, 1._dp])
       call relative_permeability_case(test, 0.2_dp, linear, [1._dp / 7._dp, 1._dp])
       call relative_permeability_case(test, 0.5_dp, linear, [4._dp / 7._dp, 4._dp / 9._dp])
       call relative_permeability_case(test, 0.9_dp, linear, [1._dp, 0._dp])

    end if

    call linear%destroy()

  end subroutine test_relative_permeability_linear

!------------------------------------------------------------------------

  subroutine test_relative_permeability_pickens(test)

    ! Pickens relative permeability functions

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(relative_permeability_pickens_type) :: pickens
    type(fson_value), pointer :: json
    character(100), parameter :: json_str = &
         '{"type": "Pickens", "power": 2.0}'
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse(str = json_str)
    call pickens%init(json)
    call fson_destroy(json)

    if (rank == 0) then

       call test%assert("Pickens", pickens%name, "Name")
       call relative_permeability_case(test, 0.01_dp, pickens, [1.e-4_dp, 1._dp])
       call relative_permeability_case(test, 0.5_dp, pickens, [0.25_dp, 1._dp])
       call relative_permeability_case(test, 0.9_dp, pickens, [0.81_dp, 1._dp])

    end if

    call pickens%destroy()

  end subroutine test_relative_permeability_pickens

!------------------------------------------------------------------------

  subroutine test_relative_permeability_corey(test)

    ! Corey relative permeability functions

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(relative_permeability_corey_type) :: corey
    type(fson_value), pointer :: json
    character(100), parameter :: json_str = &
    '{"type": "Corey", "slr": 0.3, "ssr": 0.1}'
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse(str = json_str)
    call corey%init(json)
    call fson_destroy(json)

    if (rank == 0) then

       call test%assert("Corey", corey%name, "Name")
       call relative_permeability_case(test, 0.01_dp, corey, [0._dp, 1._dp])
       call relative_permeability_case(test, 0.5_dp, corey, [1._dp / 81._dp, 32._dp/ 81._dp])
       call relative_permeability_case(test, 0.95_dp, corey, [1._dp, 0._dp])

    end if

    call corey%destroy()

  end subroutine test_relative_permeability_corey

!------------------------------------------------------------------------

  subroutine test_relative_permeability_grant(test)

    ! Grant relative permeability functions

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(relative_permeability_grant_type) :: grant
    type(fson_value), pointer :: json
    character(100), parameter :: json_str = &
    '{"type": "Grant", "slr": 0.3, "ssr": 0.1}'
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse(str = json_str)
    call grant%init(json)
    call fson_destroy(json)

    if (rank == 0) then

       call test%assert("Grant", grant%name, "Name")
       call relative_permeability_case(test, 0.01_dp, grant, [0._dp, 1._dp])
       call relative_permeability_case(test, 0.5_dp, grant, [1._dp / 81._dp, 80._dp/ 81._dp])
       call relative_permeability_case(test, 0.95_dp, grant, [1._dp, 0._dp])

    end if

    call grant%destroy()

  end subroutine test_relative_permeability_grant

!------------------------------------------------------------------------

  subroutine test_relative_permeability_van_genuchten(test)

    ! van Genuchten relative permeability functions

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(relative_permeability_van_genuchten_type) :: vg
    type(fson_value), pointer :: json
    character(100), parameter :: json_str = &
    '{"type": "van Genuchten", "slr": 0.1, "sls": 0.8, "lambda": 0.5}'
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse(str = json_str)
    call vg%init(json)
    call fson_destroy(json)

    if (rank == 0) then

       call test%assert("van Genuchten", vg%name, "Name")
       call relative_permeability_case(test, 0.01_dp, vg, [0._dp, 1._dp])
       call relative_permeability_case(test, 0.25_dp, vg, &
            [0.00024977947758877213_dp, 0.9997502205224112_dp])
       call relative_permeability_case(test, 0.5_dp, vg, &
            [0.024315039984298164_dp, 0.9756849600157018_dp])
       call relative_permeability_case(test, 0.75_dp, vg, &
            [0.38106285486468433_dp, 0.6189371451353156_dp])
       call relative_permeability_case(test, 0.95_dp, vg, [1._dp, 0._dp])

    end if

    call vg%destroy()

  end subroutine test_relative_permeability_van_genuchten

!------------------------------------------------------------------------

  subroutine test_relative_permeability_table(test)

    ! Table relative permeability functions

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(relative_permeability_table_type) :: rp
    type(fson_value), pointer :: json
    character(160) :: json_str
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json_str = '{"type": "table", ' // &
         '"liquid": [[0,0], [0.7, 0.01], [0.95, 0.99], [1,1]], ' // &
         '"vapour": [[0,0], [0.05, 0.01], [0.3, 0.99], [1,1]]}'

    json => fson_parse(str = json_str)
    call rp%init(json)
    call fson_destroy(json)

    if (rank == 0) then

       call test%assert("table", rp%name, "Name")
       call relative_permeability_case(test, 0._dp, rp, [0._dp, 1._dp])
       call relative_permeability_case(test, 0.3_dp, rp, &
            [0.01_dp * 3._dp / 7._dp, (4._dp + 3._dp * 0.99_dp) / 7._dp])
       call relative_permeability_case(test, 0.7_dp, rp, [0.01_dp, 0.99_dp])
       call relative_permeability_case(test, 0.9_dp, rp, &
            [0.2_dp * 0.01_dp + 0.8_dp * 0.99_dp, 0.8_dp * 0.01_dp + 0.2_dp * 0.99_dp])
       call relative_permeability_case(test, 1._dp, rp, [1._dp, 0._dp])

    end if

    call rp%destroy()

  end subroutine test_relative_permeability_table

!------------------------------------------------------------------------

  subroutine test_relative_permeability_fully_mobile(test)

    ! Fully mobile relative permeability functions

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(relative_permeability_fully_mobile_type) :: mobile
    type(fson_value), pointer :: json
    character(100), parameter :: json_str = &
    '{"type": "Fully mobile"}'
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse(str = json_str)
    call mobile%init(json)
    call fson_destroy(json)

    if (rank == 0) then

       call test%assert("Fully mobile", mobile%name, "Name")
       call relative_permeability_case(test, 0.2_dp, mobile, [1._dp, 1._dp])
       call relative_permeability_case(test, 0.9_dp, mobile, [1._dp, 1._dp])

    end if

  end subroutine test_relative_permeability_fully_mobile
  
!-----------------------------------------------------------------------

end module relative_permeability_test
