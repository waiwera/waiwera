module relative_permeability_test

  ! Tests for relative permeability module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use fruit
  use relative_permeability_module
  use fson

  implicit none
  private

  PetscReal, parameter :: tol = 1.e-6_dp

public :: test_relative_permeability_linear, &
     test_relative_permeability_pickens, &
     test_relative_permeability_corey, &
     test_relative_permeability_grant, &
     test_relative_permeability_van_genuchten, &
     test_relative_permeability_fully_mobile

contains

!------------------------------------------------------------------------

  subroutine relative_permeability_case(sl, relp, expected)

    ! Test one case

    PetscReal, intent(in) :: sl
    class(relative_permeability_type), intent(in) :: relp
    PetscReal, intent(in) :: expected(2)
    ! Locals:
    PetscReal :: rp(2)
    character(12) :: msg

    rp = relp%values(sl)
    write(msg, '(a, f4.2)') ' sl = ', sl
    call assert_equals(expected(1), rp(1), tol, "Liquid" // trim(msg))
    call assert_equals(expected(2), rp(2), tol, "Vapour" // trim(msg))

  end subroutine relative_permeability_case

!------------------------------------------------------------------------

  subroutine test_relative_permeability_linear

    ! Linear relative permeability functions

    type(relative_permeability_linear_type) :: linear
    ! Locals:
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

       call assert_equals("Linear", linear%name, "Name")

       call relative_permeability_case(0.01_dp, linear, [0._dp, 1._dp])
       call relative_permeability_case(0.2_dp, linear, [1._dp / 7._dp, 1._dp])
       call relative_permeability_case(0.5_dp, linear, [4._dp / 7._dp, 4._dp / 9._dp])
       call relative_permeability_case(0.9_dp, linear, [1._dp, 0._dp])

    end if

  end subroutine test_relative_permeability_linear

!------------------------------------------------------------------------

  subroutine test_relative_permeability_pickens

    ! Pickens relative permeability functions

    type(relative_permeability_pickens_type) :: pickens
    ! Locals:
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

       call assert_equals("Pickens", pickens%name, "Name")

       call relative_permeability_case(0.01_dp, pickens, [1.e-4_dp, 1._dp])
       call relative_permeability_case(0.5_dp, pickens, [0.25_dp, 1._dp])
       call relative_permeability_case(0.9_dp, pickens, [0.81_dp, 1._dp])

    end if

  end subroutine test_relative_permeability_pickens

!------------------------------------------------------------------------

  subroutine test_relative_permeability_corey

    ! Corey's relative permeability functions

    type(relative_permeability_corey_type) :: corey
    ! Locals:
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

       call assert_equals("Corey", corey%name, "Name")

       call relative_permeability_case(0.01_dp, corey, [0._dp, 1._dp])
       call relative_permeability_case(0.5_dp, corey, [1._dp / 81._dp, 32._dp/ 81._dp])
       call relative_permeability_case(0.95_dp, corey, [1._dp, 0._dp])

    end if

  end subroutine test_relative_permeability_corey

!------------------------------------------------------------------------

  subroutine test_relative_permeability_grant

    ! Grant's relative permeability functions

    type(relative_permeability_grant_type) :: grant
    ! Locals:
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

       call assert_equals("Grant", grant%name, "Name")

       call relative_permeability_case(0.01_dp, grant, [0._dp, 1._dp])
       call relative_permeability_case(0.5_dp, grant, [1._dp / 81._dp, 80._dp/ 81._dp])
       call relative_permeability_case(0.95_dp, grant, [1._dp, 0._dp])

    end if

  end subroutine test_relative_permeability_grant

!------------------------------------------------------------------------

  subroutine test_relative_permeability_van_genuchten

    ! van Genuchten relative permeability functions

    type(relative_permeability_van_genuchten_type) :: vg
    ! Locals:
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

       call assert_equals("van Genuchten", vg%name, "Name")

       call relative_permeability_case(0.01_dp, vg, [0._dp, 1._dp])
       call relative_permeability_case(0.25_dp, vg, &
            [0.00024977947758877213_dp, 0.9997502205224112_dp])
       call relative_permeability_case(0.5_dp, vg, &
            [0.024315039984298164_dp, 0.9756849600157018_dp])
       call relative_permeability_case(0.75_dp, vg, &
            [0.38106285486468433_dp, 0.6189371451353156_dp])
       call relative_permeability_case(0.95_dp, vg, [1._dp, 0._dp])

    end if

  end subroutine test_relative_permeability_van_genuchten

!------------------------------------------------------------------------

  subroutine test_relative_permeability_fully_mobile

    ! Fully mobile relative permeability functions

    type(relative_permeability_fully_mobile_type) :: mobile
    ! Locals:
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

       call assert_equals("Fully mobile", mobile%name, "Name")

       call relative_permeability_case(0.2_dp, mobile, [1._dp, 1._dp])
       call relative_permeability_case(0.9_dp, mobile, [1._dp, 1._dp])

    end if

  end subroutine test_relative_permeability_fully_mobile
  
!-----------------------------------------------------------------------

end module relative_permeability_test
