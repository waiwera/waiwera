module relative_permeability_test

  ! Tests for relative permeability module

  use kinds_module
  use mpi_module
  use fruit
  use relative_permeability_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>

  PetscReal, parameter :: tol = 1.e-6_dp

public :: test_relative_permeability_linear, &
     test_relative_permeability_corey, &
     test_relative_permeability_grant

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
    PetscReal, dimension(2), parameter :: &
         liquid_limits = [0.1_dp, 0.8_dp], &
         vapour_limits = [0.3_dp, 0.75_dp]

    call linear%init(liquid_limits, vapour_limits)

    if (mpi%rank == mpi%output_rank) then

       call assert_equals("Linear", linear%name, "Name")

       call relative_permeability_case(0.01_dp, linear, [0._dp, 1._dp])
       call relative_permeability_case(0.2_dp, linear, [1._dp / 7._dp, 1._dp])
       call relative_permeability_case(0.5_dp, linear, [4._dp / 7._dp, 4._dp / 9._dp])
       call relative_permeability_case(0.9_dp, linear, [1._dp, 0._dp])

    end if

  end subroutine test_relative_permeability_linear

!------------------------------------------------------------------------

  subroutine test_relative_permeability_corey

    ! Corey's relative permeability functions

    type(relative_permeability_corey_type) :: corey
    ! Locals:
    PetscReal, parameter :: slr = 0.3_dp, ssr = 0.1_dp

    call corey%init(slr, ssr)

    if (mpi%rank == mpi%output_rank) then

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
    PetscReal, parameter :: slr = 0.3_dp, ssr = 0.1_dp

    call grant%init(slr, ssr)

    if (mpi%rank == mpi%output_rank) then

       call assert_equals("Grant", grant%name, "Name")

       call relative_permeability_case(0.01_dp, grant, [0._dp, 1._dp])
       call relative_permeability_case(0.5_dp, grant, [1._dp / 81._dp, 80._dp/ 81._dp])
       call relative_permeability_case(0.95_dp, grant, [1._dp, 0._dp])

    end if

  end subroutine test_relative_permeability_grant

!-----------------------------------------------------------------------

end module relative_permeability_test
