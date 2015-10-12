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

  subroutine test_relative_permeability_linear

    ! Linear relative permeability functions

    type(relative_permeability_linear_type) :: linear
    ! Locals:
    PetscReal :: sl, rp(2)
    PetscReal, dimension(2), parameter :: &
         liquid_limits = [0.1_dp, 0.8_dp], &
         vapour_limits = [0.3_dp, 0.75_dp]

    call linear%init(liquid_limits, vapour_limits)

    if (mpi%rank == mpi%output_rank) then

       call assert_equals("Linear", linear%name, "Name")

       sl = 0.01_dp
       rp = linear%values(sl)
       call assert_equals(0._dp, rp(1), tol, "Liquid sl = 0.01")
       call assert_equals(1._dp, rp(2), tol, "Vapour sl = 0.01")

       sl = 0.2_dp
       rp = linear%values(sl)
       call assert_equals(1._dp / 7._dp, rp(1), tol, "Liquid sl = 0.2")
       call assert_equals(1._dp, rp(2), tol, "Vapour sl = 0.2")

       sl = 0.5_dp
       rp = linear%values(sl)
       call assert_equals(4._dp / 7._dp, rp(1), tol, "Liquid sl = 0.5")
       call assert_equals(4._dp / 9._dp, rp(2), tol, "Vapour sl = 0.5")

       sl = 0.9_dp
       rp = linear%values(sl)
       call assert_equals(1._dp, rp(1), tol, "Liquid sl = 0.9")
       call assert_equals(0._dp, rp(2), tol, "Vapour sl = 0.9")

    end if

  end subroutine test_relative_permeability_linear

!------------------------------------------------------------------------

  subroutine test_relative_permeability_corey

    ! Corey's relative permeability functions

    type(relative_permeability_corey_type) :: corey
    ! Locals:
    PetscReal, parameter :: slr = 0.3_dp, ssr = 0.1_dp
    PetscReal :: sl, rp(2)

    call corey%init(slr, ssr)

    if (mpi%rank == mpi%output_rank) then

       call assert_equals("Corey", corey%name, "Name")

       sl = 0.01_dp
       rp = corey%values(sl)
       call assert_equals(0._dp, rp(1), tol, "Liquid sl = 0.01")
       call assert_equals(1._dp, rp(2), tol, "Vapour sl = 0.01")

       sl = 0.5_dp
       rp = corey%values(sl)
       call assert_equals(1._dp / 81._dp, rp(1), tol, "Liquid sl = 0.5")
       call assert_equals(32._dp/ 81._dp, rp(2), tol, "Vapour sl = 0.5")

       sl = 0.95_dp
       rp = corey%values(sl)
       call assert_equals(1._dp, rp(1), tol, "Liquid sl = 0.5")
       call assert_equals(0._dp, rp(2), tol, "Vapour sl = 0.5")

    end if

  end subroutine test_relative_permeability_corey

!------------------------------------------------------------------------

  subroutine test_relative_permeability_grant

    ! Grant's relative permeability functions

    type(relative_permeability_grant_type) :: grant
    ! Locals:
    PetscReal, parameter :: slr = 0.3_dp, ssr = 0.1_dp
    PetscReal :: sl, rp(2)

    call grant%init(slr, ssr)

    if (mpi%rank == mpi%output_rank) then

       call assert_equals("Grant", grant%name, "Name")

       sl = 0.01_dp
       rp = grant%values(sl)
       call assert_equals(0._dp, rp(1), tol, "Liquid sl = 0.01")
       call assert_equals(1._dp, rp(2), tol, "Vapour sl = 0.01")

       sl = 0.5_dp
       rp = grant%values(sl)
       call assert_equals(1._dp / 81._dp, rp(1), tol, "Liquid sl = 0.5")
       call assert_equals(80._dp/ 81._dp, rp(2), tol, "Vapour sl = 0.5")

       sl = 0.95_dp
       rp = grant%values(sl)
       call assert_equals(1._dp, rp(1), tol, "Liquid sl = 0.5")
       call assert_equals(0._dp, rp(2), tol, "Vapour sl = 0.5")

    end if

  end subroutine test_relative_permeability_grant

!-----------------------------------------------------------------------

end module relative_permeability_test
