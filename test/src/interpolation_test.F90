module interpolation_test

  ! Tests for interpolation module

  use kinds_module
  use mpi_module
  use fruit
  use interpolation_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>

  PetscReal, dimension(5,2), parameter :: data5 = reshape([&
       0._dp, 2.1_dp, 3.7_dp,  6.3_dp,  8.9_dp, &
       1._dp, 2.0_dp, 0.5_dp, -1.1_dp, -0.1_dp], &
       [5,2])
  PetscReal, parameter :: tol = 1.e-9_dp

  public :: test_interpolation_linear, test_interpolation_step, &
       test_interpolation_step_average

contains

!------------------------------------------------------------------------

  subroutine test_interpolation_linear

    ! Test interpolation_table_type with linear interpolation

    type(interpolation_table_type) :: table

    if (mpi%rank == mpi%output_rank) then

       call table%init(data5) ! default linear interpolation

       call assert_equals(1._dp, table%interpolate(-0.5_dp), tol, "-0.5")
       call assert_equals(1, table%index, "-0.5 index")

       call assert_equals(1._dp, table%interpolate(0.0_dp), tol, "0.0")
       call assert_equals(1, table%index, "0.0 index")

       call assert_equals(1.4761904761904763_dp, table%interpolate(1.0_dp), tol, "1.0")
       call assert_equals(1, table%index, "1.0 index")

       call assert_equals(0.007692307692307665_dp, table%interpolate(4.5_dp), tol, "4.5")
       call assert_equals(3, table%index, "4.5 index")

       call assert_equals(0.59375_dp, table%interpolate(3.6_dp), tol, "3.6")
       call assert_equals(2, table%index, "3.6 index")

       call assert_equals(-1.1_dp, table%interpolate(6.3_dp), tol, "6.3")
       call assert_equals(4, table%index, "6.3 index")

       call assert_equals(-0.1_dp, table%interpolate(10.0_dp), tol, "10.0")
       call assert_equals(5, table%index, "10.0 index")

       call table%destroy()

    end if

  end subroutine test_interpolation_linear

!------------------------------------------------------------------------

  subroutine test_interpolation_step

    ! Test step interpolation

    type(interpolation_table_type) :: table

    if (mpi%rank == mpi%output_rank) then

       call table%init(data5, INTERP_STEP)

       call assert_equals(1._dp, table%interpolate(-0.5_dp), tol, "-0.5 step")

       call assert_equals(1._dp, table%interpolate(0.0_dp), tol, "0.0 step")

       call assert_equals(1._dp, table%interpolate(1.0_dp), tol, "1.0 step")

       call assert_equals(0.5_dp, table%interpolate(4.5_dp), tol, "4.5 step")

       call assert_equals(2.0_dp, table%interpolate(3.6_dp), tol, "3.6 step")

       call assert_equals(-1.1_dp, table%interpolate(6.3_dp), tol, "6.3 step")

       call assert_equals(-0.1_dp, table%interpolate(10.0_dp), tol, "10.0 step")

       call table%destroy()

    end if

  end subroutine test_interpolation_step

!------------------------------------------------------------------------

  subroutine test_interpolation_step_average

    ! Test step average interpolation

    type(interpolation_table_type) :: table

    if (mpi%rank == mpi%output_rank) then

       call table%init(data5, INTERP_STEP_AVERAGE)

       call assert_equals(1._dp, table%interpolate(-0.5_dp), tol, "-0.5 step")

       call assert_equals(1._dp, table%interpolate(0.0_dp), tol, "0.0 step")

       call assert_equals(1.5_dp, table%interpolate(1.0_dp), tol, "1.0 step")

       call assert_equals(-0.3_dp, table%interpolate(4.5_dp), tol, "4.5 step")

       call assert_equals(1.25_dp, table%interpolate(3.6_dp), tol, "3.6 step")

       call assert_equals(-0.6_dp, table%interpolate(6.3_dp), tol, "6.3 step")

       call assert_equals(-0.1_dp, table%interpolate(10.0_dp), tol, "10.0 step")

       call table%destroy()

    end if

  end subroutine test_interpolation_step_average

!------------------------------------------------------------------------

end module interpolation_test
