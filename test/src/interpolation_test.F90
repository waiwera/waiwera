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

  public :: test_interpolation_linear, test_interpolation_single, &
       test_interpolation_step, test_average_linear, test_average_step, &
       test_average_linear_integration

contains

!------------------------------------------------------------------------

  subroutine test_interpolation_linear

    ! Linear interpolation

    type(interpolation_table_type) :: table

    if (mpi%rank == mpi%output_rank) then

       call table%init(data5) ! default linear interpolation

       call assert_equals(1._dp, table%interpolate(-0.5_dp), tol, "-0.5")
       call assert_equals(0, table%index, "-0.5 index")

       call assert_equals(1._dp, table%interpolate(0.0_dp), tol, "0.0")
       call assert_equals(0, table%index, "0.0 index")

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

  subroutine test_interpolation_single

    ! Single-value data table

    type(interpolation_table_type) :: table
    ! Locals:
    PetscReal, dimension(1,2), parameter :: data1 = reshape(&
         [0._dp, 2._dp], [1,2])

    if (mpi%rank == mpi%output_rank) then

       call table%init(data1)

       call assert_equals(2._dp, table%interpolate(-0.5_dp), tol, "-0.5")

       call assert_equals(2._dp, table%interpolate(0.0_dp), tol, "0.0")

       call assert_equals(2._dp, table%interpolate(1.1_dp), tol, "1.1")

       call table%destroy()
    end if

  end subroutine test_interpolation_single

!------------------------------------------------------------------------

  subroutine test_interpolation_step

    ! Step interpolation

    type(interpolation_table_type) :: table

    if (mpi%rank == mpi%output_rank) then

       call table%init(data5, INTERP_STEP)

       call assert_equals(1._dp, table%interpolate(-0.5_dp), tol, "-0.5")

       call assert_equals(1._dp, table%interpolate(0.0_dp), tol, "0.0")

       call assert_equals(1._dp, table%interpolate(1.0_dp), tol, "1.0")

       call assert_equals(0.5_dp, table%interpolate(4.5_dp), tol, "4.5")

       call assert_equals(2.0_dp, table%interpolate(3.6_dp), tol, "3.6")

       call assert_equals(-1.1_dp, table%interpolate(6.3_dp), tol, "6.3")

       call assert_equals(-0.1_dp, table%interpolate(10.0_dp), tol, "10.0")

       call table%destroy()

    end if

  end subroutine test_interpolation_step

!------------------------------------------------------------------------

  subroutine test_average_linear

    ! Linear interpolation average

    type(interpolation_table_type) :: table

    if (mpi%rank == mpi%output_rank) then

       call table%init(data5)

       call assert_equals(1._dp, table%average([-0.5_dp, -0.1_dp]), tol, "[-0.5, -0.1]")

       call assert_equals(1.0238095238095237_dp, table%average([-0.5_dp, 0.1_dp]), &
            tol, "[-0.5, 0.1]")

       call assert_equals(1.5_dp, table%average([0.1_dp, 2._dp]), tol, "[0.1, 2.]")

       call assert_equals(1.1019345238095237_dp, table%average([0.1_dp, 3._dp]), tol, "[0.1, 3.]")

       call assert_equals(0.11586538461538454_dp, table%average([3.1_dp, 7._dp]), tol, "[3.1, 7.]")

       call assert_equals(-0.27307692307692316_dp, table%average([8._dp, 12._dp]), tol, "[8., 12.]")

       call assert_equals(1.4761904761904763_dp, table%average([1._dp, 1._dp]), tol, "[1., 1.]")

       call table%destroy()

    end if

  end subroutine test_average_linear

!------------------------------------------------------------------------

  subroutine test_average_step

    ! Step interpolation average

    type(interpolation_table_type) :: table

    if (mpi%rank == mpi%output_rank) then

       call table%init(data5, INTERP_STEP)

       call assert_equals(1._dp, table%average([-0.5_dp, -0.1_dp]), tol, "[-0.5, -0.1]")

       call assert_equals(1._dp, table%average([-0.5_dp, 0.1_dp]), tol, "[-0.5, 0.1]")

       call assert_equals(1._dp, table%average([0.1_dp, 2._dp]), tol, "[0.1, 2.]")

       call assert_equals(1.5_dp, table%average([0.1_dp, 3._dp]), tol, "[0.1, 3.]")

       call assert_equals(0.45_dp, table%average([3.1_dp, 7._dp]), tol, "[3.1, 7.]")

       call assert_equals(-0.6_dp, table%average([8._dp, 12._dp]), tol, "[8., 12.]")

       call assert_equals(1._dp, table%average([1._dp, 1._dp]), tol, "[1., 1.]")

       call table%destroy()

    end if

  end subroutine test_average_step

!------------------------------------------------------------------------

  subroutine test_average_linear_integration

    ! Linear interpolation integration averaging

    type(interpolation_table_type) :: table

    if (mpi%rank == mpi%output_rank) then

       call table%init(data5, INTERP_LINEAR, INTERP_AVERAGING_INTEGRATE)

       call assert_equals(1._dp, table%average([-0.5_dp, -0.1_dp]), tol, "[-0.5, -0.1]")

       call assert_equals(1.003968253968254_dp, table%average([-0.5_dp, 0.1_dp]), tol, "[-0.5, 0.1]")

       call assert_equals(1.5_dp, table%average([0.1_dp, 2._dp]), tol, "[0.1, 2.]")

       call assert_equals(1.5406660509031198_dp, table%average([0.1_dp, 3._dp]), tol, "[0.1, 3.]")

       call assert_equals(-0.2530818540433925_dp, table%average([3.1_dp, 7._dp]), tol, "[3.1, 7.]")

       call assert_equals(-0.1389423076923077_dp, table%average([8._dp, 12._dp]), tol, "[8., 12.]")

       call assert_equals(-0.1_dp, table%average([9._dp, 12._dp]), tol, "[9., 12.]")

       call assert_equals(1.4761904761904763_dp, table%average([1._dp, 1._dp]), tol, "[1., 1.]")

       call table%destroy()

    end if

  end subroutine test_average_linear_integration

!------------------------------------------------------------------------

end module interpolation_test
