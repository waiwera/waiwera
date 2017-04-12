module interpolation_test

  ! Tests for interpolation module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use fruit
  use interpolation_module

  implicit none
  private

  PetscReal, dimension(5,2), parameter :: data5 = reshape([&
       0._dp, 2.1_dp, 3.7_dp,  6.3_dp,  8.9_dp, &
       1._dp, 2.0_dp, 0.5_dp, -1.1_dp, -0.1_dp], &
       [5,2])
  PetscReal, dimension(3, 4), parameter :: data3 = reshape([&
       0._dp, 1._dp, 2._dp, &
       1._dp, 2._dp, 3._dp, &
       2._dp, 3._dp, 4._dp, &
       3._dp, 4._dp, 5._dp], &
       [3, 4])
  PetscReal, parameter :: tol = 1.e-9_dp

  public :: test_interpolation_linear, test_interpolation_single, &
       test_interpolation_step, test_average_linear, test_average_step, &
       test_average_linear_integration, test_array_interpolator, &
       test_ramp_interpolate,  test_interpolation_linear_array

contains

!------------------------------------------------------------------------

  subroutine test_interpolation_linear

    ! Linear interpolation

    type(interpolation_table_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call table%init(data5) ! default linear interpolation
       call assert_equals(1, table%dim, "dim")

       call assert_equals(1._dp, table%interpolate(-0.5_dp, 1), tol, "-0.5")
       call assert_equals(0, table%coord%index, "-0.5 index")

       call assert_equals(1._dp, table%interpolate(0.0_dp, 1), tol, "0.0")
       call assert_equals(0, table%coord%index, "0.0 index")

       call assert_equals(1.4761904761904763_dp, table%interpolate(1.0_dp, 1), tol, "1.0")
       call assert_equals(1, table%coord%index, "1.0 index")

       call assert_equals(0.007692307692307665_dp, table%interpolate(4.5_dp, 1), tol, "4.5")
       call assert_equals(3, table%coord%index, "4.5 index")

       call assert_equals(0.59375_dp, table%interpolate(3.6_dp, 1), tol, "3.6")
       call assert_equals(2, table%coord%index, "3.6 index")

       call assert_equals(-1.1_dp, table%interpolate(6.3_dp, 1), tol, "6.3")
       call assert_equals(4, table%coord%index, "6.3 index")

       call assert_equals(-0.1_dp, table%interpolate(10.0_dp, 1), tol, "10.0")
       call assert_equals(5, table%coord%index, "10.0 index")

       call table%destroy()

    end if

  end subroutine test_interpolation_linear

!------------------------------------------------------------------------

  subroutine test_interpolation_single

    ! Single-value data table

    type(interpolation_table_type) :: table
    PetscReal, dimension(1,2), parameter :: data1 = reshape(&
         [0._dp, 2._dp], [1,2])
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call table%init(data1)

       call assert_equals(2._dp, table%interpolate(-0.5_dp, 1), tol, "-0.5")

       call assert_equals(2._dp, table%interpolate(0.0_dp, 1), tol, "0.0")

       call assert_equals(2._dp, table%interpolate(1.1_dp, 1), tol, "1.1")

       call table%destroy()
    end if

  end subroutine test_interpolation_single

!------------------------------------------------------------------------

  subroutine test_interpolation_step

    ! Step interpolation

    type(interpolation_table_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call table%init(data5, INTERP_STEP)

       call assert_equals(1._dp, table%interpolate(-0.5_dp, 1), tol, "-0.5")

       call assert_equals(1._dp, table%interpolate(0.0_dp, 1), tol, "0.0")

       call assert_equals(1._dp, table%interpolate(1.0_dp, 1), tol, "1.0")

       call assert_equals(0.5_dp, table%interpolate(4.5_dp, 1), tol, "4.5")

       call assert_equals(2.0_dp, table%interpolate(3.6_dp, 1), tol, "3.6")

       call assert_equals(-1.1_dp, table%interpolate(6.3_dp, 1), tol, "6.3")

       call assert_equals(-0.1_dp, table%interpolate(10.0_dp, 1), tol, "10.0")

       call table%destroy()

    end if

  end subroutine test_interpolation_step

!------------------------------------------------------------------------

  subroutine test_average_linear

    ! Linear interpolation average

    type(interpolation_table_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       call table%init(data5)

       call assert_equals(1._dp, table%average([-0.5_dp, -0.1_dp], 1), tol, "[-0.5, -0.1]")

       call assert_equals(1.0238095238095237_dp, table%average([-0.5_dp, 0.1_dp], 1), &
            tol, "[-0.5, 0.1]")

       call assert_equals(1.5_dp, table%average([0.1_dp, 2._dp], 1), tol, "[0.1, 2.]")

       call assert_equals(1.1019345238095237_dp, table%average([0.1_dp, 3._dp], 1), tol, "[0.1, 3.]")

       call assert_equals(0.11586538461538454_dp, table%average([3.1_dp, 7._dp], 1), tol, "[3.1, 7.]")

       call assert_equals(-0.27307692307692316_dp, table%average([8._dp, 12._dp], 1), tol, "[8., 12.]")

       call assert_equals(1.4761904761904763_dp, table%average([1._dp, 1._dp], 1), tol, "[1., 1.]")

       call table%destroy()

    end if

  end subroutine test_average_linear

!------------------------------------------------------------------------

  subroutine test_average_step

    ! Step interpolation average

    type(interpolation_table_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       call table%init(data5, INTERP_STEP)

       call assert_equals(1._dp, table%average([-0.5_dp, -0.1_dp], 1), tol, "[-0.5, -0.1]")

       call assert_equals(1._dp, table%average([-0.5_dp, 0.1_dp], 1), tol, "[-0.5, 0.1]")

       call assert_equals(1._dp, table%average([0.1_dp, 2._dp], 1), tol, "[0.1, 2.]")

       call assert_equals(1.5_dp, table%average([0.1_dp, 3._dp], 1), tol, "[0.1, 3.]")

       call assert_equals(0.45_dp, table%average([3.1_dp, 7._dp], 1), tol, "[3.1, 7.]")

       call assert_equals(-0.6_dp, table%average([8._dp, 12._dp], 1), tol, "[8., 12.]")

       call assert_equals(1._dp, table%average([1._dp, 1._dp], 1), tol, "[1., 1.]")

       call table%destroy()

    end if

  end subroutine test_average_step

!------------------------------------------------------------------------

  subroutine test_average_linear_integration

    ! Linear interpolation integration averaging

    type(interpolation_table_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call table%init(data5, INTERP_LINEAR, INTERP_AVERAGING_INTEGRATE)

       call assert_equals(1._dp, table%average([-0.5_dp, -0.1_dp], 1), tol, "[-0.5, -0.1]")

       call assert_equals(1.003968253968254_dp, table%average([-0.5_dp, 0.1_dp], &
            1), tol, "[-0.5, 0.1]")

       call assert_equals(1.5_dp, table%average([0.1_dp, 2._dp], 1), tol, "[0.1, 2.]")

       call assert_equals(1.5406660509031198_dp, table%average([0.1_dp, 3._dp], 1), tol, "[0.1, 3.]")

       call assert_equals(-0.2530818540433925_dp, table%average([3.1_dp, 7._dp], 1), tol, "[3.1, 7.]")

       call assert_equals(-0.1389423076923077_dp, table%average([8._dp, 12._dp], 1), tol, "[8., 12.]")

       call assert_equals(-0.1_dp, table%average([9._dp, 12._dp], 1), tol, "[9., 12.]")

       call assert_equals(1.4761904761904763_dp, table%average([1._dp, 1._dp], 1), tol, "[1., 1.]")

       call table%destroy()

    end if

  end subroutine test_average_linear_integration

!------------------------------------------------------------------------

  subroutine test_array_interpolator
    ! Array interpolator

    type(array_interpolator_type) :: interp
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscInt, parameter :: size = 3
    PetscReal, parameter :: start(size) = [1._dp, 2._dp, 3._dp], &
         end(size - 1) = [4._dp, 5._dp]
    PetscReal :: v(size), xi
    PetscErrorCode :: err

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       call interp%init(start, end)
       call assert_equals(size, interp%size, "size")
       v = interp%interpolate(0._dp)
       call assert_equals(start, v, size, tol, "xi = 0")
       v = interp%interpolate(1._dp)
       call assert_equals([end, 0._dp], v, size, tol, "xi = 1")
       v = interp%interpolate(0.25_dp)
       call assert_equals([1.75_dp, 2.75_dp, 2.25_dp], v, &
            size, tol, "xi = 0.25")
       call interp%find(2, 3._dp, xi, err)
       call assert_true(err == 0, "find error")
       call assert_equals(1._dp / 3._dp, xi, tol, "find xi")
       call interp%destroy()
    end if

  end subroutine test_array_interpolator

!------------------------------------------------------------------------

  subroutine test_ramp_interpolate
    ! Ramp interpolation

    PetscReal :: x(2), y(2)
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       x = [0._dp, 1._dp]
       y = [0._dp, 1._dp]
       call assert_equals(0.3_dp, &
            ramp_interpolate(0.3_dp, x, y), tol, "case 1a")
       call assert_equals(0.0_dp, &
            ramp_interpolate(-0.3_dp, x, y), tol, "case 1b")
       call assert_equals(1.0_dp, &
            ramp_interpolate(2.3_dp, x, y), tol, "case 1c")

       x = [-1._dp, 2._dp]
       y = [2._dp, -3._dp]
       call assert_equals(1._dp / 3._dp, &
            ramp_interpolate(0._dp, x, y), tol, "case 2a")
       call assert_equals(2._dp, &
            ramp_interpolate(-1.3_dp, x, y), tol, "case 2b")
       call assert_equals(-3._dp, &
            ramp_interpolate(3._dp, x, y), tol, "case 2c")

    end if

  end subroutine test_ramp_interpolate

!------------------------------------------------------------------------

  subroutine test_interpolation_linear_array
    ! Array interpolation

    type(interpolation_table_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call table%init(data3)
       call assert_equals(3, table%dim, "dim")

       call assert_equals([1._dp, 2._dp, 3._dp], table%interpolate(0._dp), 3, tol, "0")
       call assert_equals(0, table%coord%index, "0 index")

       call assert_equals([1.5_dp, 2.5_dp, 3.5_dp], table%interpolate(0.5_dp), 3, tol, "0.5")
       call assert_equals(1, table%coord%index, "0.5 index")

       call assert_equals([2.5_dp, 3.5_dp, 4.5_dp], table%interpolate(1.5_dp), 3, tol, "1.5")
       call assert_equals(2, table%coord%index, "1.5 index")

       call table%destroy()

    end if

  end subroutine test_interpolation_linear_array

!------------------------------------------------------------------------

end module interpolation_test
