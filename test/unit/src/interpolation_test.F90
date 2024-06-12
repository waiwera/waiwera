module interpolation_test

  ! Tests for interpolation module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use zofu
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

  public :: setup, teardown, setup_test
  public :: test_interpolation_linear, test_interpolation_single, &
       test_interpolation_step, test_interpolation_pchip, &
       test_average_linear, test_average_step, &
       test_average_linear_integration, test_average_step_integration, &
       test_interpolation_linear_array, test_find, test_unsorted, &
       test_duplicate

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

  subroutine setup_test(test)

    class(unit_test_type), intent(in out) :: test

    test%tolerance = 1.e-9

  end subroutine setup_test

!------------------------------------------------------------------------

  subroutine test_interpolation_linear(test)

    ! Linear interpolation

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(interpolation_table_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call table%init(data5) ! default linear interpolation

       call test%assert(1, table%dim, "dim")
       call test%assert(table%continuous, "continuous")

       call test%assert(1._dp, table%interpolate(-0.5_dp, 1), "-0.5")
       call test%assert(0, table%coord%index, "-0.5 index")

       call test%assert(1._dp, table%interpolate(0.0_dp, 1), "0.0")
       call test%assert(0, table%coord%index, "0.0 index")

       call test%assert(1.4761904761904763_dp, table%interpolate(1.0_dp, 1), "1.0")
       call test%assert(1, table%coord%index, "1.0 index")

       call test%assert(0.007692307692307665_dp, table%interpolate(4.5_dp, 1), "4.5")
       call test%assert(3, table%coord%index, "4.5 index")

       call test%assert(0.59375_dp, table%interpolate(3.6_dp, 1), "3.6")
       call test%assert(2, table%coord%index, "3.6 index")

       call test%assert(-1.1_dp, table%interpolate(6.3_dp, 1), "6.3")
       call test%assert(4, table%coord%index, "6.3 index")

       call test%assert(-0.1_dp, table%interpolate(10.0_dp, 1), "10.0")
       call test%assert(5, table%coord%index, "10.0 index")

       call table%destroy()

    end if

  end subroutine test_interpolation_linear

!------------------------------------------------------------------------

  subroutine test_interpolation_single(test)

    ! Single-value data table

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(interpolation_table_type) :: table
    PetscReal, dimension(1,2), parameter :: data1 = reshape(&
         [0._dp, 2._dp], [1,2])
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call table%init(data1)

       call test%assert(2._dp, table%interpolate(-0.5_dp, 1), "-0.5")

       call test%assert(2._dp, table%interpolate(0.0_dp, 1), "0.0")

       call test%assert(2._dp, table%interpolate(1.1_dp, 1), "1.1")

       call table%destroy()
    end if

  end subroutine test_interpolation_single

!------------------------------------------------------------------------

  subroutine test_interpolation_step(test)

    ! Step interpolation

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(interpolation_table_step_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call table%init(data5)

       call test%assert(1._dp, table%interpolate(-0.5_dp, 1), "-0.5")
       call test%assert(.not. table%continuous, "continuous")

       call test%assert(1._dp, table%interpolate(0.0_dp, 1), "0.0")

       call test%assert(1._dp, table%interpolate(1.0_dp, 1), "1.0")

       call test%assert(0.5_dp, table%interpolate(4.5_dp, 1), "4.5")

       call test%assert(2.0_dp, table%interpolate(3.6_dp, 1), "3.6")

       call test%assert(-1.1_dp, table%interpolate(6.3_dp, 1), "6.3")

       call test%assert(-0.1_dp, table%interpolate(10.0_dp, 1), "10.0")

       call table%destroy()

    end if

  end subroutine test_interpolation_step

!------------------------------------------------------------------------

  subroutine test_interpolation_pchip(test)

    ! PCHIP interpolation

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(interpolation_table_pchip_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call table%init(data5)

       call test%assert(1._dp, table%interpolate(-0.5_dp, 1), "-0.5")

       call test%assert(1._dp, table%interpolate(0.0_dp, 1), "0.0")

       call test%assert(1.8151181429242655_dp, table%interpolate(1.0_dp, 1), "1.0")

       call test%assert(-0.15089163120635768_dp, table%interpolate(4.5_dp, 1), "4.5")

       call test%assert(0.5832445364416781_dp, table%interpolate(3.6_dp, 1), "3.6")

       call test%assert(-1.1_dp, table%interpolate(6.3_dp, 1), "6.3")

       call test%assert(-0.1_dp, table%interpolate(10.0_dp, 1), "10.0")

       call table%destroy()

    end if

  end subroutine test_interpolation_pchip

!------------------------------------------------------------------------

  subroutine test_average_linear(test)

    ! Linear interpolation average

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(interpolation_table_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       call table%init(data5)

       call test%assert(1._dp, table%average([-0.5_dp, -0.1_dp], 1), "[-0.5, -0.1]")

       call test%assert(1.0238095238095237_dp, table%average([-0.5_dp, 0.1_dp], 1), &
            "[-0.5, 0.1]")

       call test%assert(1.5_dp, table%average([0.1_dp, 2._dp], 1), "[0.1, 2.]")

       call test%assert(1.1019345238095237_dp, table%average([0.1_dp, 3._dp], 1), "[0.1, 3.]")

       call test%assert(0.11586538461538454_dp, table%average([3.1_dp, 7._dp], 1), "[3.1, 7.]")

       call test%assert(-0.27307692307692316_dp, table%average([8._dp, 12._dp], 1), "[8., 12.]")

       call test%assert(1.4761904761904763_dp, table%average([1._dp, 1._dp], 1), "[1., 1.]")

       call table%destroy()

    end if

  end subroutine test_average_linear

!------------------------------------------------------------------------

  subroutine test_average_step(test)

    ! Step interpolation average

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(interpolation_table_step_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       call table%init(data5)

       call test%assert(1._dp, table%average([-0.5_dp, -0.1_dp], 1), "[-0.5, -0.1]")

       call test%assert(1._dp, table%average([-0.5_dp, 0.1_dp], 1), "[-0.5, 0.1]")

       call test%assert(1._dp, table%average([0.1_dp, 2._dp], 1), "[0.1, 2.]")

       call test%assert(1.5_dp, table%average([0.1_dp, 3._dp], 1), "[0.1, 3.]")

       call test%assert(0.45_dp, table%average([3.1_dp, 7._dp], 1), "[3.1, 7.]")

       call test%assert(-0.6_dp, table%average([8._dp, 12._dp], 1), "[8., 12.]")

       call test%assert(1._dp, table%average([1._dp, 1._dp], 1), "[1., 1.]")

       call table%destroy()

    end if

  end subroutine test_average_step

!------------------------------------------------------------------------

  subroutine test_average_linear_integration(test)

    ! Linear interpolation integration averaging

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(interpolation_table_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call table%init(data5, INTERP_AVERAGING_INTEGRATE)

       call test%assert(1._dp, table%average([-0.5_dp, -0.1_dp], 1), "[-0.5, -0.1]")

       call test%assert(1.003968253968254_dp, table%average([-0.5_dp, 0.1_dp], &
            1), "[-0.5, 0.1]")

       call test%assert(1.5_dp, table%average([0.1_dp, 2._dp], 1), "[0.1, 2.]")

       call test%assert(1.5406660509031198_dp, table%average([0.1_dp, 3._dp], 1), "[0.1, 3.]")

       call test%assert(-0.2530818540433925_dp, table%average([3.1_dp, 7._dp], 1), "[3.1, 7.]")

       call test%assert(-0.1389423076923077_dp, table%average([8._dp, 12._dp], 1), "[8., 12.]")

       call test%assert(-0.1_dp, table%average([9._dp, 12._dp], 1), "[9., 12.]")

       call test%assert(1.4761904761904763_dp, table%average([1._dp, 1._dp], 1), "[1., 1.]")

       call table%destroy()

    end if

  end subroutine test_average_linear_integration

!------------------------------------------------------------------------

  subroutine test_average_step_integration(test)
    ! Step interpolation integration averaging

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(interpolation_table_step_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    if (rank == 0) then

       call table%init(data5, INTERP_AVERAGING_INTEGRATE)

       call test%assert(1._dp, table%average([-0.5_dp, -0.1_dp], 1), "[-0.5, -0.1]")

       call test%assert(1._dp, table%average([-0.5_dp, 0.1_dp], 1), "[-0.5, 0.1]")

       call test%assert(1._dp, table%average([0.1_dp, 2._dp], 1), "[0.1, 2.]")

       call test%assert(3.8_dp / 2.9_dp, table%average([0.1_dp, 3._dp], 1), "[0.1, 3.]")

       call test%assert(1.73_dp / 3.9_dp, table%average([3.1_dp, 7._dp], 1), "[3.1, 7.]")

       call test%assert(-0.325_dp, table%average([8._dp, 12._dp], 1), "[8., 12.]")

       call test%assert(1._dp, table%average([1._dp, 1._dp], 1), "[1., 1.]")

       call table%destroy()

    end if

  end subroutine test_average_step_integration

!------------------------------------------------------------------------

  subroutine test_interpolation_linear_array(test)
    ! Array interpolation

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(interpolation_table_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call table%init(data3)
       call test%assert(3, table%dim, "dim")

       call test%assert([1._dp, 2._dp, 3._dp], table%interpolate(0._dp), "0")
       call test%assert(0, table%coord%index, "0 index")

       call test%assert([1.5_dp, 2.5_dp, 3.5_dp], table%interpolate(0.5_dp), "0.5")
       call test%assert(1, table%coord%index, "0.5 index")

       call test%assert([2.5_dp, 3.5_dp, 4.5_dp], table%interpolate(1.5_dp), "1.5")
       call test%assert(2, table%coord%index, "1.5 index")

       call table%destroy()

    end if

  end subroutine test_interpolation_linear_array

!------------------------------------------------------------------------

  subroutine test_find(test)
    ! Find_component_at_index()

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(interpolation_table_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal :: x
    PetscErrorCode :: err
    PetscReal, dimension(2, 2), parameter :: data_const = reshape([&
         0._dp, 2.0_dp, &
         3._dp, 3.0_dp], &
         [2,2])

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call table%init(data3)

       call table%find_component_at_index(1.5_dp, 1, x, err)
       call test%assert(0.5_dp, x, "component 1 = 1.5")
       call test%assert(0, err, "component 1 = 1.5 error")

       call table%find_component_at_index(3.75_dp, 3, x, err)
       call test%assert(0.75_dp, x, "component 3 = 3.75")
       call test%assert(0, err, "component 3 = 3.75 error")

       call table%find_component_at_index(0._dp, 1, x, err)
       call test%assert(0, err, "component 1 = 0 error")
       
       call table%destroy()

       call table%init(data_const)
       call table%find_component_at_index(3._dp, 1, x, err)
       call test%assert(1, err, "constant data error")
       call table%destroy()

    end if
    

  end subroutine test_find

!------------------------------------------------------------------------

  subroutine test_unsorted(test)
    ! Unsorted table

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(interpolation_table_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, dimension(3,2), parameter :: data = reshape([&
         0._dp, 2.1_dp, 2.0_dp, &
         1._dp, 2.0_dp, 0.5_dp], &
         [3,2])

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call table%init(data)
       call test%assert(0.75_dp, table%interpolate(1._dp, 1), "1")
       call test%assert(1, table%coord%index, "1 index")
       call test%assert(1.25_dp, table%interpolate(2.05_dp, 1), "2.05")
       call test%assert(2, table%coord%index, "2.05 index")
       call test%assert(2._dp, table%interpolate(3._dp, 1), "3")
       call test%assert(3, table%coord%index, "3 index")
       call table%destroy()

    end if

  end subroutine test_unsorted

!------------------------------------------------------------------------

  subroutine test_duplicate(test)
    ! Duplicate coordinates

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(interpolation_table_type) :: table
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, dimension(4,2), parameter :: data = reshape([&
         0._dp, 1.0_dp, 1.0_dp, 2.0_dp, &
         1._dp, 2.0_dp, 0.5_dp, 0.1_dp], &
         [4,2])

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call table%init(data)
       call test%assert(0.3_dp, table%interpolate(1.5_dp, 1), "1.5")
       call test%assert(3, table%coord%index, "1.5 index")
       call test%assert(0.5_dp, table%interpolate(1.0_dp, 1), "1.0")
       call test%assert(3, table%coord%index, "1.0 index")
       call table%destroy()

    end if

  end subroutine test_duplicate

!------------------------------------------------------------------------

end module interpolation_test
