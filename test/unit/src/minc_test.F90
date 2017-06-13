module minc_test

  ! Tests for minc module

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscdm.h>

  use petscsys
  use petscdm
  use fruit
  use kinds_module
  use minc_module

  implicit none
  private 

  public :: test_proximity, test_proximity_derivative, &
       test_inner_connection_distance, test_geometry

contains

!------------------------------------------------------------------------

  subroutine test_proximity
    ! MINC proximity functions

    use fson
    use fson_mpi_module

    type(minc_type) :: minc
    DM :: dm
    type(fson_value), pointer :: json
    PetscMPIInt :: rank
    PetscErrorCode :: ierr, err
    PetscReal, parameter :: tol = 1.e-9_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '{"volume_fractions": [0.1, 0.9],' // &
         '"fracture": {"planes": 1, "spacing": 50.}}')
    call minc%init(json, dm, 0, "minc", err = err)
    if (rank == 0) then
       call assert_equals(1, minc%num_fracture_planes, '1 set of fracture planes')
       call assert_equals(0._dp, minc%proximity(0._dp), tol, '1 plane x = 0')
       call assert_equals(0.4_dp, minc%proximity(10._dp), tol, '1 plane x = 10')
       call assert_equals(0.8_dp, minc%proximity(20._dp), tol, '1 plane x = 20')
       call assert_equals(1._dp, minc%proximity(25._dp), tol, '1 plane x = 25')
       call assert_equals(1._dp, minc%proximity(30._dp), tol, '1 plane x = 30')
    end if
    call minc%destroy()
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"volume_fractions": [0.1, 0.9],' // &
         '"fracture": {"planes": 2, "spacing": [50, 80]}}')
    call minc%init(json, dm, 0, "minc", err = err)
    if (rank == 0) then
       call assert_equals(2, minc%num_fracture_planes, '2 sets of fracture planes')
       call assert_equals(0._dp, minc%proximity(0._dp), tol, '2 planes x = 0')
       call assert_equals(0.55_dp, minc%proximity(10._dp), tol, '2 planes x = 10')
       call assert_equals(0.9_dp, minc%proximity(20._dp), tol, '2 planes x = 20')
       call assert_equals(1._dp, minc%proximity(25._dp), tol, '2 planes x = 25')
       call assert_equals(1._dp, minc%proximity(30._dp), tol, '2 planes x = 30')
       call assert_equals(1._dp, minc%proximity(45._dp), tol, '2 planes x = 45')
    end if
    call minc%destroy()
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"volume_fractions": [0.1, 0.9],' // &
         '"fracture": {"planes": 3, "spacing": [50, 80, 60]}}')
    call minc%init(json, dm, 0, "minc", err = err)
    if (rank == 0) then
       call assert_equals(3, minc%num_fracture_planes, '3 sets of fracture planes')
       call assert_equals(0._dp, minc%proximity(0._dp), tol, '3 planes x = 0')
       call assert_equals(0.7_dp, minc%proximity(10._dp), tol, '3 planes x = 10')
       call assert_equals(29._dp / 30._dp, minc%proximity(20._dp), tol, '3 planes x = 20')
       call assert_equals(1._dp, minc%proximity(25._dp), tol, '3 planes x = 25')
       call assert_equals(1._dp, minc%proximity(30._dp), tol, '3 planes x = 30')
       call assert_equals(1._dp, minc%proximity(45._dp), tol, '3 planes x = 45')
    end if
    call minc%destroy()
    call fson_destroy_mpi(json)

  end subroutine test_proximity

!------------------------------------------------------------------------

  subroutine test_proximity_derivative
    ! MINC proximity function derivatives

    use fson
    use fson_mpi_module

    type(minc_type) :: minc
    DM :: dm
    type(fson_value), pointer :: json
    PetscMPIInt :: rank
    PetscErrorCode :: ierr, err
    PetscReal, parameter :: tol = 1.e-9_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '{"volume_fractions": [0.1, 0.9],' // &
         '"fracture": {"planes": 1, "spacing": 50.}}')
    call minc%init(json, dm, 0, "minc", err = err)
    if (rank == 0) then
       call assert_equals(0.04_dp, &
            minc%proximity_derivative(0._dp), tol, '1 plane x = 0')
       call assert_equals(0.04_dp, &
            minc%proximity_derivative(10._dp), tol, '1 plane x = 10')
       call assert_equals(0.04_dp, &
            minc%proximity_derivative(20._dp), tol, '1 plane x = 20')
    end if
    call minc%destroy()
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"volume_fractions": [0.1, 0.9],' // &
         '"fracture": {"planes": 2, "spacing": [50, 80]}}')
    call minc%init(json, dm, 0, "minc", err = err)
    if (rank == 0) then
       call assert_equals(0.065_dp, &
            minc%proximity_derivative(0._dp), tol, '2 planes x = 0')
       call assert_equals(0.045_dp, &
            minc%proximity_derivative(10._dp), tol, '2 planes x = 10')
       call assert_equals(0.025_dp, &
            minc%proximity_derivative(20._dp), tol, '2 planes x = 20')
    end if
    call minc%destroy()
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"volume_fractions": [0.1, 0.9],' // &
         '"fracture": {"planes": 3, "spacing": [50, 80, 60]}}')
    call minc%init(json, dm, 0, "minc", err = err)
    if (rank == 0) then
       call assert_equals(0.0983333333_dp, &
            minc%proximity_derivative(0._dp), tol, '3 planes x = 0')
       call assert_equals(0.045_dp, &
            minc%proximity_derivative(10._dp), tol, '3 planes x = 10')
       call assert_equals(0.0116666667_dp, &
            minc%proximity_derivative(20._dp), tol, '3 planes x = 20')
    end if
    call minc%destroy()
    call fson_destroy_mpi(json)

  end subroutine test_proximity_derivative

!------------------------------------------------------------------------

  subroutine test_inner_connection_distance
    ! MINC inner connection distance

    use fson
    use fson_mpi_module

    type(minc_type) :: minc
    DM :: dm
    type(fson_value), pointer :: json
    PetscMPIInt :: rank
    PetscErrorCode :: ierr, err
    PetscReal, parameter :: tol = 1.e-9_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '{"volume_fractions": [0.1, 0.9],' // &
         '"fracture": {"planes": 1, "spacing": 50.}}')
    call minc%init(json, dm, 0, "minc", err = err)
    if (rank == 0) then
       call assert_equals(25._dp / 3._dp, &
            minc%inner_connection_distance(0._dp), tol, '1 plane x = 0')
       call assert_equals(5._dp, &
            minc%inner_connection_distance(10._dp), tol, '1 plane x = 10')
       call assert_equals(5._dp / 3._dp, &
            minc%inner_connection_distance(20._dp), tol, '1 plane x = 20')
    end if
    call minc%destroy()
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"volume_fractions": [0.1, 0.9],' // &
         '"fracture": {"planes": 2, "spacing": [50, 80]}}')
    call minc%init(json, dm, 0, "minc", err = err)
    if (rank == 0) then
       call assert_equals(100._dp / 13._dp, &
            minc%inner_connection_distance(0._dp), tol, '2 planes x = 0')
       call assert_equals(5._dp, &
            minc%inner_connection_distance(10._dp), tol, '2 planes x = 10')
       call assert_equals(2._dp, &
            minc%inner_connection_distance(20._dp), tol, '2 planes x = 20')
    end if
    call minc%destroy()
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"volume_fractions": [0.1, 0.9],' // &
         '"fracture": {"planes": 3, "spacing": [50, 80, 60]}}')
    call minc%init(json, dm, 0, "minc", err = err)
    if (rank == 0) then
       call assert_equals(360._dp / 59._dp, &
            minc%inner_connection_distance(0._dp), tol, '3 planes x = 0')
       call assert_equals(4._dp, &
            minc%inner_connection_distance(10._dp), tol, '3 planes x = 10')
       call assert_equals(12._dp / 7._dp, &
            minc%inner_connection_distance(20._dp), tol, '3 planes x = 20')
    end if
    call minc%destroy()
    call fson_destroy_mpi(json)

  end subroutine test_inner_connection_distance

!------------------------------------------------------------------------

  subroutine test_geometry
    ! MINC geometry

    use fson
    use fson_mpi_module

    type(minc_type) :: minc
    DM :: dm
    type(fson_value), pointer :: json
    PetscMPIInt :: rank
    PetscErrorCode :: ierr, err
    PetscReal, parameter :: tol = 1.e-9_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '{"volume_fractions": [0.1, 0.9],' // &
         '"fracture": {"planes": 1, "spacing": 50.}}')
    call minc%init(json, dm, 0, "minc", err = err)
    if (rank == 0) then
       call assert_equals(0, err, '1 plane 1 level error')
       call assert_equals([0.1_dp, 0.9_dp], minc%volume, 2, tol, &
            '1 plane 1 level volume fractions')
       call assert_equals([0.036_dp, 0.036_dp], &
            minc%connection_area, 2, tol, '1 plane 1 level connection areas')
       call assert_equals([0._dp, 25._dp / 3._dp], &
            minc%connection_distance, 2, tol, '1 plane 1 level connection distances')
    end if
    call minc%destroy()
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"volume_fractions": [0.1, 0.3, 0.6],' // &
         '"fracture": {"planes": 1, "spacing": 100.}}')
    call minc%init(json, dm, 0, "minc", err = err)
    if (rank == 0) then
       call assert_equals(0, err, '1 plane 2 levels error')
       call assert_equals([0.1_dp, 0.3_dp, 0.6_dp], minc%volume, 3, tol, &
            '1 plane 2 levels volume fractions')
       call assert_equals([0.018_dp, 0.018_dp, 0.018_dp], &
            minc%connection_area, 3, tol, '1 plane 2 levels connection areas')
       call assert_equals([0._dp, 25._dp / 3._dp, 100._dp / 9._dp], &
            minc%connection_distance, 3, tol, '1 plane 2 levels connection distances')
    end if
    call minc%destroy()
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"volume_fractions": [10, 20, 30, 40],' // &
         '"fracture": {"planes": 1, "spacing": 100.}}')
    call minc%init(json, dm, 0, "minc", err = err)
    if (rank == 0) then
       call assert_equals(0, err, '1 plane 3 levels error')
       call assert_equals([0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp], minc%volume, 4, tol, &
            '1 plane 3 levels volume fractions')
       call assert_equals([0.018_dp, 0.018_dp, 0.018_dp, 0.018_dp], &
            minc%connection_area, 4, tol, '1 plane 3 levels connection areas')
       call assert_equals([0._dp, 50._dp / 9._dp, 25._dp / 3._dp, 7400._dp / 999._dp], &
            minc%connection_distance, 4, tol, '1 plane 3 levels connection distances')
    end if
    call minc%destroy()
    call fson_destroy_mpi(json)

  end subroutine test_geometry

!------------------------------------------------------------------------

end module minc_test
