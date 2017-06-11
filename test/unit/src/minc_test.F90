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

  public :: test_proximity

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
    PetscErrorCode :: ierr
    PetscReal, parameter :: tol = 1.e-9_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    call DMPlexCreateFromFile(PETSC_COMM_WORLD, "data/mesh/block3.exo", &
         PETSC_TRUE, dm, ierr); CHKERRQ(ierr)

    json => fson_parse_mpi(str = '{"volume_fractions": [0.1, 0.9],' // &
         '"fracture": {"planes": 1, "spacing": 50.}}')
    call minc%init(json, dm, 0, "minc")
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
    call minc%init(json, dm, 0, "minc")
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
    call minc%init(json, dm, 0, "minc")
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

    call DMDestroy(dm, ierr)

  end subroutine test_proximity

!------------------------------------------------------------------------

end module minc_test
