module face_test

  ! Test for face module

  use kinds_module
  use mpi_module
  use fruit
  use face_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>

public :: test_face_assign, test_face_permeability_direction

PetscReal, parameter :: tol = 1.e-6_dp

contains
  
!------------------------------------------------------------------------

  subroutine test_face_assign

    ! Face assign() test

    type(face_type) :: face
    PetscReal, parameter :: area = 300._dp
    PetscReal, parameter :: distance(2) = [20._dp, 30._dp]
    PetscReal, parameter :: normal(3) = [0.5_dp, -0.25_dp, 0.75_dp]
    PetscReal, parameter :: centroid(3) = [-1250._dp, 3560._dp, -2530._dp]
    PetscReal, parameter :: permeability_direction = dble(2)
    PetscInt, parameter :: offset = 6
    PetscReal :: offset_padding(offset-1) = 0._dp
    PetscReal, allocatable :: face_data(:)

    if (mpi%rank == mpi%output_rank) then

       face_data = [offset_padding, area, distance, normal, centroid, &
            permeability_direction]

       call assert_equals(face%dof(), size(face_data) - (offset-1), "face dof")

       call face%assign(face_data, offset)

       call assert_equals(area, face%area, tol, "area")
       call assert_equals(0._dp, norm2(face%distance - distance), tol, "distances")
       call assert_equals(0._dp, norm2(face%normal - normal), tol, "normal")
       call assert_equals(0._dp, norm2(face%centroid - centroid), tol, "centroid")
       call assert_equals(permeability_direction, face%permeability_direction, &
            tol, "permeability direction")

       call face%destroy()
       deallocate(face_data)

    end if

  end subroutine test_face_assign

!------------------------------------------------------------------------

  subroutine test_face_permeability_direction

    ! Face permeability_direction() test

    type(face_type) :: face
    PetscReal, parameter :: area = 10._dp
    PetscReal, parameter :: distance(2) = [10._dp, 10._dp]
    PetscReal, parameter :: centroid(3) = [0._dp, 0._dp, 0._dp]
    PetscReal, parameter :: initial_permeability_direction = dble(0)
    PetscInt,  parameter :: num_tests = 3
    PetscReal, parameter :: normal(3, num_tests) = reshape( &
         [  1._dp, -0.25_dp,  0.3_dp, &
          -0.1_dp,   2.1_dp,  0.5_dp,&
            1._dp,  -1.4_dp, -1.6_dp], [3, num_tests])
    PetscReal, parameter :: expected_permeability_direction(num_tests) = &
         [dble(1), dble(2), dble(3)]
    PetscReal, allocatable :: face_data(:)
    PetscInt :: i
    PetscInt :: offset = 1
    character(len = 32) :: msg

    if (mpi%rank == mpi%output_rank) then

       do i = 1, num_tests

          face_data = [area, distance, normal(:,i), centroid, &
               initial_permeability_direction]
          call face%assign(face_data, offset)
          call face%calculate_permeability_direction()

          write(msg, '(a, i2)') "Permeability direction test ", i
          call assert_equals(expected_permeability_direction(i), &
               face%permeability_direction, tol, trim(msg))

          call face%destroy()
          deallocate(face_data)

       end do

    end if

  end subroutine test_face_permeability_direction

!------------------------------------------------------------------------

end module face_test
