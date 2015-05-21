module face_test

  ! Test for face module

  use kinds_module
  use mpi_module
  use fruit
  use face_module

  implicit none
  private

#include <petsc-finclude/petscdef.h>

public :: test_face

PetscReal, parameter :: tol = 1.e-6_dp

contains
  
!------------------------------------------------------------------------

  subroutine test_face

    ! Face test

    type(face_type) :: face
    PetscReal, parameter :: area = 300._dp
    PetscReal, parameter :: distance(2) = [20._dp, 30._dp]
    PetscReal, parameter :: normal(3) = [0.5_dp, -0.25_dp, 0.75_dp]
    PetscReal, parameter :: centroid(3) = [-1250._dp, 3560._dp, -2530._dp]
    PetscInt, parameter :: offset = 6
    PetscReal :: offset_padding(offset-1) = 0._dp
    PetscReal, allocatable :: face_data(:)

    if (mpi%rank == mpi%output_rank) then

       face_data = [offset_padding, area, distance, normal, centroid]

       call assert_equals(face%dof(), size(face_data) - (offset-1), "face dof")

       call face%assign(face_data, offset)

       call assert_equals(area, face%area, tol, "area")
       call assert_equals(0._dp, norm2(face%distance - distance), tol, "distances")
       call assert_equals(0._dp, norm2(face%normal - normal), tol, "normal")
       call assert_equals(0._dp, norm2(face%centroid - centroid), tol, "centroid")

       call face%destroy()
       deallocate(face_data)

    end if

  end subroutine test_face

!------------------------------------------------------------------------

end module face_test
