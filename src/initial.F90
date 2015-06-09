module initial_module
  !! Module for setting up initial conditions.

  use kinds_module
  use fson
  use fson_mpi_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  public :: setup_initial

contains

!------------------------------------------------------------------------

  subroutine setup_initial(json, t, v)
    !! Initializes time t and a Vec v with initial conditions read
    !! from JSON input 'initial'.  Conditions may be specified as a
    !! constant value or as an array. The array may contain a complete
    !! set of initial conditions for all cells, or if a shorter array is
    !! given, this is repeated over initial conditions vector.

    use fson_value_m, only : TYPE_REAL, TYPE_INTEGER, TYPE_ARRAY

    type(fson_value), pointer, intent(in) :: json
    PetscReal, intent(out) :: t
    Vec, intent(out) :: v
    ! Locals:
    PetscErrorCode :: ierr
    PetscReal :: const_initial_value
    PetscInt :: int_const_initial_value
    PetscInt, allocatable :: indices(:)
    PetscReal, allocatable :: initial_input(:), initial_data(:)
    PetscInt :: i, np, count
    PetscBool :: const
    PetscReal, parameter :: default_start_time = 0.0_dp
    PetscReal, parameter :: default_initial_value = 0.0_dp

    call fson_get_mpi(json, "time.start", default_start_time, t)

    call VecGetSize(v, count, ierr); CHKERRQ(ierr)
    const = .true.

    if (fson_has_mpi(json, "initial")) then

       select case (fson_type_mpi(json, "initial"))
          case (TYPE_REAL)
             call fson_get_mpi(json, "initial", val = const_initial_value)
          case (TYPE_INTEGER)
             call fson_get_mpi(json, "initial", val = int_const_initial_value)
             const_initial_value = real(int_const_initial_value)
          case (TYPE_ARRAY)
             const = .false.
             call fson_get_mpi(json, "initial", val = initial_input)
             np = size(initial_input)
             if (np >= count) then
                initial_data = initial_input(1:count)
             else ! repeat input over array:
                do i = 1, np
                   initial_data(i:count:np) = initial_input(i)
                end do
             end if
             deallocate(initial_input)
       end select
    else
       const_initial_value = default_initial_value
    end if

    if (const) then
       call VecSet(v, const_initial_value, ierr); CHKERRQ(ierr)
    else
       allocate(indices(count))
       do i = 1, count
          indices(i) = i-1
       end do
       call VecSetValues(v, count, indices, &
            initial_data, INSERT_VALUES, ierr); CHKERRQ(ierr)
       call VecAssemblyBegin(v, ierr); CHKERRQ(ierr)
       call VecAssemblyEnd(v, ierr); CHKERRQ(ierr)
       deallocate(indices, initial_data)
    end if

  end subroutine setup_initial

!------------------------------------------------------------------------

end module
