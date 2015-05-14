module initial_module
  !! Module for setting up initial conditions.

  use mpi_module
  use kinds_module
  use fson
  use fson_mpi_module

  implicit none
  private

#include <petsc-finclude/petsc.h90>

  public :: setup_initial

contains

!------------------------------------------------------------------------

  subroutine setup_initial(json, dm, initial)
    !! Reads initial conditions from JSON input and creates a Vec
    !! 'initial'.  Conditions may be specified as a constant value or
    !! as an array. The array may contain a complete of initial
    !! conditions for all cells, or if a shorter array is given, this
    !! is repeated over initial conditions vector.

    use fson_value_m, only : TYPE_REAL, TYPE_INTEGER, TYPE_ARRAY

    type(fson_value), pointer, intent(in) :: json
    DM, intent(in) :: dm
    Vec, intent(out) :: initial
    ! Locals:
    PetscErrorCode :: ierr
    real(dp) :: const_initial_value
    integer :: int_const_initial_value
    integer, allocatable :: indices(:)
    real(dp), allocatable :: initial_input(:), initial_data(:)
    integer :: i, np, count
    logical :: const
    real(dp), parameter :: default_initial_value = 0.0_dp

    call DMCreateGlobalVector(dm, initial, ierr); CHKERRQ(ierr)
    call VecGetSize(initial, count, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(initial, "initial", ierr); CHKERRQ(ierr)

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
       call VecSet(initial, const_initial_value, ierr); CHKERRQ(ierr)
    else
       allocate(indices(count))
       do i = 1, count
          indices(i) = i-1
       end do
       call VecSetValues(initial, count, indices, &
            initial_data, INSERT_VALUES, ierr); CHKERRQ(ierr)
       call VecAssemblyBegin(initial, ierr); CHKERRQ(ierr)
       call VecAssemblyEnd(initial, ierr); CHKERRQ(ierr)
       deallocate(indices, initial_data)
    end if

  end subroutine setup_initial

!------------------------------------------------------------------------

end module
