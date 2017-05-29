module list_test

  ! Test for list module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use fruit
  use list_module

  implicit none
  private

  type :: thing_type
     character(len = 5) :: name
     PetscReal, allocatable :: x(:)
   contains
     procedure :: destroy => thing_destroy
  end type thing_type

public :: test_list

contains

!------------------------------------------------------------------------

  subroutine thing_destroy(self)
    !! Destroys a thing.

    class(thing_type), intent(in out) :: self

    deallocate(self%x)

  end subroutine thing_destroy

!------------------------------------------------------------------------

  subroutine test_list
    ! Test list operations
    
    type(list_type) :: list, list2
    PetscInt :: i, j
    PetscReal :: x, y
    character(3) :: str, str2
    type(thing_type) :: thing
    type(list_node_type), pointer :: node
    PetscReal, parameter :: tol = 1.e-8_dp
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       i = 1
       x = -5.67_dp
       str = 'bob'
       thing%name = 'clive'
       allocate(thing%x(10))

       call list%init()
       call list%append(i, 'int')
       call list%prepend(x, 'real')
       call list%append(str)
       call list%append(thing)

       call assert_equals(4, list%count, 'initial list count')
       call assert_true(associated(list%head), 'initial head associated')
       call assert_true(associated(list%tail), 'initial tail associated')

       call list%delete('int')
       call assert_equals(3, list%count, 'list count after deletion')
       i = 3
       call list%append(i, 'number')
       call assert_equals(4, list%count, 'list count after append')

       call list%traverse(list_node_adjust, .true.)
       call assert_equals(4, list%count, 'list count after adjust')

       node => list%find('number')
       call assert_true(associated(node), 'find tag')
       if (associated(node)) then
          select type (d => node%data)
          type is (PetscInt)
             call assert_equals(4, d, 'tagged value')
          end select
       end if

       node => list%find('brian')
       call assert_false(associated(node), 'find missing tag')

       node => list%get(0)
       call assert_true(associated(node), 'get(0)')
       if (associated(node)) then
          select type (d => node%data)
          type is (PetscReal)
             call assert_equals(-8.505_dp, d, tol, 'get(0) value')
          end select
       end if

       i = 2
       node => list%get(-1)
       call assert_true(associated(node), 'get(-1)')
       if (associated(node)) then
          select type (d => node%data)
          type is (PetscInt)
             call assert_equals(i, d, 'get(-1) value')
          end select
       end if

       node => list%get(-2)
       call assert_true(associated(node), 'get(-2)')
       if (associated(node)) then
          select type (d => node%data)
          type is (thing_type)
             call assert_equals('boris', d%name, 'get(-2) value')
          end select
       end if

       node => list%get(10)
       call assert_false(associated(node), 'get(10)')
       node => list%get(-8)
       call assert_false(associated(node), 'get(-8)')

       j = 10
       str2 = 'foo'
       y = 101._dp
       call list2%init()
       call list2%append(j, 'int2')
       call list2%append(str2)
       call list2%append(y, 'real2')

       call list%add(list2)
       call list2%destroy()
       call assert_equals(7, list%count, 'list count after add')

       node => list%get(-1)
       call assert_true(associated(node), 'get(-1) after add')
       if (associated(node)) then
          select type (d => node%data)
          type is (PetscReal)
             call assert_equals(y, d, tol, 'get(-1) value after add')
          end select
       end if

       call list%destroy(list_node_destroy_data)
       call assert_equals(0, list%count, 'list count after destroy')
       call assert_false(associated(list%head), 'head associated after destroy')
       call assert_false(associated(list%tail), 'tail associated after destroy')

    end if

  contains

    subroutine list_node_adjust(node, stopped)
      ! Iterator function to apply to nodes.
      implicit none
      type(list_node_type), pointer, intent(in out) :: node
      logical, intent(out) :: stopped

      select type(d => node%data)
      type is (PetscInt)
         d = d + 1
      type is (PetscReal)
         d = d * 1.5_dp
      type is (thing_type)
         d%name = 'boris'
      end select

      stopped = .false.

    end subroutine list_node_adjust

    subroutine list_node_destroy_data(node)
      ! Destroy data in list node.
      implicit none
      type(list_node_type), pointer, intent(in out) :: node

      select type(d => node%data)
      type is (thing_type)
         call d%destroy()
      end select

    end subroutine list_node_destroy_data

  end subroutine test_list

!------------------------------------------------------------------------

end module list_test

