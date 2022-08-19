module dictionary_test

  ! Tests for dictionaries.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use zofu
  use dictionary_module
  use list_module, only: list_node_type, list_type

  implicit none
  private

  type :: thing_type
     PetscReal, allocatable :: x(:)
   contains
     procedure :: destroy => thing_destroy
  end type thing_type

  public :: setup, teardown
  public :: test_integer, test_thing, test_dict_list, test_dict_no_data, &
       test_long_key

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

  subroutine thing_destroy(self)
    !! Destroys a thing.

    class(thing_type), intent(in out) :: self

    deallocate(self%x)

  end subroutine thing_destroy

!------------------------------------------------------------------------

  subroutine test_integer(test)
    ! dictionary of integers

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(dictionary_type) :: dict
    PetscInt :: i, count
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    call dict%init()

    do i = ichar("a"), ichar("z")
       call dict%add(char(i), i)
    end do
    call dict%add("fred", 100)
    call dict%add("egeszegedre", 101)
    call dict%add("zone100", 102)
    call dict%add("ZZZZZZ", 103)
    
    if (rank == 0) then
       call test%assert(dict%has("fred"), 'fred present')
    end if
    
    count = dict%count()
    call dict%delete("fred")
    if (rank == 0) then
       call test%assert(.not. dict%has("fred"), 'fred gone')
       call test%assert(count - 1, dict%count(), 'count with fred gone')
    end if

    do i = ichar("a"), ichar("z")
       call dict_test(char(i), i)
    end do
    call dict_test("egeszegedre", 101)
    call dict_test("zone100", 102)
    call dict_test("ZZZZZZ", 103)

    call dict%destroy()

  contains

    subroutine dict_test(key, expected)

      character(*), intent(in) :: key
      PetscInt, intent(in) :: expected
      ! Locals:
      type(list_node_type), pointer :: node

      node => dict%get(key)
      if (rank == 0) then
         select type (data => node%data)
         type is (PetscInt)
            call test%assert(expected, data, key)
         end select
      end if

    end subroutine dict_test

  end subroutine test_integer

!------------------------------------------------------------------------

  subroutine test_thing(test)
    ! dictionary of objects

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(thing_type), pointer :: thing
    type(dictionary_type) :: dict
    type(list_node_type), pointer :: node
    
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    call dict%init(owner = PETSC_TRUE)

    allocate(thing_type :: thing)
    thing%x = [1._dp, 2._dp, 3._dp]

    call dict%add("thing1", thing)

    if (rank == 0) then
       call test%assert(dict%has("thing1"), "has thing1")
    end if

    node => dict%get("another thing")
    if (rank == 0) then
       call test%assert(.not. associated(node), "get thing2")
    end if

    node => dict%get("thing1")
    select type (data => node%data)
    type is (thing_type)
       if (rank == 0) then
          call test%assert([1._dp, 2._dp, 3._dp], data%x, "thing x")
       end if
    end select

    call dict%destroy(list_node_destroy_data)

  contains

    subroutine list_node_destroy_data(node)
      ! Destroy data in list node.
      implicit none
      type(list_node_type), pointer, intent(in out) :: node

      select type(d => node%data)
      type is (thing_type)
         call d%destroy()
      end select

    end subroutine list_node_destroy_data

  end subroutine test_thing

!------------------------------------------------------------------------

  subroutine test_dict_list(test)
    ! list to dict

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(list_type) :: list
    type(dictionary_type) :: dict
    PetscInt :: i, j
    PetscReal :: x
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    i = 5
    j = 9
    x = -5.67_dp

    call list%init()
    call list%append(i, "i")
    call list%append(j, "j")
    call list%append(x, "x")

    call dict%init(list)

    if (rank == 0) then
       call test%assert(3, dict%count(), "count")
       call test%assert(dict%has("i"), "has i")
       call test%assert(dict%has("j"), "has j")
       call test%assert(dict%has("x"), "has x")
    end if

    call dict%destroy()
    call list%destroy()

  end subroutine test_dict_list

!------------------------------------------------------------------------

  subroutine test_dict_no_data(test)
    ! dict with no data

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(dictionary_type) :: dict
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    call dict%init()

    call dict%add("foo")
    call dict%add("bar")

    call test%assert(2, dict%count(), "count")

    call test%assert(dict%has("foo"), "has foo")
    call test%assert(dict%has("bar"), "has bar")
    call test%assert(.not. dict%has("bat"), "has no bat")

    call dict%destroy()

  end subroutine test_dict_no_data

!------------------------------------------------------------------------

  subroutine test_long_key(test)
    ! dict with long keys

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(dictionary_type) :: dict
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    call dict%init()

    call dict%add("this is a rather long name 1", 101)
    call dict%add("this is a rather long name 2", 102)

    call test%assert(2, dict%count(), "count")

    call test%assert(dict%has("this is a rather long name 1"), "has 1")
    call test%assert(dict%has("this is a rather long name 2"), "has 2")
    call test%assert(.not. dict%has("bar"), "has no bar")

    call dict_test("this is a rather long name 1", 101)
    call dict_test("this is a rather long name 2", 102)

    call dict%destroy()

  contains

    subroutine dict_test(key, expected)

      character(*), intent(in) :: key
      PetscInt, intent(in) :: expected
      ! Locals:
      type(list_node_type), pointer :: node

      node => dict%get(key)
      if (rank == 0) then
         select type (data => node%data)
         type is (PetscInt)
            call test%assert(expected, data, key)
         end select
      end if

    end subroutine dict_test

  end subroutine test_long_key

!------------------------------------------------------------------------

end module dictionary_test
