module list_module
  !! Module for linked lists with arbitrary data and optional character tags.

  implicit none
  private

#include <petsc/finclude/petscdef.h>
#include <petsc/finclude/petscsys.h>

  type, public :: list_node_type
     !! List node type.
     private
     character(:), allocatable, public :: tag !! Optional character tag
     class(*), pointer, public :: data => null() !! Node data
     type(list_node_type), pointer :: next => null() !! Pointer to next node
     type(list_node_type), pointer :: previous => null() !! Pointer to previous node
  end type list_node_type

  type, public :: list_type
     !! Linked list type.
     private
     PetscBool :: delete_deallocates !! Whether deleting a node also deallocates its data
     PetscInt, public :: count  !! Number of nodes in the list
     type(list_node_type), pointer, public :: head !! Node at start of list
     type(list_node_type), pointer, public :: tail !! Node at end of list
   contains
     private
     procedure :: list_init
     procedure :: list_init_default
     generic, public :: init => list_init, list_init_default
     procedure :: list_append
     procedure :: list_append_with_tag
     generic, public :: append => list_append, list_append_with_tag
     procedure :: list_prepend
     procedure :: list_prepend_with_tag
     generic, public :: prepend => list_prepend, list_prepend_with_tag
     procedure :: list_delete
     procedure :: list_delete_tag
     generic, public :: delete => list_delete, list_delete_tag
     procedure :: list_traverse
     procedure :: list_traverse_default
     generic, public :: traverse => list_traverse, list_traverse_default
     procedure :: list_find
     procedure :: list_find_default
     generic, public :: find => list_find, list_find_default
     procedure, public :: get => list_get
     procedure, public :: destroy => list_destroy
  end type list_type

  abstract interface
     subroutine list_iterator(node, stopped)
       import :: list_node_type
       type(list_node_type), pointer, intent(in out)  :: node
       PetscBool, intent(out) :: stopped
     end subroutine list_iterator
  end interface

contains

!------------------------------------------------------------------------

  subroutine list_init(self, delete_deallocates)
    !! Initialises list, specifying whether deletion of nodes should
    !! destroy (i.e. deallocate) their data.
    
    class(list_type), intent(in out) :: self
    PetscBool, intent(in) :: delete_deallocates

    self%delete_deallocates = delete_deallocates
    self%count = 0
    self%head => null()
    self%tail => null()

  end subroutine list_init

  subroutine list_init_default(self)
    !! Initialises list with default delete_deallocates = false.

    class(list_type), intent(in out) :: self

    call self%init(PETSC_FALSE)

  end subroutine list_init_default

!------------------------------------------------------------------------

  subroutine list_append(self, data)
    !! Adds node with specified data to end of list.

    class(list_type), intent(in out) :: self
    class(*), target, intent(in) :: data
    ! Locals:
    type(list_node_type), pointer :: old_tail

    if (associated(self%tail)) then
       old_tail => self%tail
       allocate(self%tail%next)
       self%tail => self%tail%next
       self%tail%previous => old_tail
    else ! empty:
       allocate(self%tail)
       self%head => self%tail
    end if

    self%tail%data => data
    self%count = self%count + 1

  end subroutine list_append

  subroutine list_append_with_tag(self, data, tag)
    !! Adds node with specified data to end of list, with a tag.

    class(list_type), intent(in out) :: self
    class(*), target, intent(in) :: data
    character(len = *), intent(in) :: tag

    call self%append(data)
    allocate(self%tail%tag, source = tag)

  end subroutine list_append_with_tag

!------------------------------------------------------------------------

  subroutine list_prepend(self, data)
    !! Adds node with specified data to start of list.

    class(list_type), intent(in out) :: self
    class(*), target, intent(in) :: data
    ! Locals:
    type(list_node_type), pointer :: old_head

    if (associated(self%head)) then
       old_head => self%head
       allocate(self%head%previous)
       self%head => self%head%previous
       self%head%next => old_head
    else ! empty:
       allocate(self%head)
       self%tail => self%head
    end if

    self%head%data => data
    self%count = self%count + 1

  end subroutine list_prepend

  subroutine list_prepend_with_tag(self, data, tag)
    !! Adds node with specified data to start of list, with a tag.

    class(list_type), intent(in out) :: self
    class(*), target, intent(in) :: data
    character(len = *), intent(in) :: tag

    call self%prepend(data)
    allocate(self%head%tag, source = tag)

  end subroutine list_prepend_with_tag

!------------------------------------------------------------------------

  subroutine list_delete(self, node)
    !! Deletes specified list node. If self%delete_deallocates is true,
    !! this will also deallocate the data in the node.

    class(list_type), intent(in out) :: self
    type(list_node_type), pointer :: node
    ! Locals:
    type(list_node_type), pointer :: previous, next

    if (associated(node)) then

       previous => node%previous
       next => node%next

       if (self%delete_deallocates) then
          deallocate(node%data)
       end if
       node%data => null()
       node%next => null()
       node%previous => null()
       deallocate(node)
       node => null()

       if (associated(previous)) then
          previous%next => next
       else
          self%head => next
       end if

       if (associated(next)) then
          next%previous => previous
       else
          self%tail => previous
       end if

       self%count = self%count - 1

    end if

  end subroutine list_delete

  subroutine list_delete_tag(self, tag)
    !! Deletes first list node with specified tag, if one exists.

    class(list_type), intent(in out) :: self
    character(len = *), intent(in) :: tag
    ! Locals:
    type(list_node_type), pointer :: node

    node => self%find(tag)
    call self%delete(node)
    
  end subroutine list_delete_tag

!------------------------------------------------------------------------

  subroutine list_traverse(self, iterator, backwards)
    !! Traverses the list, applying the iterator function to each
    !! node. If backwards is true, the list is traversed from tail to
    !! head, otherwise from head to tail.

    class(list_type), intent(in out) :: self
    procedure(list_iterator)  :: iterator
    PetscBool, intent(in) :: backwards
    ! Locals:
    type(list_node_type), pointer :: node
    PetscBool :: stopped

    stopped = PETSC_FALSE

    if (backwards) then
       node => self%tail
    else
       node => self%head
    end if
    
    if (associated(node)) then
       do
          call iterator(node, stopped)
          if (backwards) then
             node => node%previous
          else
             node => node%next
          end if
          if (stopped .or. .not. (associated(node))) then
             exit
          end if
       end do
    end if

  end subroutine list_traverse

  subroutine list_traverse_default(self, iterator)
    !! Traverses the list in the default forward direction.

    class(list_type), intent(in out) :: self
    procedure(list_iterator)  :: iterator

    call self%traverse(iterator, PETSC_FALSE)

  end subroutine list_traverse_default

!------------------------------------------------------------------------

  function list_find(self, tag, backwards) result(found_node)
    !! Returns pointer to first list node with specified tag, if one
    !! exists, or a null pointer otherwise. If backwards is true, the
    !! list is searched from tail to head, otherwise from head to
    !! tail.

    class(list_type), intent(in out) :: self
    character(len = *), intent(in) :: tag
    PetscBool, intent(in) :: backwards
    type(list_node_type), pointer :: found_node

    found_node => null()
    call self%traverse(tag_find_iterator, backwards)

  contains

    subroutine tag_find_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      if (allocated(node%tag) .and. (node%tag == tag)) then
         stopped = PETSC_TRUE
         found_node => node
      else
         stopped = PETSC_FALSE
      end if
    end subroutine tag_find_iterator

  end function list_find

  function list_find_default(self, tag) result(found_node)
    !! Returns pointer to first list node with specified tag,
    !! searching in the default forward direction.

    class(list_type), intent(in out) :: self
    character(len = *), intent(in) :: tag
    type(list_node_type), pointer :: found_node

    found_node => self%find(tag, PETSC_FALSE)

  end function list_find_default

!------------------------------------------------------------------------

  function list_get(self, index) result(found_node)
    !! Returns pointer to list node with specified (zero-based) index,
    !! if it exists, or a null pointer otherwise.  index = 0
    !! corresponds to the head node, index = 1 the second node, index
    !! = -1 the tail node, index = -2 the second to last node, etc.

    class(list_type), intent(in out) :: self
    PetscInt, intent(in) :: index
    type(list_node_type), pointer :: found_node
    ! Locals:
    PetscInt :: count, target_count
    PetscBool :: backwards

    found_node => null()
    count = 0
    if (index >= 0) then
       backwards = PETSC_FALSE
       target_count = index
    else
       backwards = PETSC_TRUE
       target_count = -index - 1
    end if

    call self%traverse(get_iterator, backwards)

  contains

    subroutine get_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      if (count == target_count) then
         stopped = PETSC_TRUE
         found_node => node
      else
         count = count + 1
      end if
    end subroutine get_iterator

  end function list_get

!------------------------------------------------------------------------

  subroutine list_destroy(self)
    !! Destroys list.

    class(list_type), intent(in out) :: self

    do while (associated(self%head))
       call self%delete(self%head)
    end do

  end subroutine list_destroy

!------------------------------------------------------------------------

end module list_module
