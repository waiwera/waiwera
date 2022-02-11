!   Copyright 2017 University of Auckland.

!   This file is part of Waiwera.

!   Waiwera is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.

!   Waiwera is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.

!   You should have received a copy of the GNU Lesser General Public License
!   along with Waiwera.  If not, see <http://www.gnu.org/licenses/>.

module dictionary_module
  !! Module for dictionary, associating a character key with a data node.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use list_module

  implicit none
  private

  PetscInt, parameter :: max_key_size  = 5 !! Maximum key length for generating hash
  PetscInt, parameter :: hash_size  = 4993 !! Number of hash lists
  
  type, public :: dictionary_type
     !! Dictionary type.
     private
     type(list_type), allocatable :: list(:) !! Array of hash lists
     PetscBool :: owner !! Whether to delete data on destroy()
   contains
     private
     procedure :: hash_key => dictionary_hash_key
     procedure :: init_default => dictionary_init_default
     procedure :: init_list => dictionary_init_list
     generic, public :: init => init_default, init_list
     procedure :: destroy_default => dictionary_destroy_default
     procedure :: destroy_proc => dictionary_destroy_proc
     generic, public :: destroy => destroy_default, destroy_proc
     procedure, public :: add => dictionary_add
     procedure :: delete_default => dictionary_delete_default
     procedure :: delete_destroy => dictionary_delete_destroy
     generic, public :: delete => delete_default, delete_destroy
     procedure, public :: has => dictionary_has
     procedure, public :: get => dictionary_get
     procedure, public :: count => dictionary_count
  end type dictionary_type

contains

!------------------------------------------------------------------------

  PetscInt function dictionary_hash_key(self, key) result(hash_key)
    !! Returns hash key for specified character key. The hash function
    !! is from Kernighan and Pike (1999), "The Practice of
    !! Programming".

    class(dictionary_type), intent(in) :: self
    character(*), intent(in) :: key !! Character key
    ! Locals:
    PetscInt :: i
    PetscInt, parameter :: hash_multiplier = 31

    hash_key = 0
    do i = 1, min(len_trim(key), max_key_size)
       hash_key = hash_multiplier * hash_key + ichar(key(i:i))
    end do
    hash_key = 1 + mod(hash_key - 1, hash_size)

  end function dictionary_hash_key

!------------------------------------------------------------------------

  subroutine dictionary_init_default(self, owner)
    !! Initialise dictionary.

    class(dictionary_type), intent(in out) :: self
    PetscBool, intent(in), optional :: owner !! Whether to delete data on destroy()
    ! Locals:
    PetscInt :: i
    PetscBool, parameter :: default_owner = PETSC_FALSE

    if (present(owner)) then
       self%owner = owner
    else
       self%owner = default_owner
    end if

    allocate(self%list(hash_size))
    do i = 1, hash_size
       call self%list(i)%init(owner = self%owner)
    end do

  end subroutine dictionary_init_default

!------------------------------------------------------------------------

  subroutine dictionary_init_list(self, list, transfer_owner)
    !! Initialise dictionary using a tagged list, with the (non-empty)
    !! list node tags becoming the dictionary keys.

    class(dictionary_type), intent(in out) :: self
    type(list_type), intent(in out) :: list !! List to build dictionary on
    PetscBool, intent(in), optional :: transfer_owner !! Whether to transfer data ownership from the list to the dictionary
    ! Locals:
    PetscBool :: transfer, dict_owner
    PetscBool, parameter :: default_transfer_owner = PETSC_FALSE

    if (present(transfer_owner)) then
       transfer = transfer_owner
    else
       transfer = default_transfer_owner
    end if

    if (transfer .and. list%owner) then
       list%owner = PETSC_FALSE
       dict_owner = PETSC_TRUE
    else
       dict_owner = PETSC_FALSE
    end if

    call self%init(dict_owner)
    call list%traverse(add_node_to_dict)

  contains

    subroutine add_node_to_dict(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      if (node%tag /= "") then
         if (associated(node%data)) then
            call self%add(node%tag, node%data)
         else
            call self%add(node%tag)
         end if
      end if
      
    end subroutine add_node_to_dict

  end subroutine dictionary_init_list

!------------------------------------------------------------------------

  subroutine dictionary_destroy_default(self)
    !! Destroy dictionary (without destroying any data objects).

    class(dictionary_type), intent(in out) :: self
    ! Locals:
    PetscInt :: i

    do i = 1, size(self%list)
       call self%list(i)%destroy()
    end do
    deallocate(self%list)

  end subroutine dictionary_destroy_default

!------------------------------------------------------------------------

  subroutine dictionary_destroy_proc(self, node_data_destroy_procedure)
    !! Destroy dictionary, applying specified procedure to destroy the
    !! data in the given node.

    class(dictionary_type), intent(in out) :: self
    procedure(list_node_data_destroy_procedure) :: node_data_destroy_procedure
    ! Locals:
    PetscInt :: i

    do i = 1, size(self%list)
       call self%list(i)%destroy(node_data_destroy_procedure)
    end do
    deallocate(self%list)

  end subroutine dictionary_destroy_proc

!------------------------------------------------------------------------

  function dictionary_get(self, key) result(node)
    !! Returns pointer to dictionary list node for specified key.

    class(dictionary_type), intent(in) :: self
    character(*), intent(in) :: key !! Key to find
    type(list_node_type), pointer :: node
    ! Locals:
    PetscInt :: hash
    type(list_type) :: list

    hash = self%hash_key(key)
    list = self%list(hash)
    node => list%find(key)
    
  end function dictionary_get

!------------------------------------------------------------------------

  subroutine dictionary_add(self, key, data)
    !! Add key, data pair to dictionary. If data is not present, a
    !! node is still added. This can be used to create a dictionary
    !! without data, which can be used essentially as a set of
    !! strings.

    class(dictionary_type), intent(in out) :: self
    character(*), intent(in) :: key !! Key to add
    class(*), target, intent(in), optional :: data !! Data to add
    ! Locals:
    type(list_node_type), pointer :: node
    PetscInt :: hash

    node => self%get(key)

    if (associated(node)) then
       if (present(data)) then
          node%data => data
       end if
    else
       hash = self%hash_key(key)
       if (present(data)) then
          call self%list(hash)%append(data, key)
       else
          call self%list(hash)%append(tag = key)
       end if
    end if
    
  end subroutine dictionary_add

!------------------------------------------------------------------------

  subroutine dictionary_delete_default(self, key)
    !! Delete key from dictionary, without destroying node data.

    class(dictionary_type), intent(in out) :: self
    character(*), intent(in) :: key !! Key to delete
    ! Locals:
    type(list_node_type), pointer :: node
    PetscInt :: hash
    
    node => self%get(key)

    if (associated(node)) then
       hash = self%hash_key(key)
       call self%list(hash)%delete(node)
    end if

  end subroutine dictionary_delete_default

!------------------------------------------------------------------------

  subroutine dictionary_delete_destroy(self, key, node_data_destroy_procedure)
    !! Delete key from dictionary, and calls given procedure to
    !! destroy its data.

    class(dictionary_type), intent(in out) :: self
    character(*), intent(in) :: key !! Key to delete
    procedure(list_node_data_destroy_procedure) :: node_data_destroy_procedure
    ! Locals:
    type(list_node_type), pointer :: node
    PetscInt :: hash

    node => self%get(key)

    if ((self%owner) .and. (associated(node))) then
       call node_data_destroy_procedure(node)
       hash = self%hash_key(key)
       call self%list(hash)%delete(node)
    end if

  end subroutine dictionary_delete_destroy

!------------------------------------------------------------------------

  PetscBool function dictionary_has(self, key) result(has)
    !! Returns true if dictionary contains the specified key.

    class(dictionary_type), intent(in out) :: self
    character(*), intent(in) :: key !! Key to find
    ! Locals:
    type(list_node_type), pointer :: node

    node => self%get(key)
    has = associated(node)

  end function dictionary_has

!------------------------------------------------------------------------

  PetscInt function dictionary_count(self) result(count)
    !! Returns number of items in the dictionary.

    class(dictionary_type), intent(in out) :: self

    count = sum(self%list%count)

  end function dictionary_count

!------------------------------------------------------------------------

end module dictionary_module
