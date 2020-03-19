!   Copyright 2020 University of Auckland.

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

module kdtree_module

  !! Module for k-d tree data structure. Implementation based on
  !! Kennel (2004), "KDTREE2: Fortran 95 and C++ software to
  !! efficiently search for near neighbors in a multi-dimensional
  !! Euclidean space", https://arxiv.org/abs/physics/0408067.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module

  implicit none
  private

  type interval_type
     !! real interval
      PetscReal :: lower, upper
  end type interval_type

  type, public :: kdtree_result_type
     !! search result
      PetscReal :: priority
      PetscInt :: index
  end type kdtree_result_type

  type :: priority_queue_type
     !! priority queue type
     private
     PetscInt :: heap_size = 0
     type(kdtree_result_type), pointer :: elems(:)
   contains
     private
     procedure :: heapify => priority_queue_heapify
     procedure, public :: init => priority_queue_init
     procedure, public :: insert => priority_queue_insert
     procedure, public :: replace_max => priority_queue_replace_max
  end type priority_queue_type

  type :: kdtree_node_type
      !! k-d tree node type
      private
      PetscInt :: cut_dim
      PetscReal :: cut_val, cut_val_left, cut_val_right
      PetscInt :: l, u
      type (kdtree_node_type), pointer :: left, right
      type(interval_type), pointer :: box(:) => null()
    contains
      private
      procedure :: process_terminal_node => kdtree_node_process_terminal_node
      procedure :: process_terminal_node_fixed_ball => kdtree_node_process_terminal_node_fixed_ball
      procedure, public :: search => kdtree_node_search
   end type kdtree_node_type

  type, public :: kdtree_type
     !! k-d tree type
      PetscInt :: dim, num_data
      PetscReal, pointer :: data(:,:) => null()
      PetscInt, pointer :: ind(:) => null()
      type (kdtree_node_type), pointer :: root => null()
    contains
      private
      procedure :: build => kdtree_build
      procedure :: build_for_range => kdtree_build_for_range
      procedure :: coordinate_spread => kdtree_coordinate_spread
      procedure, public :: init => kdtree_init
      procedure, public :: destroy => kdtree_destroy
      procedure :: kdtree_search_n
      procedure :: kdtree_search_nb
      procedure :: kdtree_search_1
      procedure :: kdtree_search_1b
      generic, public :: search => kdtree_search_n, kdtree_search_nb, &
           kdtree_search_1, kdtree_search_1b
   end type kdtree_type

  type :: kdtree_search_type
     !! k-d tree search data
      private
      PetscInt :: dim
      PetscInt :: nn, nfound
      PetscReal :: ballsize
      PetscInt :: centeridx = 999, correltime = 9999
      PetscInt :: nalloc
      PetscBool :: overflow
      PetscReal, pointer :: qv(:)
      type(kdtree_result_type), pointer :: results(:)
      type(priority_queue_type) :: pq
      PetscReal, pointer :: data(:,:)
      PetscInt, pointer :: ind(:)
   end type kdtree_search_type

contains

!------------------------------------------------------------------------
! Utilities
!------------------------------------------------------------------------

  PetscInt function select_on_coordinate_value(v, ind, c, alpha, li, ui) &
       result(last_index)

    !! Move elements of ind around between l and u, so that all points
    !! <= than alpha (in c cooordinate) are first, and then all points
    !! > alpha are second.

    PetscInt, intent (in) :: c, li, ui
    PetscReal, intent(in) :: alpha
    PetscReal, intent(in) :: v(1:,1:)
    PetscInt, intent(in out) :: ind(1:)
    ! Locals:
    PetscInt :: tmp
    PetscInt :: lb, rb

    lb = li; rb = ui

    do while (lb < rb)
       if ( v(c,ind(lb)) <= alpha ) then
          lb = lb + 1
       else
          tmp = ind(lb); ind(lb) = ind(rb); ind(rb) = tmp
          rb = rb - 1
       end if
    end do

    if (v(c, ind(lb)) <= alpha) then
       last_index = lb
    else
       last_index = lb - 1
    end if

  end function select_on_coordinate_value

!------------------------------------------------------------------------

  PetscReal function bound_square_distance(x, amin, amax) result (d2)
    !! Return squared distance from x to bounds [amin, amax].

    PetscReal, intent(in) :: x, amin, amax

    if (x > amax) then
       d2 = (x - amax)**2
    else
       if (x < amin) then
          d2 = (amin - x)**2
       else
          d2 = 0.0_dp
       end if
    end if

  end function bound_square_distance

!------------------------------------------------------------------------
! Priority queue
!------------------------------------------------------------------------

  subroutine priority_queue_init(self, results, err)

    !! Initialise priority queue from existing results.

    class(priority_queue_type), intent(in out) :: self
    type(kdtree_result_type), target :: results(:)
    PetscErrorCode, intent(out) :: err

    if (size(results, 1) > 0) then
       self%elems => results
       self%heap_size = 0
       err = 0
    else
       err = 1
    end if

  end subroutine priority_queue_init

!------------------------------------------------------------------------

  subroutine priority_queue_heapify(self, iroot)

    !! Takes a heap rooted at 'iroot' and force it to be in the
    !! heap canonical form.

    class(priority_queue_type), intent(in out) :: self
    PetscInt, intent(in) :: iroot
    ! Locals:
    PetscInt :: i, l, r, largest
    PetscReal :: priority_i, priority_l, priority_r, priority_largest
    type(kdtree_result_type) :: temp

    i = iroot

    do
       l = 2 * i
       r = l + 1
       if (l > self%heap_size) then
          exit
       else
          priority_i = self%elems(i)%priority
          priority_l = self%elems(l)%priority
          if (priority_l > priority_i) then
             largest = l
             priority_largest = priority_l
          else
             largest = i
             priority_largest = priority_i
          end if
          if (r <= self%heap_size) then
             priority_r = self%elems(r)%priority
             if (priority_r > priority_largest) then
                largest = r
             end if
          end if
       end if
       if (i /= largest) then
          temp = self%elems(i)
          self%elems(i) = self%elems(largest)
          self%elems(largest) = temp
          i = largest
          cycle
       else
          exit
       end if
    end do

  end subroutine priority_queue_heapify

!------------------------------------------------------------------------

  PetscReal function priority_queue_insert(self, priority, index) &
       result(max_priority)
    !
    ! Insert a new element and return the new maximum priority,
    ! which may or may not be the same as the old maximum priority.
    !
    class(priority_queue_type), intent(in out)  :: self
    PetscReal, intent(in) :: priority
    PetscInt, intent(in) :: index
    ! Locals:
    PetscInt :: i, parent_i
    PetscReal :: parent_priority

    self%heap_size = self%heap_size + 1
    i = self%heap_size

    do while (i > 1)
       parent_i = int(i / 2)
       parent_priority = self%elems(parent_i)%priority
       if (priority > parent_priority) then
          self%elems(i) = kdtree_result_type(parent_priority, &
               self%elems(parent_i)%index)
          i = parent_i
       else
          exit
       end if
    end do

    self%elems(i) = kdtree_result_type(priority, index)
    max_priority = self%elems(1)%priority

  end function priority_queue_insert

!------------------------------------------------------------------------

  PetscReal function priority_queue_replace_max(self, priority, index) &
       result(max_priority)

    !! Replace the maximum priority element in the queue with
    !! (priority, index).  Return the new maximum priority, which may
    !! be larger or smaller than the old one.

    class(priority_queue_type), intent(in out) :: self
    PetscReal, intent(in) :: priority
    PetscInt, intent(in) :: index
    PetscInt :: parent, child, n
    PetscReal :: child_priority, child1_priority

    n = self%heap_size
    if (n >= 1) then
       parent = 1
       child = 2
       do while (child <= n)
          child_priority = self%elems(child)%priority
          if (child < n) then
             child1_priority = self%elems(child + 1)%priority
             if (child_priority < child1_priority) then
                child = child + 1
                child_priority = child1_priority
             end if
          end if
          if (priority >= child_priority) then
             exit
          else
             self%elems(parent) = self%elems(child)
             parent = child
             child = parent * 2
          end if
       end do
       self%elems(parent) = kdtree_result_type(priority, index)
       max_priority = self%elems(1)%priority
    else
       self%elems(1) = kdtree_result_type(priority, index)
       max_priority = priority
    end if

  end function priority_queue_replace_max

!------------------------------------------------------------------------
! k-d tree node
!------------------------------------------------------------------------

  recursive subroutine kdtree_node_search(self, search)

    !! Search k-d tree node

    class(kdtree_node_type), intent(in out) :: self
    type(kdtree_search_type), target :: search
    ! Locals:
    type(kdtree_node_type), pointer :: ncloser, nfarther
    PetscInt :: cut_dim, i
    PetscReal :: qval, d2
    PetscReal :: ballsize
    PetscReal, pointer :: qv(:)
    type(interval_type), pointer :: box(:)

    if (associated(self%left) .and. associated(self%right)) then

       qv => search%qv(1:)
       cut_dim = self%cut_dim
       qval = qv(cut_dim)

       if (qval < self%cut_val) then
          ncloser => self%left
          nfarther => self%right
          d2 = (self%cut_val_right - qval)**2
       else
          ncloser => self%right
          nfarther => self%left
          d2 = (self%cut_val_left - qval)**2
       end if

       if (associated(ncloser)) call ncloser%search(search)

       if (associated(nfarther)) then
          ballsize = search%ballsize
          if (d2 <= ballsize) then
             box => self%box(1:)
             do i = 1, search%dim
                if (i /= cut_dim) then
                   d2 = d2 + bound_square_distance(qv(i), box(i)%lower, box(i)%upper)
                   if (d2 > ballsize) then
                      return
                   end if
                end if
             end do
             call nfarther%search(search)
          end if
       end if

    else
       if (search%nn == 0) then
          call self%process_terminal_node_fixed_ball(search)
       else
          call self%process_terminal_node(search)
       end if
    end if

  end subroutine kdtree_node_search

!------------------------------------------------------------------------

  subroutine kdtree_node_process_terminal_node(self, search)

    !! Look for actual near neighbors, and update the search results
    !! on the search data structure.

    class(kdtree_node_type), intent(in out) :: self
    type(kdtree_search_type), target :: search
    ! Locals:
    PetscReal, pointer :: qv(:)
    PetscInt, pointer :: ind(:)
    PetscReal, pointer :: data(:,:)
    PetscInt :: dim, i, i_index, k, centeridx, correltime
    PetscReal :: ballsize, sd, newpri
    type(priority_queue_type), pointer :: pqp

    qv => search%qv(1:)
    pqp => search%pq
    dim = search%dim
    ballsize = search%ballsize
    ind => search%ind(1:)
    data => search%data(1:, 1:)
    centeridx = search%centeridx
    correltime = search%correltime

    mainloop: do i = self%l, self%u
       i_index = ind(i)
       sd = 0.0_dp
       do k = 1,dim
          sd = sd + (data(k, i_index) - qv(k))**2
          if (sd > ballsize) cycle mainloop
       end do

       if (centeridx > 0) then
          if (abs(i_index - centeridx) < correltime) cycle mainloop
       end if

       if (search%nfound < search%nn) then
          search%nfound = search%nfound + 1
          newpri = pqp%insert(sd, i_index)
          if (search%nfound == search%nn) ballsize = newpri
       else
          ballsize = pqp%replace_max(sd, i_index)
       end if

    end do mainloop

    search%ballsize = ballsize

  end subroutine kdtree_node_process_terminal_node

!------------------------------------------------------------------------

  subroutine kdtree_node_process_terminal_node_fixed_ball(self, search)

    !! Look for actual near neighbors, and update the search results on
    !! the search data structure, i.e. save all within a fixed ball.

    class(kdtree_node_type), intent(in out) :: self
    type(kdtree_search_type), target :: search
    ! Locals:
    PetscReal, pointer :: qv(:)
    PetscInt, pointer :: ind(:)
    PetscReal, pointer :: data(:,:)
    PetscInt :: nfound
    PetscInt :: dim, i, i_index, k
    PetscInt :: centeridx, correltime, nn
    PetscReal :: ballsize, sd

    qv => search%qv(1:)
    dim = search%dim
    ballsize = search%ballsize
    ind => search%ind(1:)
    data => search%data(1:,1:)
    centeridx = search%centeridx
    correltime = search%correltime
    nn = search%nn
    nfound = search%nfound

    mainloop: do i = self%l, self%u

       i_index = ind(i)
       sd = 0.0_dp
       do k = 1, dim
          sd = sd + (data(k, i_index) - qv(k))**2
          if (sd > ballsize) cycle mainloop
       end do

       if (centeridx > 0) then
          if (abs(i_index - centeridx) < correltime) cycle
       end if

       nfound = nfound + 1
       if (nfound > search%nalloc) then
          search%overflow = .true.
       else
          search%results(nfound)%priority = sd
          search%results(nfound)%index = i_index
       end if

    end do mainloop

    search%nfound = nfound

  end subroutine kdtree_node_process_terminal_node_fixed_ball

!------------------------------------------------------------------------
! k-d tree
!------------------------------------------------------------------------

  subroutine kdtree_init(self, data)

    !! Initialise k-d tree

    class(kdtree_type), intent(in out) :: self
    PetscReal, target :: data(:,:)

    self%data => data
    self%dim = size(data, 1)
    self%num_data = size(data, 2)
    call self%build()

  end subroutine kdtree_init

!------------------------------------------------------------------------

  subroutine kdtree_destroy(self)

    !! Destroy k-d tree

    class(kdtree_type), intent(in out) :: self

    call destroy_node(self%root)
    deallocate (self%ind)
    nullify (self%ind)
    return

  contains

    recursive subroutine destroy_node(np)

      !! Destroy k-d tree node

      type (kdtree_node_type), pointer :: np

      if (associated(np%left)) then
         call destroy_node(np%left)
         nullify (np%left)
      end if

      if (associated(np%right)) then
         call destroy_node(np%right)
         nullify (np%right)
      end if

      if (associated(np%box)) deallocate(np%box)
      deallocate(np)

    end subroutine destroy_node

  end subroutine kdtree_destroy

!------------------------------------------------------------------------

  subroutine kdtree_build(self)

    !! Build k-d tree

    class(kdtree_type), intent(in out) :: self
    ! Locals:
    PetscInt :: j
    type(kdtree_node_type), pointer :: dummy => null()

    allocate(self%ind(self%num_data))
    forall (j = 1: self%num_data)
       self%ind(j) = j
    end forall
    self%root => self%build_for_range(1, self%num_data, dummy)

  end subroutine kdtree_build

!------------------------------------------------------------------------

  recursive function kdtree_build_for_range(self, l, u, parent) result (node)

    !! Build tree for specified range

    class(kdtree_type), intent(in out) :: self
    type(kdtree_node_type), pointer :: node
    type(kdtree_node_type), pointer :: parent
    PetscInt, intent(in) :: l, u
    ! Locals:
    PetscInt :: i, c, m, dim
    PetscBool :: recompute
    PetscReal :: average
    PetscInt, parameter :: bucket_size = 12 !! The maximum number of points to keep in a terminal node.

    dim = self%dim
    allocate (node)
    allocate(node%box(dim))

    if (u < l) then
       nullify(node)
       return
    end if

    if ((u - l) <= bucket_size) then
       do i = 1, dim
          node%box(i) = self%coordinate_spread(i, l, u)
       end do
       node%cut_dim = 0
       node%cut_val = 0.0_dp
       node%l = l
       node%u = u
       node%left => null()
       node%right => null()
    else

       do i = 1, dim
          recompute = .true.
          if (associated(parent)) then
             if (i /= parent%cut_dim) then
                recompute=.false.
             end if
          end if
          if (recompute) then
             node%box(i) = self%coordinate_spread(i, l, u)
          else
             node%box(i) = parent%box(i)
          end if
       end do

       c = maxloc(node%box(1:dim)%upper - node%box(1:dim)%lower, 1)

       average = sum(self%data(c, self%ind(l: u))) / real(u - l + 1, dp)
       node%cut_val = average
       m = select_on_coordinate_value(self%data, self%ind, c, average, l, u)

       node%cut_dim = c
       node%l = l
       node%u = u

       node%left => self%build_for_range(l, m, node)
       node%right => self%build_for_range(m + 1, u, node)

       if (.not.associated(node%right)) then
          node%box = node%left%box
          node%cut_val_left = node%left%box(c)%upper
          node%cut_val = node%cut_val_left
       elseif (.not.associated(node%left)) then
          node%box = node%right%box
          node%cut_val_right = node%right%box(c)%lower
          node%cut_val = node%cut_val_right
       else
          node%cut_val_right = node%right%box(c)%lower
          node%cut_val_left = node%left%box(c)%upper
          node%cut_val = (node%cut_val_left + node%cut_val_right) / 2
          node%box%upper = max(node%left%box%upper, node%right%box%upper)
          node%box%lower = min(node%left%box%lower, node%right%box%lower)
       end if

    end if

  end function kdtree_build_for_range

!------------------------------------------------------------------------

  type(interval_type) function kdtree_coordinate_spread(self, c, l, u) &
       result(interv)

    ! Returns the spread in coordinate c, between l and u.

    class(kdtree_type), intent(in out) :: self
    PetscInt, intent (in) :: c, l, u
    ! Locals:
    PetscReal :: last, lmax, lmin, t, smin,smax
    PetscInt :: i, ulocal
    PetscReal, pointer :: v(:,:)
    PetscInt, pointer :: ind(:)

    v => self%data(1:,1:)
    ind => self%ind(1:)
    smin = v(c,ind(l))
    smax = smin
    ulocal = u

    do i = l + 2, ulocal, 2
       lmin = v(c, ind(i - 1))
       lmax = v(c, ind(i))
       if (lmin > lmax) then
          t = lmin
          lmin = lmax
          lmax = t
       end if
       if (smin > lmin) smin = lmin
       if (smax < lmax) smax = lmax
    end do
    if (i == ulocal + 1) then
       last = v(c, ind(ulocal))
       if (smin > last) smin = last
       if (smax < last) smax = last
    end if

    interv%lower = smin
    interv%upper = smax

  end function kdtree_coordinate_spread

!------------------------------------------------------------------------

  subroutine kdtree_search_n(self, x, n, results, err)

    ! Find the n vectors in the tree nearest to x in Euclidean norm,
    ! returning results in an array of kdtree_result_type, allocated
    ! in the calling program.

    class(kdtree_type), intent(in out) :: self
    PetscReal, target, intent(in) :: x(:)
    PetscInt, intent (in) :: n
    type(kdtree_result_type), target :: results(:)
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(kdtree_search_type) :: search

    err = 0

    search%ballsize = huge(1.0)
    search%qv => x
    search%nn = n
    search%nfound = 0
    search%centeridx = -1
    search%correltime = 0
    search%overflow = .false.
    search%results => results
    search%nalloc = n
    search%ind => self%ind
    search%data => self%data
    search%dim = self%dim

    if (size(search%results, 1) >= n) then
       call search%pq%init(results, err)
       if (err == 0) then
          call self%root%search(search)
       end if
    else
       err = 1
    end if

  end subroutine kdtree_search_n

!------------------------------------------------------------------------

  subroutine kdtree_search_nb(self, x, n, index, distance2, err)

    !! Find n vectors nearest to x, returning results in separate
    !! variables for index and (optional) square distance.

    class(kdtree_type), intent(in out) :: self
    PetscReal, target, intent(in) :: x(:)
    PetscInt, intent (in) :: n
    PetscInt, intent(out) :: index(:)
    PetscReal, optional, intent(out) :: distance2(:)
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(kdtree_result_type) :: results(size(index))

    call self%search(x, n, results, err)
    if (err == 0) then
       index = results%index
       if (present(distance2)) then
          distance2 = sqrt(results%priority)
       end if
    end if

  end subroutine kdtree_search_nb

!------------------------------------------------------------------------

  subroutine kdtree_search_1(self, x, result, err)

    !! Search and return only nearest result.

    class(kdtree_type), intent(in out) :: self
    PetscReal, target, intent(in) :: x(:)
    type(kdtree_result_type), intent(out) :: result
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(kdtree_result_type) :: results(1)

    call self%search(x, 1, results, err)
    if (err == 0) then
       result = results(1)
    end if

  end subroutine kdtree_search_1

!------------------------------------------------------------------------

  subroutine kdtree_search_1b(self, x, index, distance2, err)

    !! Search and return only nearest result, in separate variables
    !! for index and (optional) square distance.

    class(kdtree_type), intent(in out) :: self
    PetscReal, target, intent(in) :: x(:)
    PetscInt, intent(out) :: index
    PetscReal, optional, intent(out) :: distance2
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(kdtree_result_type) :: result

    call self%search(x, result, err)
    if (err == 0) then
       index = result%index
       if (present(distance2)) then
          distance2 = result%priority
       end if
    end if

  end subroutine kdtree_search_1b

!------------------------------------------------------------------------

end module kdtree_module
