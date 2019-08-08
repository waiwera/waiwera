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

module dag_module

  !! Module for directed acyclic graphs (DAG).
  !! Based on code made available by Jacob Williams: 
  !! http://degenerateconic.com/topological-sorting/

#include <petsc/finclude/petscsys.h>

  use petscsys

  implicit none
  private

  type :: dag_vertex_type
     !! Vertex of a DAG.
     private
     PetscInt, allocatable, public :: edges(:) !! Other vertices this one depends on
     PetscInt, public :: index = 0 !! Vertex index in DAG (zero-based)
     PetscBool, public :: checking !! Whether vertex is being checked during sort
     PetscBool, public :: done !! Whether vertex has been marked as done during sort
   contains
     private
     procedure, public :: set_edges => vertex_set_edges
  end type dag_vertex_type

  type, public :: dag_type
     !! Directed acyclic graph (DAG).
     private
     type(dag_vertex_type), allocatable, public :: vertices(:) !! DAG vertices
     PetscInt, public :: count !! Number of vertices
   contains
     private
     procedure, public :: init => dag_init
     procedure, public :: set_edges => dag_set_edges
     procedure, public :: sort => dag_sort
  end type dag_type

contains

!------------------------------------------------------------------------
! DAG vertex methods
!------------------------------------------------------------------------

  subroutine vertex_set_edges(self, edges)
    !! Specify the edge indices for this vertex.

    class(dag_vertex_type), intent(in out) :: self
    PetscInt, intent(in) :: edges(:) !! Edge indices array for this vertex

    self%edges = edges

  end subroutine vertex_set_edges

!------------------------------------------------------------------------
! DAG methods
!------------------------------------------------------------------------

  subroutine dag_init(self, count)
    !! Sets the number of vertices in the DAG.

    class(dag_type), intent(in out) :: self
    PetscInt, intent(in) :: count !! Number of vertices
    ! Locals:
    PetscInt :: i

    self%count = count
    allocate(self%vertices(0: count - 1))
    self%vertices%index = [(i, i = 0, count - 1)]

  end subroutine dag_init

!------------------------------------------------------------------------

  subroutine dag_set_edges(self, index, edges)
    !! Set the edges for a specified vertex.

    class(dag_type), intent(in out) :: self
    PetscInt, intent(in) :: index !! Vertex index
    PetscInt, intent(in) :: edges(:) !! Edge indices array for the vertex

    call self%vertices(index)%set_edges(edges)

  end subroutine dag_set_edges

!------------------------------------------------------------------------

  subroutine dag_sort(self, order, err)
    !! Return an array of indices corresponding to the topological
    !! sort order of the DAG. The error flag returns 1 if a circular
    !! dependency is encountered, zero otherwise.

    class(dag_type), intent(in out) :: self
    PetscInt, allocatable, intent(out) :: order(:) !! Returned order array
    PetscInt, intent(out) :: err !! Error flag
    ! Locals:
    PetscInt :: i, iorder

    err = 0
    allocate(order(0: self%count - 1))
    order = 0
    iorder = 0  ! index in order array

    self%vertices%checking = .false.
    self%vertices%done = .false.

    do i = 0, self%count - 1
       if (.not. self%vertices(i)%done) then
          call traverse(self%vertices(i), err)
          if (err > 0) exit
       end if
    end do
 
  contains

    recursive subroutine traverse(vertex, err)
      !! Depth-first graph traversal. The error flag returns 1 if a
      !! vertex is encountered which is already being checked.

      type(dag_vertex_type), intent(in out) :: vertex
      PetscInt, intent(out) :: err
      ! Locals:
      PetscInt :: j

      if (vertex%checking) then
         err = 1
      else
         err = 0
         if (.not. vertex%done) then
            vertex%checking = .true.
            if (allocated(vertex%edges)) then
               do j = 1, size(vertex%edges)
                  call traverse(self%vertices(vertex%edges(j)), err)
                  if (err > 0) exit
               end do
            end if
            if (err == 0) then
               vertex%checking = .false.
               vertex%done = .true.
               order(iorder) = vertex%index
               iorder = iorder + 1
            end if
         end if
      end if

    end subroutine traverse

  end subroutine dag_sort

!------------------------------------------------------------------------

end module dag_module
