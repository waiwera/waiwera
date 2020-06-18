!   Copyright 2016 University of Auckland.

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

module ode_module
  !! Abstract base class for ordinary differential equations defined over a mesh,
  !! to be solved by timestepper class.

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

  use petscsys
  use petscvec
  use mesh_module
  use logfile_module

  implicit none

  private

  type, public, abstract :: ode_type
     private
     PetscReal, public :: time
     Vec, public :: solution
     Vec, public :: aux_solution
     type(mesh_type),  public :: mesh
     type(logfile_type), public :: logfile
     PetscBool, public :: auxiliary
   contains
     private
     procedure(lhs_function), public, deferred :: lhs
     procedure(rhs_function), public, deferred :: rhs
     procedure, public :: pre_solve => ode_pre_eval
     procedure, public :: pre_iteration => ode_pre_iteration
     procedure, public :: pre_eval => ode_pre_eval
     procedure, public :: pre_timestep => ode_pre_timestep
     procedure, public :: post_timestep => ode_post_timestep
     procedure, public :: pre_retry_timestep => ode_pre_retry_timestep
     procedure, public :: post_linesearch => ode_post_linesearch
     procedure, public :: boundary_residuals => ode_boundary_residuals
     procedure(ode_output_procedure), public, deferred :: output
  end type ode_type

  abstract interface

     subroutine lhs_function(self, t, interval, y, lhs, err)
       !! LHS function lhs = L(t, y)
       use petscvec
       import :: ode_type
       class(ode_type), intent(in out) :: self
       PetscReal, intent(in) :: t, interval(2)
       Vec, intent(in) :: y
       Vec, intent(in out) :: lhs
       PetscErrorCode, intent(out) :: err
     end subroutine lhs_function

     subroutine rhs_function(self, t, interval, y, rhs, err)
       !! RHS function rhs = R(t, y)
       use petscvec
       import :: ode_type
       class(ode_type), intent(in out) :: self
       PetscReal, intent(in) :: t, interval(2)
       Vec, intent(in) :: y
       Vec, intent(in out) :: rhs
       PetscErrorCode, intent(out) :: err
     end subroutine rhs_function

     subroutine ode_output_procedure(self, time_index, time)
       !! Routine for output of ode solution.
       import :: ode_type
       class(ode_type), intent(in out) :: self
       PetscInt, intent(in) :: time_index
       PetscReal, intent(in) :: time
     end subroutine ode_output_procedure

  end interface

  public :: lhs_function, rhs_function
  public :: ode_output_procedure

!------------------------------------------------------------------------

contains

!------------------------------------------------------------------------

  subroutine ode_pre_iteration(self, y, err)
    !! Default routine to be called before each nonlinear solver
    !! iteration during solution at each time step.

    class(ode_type), intent(in out) :: self
    Vec, intent(in out) :: y
    PetscErrorCode, intent(out) :: err

    ! Do nothing
    err = 0

  end subroutine ode_pre_iteration

!------------------------------------------------------------------------

  subroutine ode_pre_eval(self, t, y, perturbed_columns, err)
    !! Default routine to be called before each evaluation
    !! of LHS and RHS functions.

    class(ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    PetscInt, intent(in), optional :: perturbed_columns(:)
    PetscErrorCode, intent(out) :: err

    ! Do nothing
    err = 0

  end subroutine ode_pre_eval

!------------------------------------------------------------------------

  subroutine ode_pre_timestep(self)
    !! Default routine to be called before starting each timestep.

    class(ode_type), intent(in out) :: self

    ! Do nothing

  end subroutine ode_pre_timestep

!------------------------------------------------------------------------

  subroutine ode_post_timestep(self)
    !! Default routine to be called after completing each timestep.

    class(ode_type), intent(in out) :: self

    ! Do nothing

  end subroutine ode_post_timestep

!------------------------------------------------------------------------

  subroutine ode_pre_retry_timestep(self)
    !! Default routine to be called before re-trying a timestep with a
    !! reduced step size.

    class(ode_type), intent(in out) :: self

    ! Do nothing

  end subroutine ode_pre_retry_timestep

!------------------------------------------------------------------------

  subroutine ode_post_linesearch(self, y_old, search, y, &
       changed_search, changed_y, err)
    !! Default routine to be called after nonlinear solver line search
    !! at each time step.

    class(ode_type), intent(in out) :: self
    Vec, intent(in) :: y_old
    Vec, intent(in out) :: search, y
    PetscBool, intent(out) :: changed_search, changed_y
    PetscErrorCode, intent(out) :: err

    ! Do nothing
    err = 0

  end subroutine ode_post_linesearch

!------------------------------------------------------------------------

  subroutine ode_boundary_residuals(self, y, lhs, residual, err)
    !! Computes residual terms for boundary ghost cells and adds them
    !! in to the specified residual vector. The residuals are computed
    !! from the difference between the LHS term for the specified
    !! primary variables y, and those in the specified initial lhs
    !! vector. This routine is used for direct steady state runs.

    class(ode_type), intent(in out) :: self
    Vec, intent(in) :: y !! primary variables
    Vec, intent(in) :: lhs !! initial LHS vector
    Vec, intent(in out) :: residual !! residual vector
    PetscErrorCode, intent(out) :: err !! error code

    ! Do nothing
    err = 0

  end subroutine ode_boundary_residuals

!------------------------------------------------------------------------

end module ode_module
