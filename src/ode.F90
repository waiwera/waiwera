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
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscdm.h>

  use petscsys
  use petscvec
  use petscmat
  use petscdm
  use mesh_module
  use logfile_module
  use kinds_module

  implicit none

  private

  type, public, abstract :: ode_type
     !! ODE type, of the form \(\frac{d}{dt} L(t,y) = R(t,y)\). The
     !! main ODE is assumed to be non-linear. An auxiliary linear
     !! problem can also optionally be solved, in which L = Al * y and
     !! R = Ar * y + br, where Al and Ar are matrices (Al is assumed
     !! diagonal) and br is a vector.
     private
     PetscReal, public :: time
     Vec, public :: solution
     Vec, public :: aux_solution
     type(mesh_type),  public :: mesh
     Mat, public :: jacobian
     Mat, public :: A_aux !! Left-hand side matrix for auxiliary linear problem
     Vec, public :: b_aux !! Right-hand side vector for auxiliary linear problem
     type(logfile_type), public :: logfile
     PetscBool, public :: auxiliary
   contains
     private
     procedure(lhs_function), public, deferred :: lhs
     procedure(rhs_function), public, deferred :: rhs
     procedure, public :: aux_lhs => ode_aux_lhs
     procedure, public :: aux_rhs => ode_aux_rhs
     procedure, public :: aux_pre_solve => ode_aux_pre_solve
     procedure, public :: pre_solve => ode_pre_eval
     procedure, public :: pre_iteration => ode_pre_iteration
     procedure, public :: pre_eval => ode_pre_eval
     procedure, public :: pre_timestep => ode_pre_timestep
     procedure, public :: post_timestep => ode_post_timestep
     procedure, public :: pre_try_timestep => ode_pre_try_timestep
     procedure, public :: pre_retry_timestep => ode_pre_retry_timestep
     procedure, public :: post_linesearch => ode_post_linesearch
     procedure, public :: setup_jacobian => ode_setup_jacobian
     procedure, public :: modify_jacobian => ode_modify_jacobian
     procedure, public :: setup_auxiliary => ode_setup_auxiliary
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

  subroutine ode_pre_try_timestep(self, t)
    !! Default routine to be called before trying a timestep (on first
    !! attempt or a re-try).

    class(ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time

    ! Do nothing

  end subroutine ode_pre_try_timestep

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

  subroutine ode_aux_lhs(self, t, interval, Al, err)
    !! Default routine returning diagonal left-hand side matrix Al
    !! where lhs = Al * X. Here Al = I.

    class(ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Vec, intent(in out) :: Al
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr

    err = 0
    call VecSet(Al, 1._dp, ierr); CHKERRQ(ierr)

  end subroutine ode_aux_lhs

!------------------------------------------------------------------------

  subroutine ode_aux_rhs(self, t, interval, Ar, br, err)
    !! Default routine returning right-hand side matrix Al and vector
    !! br where rhs = Ar * X + br. Here Ar and br are zero.

    class(ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Mat, intent(in out) :: Ar
    Vec, intent(in out) :: br
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr

    err = 0
    call MatZeroEntries(Ar, ierr); CHKERRQ(ierr)
    call VecSet(br, 0._dp, ierr); CHKERRQ(ierr)

  end subroutine ode_aux_rhs

!------------------------------------------------------------------------

  subroutine ode_aux_pre_solve(self)
    !! Default routine for modifying auxiliary linear system Ax = b
    !! before solving it.

    class(ode_type), intent(in out) :: self

    ! Do nothing

  end subroutine ode_aux_pre_solve

!------------------------------------------------------------------------

  subroutine ode_setup_jacobian(self)
    !! Set up Jacobian matrix.

    class(ode_type), intent(in out) :: self
    ! Locals:
    PetscInt :: blocksize
    MatType :: mat_type
    PetscErrorCode :: ierr

    call VecGetBlockSize(self%solution, blocksize, ierr); CHKERRQ(ierr)
    if (blocksize == 1) then
       mat_type = MATAIJ
    else
       mat_type = MATBAIJ
    end if
    call DMSetMatType(self%mesh%interior_dm, mat_type, ierr); CHKERRQ(ierr)
    call DMCreateMatrix(self%mesh%interior_dm, self%jacobian, ierr)
    CHKERRQ(ierr)
    call self%modify_jacobian()
    call MatSetFromOptions(self%jacobian, ierr); CHKERRQ(ierr)

  end subroutine ode_setup_jacobian

!------------------------------------------------------------------------

  subroutine ode_modify_jacobian(self)
    !! Carries out any required modification of the Jacobian matrix.

    class(ode_type), intent(in out) :: self

    ! Do nothing

  end subroutine ode_modify_jacobian

!------------------------------------------------------------------------

  subroutine ode_setup_auxiliary(self)
    !! Sets up linear system for auxiliary problem.

    class(ode_type), intent(in out) :: self
    ! Locals:
    DM :: dm_aux
    PetscErrorCode :: ierr

    if (self%auxiliary) then
       call VecGetDM(self%aux_solution, dm_aux, ierr); CHKERRQ(ierr)
       call DMCreateMatrix(dm_aux, self%A_aux, ierr); CHKERRQ(ierr)
       call MatSetOption(self%A_aux, MAT_KEEP_NONZERO_PATTERN, &
            PETSC_TRUE, ierr); CHKERRQ(ierr)
       call MatSetFromOptions(self%A_aux, ierr); CHKERRQ(ierr)
       call VecDuplicate(self%aux_solution, self%b_aux, ierr); CHKERRQ(ierr)
    end if

  end subroutine ode_setup_auxiliary

!------------------------------------------------------------------------

end module ode_module
