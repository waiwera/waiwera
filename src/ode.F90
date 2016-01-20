module ode_module
  !! Abstract base class for ordinary differential equations defined over a mesh,
  !! to be solved by timestepper class.

  use mesh_module
  use logfile_module

  implicit none

  private

#include <petsc/finclude/petsc.h90>

  type, public, abstract :: ode_type
     private
     PetscReal, public :: time
     Vec, public :: solution
     type(mesh_type),  public :: mesh
     type(logfile_type), allocatable, public :: logfile
   contains
     private
     procedure(lhs_function), public, deferred :: lhs
     procedure(rhs_function), public, deferred :: rhs
     procedure, public :: pre_solve => ode_pre_eval
     procedure, public :: pre_iteration => ode_pre_iteration
     procedure, public :: pre_eval => ode_pre_eval
     procedure, public :: pre_timestep => ode_pre_timestep
     procedure, public :: pre_retry_timestep => ode_pre_retry_timestep
     procedure, public :: post_linesearch => ode_post_linesearch
     procedure(ode_output_procedure), public, deferred :: output
  end type ode_type

  abstract interface

     subroutine lhs_function(self, t, y, lhs, err)
       !! LHS function lhs = L(t, y)
       import :: ode_type
       class(ode_type), intent(in out) :: self
       PetscReal, intent(in) :: t
       Vec, intent(in) :: y
       Vec, intent(out) :: lhs
       PetscErrorCode, intent(out) :: err
     end subroutine lhs_function

     subroutine rhs_function(self, t, y, rhs, err)
       !! RHS function rhs = R(t, y)
       import :: ode_type
       class(ode_type), intent(in out) :: self
       PetscReal, intent(in) :: t
       Vec, intent(in) :: y
       Vec, intent(out) :: rhs
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

  subroutine ode_pre_eval(self, t, y, err)
    !! Default routine to be called before each evaluation
    !! of LHS and RHS functions.

    class(ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
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

end module ode_module
