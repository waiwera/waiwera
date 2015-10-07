module ode_module
  !! Abstract base class for ordinary differential equations defined over a mesh,
  !! to be solved by timestepper class.

  use mesh_module

  implicit none

  private

#include <petsc/finclude/petsc.h90>

  type, public, abstract :: ode_type
     private
     PetscReal, public :: time
     Vec, public :: solution
     type(mesh_type),  public :: mesh
   contains
     private
     procedure(lhs_function), public, deferred :: lhs
     procedure(rhs_function), public, deferred :: rhs
     procedure, public :: pre_eval => ode_pre_eval
     procedure, public :: pre_iteration => ode_pre_iteration
     procedure(ode_output_procedure), public, deferred :: output
  end type ode_type

  abstract interface

     subroutine lhs_function(self, t, y, lhs)
       !! LHS function lhs = L(t, y)
       import :: ode_type
       class(ode_type), intent(in out) :: self
       PetscReal, intent(in) :: t
       Vec, intent(in) :: y
       Vec, intent(out) :: lhs
     end subroutine lhs_function

     subroutine rhs_function(self, t, y, rhs)
       !! RHS function rhs = R(t, y)
       import :: ode_type
       class(ode_type), intent(in out) :: self
       PetscReal, intent(in) :: t
       Vec, intent(in) :: y
       Vec, intent(out) :: rhs
     end subroutine rhs_function

     subroutine ode_output_procedure(self)
       !! Routine for output of ode solution.
       import :: ode_type
       class(ode_type), intent(in out) :: self
     end subroutine ode_output_procedure

  end interface

  public :: lhs_function, rhs_function
  public :: ode_output_procedure

!------------------------------------------------------------------------

contains

!------------------------------------------------------------------------

  subroutine ode_pre_eval(self, t, y)
    !! Default routine to be called before each evaluation
    !! of LHS and RHS functions.

    class(ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y

    ! Do nothing

  end subroutine ode_pre_eval

!------------------------------------------------------------------------

  subroutine ode_pre_iteration(self, ierr)
    !! Default routine to be called before each nonlinear solver
    !! iteration during solution at each time step.

    class(ode_type), intent(in out) :: self
    PetscErrorCode, intent(out) :: ierr

    ierr = 0

  end subroutine ode_pre_iteration

!------------------------------------------------------------------------

end module ode_module
