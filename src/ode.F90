module ode_module
  !! Abstract base class for ordinary differential equations defined over a mesh,
  !! to be solved by timestepper class.

  use mesh_module

  implicit none

  private

#include <petsc-finclude/petsc.h90>

  type, public, abstract :: ode_type
     private
     Vec, public :: initial
     type(mesh_type),  public :: mesh
   contains
     private
     procedure(lhs_function), public, deferred :: lhs
     procedure(rhs_function), public, deferred :: rhs
     procedure(pre_eval_procedure), public, deferred :: pre_eval
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

     subroutine pre_eval_procedure(self, t, y)
       !! Optional routine to be called before each evaluation
       !! of LHS and RHS functions.
       import :: ode_type
       class(ode_type), intent(in out) :: self
       PetscReal, intent(in) :: t
       Vec, intent(in) :: y
     end subroutine pre_eval_procedure

  end interface

  public :: lhs_function, rhs_function, lhs_identity
  public :: pre_eval_procedure

contains

!------------------------------------------------------------------------

  subroutine lhs_identity(t, y, lhs)
    !! Default identity LHS function L(t,y) = y.

    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    Vec, intent(out) :: lhs
    ! Locals:
    PetscErrorCode :: ierr

    call VecCopy(y, lhs, ierr); CHKERRQ(ierr)

  end subroutine lhs_identity

!------------------------------------------------------------------------

end module ode_module
