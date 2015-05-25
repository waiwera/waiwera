module ode_module
  !! Module for ordinary differential equations, to be solved by timestepper class.

  implicit none

  private

#include <petsc-finclude/petsc.h90>

  PetscInt, parameter, public :: max_title_length = 120

  type, public, abstract :: ode_type
     private
     character(max_title_length), public :: title
     Vec, public :: initial
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

end module ode_module
