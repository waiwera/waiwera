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

module timestepper_module
  !! Time stepping methods for solving \(\frac{d}{dt} L(t,y) = R(t,y)\), where y is a parallel vector.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use ode_module
  use logfile_module

  implicit none

  private

  ! Timestepping methods
     PetscInt, parameter, public :: TS_BEULER = 0, TS_BDF2 = 1, TS_DIRECTSS = 2

  ! Timestepping adaptor methods
     PetscInt, parameter, public :: TS_ADAPT_CHANGE = 0, TS_ADAPT_ITERATION = 1
     PetscInt, parameter, public :: timestepper_max_monitor_name_length = 12

  ! Timestep status
     PetscInt, parameter, public :: TIMESTEP_OK = 0, TIMESTEP_NOT_CONVERGED = 1, &
     TIMESTEP_TOO_SMALL = 2, TIMESTEP_TOO_BIG = 3, TIMESTEP_ABORTED = 4, &
     TIMESTEP_FINAL = 5

     PetscInt, parameter :: max_ksp_type_str_len = 8, max_pc_type_str_len = 8

  type timestepper_step_type
     !! Results of a time step.
     private
     PetscReal, public :: time !! End time for step
     PetscReal, public :: stepsize !! Step size
     Vec, public :: solution !! Global solution vector
     Vec, public :: lhs !! Global left-hand side vector
     Vec, public :: rhs !! Global right-hand side vector
     PetscInt, public :: num_tries !! Number of attempts at carrying out the step (with step size reductions after first attempt)
     PetscInt, public :: num_iterations !! Number of non-linear solver iterations
     PetscInt, public :: status = TIMESTEP_OK !! Step status
     PetscReal, public :: max_residual !! Maximum residual over the mesh cells
     PetscInt, public :: max_residual_cell !! Global cell index at which maximum residual occurs
     PetscInt, public :: max_residual_equation !! Equation number at which maximum residual occurs
   contains
     private
     procedure, public :: init => timestepper_step_init
     procedure, public :: destroy => timestepper_step_destroy
     procedure, public :: print => timestepper_step_print
     procedure, public :: status_str => timestepper_step_status_str
  end type timestepper_step_type

  type ptimestepper_step_type
     !! Pointer to timestepper_step_type.
     private
     type(timestepper_step_type), pointer, public :: p
  end type ptimestepper_step_type

  type timestep_adaptor_type
     !! For adaptive time step size control. Stepsize is increased if monitor value is
     !! below monitor_min, decreased if it is above monitor_max.
     private
     PetscBool, public :: on !! Whether the step size adaptor is active
     character(timestepper_max_monitor_name_length), public :: name !! Name of adaptor type
     procedure(monitor_function), pointer, nopass, public :: monitor => relative_change_monitor !! Monitor function
     PetscReal, public :: monitor_min !! Minimum value of monitor function (below which step size is increased)
     PetscReal, public :: monitor_max !! Maximum value of monitor function (above which step size is decreased)
     PetscReal, public :: reduction !! Factor by which step size is reduced
     PetscReal, public :: amplification !! Factor by which step size is increased
     PetscReal, public :: max_stepsize !! Maximum allowable step size
   contains
     private
     procedure :: increase => timestep_adaptor_increase
     procedure :: reduce => timestep_adaptor_reduce
  end type timestep_adaptor_type

  type, public :: timestepper_checkpoints_type
     !! Checkpoints for output at specified times.
     private
     PetscReal :: tolerance !! Relative tolerance for deciding in checkpoint has been hit
     PetscReal :: repeat_shift !! Time shift for each repeat cycle
     PetscInt, public :: index !! Index in time array
     PetscInt, public :: repeat_index !! Index of repeat
     PetscReal, allocatable, public :: time(:) !! Checkpoint times
     PetscReal, public :: restore_stepsize !! Time step size to restore after a checkpoint is hit
     PetscReal, public :: next_time
     PetscInt, public :: repeat !! Number of repeat cycles
     PetscBool, public :: hit !! Whether a checkpoint has been hit on the current time step
     PetscBool, public :: done !! Whether all checkpoints have been processed
   contains
     procedure, public :: init => timestepper_checkpoints_init
     procedure, public :: destroy => timestepper_checkpoints_destroy
     procedure, public :: update => timestepper_checkpoints_update
     procedure, public :: check => timestepper_checkpoints_check
  end type timestepper_checkpoints_type

  type, public :: timestepper_steps_type
     !! Storage for current and immediate past steps in timestepper.
     private
     type(timestepper_step_type), allocatable :: store(:)
     type(ptimestepper_step_type), pointer :: pstore(:)
     type(timestepper_step_type), pointer, public :: current !! Current time step
     type(timestepper_step_type), pointer, public :: last !! Previous time step
     PetscInt, public :: num_stored !! Number of steps stored
     PetscInt, public :: taken !! Number of steps taken so far
     PetscInt, public :: max_num !! Maximum allowable number of steps
     PetscInt, public :: max_num_tries !! Maximum allowable number of tries per step
     PetscInt, public :: fixed_step_index !! Current index in array of fixed step sizes
     PetscReal, public :: next_stepsize !! Step size to be used for next step
     PetscReal, public :: stop_time !! Maximum allowable time
     PetscReal, public :: stop_min_stepsize !! Minimum allowable step size
     PetscReal, public :: stop_max_stepsize !! Maximum allowable step size
     PetscReal, public :: nonlinear_solver_relative_tol !! Relative tolerance for non-linear solver function value
     PetscReal, public :: nonlinear_solver_abs_tol !! Absolute tolerance for non-linear solver function value
     PetscReal, public :: nonlinear_solver_update_relative_tol !! Relative tolerance for non-linear solver update
     PetscReal, public :: nonlinear_solver_update_abs_tol !! Absolute tolerance for non-linear solver update
     PetscInt, public :: nonlinear_solver_minimum_iterations !! Minimum number of non-linear solver iterations to do before checking for convergence
     PetscReal, public :: termination_tol = 1.e-3_dp !! Tolerance for detecting stop time reached
     PetscReal, allocatable, public :: sizes(:) !! Pre-specified time step sizes
     type(timestep_adaptor_type), public :: adaptor !! Time step size adaptor
     type(timestepper_checkpoints_type), public :: checkpoints !! Checkpoints
     PetscBool, public :: stop_time_specified !! Whether stop time is specified or not
     PetscBool, public :: fixed !! If fixed steps sizes are specified
     PetscBool :: finished !! Whether simulation has yet finished
     PetscBool :: steady_state !! If simulation is direct steady state or transient
   contains
     private
     procedure :: set_aliases => timestepper_steps_set_aliases
     procedure :: set_pstore => timestepper_steps_set_pstore
     procedure :: update => timestepper_steps_update
     procedure :: initialize_try => timestepper_steps_initialize_try
     procedure :: rotate => timestepper_steps_rotate
     procedure, public :: init => timestepper_steps_init
     procedure, public :: destroy => timestepper_steps_destroy
     procedure, public :: check_finished => timestepper_steps_check_finished
     procedure, public :: check_checkpoints => timestepper_steps_check_checkpoints
     procedure, public :: set_current_status => &
          timestepper_steps_set_current_status
     procedure, public :: adapt => timestepper_steps_adapt
     procedure, public :: get_next_fixed_stepsize => &
          timestepper_steps_get_next_fixed_stepsize
     procedure, public :: set_next_stepsize => timestepper_steps_set_next_stepsize
  end type timestepper_steps_type

  type timestepper_method_type
     !! Timestepping method.
     private
     character(20), public :: name
     procedure(method_residual), pointer, nopass, public :: residual
     PetscInt, public :: num_steps, num_stored_steps
   contains
     private
     procedure, public :: init => timestepper_method_init
  end type timestepper_method_type

  type timestepper_solver_context_type
     !! Context for SNES solver.
     private
     class(ode_type), pointer, public :: ode
     type(timestepper_steps_type), pointer, public :: steps
     procedure(method_residual), pointer, nopass, public :: residual
     MatFDColoring, public :: fd_coloring !! Finite difference colouring for Jacobian matrix
  end type timestepper_solver_context_type

  type, public :: timestepper_type
     !! Timestepper class.
     private
     SNES, public :: solver
     Vec :: residual
     Mat :: jacobian
     type(timestepper_solver_context_type) :: context
     class(ode_type), pointer, public :: ode !! ODE to be solved
     type(timestepper_steps_type), public :: steps !! Time steps object
     type(timestepper_method_type), public :: method !! Time stepping method
     procedure(step_output_routine), pointer, public :: &
          before_step_output => before_step_output_default !! Output function to be called before each step
     procedure(step_output_routine), pointer, public :: &
          after_step_output => after_step_output_default !! Output function to be called after each step
     PetscInt, public :: output_frequency !! Time step frequency for main output of results
     PetscInt, public :: output_index !! How many sets of results have been output
     PetscBool, public :: output_initial !! Whether to output initial conditions 
     PetscBool, public :: output_final !! Whether to output final results
   contains
     private
     procedure :: setup_jacobian => timestepper_setup_jacobian
     procedure :: setup_nonlinear_solver => timestepper_setup_nonlinear_solver
     procedure :: setup_linear_solver => timestepper_setup_linear_solver
     procedure :: setup_linear_solver_options => timestepper_setup_linear_solver_options
     procedure :: setup_linear_sub_solver => timestepper_setup_linear_sub_solver
     procedure :: step => timestepper_step
     procedure :: log_step_status => timestepper_log_step_status
     procedure, public :: init => timestepper_init
     procedure, public :: destroy => timestepper_destroy
     procedure, public :: initial_function_calls => timestepper_initial_function_calls
     procedure, public :: run => timestepper_run
     procedure, public :: input_summary => timestepper_input_summary
  end type timestepper_type

  interface

     subroutine method_residual(solver, y, residual, context, err)
       !! Residual routine to be minimised by nonlinear solver.
       use petscsnes
       use petscvec
       import :: timestepper_solver_context_type
       SNES, intent(in) :: solver
       Vec, intent(in) :: y
       Vec, intent(in out) :: residual
       type(timestepper_solver_context_type), intent(in out) :: context
       PetscErrorCode, intent(out) :: err
     end subroutine method_residual

     PetscReal function monitor_function(current, last)
       !! Function for monitoring timestep acceptability.
       import :: timestepper_step_type
       type(timestepper_step_type), intent(in) :: current, last
     end function monitor_function

     subroutine step_output_routine(self)
       !! Routine for producing output before or after each time step.
       import :: timestepper_type
       class(timestepper_type), intent(in out) :: self
     end subroutine step_output_routine

     subroutine SNESGetApplicationContext(solver, context, ierr)
       !! Interface for getting context from SNES solver- to cast it
       !! as the correct type.
       use petscsnes
       import :: timestepper_solver_context_type
       SNES, intent(in) :: solver
       type(timestepper_solver_context_type), pointer, &
            intent(out) :: context
       PetscErrorCode, intent(out) :: ierr
     end subroutine SNESGetApplicationContext

  end interface

  ! Subroutines to be available outside this module:
  public :: iteration_monitor, relative_change_monitor
  public :: step_output_routine, before_step_output_default, &
       after_step_output_default

contains

!------------------------------------------------------------------------
! Monitoring functions
!------------------------------------------------------------------------

  PetscReal function iteration_monitor(current, last)
    !! Monitors number of nonlinear iterations used to compute solution.

    type(timestepper_step_type), intent(in) :: current, last

    iteration_monitor = real(current%num_iterations)

  end function iteration_monitor

!------------------------------------------------------------------------

  PetscReal function relative_change_monitor(current, last)

    !! Monitors relative change in left hand side vector. The
    !! parameter eps avoids division by zero if any elements of
    !! last%lhs are zero.

    use dm_utils_module, only: vec_max_pointwise_abs_scale

    type(timestepper_step_type), intent(in) :: current, last
    ! Locals:
    PetscReal, parameter :: eps = 1.e-3_dp
    Vec :: diff
    PetscInt :: index
    PetscErrorCode :: ierr

    call VecDuplicate(current%lhs, diff, ierr); CHKERRQ(ierr)
    call VecCopy(current%lhs, diff, ierr); CHKERRQ(ierr)
    call VecAXPY(diff, -1.0_dp, last%lhs, ierr); CHKERRQ(ierr)
    call vec_max_pointwise_abs_scale(diff, last%lhs, &
         eps, relative_change_monitor, index)

  end function relative_change_monitor

!------------------------------------------------------------------------
! Step output routines
!------------------------------------------------------------------------

  subroutine before_step_output_default(self)
    !! Default routine for printing information at the start of each
    !! time step.

    class(timestepper_type), intent(in out) :: self

    call self%ode%logfile%write_blank()
    call self%ode%logfile%write(LOG_LEVEL_INFO, 'timestep', 'start', &
         ['count           '], [self%steps%taken + 1], &
         ['size            '], [self%steps%next_stepsize])

  end subroutine before_step_output_default

!------------------------------------------------------------------------

  subroutine after_step_output_default(self)
    !! Default routine for printing information at the end of each
    !! time step.

    class(timestepper_type), intent(in out) :: self

    call self%steps%current%print(self%ode%logfile)

  end subroutine after_step_output_default

!------------------------------------------------------------------------
! Residual routines
!------------------------------------------------------------------------

  subroutine backwards_Euler_residual(solver, y, residual, context, err)
    !! Residual for backwards Euler method.
    !! residual = L(1) - L(0)  - dt * R(1)

    SNES, intent(in) :: solver
    Vec, intent(in) :: y
    Vec, intent(in out) :: residual
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: t, dt, interval(2)
    PetscErrorCode :: ierr

    err = 0
    t = context%steps%current%time
    interval = [context%steps%last%time, t]
    dt = context%steps%current%stepsize

    call context%ode%lhs(t, interval, y, context%steps%current%lhs, err)
    if (err == 0) then
       call VecCopy(context%steps%current%lhs, residual, ierr); CHKERRQ(ierr)
       call VecAXPY(residual, -1.0_dp, context%steps%last%lhs, ierr); CHKERRQ(ierr)
       call context%ode%rhs(t, interval, y, context%steps%current%rhs, err)
       if (err == 0) then
          call VecAXPY(residual, -dt, context%steps%current%rhs, ierr)
          CHKERRQ(ierr)
       end if
    end if

  end subroutine backwards_Euler_residual

!------------------------------------------------------------------------

  subroutine BDF2_residual(solver, y, residual, context, err)
    !! Residual for variable-stepsize BDF2 method.
    !! residual = (1 + 2r) * L(1) - (r+1)^2 * L(0) + r^2 * L(-1) - dt * (r+1) * R(1)
    !! where r = dt / (last dt)

    SNES, intent(in) :: solver
    Vec, intent(in) :: y
    Vec, intent(in out) :: residual
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(timestepper_step_type), pointer :: last2
    PetscReal :: t, dt, dtlast, interval(2)
    PetscReal :: r, r1
    PetscErrorCode :: ierr

    if (context%steps%taken == 0) then

       ! Startup- use backwards Euler
       call backwards_Euler_residual(solver, y, residual, context, err)

    else

       err = 0
       t  = context%steps%current%time
       interval = [context%steps%last%time, t]
       dt = context%steps%current%stepsize
       dtlast = context%steps%last%stepsize
       r = dt / dtlast
       r1 = r + 1._dp

       last2 => context%steps%pstore(3)%p

       call context%ode%lhs(t, interval, y, context%steps%current%lhs, err)
       if (err == 0) then
          call VecCopy(context%steps%current%lhs, residual, ierr); CHKERRQ(ierr)
          call VecScale(residual, 1._dp + 2._dp * r, ierr); CHKERRQ(ierr)
          call VecAXPY(residual, -r1 * r1, context%steps%last%lhs, ierr)
          CHKERRQ(ierr)
          call VecAXPY(residual, r * r, last2%lhs, ierr); CHKERRQ(ierr)
          call context%ode%rhs(t, interval, y, context%steps%current%rhs, err)
          if (err == 0) then
             call VecAXPY(residual, -dt * r1, context%steps%current%rhs, ierr)
             CHKERRQ(ierr)
          end if
       end if

    end if

  end subroutine BDF2_residual

!------------------------------------------------------------------------

  subroutine direct_ss_residual(solver, y, residual, context, err)
    !! Residual for direct solution of steady state equations R(y) =
    !! 0.  Here we evaluate R() at the steps final time, and add
    !! residuals for boundary ghost cells.

    SNES, intent(in) :: solver
    Vec, intent(in) :: y
    Vec, intent(in out) :: residual
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: t, interval(2)
    PetscErrorCode :: ierr

    err = 0
    t = context%steps%stop_time
    interval = [t, t]
    call context%ode%rhs(t, interval, y, context%steps%current%rhs, err)
    if (err == 0) then
       call VecCopy(context%steps%current%rhs, residual, ierr)
       CHKERRQ(ierr)
       call context%ode%boundary_residuals(y, context%steps%current%lhs, &
            residual, err)
    end if

  end subroutine direct_ss_residual

!------------------------------------------------------------------------

  subroutine SNES_residual(solver, y, residual, context, err)
    !! Residual routine to be minimized by SNES solver. This calls the
    !! pre-evaluation routine first (if needed) before calling the
    !! timestepper method residual routine.
    !! If an error occurs a SNES function domain error is raised, but
    !! this routine still returns a zero error code (err), otherwise
    !! the SNES will stop.

    use dm_utils_module, only: mat_coloring_perturbed_columns

    SNES, intent(in) :: solver
    Vec, intent(in) :: y
    Vec, intent(in out) :: residual
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode, intent(out) :: err
    ! Locals:
    Mat :: J
    PetscErrorCode :: ierr, ferr
    PetscInt, allocatable :: perturbed_columns(:)

    err = 0; ferr = 0

    call SNESGetJacobian(solver, J, PETSC_NULL_MAT, PETSC_NULL_FUNCTION, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    perturbed_columns = mat_coloring_perturbed_columns(J, context%fd_coloring)
    call context%ode%pre_eval(context%steps%current%time, y, &
         perturbed_columns, ferr)
    deallocate(perturbed_columns)

    if (ferr == 0) then
       call context%residual(solver, y, residual, context, ferr)
    end if

    if (ferr > 0) then
       call SNESSetFunctionDomainError(solver, ierr); CHKERRQ(ierr)
    end if

  end subroutine SNES_residual

!------------------------------------------------------------------------

  subroutine SNES_Jacobian(solver, y, J, B, context, err)
    !! Jacobian calculation routine for SNES solver.

    SNES, intent(in) :: solver
    Vec, intent(in) :: y
    Mat, intent(in) :: J
    Mat, intent(in out) :: B
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode, intent(out) :: err

    call SNESComputeJacobianDefaultColor(solver, y, J, B, &
         context%fd_coloring, err); CHKERRQ(err)

  end subroutine SNES_Jacobian

!------------------------------------------------------------------------

  PetscErrorCode function SNES_pre_iteration_update(solver, step)
    !! Function to be called before each nonlinear solver iteration.

    SNES, intent(in out) :: solver
    PetscInt, intent(in) :: step
    ! Locals:
    type(timestepper_solver_context_type), pointer :: context
    PetscErrorCode :: ierr, err

    call SNESGetApplicationContext(solver, context, ierr); CHKERRQ(ierr)
    call context%ode%pre_iteration(context%steps%current%solution, err)

    if (err > 0) then
       call SNESSetFunctionDomainError(solver, ierr); CHKERRQ(ierr)
    end if
    SNES_pre_iteration_update = 0

  end function SNES_pre_iteration_update

!------------------------------------------------------------------------

  PetscErrorCode function SNES_linesearch_post_check(linesearch, &
       y_old, search, y, changed_search, changed_y, context)
    !! Function to be called after each nonlinear solver line search.

    SNESLineSearch, intent(in out) :: linesearch
    Vec, intent(in) :: y_old
    Vec, intent(in out) :: search, y
    PetscBool, intent(out) :: changed_search, changed_y
    type(timestepper_solver_context_type), intent(in out) :: context
    ! Locals:
    PetscErrorCode :: err, ierr
    SNES :: solver

    call context%ode%post_linesearch(y_old, search, y, changed_search, &
         changed_y, err)

    if (err > 0) then
       call SNESLineSearchGetSNES(linesearch, solver, ierr); CHKERRQ(ierr)
       call SNESSetFunctionDomainError(solver, ierr); CHKERRQ(ierr)
       call SNESLineSearchSetReason(linesearch, SNES_LINESEARCH_FAILED_DOMAIN, &
            ierr); CHKERRQ(ierr)
    end if
    SNES_linesearch_post_check = 0

  end function SNES_linesearch_post_check

!------------------------------------------------------------------------
! Timestepper_step procedures
!------------------------------------------------------------------------

  subroutine timestepper_step_init(self, template_vec)
    !! Initializes a timestep.

    class(timestepper_step_type), intent(in out) :: self
    Vec, intent(in) :: template_vec
    ! Locals:
    PetscErrorCode :: ierr

    call VecDuplicate(template_vec, self%solution, ierr); CHKERRQ(ierr)
    call VecDuplicate(template_vec, self%lhs, ierr); CHKERRQ(ierr)
    call VecDuplicate(template_vec, self%rhs, ierr); CHKERRQ(ierr)

  end subroutine timestepper_step_init

!------------------------------------------------------------------------

  subroutine timestepper_step_destroy(self)
    !! Destroys a timestep.

    class(timestepper_step_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call VecDestroy(self%solution,ierr); CHKERRQ(ierr)
    call VecDestroy(self%lhs, ierr); CHKERRQ(ierr)
    call VecDestroy(self%rhs, ierr); CHKERRQ(ierr)

  end subroutine timestepper_step_destroy

!------------------------------------------------------------------------

  character(16) function timestepper_step_status_str(self)
    !! Returns a string representing the timestep status.

    class(timestepper_step_type), intent(in) :: self

    select case (self%status)
       case (TIMESTEP_OK)
          timestepper_step_status_str = 'OK'
       case (TIMESTEP_NOT_CONVERGED)
          timestepper_step_status_str = 'not converged'
       case (TIMESTEP_TOO_SMALL)
          timestepper_step_status_str = 'increase'
       case (TIMESTEP_TOO_BIG)
          timestepper_step_status_str = 'reduce'
       case (TIMESTEP_ABORTED)
          timestepper_step_status_str = 'aborted'
       case (TIMESTEP_FINAL)
          timestepper_step_status_str = 'final'
       case default
          timestepper_step_status_str = 'unknown'
    end select

  end function timestepper_step_status_str

!------------------------------------------------------------------------

  subroutine timestepper_step_print(self, logfile)
    !! Prints a timestep.

    class(timestepper_step_type), intent(in) :: self
    type(logfile_type), intent(in out) :: logfile

    call logfile%write(LOG_LEVEL_INFO, 'timestep', 'end', &
         ['tries'], [self%num_tries], &
         ['size            ', 'time            '], &
         [self%stepsize, self%time], &
         str_key = 'status          ', &
         str_value = self%status_str())

  end subroutine timestepper_step_print

!------------------------------------------------------------------------
!  Timestep adaptor procedures
!------------------------------------------------------------------------

  PetscReal function timestep_adaptor_reduce(self, stepsize)
    !! Reduces stepsize.

    class(timestep_adaptor_type), intent(in) :: self
    PetscReal, intent(in) :: stepsize

    timestep_adaptor_reduce = stepsize * self%reduction

  end function timestep_adaptor_reduce

!------------------------------------------------------------------------

  PetscReal function timestep_adaptor_increase(self, stepsize)
    !! Increases stepsize.

    class(timestep_adaptor_type), intent(in) :: self
    PetscReal, intent(in) :: stepsize

    timestep_adaptor_increase = stepsize * self%amplification
    if (self%max_stepsize > 0._dp) then
       timestep_adaptor_increase  = min(timestep_adaptor_increase, &
            self%max_stepsize)
    end if

  end function timestep_adaptor_increase

!------------------------------------------------------------------------
! Timestepper_checkpoints procedures
!------------------------------------------------------------------------

  subroutine timestepper_checkpoints_init(self, time, repeat, tolerance, &
       start_time)
    !! Initialises checkpoints object.

    class(timestepper_checkpoints_type), intent(in out) :: self
    PetscReal, allocatable, intent(in) :: time(:)
    PetscInt, intent(in) :: repeat
    PetscReal, intent(in) :: tolerance
    PetscReal, intent(in) :: start_time
    ! Locals:
    PetscReal, parameter :: min_checkpoint_tol = 1.e-6_dp ! Don't allow zero tolerance

    if (allocated(time)) then

       self%time = time
       self%repeat = repeat
       self%tolerance = max(tolerance, min_checkpoint_tol)
       self%index = 1
       self%repeat_index = 1
       self%repeat_shift = self%time(size(self%time)) - start_time
       self%next_time = self%time(1)
       self%done = PETSC_FALSE

    else
       self%done = PETSC_TRUE
    end if

    self%hit = PETSC_FALSE

  end subroutine timestepper_checkpoints_init

!------------------------------------------------------------------------

  subroutine timestepper_checkpoints_destroy(self)
    !! Destroys checkpoints object.

    class(timestepper_checkpoints_type), intent(in out) :: self

    if (allocated(self%time)) then
       deallocate(self%time)
    end if

  end subroutine timestepper_checkpoints_destroy

!------------------------------------------------------------------------

  subroutine timestepper_checkpoints_update(self)
    !! Updates checkpoint status after checkpoint output.

    class(timestepper_checkpoints_type), intent(in out) :: self

    if (.not. self%done) then

       self%index = self%index + 1

       if (self%index > size(self%time)) then

          if ((self%repeat > 0) .and. &
               (self%repeat_index >= self%repeat)) then
             self%done = PETSC_TRUE
          else
             self%repeat_index = self%repeat_index + 1
             self%index = 1
          end if

       end if

       if (.not. (self%done)) then
          self%next_time = self%time(self%index) + &
               (self%repeat_index - 1) * self%repeat_shift
       end if

    end if

    self%hit = PETSC_FALSE

  end subroutine timestepper_checkpoints_update

!------------------------------------------------------------------------

  subroutine timestepper_checkpoints_check(self, time, stepsize)
    !! Checks if the next checkpoint will be hit during a timestep
    !! with the specified end time and stepsize.

    class(timestepper_checkpoints_type), intent(in out) :: self
    PetscReal, intent(in) :: time
    PetscReal, intent(in) :: stepsize

    if (self%done) then
       self%hit = PETSC_FALSE
    else
       if (time + self%tolerance * stepsize >= self%next_time) then
          self%hit = PETSC_TRUE
          self%restore_stepsize = stepsize
       else
          self%hit = PETSC_FALSE
       end if
    end if

  end subroutine timestepper_checkpoints_check

!------------------------------------------------------------------------
! Timestepper_steps procedures
!------------------------------------------------------------------------

  subroutine timestepper_steps_init(self, num_stored, &
       time, solution, stop_time_specified, &
       stop_time, max_num_steps, max_stepsize, &
       adapt_on, adapt_method, adapt_min, adapt_max, &
       adapt_reduction, adapt_amplification, step_sizes, &
       nonlinear_solver_relative_tol, nonlinear_solver_abs_tol, &
       nonlinear_solver_update_relative_tol, nonlinear_solver_update_abs_tol, &
       nonlinear_solver_minimum_iterations, &
       max_num_tries, steady_state, &
       checkpoint_time, checkpoint_repeat, checkpoint_tol, &
       stop_minimum_stepsize, stop_maximum_stepsize)
    !! Sets up timestepper steps object from specified parameters.

    class(timestepper_steps_type), intent(in out) :: self
    PetscInt, intent(in) :: num_stored
    PetscReal, intent(in) :: time
    PetscBool, intent(in) :: stop_time_specified
    Vec, intent(in) :: solution
    PetscReal, intent(in) :: stop_time
    PetscInt, intent(in) :: max_num_steps
    PetscReal, intent(in) :: max_stepsize
    PetscBool, intent(in) :: adapt_on
    PetscInt, intent(in) :: adapt_method
    PetscReal, intent(in) :: adapt_min, adapt_max
    PetscReal, intent(in) :: adapt_reduction, adapt_amplification
    PetscReal, intent(in), allocatable :: step_sizes(:)
    PetscReal, intent(in) :: nonlinear_solver_relative_tol, &
         nonlinear_solver_abs_tol
    PetscReal, intent(in) :: nonlinear_solver_update_relative_tol, &
         nonlinear_solver_update_abs_tol
    PetscInt, intent(in) :: nonlinear_solver_minimum_iterations
    PetscInt, intent(in) :: max_num_tries
    PetscBool, intent(in) :: steady_state
    PetscReal, allocatable, intent(in) :: checkpoint_time(:)
    PetscInt, intent(in) :: checkpoint_repeat
    PetscReal, intent(in) :: checkpoint_tol
    PetscReal, intent(in) :: stop_minimum_stepsize, stop_maximum_stepsize
    ! Locals:
    PetscInt :: i
    PetscErrorCode :: ierr

    ! Set up array of timesteps and pointers to them. This array stores the
    ! current step and one or more previous steps. The number of stored
    ! steps depends on the timestepping method (e.g. 2 for single-step methods,
    ! > 2 for multistep methods).
    self%num_stored = num_stored
    allocate(self%store(num_stored), self%pstore(num_stored))
    do i = 1, num_stored
       call self%store(i)%init(solution)
       call self%set_pstore(self%store, i, i)
    end do
    call self%set_aliases()

    self%taken = 0
    self%finished = PETSC_FALSE
    self%current%time = time
    call VecCopy(solution, self%current%solution, ierr); CHKERRQ(ierr)

    self%stop_time_specified = stop_time_specified
    self%stop_time = stop_time
    self%max_num = max_num_steps
    self%max_num_tries = max_num_tries
    self%stop_min_stepsize = stop_minimum_stepsize
    self%stop_max_stepsize = stop_maximum_stepsize

    select case (adapt_method)
    case (TS_ADAPT_CHANGE)
       self%adaptor%monitor => relative_change_monitor
       self%adaptor%name = "change"
    case (TS_ADAPT_ITERATION)
       self%adaptor%monitor => iteration_monitor
       self%adaptor%name = "iteration"
    case default
       self%adaptor%monitor => relative_change_monitor
       self%adaptor%name = "change"
    end select
    self%adaptor%monitor_min = adapt_min
    self%adaptor%monitor_max = adapt_max
    self%adaptor%reduction = adapt_reduction
    self%adaptor%amplification = adapt_amplification
    self%adaptor%max_stepsize = max_stepsize

    self%fixed = .not. adapt_on
    self%adaptor%on = PETSC_FALSE  ! at least until fixed steps are done
    self%fixed_step_index = 1
    self%sizes = step_sizes
    self%next_stepsize = self%sizes(self%fixed_step_index)
    self%steady_state = steady_state

    call self%checkpoints%init(checkpoint_time, checkpoint_repeat, &
            checkpoint_tol, time)

    self%nonlinear_solver_relative_tol = nonlinear_solver_relative_tol
    self%nonlinear_solver_abs_tol = nonlinear_solver_abs_tol
    self%nonlinear_solver_update_relative_tol = nonlinear_solver_update_relative_tol
    self%nonlinear_solver_update_abs_tol = nonlinear_solver_update_abs_tol
    self%nonlinear_solver_minimum_iterations = nonlinear_solver_minimum_iterations

  end subroutine timestepper_steps_init

!------------------------------------------------------------------------

  subroutine timestepper_initial_function_calls(self, err)
    !! Performs initial LHS function call at start of run (and
    !! pre-evaluation procedure if needed). If self%output_initial is
    !! true, an initial RHS function call is also carried out, to
    !! initialise e.g. source parameters which may appear in the
    !! output.

    class(timestepper_type), intent(in out) :: self
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: t, interval(2)

    err = 0
    t = self%steps%current%time
    interval = [t, t]

    call self%ode%pre_solve(t, self%steps%current%solution, err = err)

    if (err == 0) then
       call self%ode%lhs(t, interval, self%steps%current%solution, &
            self%steps%current%lhs, err)
       if ((err == 0) .and. (self%output_initial)) then
          call self%ode%rhs(t, interval, self%steps%current%solution, &
               self%steps%current%rhs, err)
       end if
    end if

  end subroutine timestepper_initial_function_calls

!------------------------------------------------------------------------

  subroutine timestepper_steps_destroy(self)
    !! Destroys store array and de-assigns pointers to them.

    class(timestepper_steps_type), intent(in out) :: self
    ! Locals:
    PetscInt :: i

    self%current => NULL()
    self%last => NULL()

    do i = 1, self%num_stored
       self%pstore(i)%p => NULL()
       call self%store(i)%destroy()
    end do

    deallocate(self%store, self%pstore)
    if (allocated(self%sizes)) then
       deallocate(self%sizes)
    end if

    call self%checkpoints%destroy()

  end subroutine timestepper_steps_destroy

!------------------------------------------------------------------------

  subroutine timestepper_steps_set_aliases(self)
    !! Sets current and last pointers to store array.

    class(timestepper_steps_type), intent(in out) :: self

    self%current => self%pstore(1)%p
    if (self%num_stored > 1) then
       self%last => self%pstore(2)%p
    else
       self%last => self%current
    end if

  end subroutine timestepper_steps_set_aliases

!------------------------------------------------------------------------

  subroutine timestepper_steps_set_pstore(self, store, i, j)
    !! Sets pointer i to element j of steps array.
    !! The store array is passed in so it can be given the target
    !! attribute.

    class(timestepper_steps_type), intent(in out) :: self
    type(timestepper_step_type), target :: store(:)
    PetscInt, intent(in) :: i, j

    self%pstore(i)%p => store(j)

  end subroutine timestepper_steps_set_pstore

!------------------------------------------------------------------------

  subroutine timestepper_steps_initialize_try(self)
    !! Initializes step try.

    class(timestepper_steps_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    if (self%num_stored > 1) then
       ! Last solution is initial guess for current solution:
       call VecCopy(self%last%solution, self%current%solution, ierr); CHKERRQ(ierr)
    end if

    self%current%stepsize = self%next_stepsize
    self%current%time = self%last%time + self%current%stepsize

  end subroutine timestepper_steps_initialize_try

!------------------------------------------------------------------------

  subroutine timestepper_steps_update(self)
    !! Updates pointers to store array, when timestep is advanced.

    class(timestepper_steps_type), intent(in out) :: self

    call self%rotate()
    call self%set_aliases()

    self%current%num_tries = 0
    self%current%status = TIMESTEP_OK

  end subroutine timestepper_steps_update

!------------------------------------------------------------------------

  subroutine timestepper_steps_rotate(self)
    !! Rotates pstore pointers to store array, so that e.g. current becomes
    !! last, etc, when starting a new time step.

    class(timestepper_steps_type), intent(in out) :: self
    ! Locals:
    PetscInt :: i
    type(ptimestepper_step_type) :: plast

    plast%p => self%pstore(self%num_stored)%p

    do i = self%num_stored, 2, -1
       self%pstore(i)%p => self%pstore(i-1)%p
    end do
    self%pstore(1)%p => plast%p
    nullify(plast%p)

  end subroutine timestepper_steps_rotate

!------------------------------------------------------------------------

  subroutine timestepper_steps_check_finished(self, logfile)
    !! Checks if any termination criteria have been met, and reduces
    !! timestep if needed.

    class(timestepper_steps_type), intent(in out) :: self
    type(logfile_type), intent(in out) :: logfile
    
    self%finished = PETSC_FALSE

    if (self%steady_state) then
       self%finished = (self%taken == 1)
    else
       if (self%stop_time_specified .and. &
            (self%current%time + self%termination_tol * &
            self%current%stepsize > self%stop_time)) then
          self%current%stepsize = self%stop_time - self%last%time
          self%current%time = self%stop_time
          self%finished = PETSC_TRUE
          call logfile%write(LOG_LEVEL_INFO, 'timestep', 'end_time_reached', &
               real_keys = ['size'], real_values = [self%current%stepsize])
       else
          if ((self%stop_min_stepsize > 0._dp) .and. &
               self%current%stepsize <= self%stop_min_stepsize) then
             self%current%stepsize = self%stop_min_stepsize
             self%finished = PETSC_TRUE
             call logfile%write(LOG_LEVEL_INFO, 'timestep', 'stop_size_minimum_reached', &
                  real_keys = ['size'], real_values = [self%current%stepsize])
          else if ((self%stop_max_stepsize > 0._dp) .and. &
               self%current%stepsize >= self%stop_max_stepsize) then
             self%current%stepsize = self%stop_max_stepsize
             self%finished = PETSC_TRUE
             call logfile%write(LOG_LEVEL_INFO, 'timestep', 'stop_size_maximum_reached', &
                  real_keys = ['size'], real_values = [self%current%stepsize])
          end if
       end if
       if (self%taken + 1 >= self%max_num) then
          self%finished = PETSC_TRUE
          call logfile%write(LOG_LEVEL_INFO, 'timestep', 'max_timesteps_reached')
       end if
    end if

  end subroutine timestepper_steps_check_finished

!------------------------------------------------------------------------

  subroutine timestepper_steps_check_checkpoints(self, logfile)
    !! Checks if a checkpoint is being passed, and reduces
    !! timestep if needed.

    class(timestepper_steps_type), intent(in out) :: self
    type(logfile_type), intent(in out) :: logfile

    if (.not. (self%steady_state)) then

       call self%checkpoints%check(self%current%time, self%current%stepsize)

       if (self%checkpoints%hit) then

          self%current%stepsize = self%checkpoints%next_time - self%last%time
          self%current%time = self%checkpoints%next_time
          call logfile%write(LOG_LEVEL_INFO, 'timestep', 'checkpoint_time_reached', &
               real_keys = ['time'], real_values = [self%current%time])

       end if

    end if

  end subroutine timestepper_steps_check_checkpoints

!------------------------------------------------------------------------

  subroutine timestepper_steps_set_current_status(self, converged_reason)
    !! Sets status of current step.

    class(timestepper_steps_type), intent(in out) :: self
    SNESConvergedReason, intent(in) :: converged_reason
    ! Locals:
    PetscReal :: eta

    if (self%steady_state) then
       if (converged_reason >= 0) then
          self%current%status = TIMESTEP_FINAL
       else
          self%current%status = TIMESTEP_ABORTED
       end if
       self%finished = PETSC_TRUE
    else
       if (converged_reason >= 0 ) then
          if ((self%finished) .and. (self%current%status /= &
               TIMESTEP_ABORTED)) then
             self%current%status = TIMESTEP_FINAL
          else
             eta = self%adaptor%monitor(self%current, self%last)
             if ((self%adaptor%on) .or. &
                  (self%fixed_step_index == size(self%sizes)) .and. &
                  (.not. self%fixed)) then
                if (eta < self%adaptor%monitor_min) then
                   self%current%status = TIMESTEP_TOO_SMALL
                else if (eta > self%adaptor%monitor_max) then
                   self%current%status = TIMESTEP_TOO_BIG
                else
                   self%current%status = TIMESTEP_OK
                end if
             else
                self%current%status = TIMESTEP_OK
             end if
          end if
       else if (self%current%num_tries >= self%max_num_tries) then
          self%current%status = TIMESTEP_ABORTED
          self%finished = PETSC_TRUE
       else
          self%current%status = TIMESTEP_NOT_CONVERGED
          self%finished = PETSC_FALSE
       end if
    end if
    
  end subroutine timestepper_steps_set_current_status

!------------------------------------------------------------------------

  subroutine timestepper_steps_get_next_fixed_stepsize(self, accepted)
    !! Gets next timestep size if using fixed time step sizes.

    class(timestepper_steps_type), intent(in out) :: self
    PetscBool, intent(out) :: accepted

    associate (fixed_index => self%fixed_step_index, &
         num_sizes => size(self%sizes))

      accepted = PETSC_TRUE
      fixed_index = fixed_index + 1
      if (fixed_index <= num_sizes) then
         self%next_stepsize = self%sizes(fixed_index)
      else
         fixed_index = num_sizes
         if (self%fixed) then
            self%next_stepsize = self%sizes(fixed_index)
         else
            self%adaptor%on = PETSC_TRUE
            call self%adapt(accepted)
         end if
      end if

    end associate

  end subroutine timestepper_steps_get_next_fixed_stepsize

!------------------------------------------------------------------------

  subroutine timestepper_steps_set_next_stepsize(self, accepted)

    class(timestepper_steps_type), intent(in out) :: self
    PetscBool, intent(out) :: accepted

    if (.not.(self%steady_state)) then

       if (self%adaptor%on) then

          call self%adapt(accepted)

          if ((self%fixed_step_index < size(self%sizes)) .or. &
               ((self%fixed_step_index >= size(self%sizes)) .and. self%fixed)) then
             if (self%next_stepsize >= self%sizes(self%fixed_step_index)) then
                ! switch back to fixed time stepping:
                self%adaptor%on = PETSC_FALSE
                self%next_stepsize = self%sizes(self%fixed_step_index)
             end if
          end if

       else

          select case (self%current%status)
          case (TIMESTEP_TOO_BIG, TIMESTEP_NOT_CONVERGED)
             ! temporarily switch to adaptive time stepping:
             self%adaptor%on = PETSC_TRUE
             call self%adapt(accepted)
          case default
             call self%get_next_fixed_stepsize(accepted)
          end select

       end if
    end if

end subroutine timestepper_steps_set_next_stepsize

!------------------------------------------------------------------------

  subroutine timestepper_steps_adapt(self, accepted)
    !! Checks current time stepsize for acceptability and adapts if 
    !! necessary.

    class(timestepper_steps_type), intent(in out) :: self
    PetscBool, intent(out) :: accepted

    select case (self%current%status)
    case (TIMESTEP_TOO_SMALL)
       accepted = PETSC_TRUE
       self%next_stepsize = self%adaptor%increase(self%current%stepsize)
    case (TIMESTEP_TOO_BIG, TIMESTEP_NOT_CONVERGED)
       accepted = PETSC_FALSE
       self%next_stepsize = self%adaptor%reduce(self%current%stepsize)
    case default
       accepted = PETSC_TRUE
       self%next_stepsize = self%current%stepsize
    end select

  end subroutine timestepper_steps_adapt

!------------------------------------------------------------------------
! Timestepper method procedures
!------------------------------------------------------------------------

  subroutine timestepper_method_init(self, method)
    !! Sets the timestepping method.

    class(timestepper_method_type), intent(in out) :: self
    PetscInt, intent(in) :: method

    select case (method)
    case (TS_BDF2)
       self%residual => BDF2_residual
       self%num_steps = 2
       self%name = "BDF2"
    case (TS_DIRECTSS)
       self%residual => direct_ss_residual
       self%num_steps = 0
       self%name = "Direct steady state"
    case default
       self%residual => backwards_Euler_residual
       self%num_steps = 1
       self%name = "Backward Euler"
    end select

    ! Store current as well as previous steps:
    self%num_stored_steps = self%num_steps + 1

  end subroutine timestepper_method_init

!------------------------------------------------------------------------
! Timestepper solver context procedures
! Note these cannot be bound to the solver context type, as the context
! variable goes through PETSc function interfaces where it is declared as
! type(*), which does not allow type-bound procedures.
! ------------------------------------------------------------------------

  subroutine timestepper_solver_context_init(self, ode, steps, residual)
    !! Sets up timestepper SNES solver context.

    class(timestepper_solver_context_type), intent(in out) :: self
    class(ode_type), intent(in), target :: ode
    type(timestepper_steps_type), target :: steps
    procedure(method_residual) :: residual

    self%ode => ode
    self%steps => steps
    self%residual => residual

  end subroutine timestepper_solver_context_init

!------------------------------------------------------------------------

  subroutine timestepper_solver_context_destroy(self)
    !! Destroys timestepper SNES solver context.

    class(timestepper_solver_context_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    nullify(self%ode)
    nullify(self%steps)
    nullify(self%residual)
    call MatFDColoringDestroy(self%fd_coloring, ierr); CHKERRQ(ierr)

  end subroutine timestepper_solver_context_destroy

!------------------------------------------------------------------------
! Timestepper procedures
!------------------------------------------------------------------------

  subroutine timestepper_setup_nonlinear_solver(self, max_iterations)
    !! Sets up SNES nonlinear solver for the timestepper.

    class(timestepper_type), intent(in out) :: self
    PetscInt, intent(in) :: max_iterations
    ! Locals:
    PetscErrorCode :: ierr
    MatColoring :: matrix_coloring
    ISColoring :: is_coloring
    SNESLineSearch :: linesearch
    ! This tolerance needs to be set very small so it doesn't override
    ! time step reduction when primary variables go out of bounds:
    PetscReal, parameter :: stol = 1.e-99_dp
    PetscReal, parameter :: dtol = 1.e8_dp
    PetscReal, parameter :: default_mat_fd_err = 1.e-8_dp
    PetscReal, parameter :: default_mat_fd_umin = 1.e-2_dp

    call timestepper_solver_context_init(self%context, self%ode, self%steps, &
         self%method%residual)

    call SNESCreate(PETSC_COMM_WORLD, self%solver, ierr); CHKERRQ(ierr)
    call SNESSetApplicationContext(self%solver, self%context, ierr); CHKERRQ(ierr)
    call SNESSetFunction(self%solver, self%residual, &
         SNES_residual, self%context, ierr); CHKERRQ(ierr)
    call SNESSetDM(self%solver, self%ode%mesh%dm, ierr); CHKERRQ(ierr)

    call MatColoringCreate(self%jacobian, matrix_coloring, ierr); CHKERRQ(ierr)
    call MatColoringSetFromOptions(matrix_coloring, ierr); CHKERRQ(ierr)
    call MatColoringApply(matrix_coloring, is_coloring, ierr); CHKERRQ(ierr)
    call MatColoringDestroy(matrix_coloring, ierr); CHKERRQ(ierr)
    call MatFDColoringCreate(self%jacobian, is_coloring, self%context%fd_coloring, &
         ierr); CHKERRQ(ierr)
    call MatFDColoringSetFunction(self%context%fd_coloring, SNES_residual, &
         self%context, ierr); CHKERRQ(ierr)
    call MatFDColoringSetType(self%context%fd_coloring, MATMFFD_DS, ierr)
    CHKERRQ(ierr)
    call MatFDColoringSetParameters(self%context%fd_coloring, &
         default_mat_fd_err, default_mat_fd_umin, ierr); CHKERRQ(ierr)
    call MatFDColoringSetFromOptions(self%context%fd_coloring, ierr); CHKERRQ(ierr)
    call MatFDColoringSetUp(self%jacobian, is_coloring, self%context%fd_coloring, &
         ierr); CHKERRQ(ierr)
    call ISColoringDestroy(is_coloring, ierr); CHKERRQ(ierr)

    call SNESSetJacobian(self%solver, self%jacobian, self%jacobian, &
         SNES_Jacobian, self%context, ierr)
    CHKERRQ(ierr)

    call SNESSetTolerances(self%solver, PETSC_DEFAULT_REAL, &
         PETSC_DEFAULT_REAL, stol, max_iterations, &
         PETSC_DEFAULT_INTEGER, ierr); CHKERRQ(ierr)
    call SNESSetDivergenceTolerance(self%solver, dtol, ierr); CHKERRQ(ierr)
    call SNESSetConvergenceTest(self%solver, SNES_convergence, self%context, &
         PETSC_NULL_FUNCTION, ierr); CHKERRQ(ierr)
    call SNESMonitorSet(self%solver, SNES_monitor, self%context, &
         PETSC_NULL_FUNCTION, ierr); CHKERRQ(ierr)
    call SNESSetFromOptions(self%solver, ierr); CHKERRQ(ierr)

    ! Set function to be called at start of each solver iteration:
    call SNESSetUpdate(self%solver, SNES_pre_iteration_update, ierr)
    CHKERRQ(ierr)

    ! Nonlinear solver line search:
    call SNESGetLineSearch(self%solver, linesearch, ierr); CHKERRQ(ierr)
    call SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC, ierr)
    CHKERRQ(ierr)
    call SNESLineSearchSetPostCheck(linesearch, SNES_linesearch_post_check, &
         self%context, ierr); CHKERRQ(ierr)

  end subroutine timestepper_setup_nonlinear_solver

!------------------------------------------------------------------------

  subroutine timestepper_setup_linear_solver(self, ksp_type, &
       relative_tolerance, max_iterations, pc_type)
    !! Sets up KSP linear solver and PC preconditioner for the timestepper.

    class(timestepper_type), intent(in out) :: self
    KSPType, intent(in) :: ksp_type
    PetscReal, intent(in) :: relative_tolerance
    PetscInt, intent(in) :: max_iterations
    PCType, intent(in) :: pc_type
    ! Locals:
    KSP :: ksp
    PC :: pc
    PetscErrorCode :: ierr

    call SNESGetKSP(self%solver, ksp, ierr); CHKERRQ(ierr)
    call KSPSetType(ksp, ksp_type, ierr); CHKERRQ(ierr)
    call KSPSetTolerances(ksp, relative_tolerance, PETSC_DEFAULT_REAL, &
         PETSC_DEFAULT_REAL, max_iterations, ierr); CHKERRQ(ierr)

    call KSPGetPC(ksp, pc, ierr); CHKERRQ(ierr)
    call PCSetType(pc, pc_type, ierr); CHKERRQ(ierr)
    call PCSetFromOptions(pc, ierr); CHKERRQ(ierr)

  end subroutine timestepper_setup_linear_solver

!------------------------------------------------------------------------

  subroutine timestepper_setup_linear_solver_options(self, json, ksp_type)
    !! Sets up KSP linear solver options.

    use fson
    use fson_mpi_module

    class(timestepper_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    KSPType, intent(in) :: ksp_type
    ! Locals:
    KSP :: ksp
    PetscErrorCode :: ierr
    type(fson_value), pointer :: opt_json
    PetscInt :: restart

    call SNESGetKSP(self%solver, ksp, ierr); CHKERRQ(ierr)

    if (((ksp_type == KSPGMRES) .or. (ksp_type == KSPLGMRES)) .and. &
         fson_has_mpi(json, 'time.step.solver.linear.options.gmres')) then
       call fson_get_mpi(json, 'time.step.solver.linear.options.gmres', opt_json)

       if (fson_has_mpi(opt_json, 'restart')) then
          call fson_get_mpi(opt_json, 'restart', val = restart)
          call KSPGMRESSetRestart(ksp, restart, ierr); CHKERRQ(ierr)
       end if

    end if

    call KSPSetFromOptions(ksp, ierr); CHKERRQ(ierr)

  end subroutine timestepper_setup_linear_solver_options

!------------------------------------------------------------------------

  subroutine timestepper_setup_linear_sub_solver(self, sub_pc_type, &
       factor_levels)
    !! Sets up linear sub-solver for block Jacobi and ASM preconditioners.

    class(timestepper_type), intent(in out) :: self
    PCType, intent(in) :: sub_pc_type
    PetscInt, intent(in) :: factor_levels
    ! Locals:
    KSP :: ksp
    PC :: pc, sub_pc
    PCType :: pc_type
    PetscInt :: num_local, first_local, i
    KSP, allocatable :: sub_ksp(:)
    PetscErrorCode :: ierr

    call SNESGetKSP(self%solver, ksp, ierr); CHKERRQ(ierr)
    call KSPSetUp(ksp, ierr); CHKERRQ(ierr)
    call KSPGetPC(ksp, pc, ierr); CHKERRQ(ierr)
    call PCGetType(pc, pc_type, ierr); CHKERRQ(ierr)

    num_local = 0

    select case (pc_type)
    case (PCBJACOBI)
       call PCBJacobiGetSubKSP(pc, num_local, first_local, &
            PETSC_NULL_KSP, ierr); CHKERRQ(ierr)
       allocate(sub_ksp(num_local))
       call PCBJacobiGetSubKSP(pc, num_local, first_local, &
            sub_ksp, ierr); CHKERRQ(ierr)
    case (PCASM)
       call PCASMGetSubKSP(pc, num_local, first_local, &
            PETSC_NULL_KSP, ierr); CHKERRQ(ierr)
       allocate(sub_ksp(num_local))
       call PCASMGetSubKSP(pc, num_local, first_local, &
            sub_ksp, ierr); CHKERRQ(ierr)
    end select

    if (allocated(sub_ksp)) then
       do i = 1, num_local
          call KSPGetPC(sub_ksp(i), sub_pc, ierr); CHKERRQ(ierr)
          call PCSetType(sub_pc, sub_pc_type, ierr); CHKERRQ(ierr)
          call PCFactorSetLevels(sub_pc, factor_levels, ierr)
          CHKERRQ(ierr)
          call PCSetFromOptions(sub_pc, ierr); CHKERRQ(ierr)
       end do
       deallocate(sub_ksp)
    end if

  end subroutine timestepper_setup_linear_sub_solver

!------------------------------------------------------------------------

  subroutine SNES_monitor(solver, num_iterations, fnorm, context, ierr)
    !! SNES monitor routine for output at each nonlinear iteration.

    SNES, intent(in) :: solver
    PetscInt, intent(in) :: num_iterations
    PetscReal, intent(in) :: fnorm
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode :: ierr
    ! Locals:
    PetscInt :: natural, minc_level, num_linear_solver_iterations
    character(len = 24), allocatable :: cell_keys(:)
    PetscInt, allocatable :: cell_values(:)
    KSP :: ksp

    if (num_iterations > 0) then
       call SNESGetKSP(solver, ksp, ierr); CHKERRQ(ierr)
       call KSPGetIterationNumber(ksp, num_linear_solver_iterations, ierr)
       CHKERRQ(ierr)
       associate(mesh => context%ode%mesh)
         call mesh%global_to_parent_natural( &
              context%steps%current%max_residual_cell, &
              natural, minc_level)
         call mesh%natural_cell_output_arrays( &
              natural, minc_level, cell_keys, cell_values)
         call context%ode%logfile%write(LOG_LEVEL_INFO, &
              'nonlinear_solver', 'iteration', &
              [['count                   ', 'linear_solver_iterations'], &
              cell_keys, ['equation                ']], &
              [[num_iterations, num_linear_solver_iterations], cell_values, &
              [context%steps%current%max_residual_equation]], &
              ['residual'], [context%steps%current%max_residual])
       end associate
    end if
    call context%ode%logfile%flush()

  end subroutine SNES_monitor

!------------------------------------------------------------------------

  subroutine SNES_convergence(solver, num_iterations, xnorm, pnorm, &
       fnorm, reason, context, ierr)
    !! Tests for convergence of nonlinear solver.

    use dm_utils_module, only: vec_max_pointwise_abs_scale

    SNES, intent(in) :: solver
    PetscInt, intent(in) :: num_iterations
    PetscReal, intent(in) :: xnorm, pnorm, fnorm
    SNESConvergedReason, intent(out) :: reason
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode, intent(out) :: ierr
    ! Locals:
    Vec :: residual, solution, update
    PetscReal :: max_update
    PetscInt :: index, bs

    call SNESGetFunction(solver, residual, PETSC_NULL_FUNCTION, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)

    call vec_max_pointwise_abs_scale(residual, &
         context%steps%last%lhs, &
         context%steps%nonlinear_solver_abs_tol, &
         context%steps%current%max_residual, index)
    call VecGetBlockSize(residual, bs, ierr); CHKERRQ(ierr)
    context%steps%current%max_residual_cell = int(index / bs)
    context%steps%current%max_residual_equation = 1 + index - &
         context%steps%current%max_residual_cell * bs

    call SNESConvergedDefault(solver, num_iterations, xnorm, pnorm, &
         fnorm, reason, 0, ierr); CHKERRQ(ierr)

    if (num_iterations < &
         context%steps%nonlinear_solver_minimum_iterations) then
       reason = SNES_CONVERGED_ITERATING
    else
       if (context%steps%current%max_residual < &
            context%steps%nonlinear_solver_relative_tol) then
          reason = SNES_CONVERGED_FNORM_RELATIVE
       else
          if (num_iterations > 0) then
             call SNESGetSolutionUpdate(solver, update, ierr); CHKERRQ(ierr)
             call SNESGetSolution(solver, solution, ierr); CHKERRQ(ierr)
             call vec_max_pointwise_abs_scale(update, solution, &
                  context%steps%nonlinear_solver_update_abs_tol, &
                  max_update, index)
             if (max_update <= context%steps%nonlinear_solver_update_relative_tol) then
                reason = SNES_CONVERGED_SNORM_RELATIVE
             end if
          end if
       end if
    end if

  end subroutine SNES_convergence

!------------------------------------------------------------------------

  subroutine timestepper_setup_jacobian(self)
    !! Sets up Jacobian matrix for nonlinear solver.

    class(timestepper_type), intent(in out) :: self
    ! Locals:
    PetscInt :: blocksize
    MatType :: mat_type
    PetscErrorCode :: ierr

    call VecGetBlockSize(self%ode%solution, blocksize, ierr); CHKERRQ(ierr)
    if (blocksize == 1) then
       mat_type = MATAIJ
    else
       mat_type = MATBAIJ
    end if
    call DMSetMatType(self%ode%mesh%dm, mat_type, ierr); CHKERRQ(ierr)
    call DMCreateMatrix(self%ode%mesh%dm, self%jacobian, ierr); CHKERRQ(ierr)
    call MatSetFromOptions(self%jacobian, ierr); CHKERRQ(ierr)

  end subroutine timestepper_setup_jacobian

!------------------------------------------------------------------------

  subroutine timestepper_init(self, json, ode)
    !! Initializes a timestepper.

    use utils_module, only : str_to_lower
    use fson
    use fson_mpi_module
    use fson_value_m, only : TYPE_NULL, TYPE_REAL, TYPE_INTEGER, &
         TYPE_LOGICAL, TYPE_ARRAY

    class(timestepper_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    class(ode_type), intent(in), target :: ode
    ! Locals:
    PetscInt :: max_num_steps
    PetscReal, parameter :: default_stop_time = 1.0_dp
    PetscReal, parameter :: default_stepsize = 0.1_dp
    PetscInt :: step_size_type
    PetscReal :: step_size_single
    PetscReal, allocatable :: step_sizes(:)
    PetscReal :: max_stepsize, stop_time
    PetscReal, parameter :: default_max_stepsize = 0.0_dp
    PetscReal, parameter :: default_stop_minimum_stepsize = -1.0_dp
    PetscReal, parameter :: default_stop_maximum_stepsize = -1.0_dp
    PetscReal :: stop_minimum_stepsize, stop_maximum_stepsize
    PetscInt :: method
    PetscInt, parameter :: max_method_str_len = 12
    character(max_method_str_len) :: method_str
    character(max_method_str_len), parameter :: default_method_str = "beuler"
    PetscInt, parameter :: default_max_num_steps = 100
    PetscBool, parameter :: default_adapt_on = PETSC_FALSE
    PetscBool :: adapt_on
    PetscInt, parameter :: max_adapt_method_str_len = 12
    character(max_adapt_method_str_len) :: adapt_method_str
    character(max_adapt_method_str_len), parameter :: &
         default_adapt_method_str = "iteration"
    PetscInt :: adapt_method
    PetscReal, parameter :: default_adapt_min = 5._dp, &
         default_adapt_max = 8._dp
    PetscReal :: adapt_min, adapt_max
    PetscReal, parameter :: default_adapt_reduction = 0.2_dp, &
         default_adapt_amplification = 2.0_dp
    PetscReal :: adapt_reduction, adapt_amplification
    PetscInt, parameter :: default_nonlinear_max_iterations = 8
    PetscReal, parameter :: default_nonlinear_relative_tol = 1.e-5_dp
    PetscReal, parameter :: default_nonlinear_abs_tol = 1._dp
    PetscReal, parameter :: default_nonlinear_update_relative_tol = 1.e-10_dp
    PetscReal, parameter :: default_nonlinear_update_abs_tol = 1._dp
    PetscInt, parameter :: default_nonlinear_min_iterations = 0
    PetscInt :: nonlinear_max_iterations, nonlinear_min_iterations
    PetscReal :: nonlinear_relative_tol, nonlinear_abs_tol
    PetscReal :: nonlinear_update_relative_tol, nonlinear_update_abs_tol
    PetscInt :: max_num_tries
    PetscInt, parameter :: default_max_num_tries = 10
    PetscInt, parameter :: default_output_frequency = 1
    PetscBool, parameter :: default_output_initial = PETSC_TRUE
    PetscBool, parameter :: default_output_final = PETSC_TRUE
    PetscBool :: stop_time_specified, steady_state
    character(max_ksp_type_str_len) :: ksp_type_str
    character(max_ksp_type_str_len), parameter :: default_ksp_type_str = "bcgs"
    KSPType :: ksp_type
    character(max_pc_type_str_len) :: pc_type_str, sub_pc_type_str
    character(max_pc_type_str_len), parameter :: default_pc_type_str = "asm"
    character(max_pc_type_str_len), parameter :: default_sub_pc_type_str = "ilu"
    PetscInt, parameter :: default_sub_pc_factor_levels = 0
    PCType :: pc_type, sub_pc_type
    PetscReal :: linear_relative_tol
    PetscInt :: linear_max_iterations, sub_pc_factor_levels
    PetscInt, parameter :: default_checkpoint_repeat = 1
    PetscInt :: checkpoint_repeat_type, checkpoint_repeat, i
    PetscBool :: checkpoint_do_repeat
    PetscReal, allocatable :: checkpoint_time(:), checkpoint_step(:)
    PetscReal :: checkpoint_tol
    PetscReal, parameter :: default_checkpoint_tol = 0.1_dp
    PetscErrorCode :: ierr

    self%ode => ode
    call VecDuplicate(self%ode%solution, self%residual, ierr); CHKERRQ(ierr)
    call self%setup_jacobian()
       
    call fson_get_mpi(json, "time.step.method", &
         default_method_str, method_str, self%ode%logfile)
    method = time_step_method_from_str(method_str)

    call self%method%init(method)

    call fson_get_mpi(json, &
         "time.step.solver.nonlinear.maximum.iterations", &
         default_nonlinear_max_iterations, &
         nonlinear_max_iterations, self%ode%logfile)
    call fson_get_mpi(json, &
         "time.step.solver.nonlinear.tolerance.function.relative", &
         default_nonlinear_relative_tol, &
         nonlinear_relative_tol, self%ode%logfile)
    call fson_get_mpi(json, &
         "time.step.solver.nonlinear.tolerance.function.absolute", &
         default_nonlinear_abs_tol, nonlinear_abs_tol, &
         self%ode%logfile)
    call fson_get_mpi(json, &
         "time.step.solver.nonlinear.tolerance.update.relative", &
         default_nonlinear_update_relative_tol, &
         nonlinear_update_relative_tol, self%ode%logfile)
    call fson_get_mpi(json, &
         "time.step.solver.nonlinear.tolerance.update.absolute", &
         default_nonlinear_update_abs_tol, nonlinear_update_abs_tol, &
         self%ode%logfile)
    call fson_get_mpi(json, &
         "time.step.solver.nonlinear.minimum.iterations", &
         default_nonlinear_min_iterations, &
         nonlinear_min_iterations, self%ode%logfile)
    call fson_get_mpi(json, &
         "time.step.solver.linear.type", &
         default_ksp_type_str, ksp_type_str, self%ode%logfile)
    if (fson_has_mpi(json, &
         "time.step.solver.linear.tolerance.relative")) then
       call fson_get_mpi(json, &
         "time.step.solver.linear.tolerance.relative", &
         val = linear_relative_tol)
    else
       linear_relative_tol = PETSC_DEFAULT_REAL
       call self%ode%logfile%write(LOG_LEVEL_INFO, 'input', &
            'default', str_key = &
            'time.step.solver.linear.tolerance.relative', &
            str_value = "PETSc_default")
    end if
    if (fson_has_mpi(json, &
         "time.step.solver.linear.maximum.iterations")) then
       call fson_get_mpi(json, &
         "time.step.solver.linear.maximum.iterations", &
         val = linear_max_iterations)
    else
       linear_max_iterations = PETSC_DEFAULT_INTEGER
       call self%ode%logfile%write(LOG_LEVEL_INFO, 'input', &
            'default', str_key = &
            'time.step.solver.linear.maximum.iterations', &
            str_value = "PETSc_default")
    end if
    ksp_type = ksp_type_from_str(ksp_type_str)
    call fson_get_mpi(json, &
         "time.step.solver.linear.preconditioner.type", &
         default_pc_type_str, pc_type_str, self%ode%logfile)
    pc_type = pc_type_from_str(pc_type_str)

    call self%setup_nonlinear_solver(nonlinear_max_iterations)
    call self%setup_linear_solver(ksp_type, linear_relative_tol, &
         linear_max_iterations, pc_type)
    call self%setup_linear_solver_options(json, ksp_type)

    if ((pc_type == PCBJACOBI) .or. (pc_type == PCASM)) then
       call fson_get_mpi(json, &
            "time.step.solver.linear.preconditioner.sub.preconditioner.type", &
            default_sub_pc_type_str, sub_pc_type_str, self%ode%logfile)
       sub_pc_type = pc_type_from_str(sub_pc_type_str)
       call fson_get_mpi(json, &
            "time.step.solver.linear.preconditioner.sub.preconditioner.factor.levels", &
            default_sub_pc_factor_levels, sub_pc_factor_levels, self%ode%logfile)
       call self%setup_linear_sub_solver(sub_pc_type, sub_pc_factor_levels)
    end if

    if (method == TS_DIRECTSS) then

       steady_state = PETSC_TRUE
       step_sizes = [0._dp]
       max_num_steps = 1
       max_num_tries = 1
       adapt_on = PETSC_FALSE
       stop_time = 0._dp

    else ! transient

       steady_state = PETSC_FALSE

       if (fson_has_mpi(json, "time.stop")) then
          if (fson_type_mpi(json, "time.stop") == TYPE_NULL) then
             stop_time_specified = PETSC_FALSE
          else
             stop_time_specified = PETSC_TRUE
             call fson_get_mpi(json, "time.stop", default_stop_time, &
                  stop_time, self%ode%logfile)
          end if
       else
          stop_time_specified = PETSC_FALSE
       end if

       if (.not. stop_time_specified) then
          call self%ode%logfile%write(LOG_LEVEL_WARN, 'input', &
               '"time.stop not specified"')
       end if

       if (fson_has_mpi(json, "time.step.maximum.size")) then
          if (fson_type_mpi(json, "time.step.maximum.size") == TYPE_NULL) then
             max_stepsize = default_max_stepsize
          else
             call fson_get_mpi(json, "time.step.maximum.size", &
                  val = max_stepsize)
          end if
       else
          call fson_get_mpi(json, "time.step.maximum.size", &
               default_max_stepsize, max_stepsize, self%ode%logfile)
       end if
       if (fson_has_mpi(json, "time.step.stop.size.minimum")) then
          if (fson_type_mpi(json, "time.step.stop.size.minimum") == TYPE_NULL) then
             stop_minimum_stepsize = default_stop_minimum_stepsize
          else
             call fson_get_mpi(json, "time.step.stop.size.minimum", &
                  val = stop_minimum_stepsize)
          end if
       else
          stop_minimum_stepsize = default_stop_minimum_stepsize
          call self%ode%logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  str_key = "time.step.stop.size.minimum", &
                  str_value = "null")
       end if
       if (fson_has_mpi(json, "time.step.stop.size.maximum")) then
          if (fson_type_mpi(json, "time.step.stop.size.maximum") == TYPE_NULL) then
             stop_maximum_stepsize = default_stop_maximum_stepsize
          else
             call fson_get_mpi(json, "time.step.stop.size.maximum", &
                  val = stop_maximum_stepsize)
          end if
       else
          stop_maximum_stepsize = default_stop_maximum_stepsize
          call self%ode%logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  str_key = "time.step.stop.size.maximum", &
                  str_value = "null")
       end if
       call fson_get_mpi(json, "time.step.maximum.number", &
            default_max_num_steps, max_num_steps, self%ode%logfile)
       call fson_get_mpi(json, "time.step.maximum.tries", &
            default_max_num_tries, max_num_tries, self%ode%logfile)

       call fson_get_mpi(json, "time.step.adapt.on", &
            default_adapt_on, adapt_on, self%ode%logfile)

       call fson_get_mpi(json, "time.step.adapt.method", &
            default_adapt_method_str, adapt_method_str, self%ode%logfile)
       adapt_method = adapt_method_from_str(adapt_method_str)

       call fson_get_mpi(json, "time.step.adapt.minimum", &
            default_adapt_min, adapt_min, self%ode%logfile)
       call fson_get_mpi(json, "time.step.adapt.maximum", &
            default_adapt_max, adapt_max, self%ode%logfile)

       call fson_get_mpi(json, "time.step.adapt.reduction", &
            default_adapt_reduction, adapt_reduction, self%ode%logfile)
       call fson_get_mpi(json, "time.step.adapt.amplification", &
            default_adapt_amplification, adapt_amplification, self%ode%logfile)

       if (fson_has_mpi(json, "time.step.size")) then
          step_size_type = fson_type_mpi(json, "time.step.size")
          select case (step_size_type)
          case (TYPE_REAL, TYPE_INTEGER)
             call fson_get_mpi(json, "time.step.size", val = step_size_single)
             step_sizes = [step_size_single]
          case (TYPE_ARRAY)
             call fson_get_mpi(json, "time.step.size", val = step_sizes)
          case default
             step_sizes = [default_stepsize]
          end select
       else
          step_sizes = [default_stepsize]
          call self%ode%logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
               real_keys = ["time.step.size"], real_values = [default_stepsize])
       end if

       call fson_get_mpi(json, "output.frequency", &
            default_output_frequency, self%output_frequency, self%ode%logfile)

       if (fson_has_mpi(json, "output.checkpoint")) then
          if (fson_has_mpi(json, "output.checkpoint.repeat")) then
             checkpoint_repeat_type = fson_type_mpi(json, "output.checkpoint.repeat")
             select case (checkpoint_repeat_type)
             case (TYPE_INTEGER)
                call fson_get_mpi(json, "output.checkpoint.repeat", &
                     val = checkpoint_repeat)
             case (TYPE_LOGICAL)
                call fson_get_mpi(json, "output.checkpoint.repeat", &
                     val = checkpoint_do_repeat)
                if (checkpoint_do_repeat) then
                   checkpoint_repeat = -1
                else
                   checkpoint_repeat = default_checkpoint_repeat
                end if
             case default
                call self%ode%logfile%write(LOG_LEVEL_WARN, 'input', 'unrecognised', &
                     str_key = "output.checkpoint.repeat", str_value = '...')
             end select
          else
             checkpoint_repeat = default_checkpoint_repeat
             call self%ode%logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  logical_keys = ['output.checkpoint.repeat'], &
                  logical_values = [PETSC_FALSE])
          end if
          if (fson_has_mpi(json, "output.checkpoint.time")) then
             call fson_get_mpi(json, "output.checkpoint.time", &
                  val = checkpoint_time)
          else if (fson_has_mpi(json, "output.checkpoint.step")) then
             call fson_get_mpi(json, "output.checkpoint.step", &
                  val = checkpoint_step)
             allocate(checkpoint_time(size(checkpoint_step)))
             checkpoint_time(1) = self%ode%time + checkpoint_step(1)
             do i = 2, size(checkpoint_step)
                checkpoint_time(i) = checkpoint_time(i-1) + checkpoint_step(i)
             end do
             deallocate(checkpoint_step)
          end if
          call fson_get_mpi(json, "output.checkpoint.tolerance", &
            default_checkpoint_tol, checkpoint_tol, self%ode%logfile)
       else
          call self%ode%logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
               str_key = 'output.checkpoint', str_value = 'none')
       end if

    end if

    call self%steps%init(self%method%num_stored_steps, &
         self%ode%time, self%ode%solution, &
         stop_time_specified, stop_time, max_num_steps, max_stepsize, &
         adapt_on, adapt_method, adapt_min, adapt_max, &
         adapt_reduction, adapt_amplification, step_sizes, &
         nonlinear_relative_tol, nonlinear_abs_tol, &
         nonlinear_update_relative_tol, nonlinear_update_abs_tol, &
         nonlinear_min_iterations, max_num_tries, steady_state, &
         checkpoint_time, checkpoint_repeat, checkpoint_tol, &
         stop_minimum_stepsize, stop_maximum_stepsize)

    call fson_get_mpi(json, "output.initial", &
         default_output_initial, self%output_initial, self%ode%logfile)
    call fson_get_mpi(json, "output.final", &
         default_output_final, self%output_final, self%ode%logfile)

    call self%ode%logfile%write_blank()

    if (allocated(checkpoint_time)) then
       deallocate(checkpoint_time)
    end if

  contains

    PetscInt function time_step_method_from_str(method_str) result(method)
      character(*), intent(in) :: method_str
      select case (str_to_lower(method_str))
      case ("beuler")
         method = TS_BEULER
      case ("bdf2")
         method = TS_BDF2
      case ("directss")
         method = TS_DIRECTSS
      case default
         method = TS_BEULER
      end select
    end function time_step_method_from_str

    PetscInt function adapt_method_from_str(adapt_method_str) result(method)
      character(*), intent(in) :: adapt_method_str
      select case (str_to_lower(adapt_method_str))
      case ("change")
         method = TS_ADAPT_CHANGE
      case ("iteration")
         method = TS_ADAPT_ITERATION
      case default
         method = TS_ADAPT_CHANGE
      end select
    end function adapt_method_from_str

    KSPType function ksp_type_from_str(ksp_type_str) result(ksp_type)
      character(*), intent(in) :: ksp_type_str
      select case(str_to_lower(ksp_type_str))
      case ("gmres")
         ksp_type = KSPGMRES
      case ("lgmres")
         ksp_type = KSPLGMRES
      case ("bcgs")
         ksp_type = KSPBCGS
      case ("bcgsl")
         ksp_type = KSPBCGSL
      case default
         ksp_type = KSPBCGS
      end select
    end function ksp_type_from_str

    PCType function pc_type_from_str(pc_type_str) result(pc_type)
      character(*), intent(in) :: pc_type_str
      select case(str_to_lower(pc_type_str))
      case ("none")
         pc_type = PCNONE
      case ("lu")
         pc_type = PCLU
      case ("ilu")
         pc_type = PCILU
      case ("bjacobi")
         pc_type = PCBJACOBI
      case ("asm")
         pc_type = PCASM
      case default
         pc_type = PCASM
      end select
    end function pc_type_from_str
    
  end subroutine timestepper_init

!------------------------------------------------------------------------

  subroutine timestepper_destroy(self)
    !! Destroys a timestepper.

    class(timestepper_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call timestepper_solver_context_destroy(self%context)

    call SNESDestroy(self%solver, ierr);  CHKERRQ(ierr)
    call VecDestroy(self%residual, ierr); CHKERRQ(ierr)
    call MatDestroy(self%jacobian, ierr); CHKERRQ(ierr)

    call self%steps%destroy()

    nullify(self%ode)

  end subroutine timestepper_destroy

!------------------------------------------------------------------------

  subroutine timestepper_step(self)
    !! Takes a single time step.

    class(timestepper_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    PetscBool :: accepted
    SNESConvergedReason :: converged_reason

    call self%steps%update()
    call self%ode%pre_timestep()
    accepted = PETSC_FALSE

    do while (.not. (accepted .or. (self%steps%finished)))

       call self%steps%initialize_try()
       if (self%steps%current%num_tries > 0) then
          call self%ode%pre_retry_timestep()
       end if

       call self%steps%check_checkpoints(self%ode%logfile)
       call self%steps%check_finished(self%ode%logfile)

       call SNESSolve(self%solver, PETSC_NULL_VEC, self%steps%current%solution, &
            ierr); CHKERRQ(ierr)
       call SNESGetIterationNumber(self%solver, self%steps%current%num_iterations, &
            ierr); CHKERRQ(ierr)
       call SNESGetConvergedReason(self%solver, converged_reason, ierr); CHKERRQ(ierr)

       self%steps%current%num_tries = self%steps%current%num_tries + 1
       call self%steps%set_current_status(converged_reason)
       call self%steps%set_next_stepsize(accepted)
       call self%log_step_status(accepted, converged_reason)

    end do

    self%steps%taken = self%steps%taken + 1
    self%ode%time = self%steps%current%time
    call VecCopy(self%steps%current%solution, self%ode%solution, ierr)
    CHKERRQ(ierr)
    call self%ode%post_timestep()
    call self%ode%logfile%flush()

  end subroutine timestepper_step

!------------------------------------------------------------------------

  subroutine timestepper_log_step_status(self, accepted, converged_reason)
    !! Writes logfile messages about timestep status.

    class(timestepper_type), intent(in out) :: self
    PetscBool, intent(in) :: accepted
    SNESConvergedReason, intent(in) :: converged_reason
    ! Locals:
    PetscInt, parameter :: reason_str_len = 80
    character(len = reason_str_len) :: reason_str
    KSP :: linear_solver
    PetscInt :: iterations, log_level
    PetscBool :: snes_converged
    KSPConvergedReason :: ksp_reason
    PetscErrorCode :: ierr

    if (converged_reason == SNES_DIVERGED_LINEAR_SOLVE) then
       call SNESGetKSP(self%solver, linear_solver, ierr); CHKERRQ(ierr)
       call KSPGetIterationNumber(linear_solver, iterations, ierr); CHKERRQ(ierr)
       call KSPGetConvergedReason(linear_solver, ksp_reason, ierr); CHKERRQ(ierr)
       reason_str = KSP_reason_str(ksp_reason)
       call self%ode%logfile%write(LOG_LEVEL_WARN, 'linear_solver', 'end', &
            logical_keys = ['converged'], logical_values = [ksp_reason >= 0], &
            int_keys = ['iterations'], int_values = [iterations], &
            str_key = 'reason', str_value = trim(reason_str))
    end if

    reason_str = SNES_reason_str(converged_reason)
    snes_converged = (converged_reason >= 0)
    if (snes_converged) then
       log_level = LOG_LEVEL_INFO
    else
       log_level = LOG_LEVEL_WARN
    end if

    call self%ode%logfile%write(log_level, 'nonlinear_solver', 'end', &
         logical_keys = ['converged'], logical_values = [snes_converged], &
         int_keys = ['iterations'], &
         int_values = [self%steps%current%num_iterations], &
         str_key = 'reason', str_value = trim(reason_str))

    if (self%steps%current%status == TIMESTEP_ABORTED) then
       call self%ode%logfile%write(LOG_LEVEL_WARN, 'timestep', 'aborted', &
            int_keys = ['num_tries'], &
            int_values = [self%steps%current%num_tries])
    end if

    if (.not. (accepted .or. self%steps%steady_state)) then
       call self%ode%logfile%write(LOG_LEVEL_WARN, 'timestep', 'reduction', &
            real_keys = ['new_size        '], &
            real_values = [self%steps%next_stepsize])
    end if

  contains

    function SNES_reason_str(converged_reason) result (s)
      SNESConvergedReason, intent(in) :: converged_reason
      character(len = reason_str_len) :: s
      select case (converged_reason)
      case (SNES_CONVERGED_FNORM_ABS)
         s = "function_absolute"
      case (SNES_CONVERGED_FNORM_RELATIVE)
         s = "function_relative"
      case (SNES_CONVERGED_SNORM_RELATIVE)
         s = "update_relative"
      case (SNES_CONVERGED_ITS)
         s = "iterations"
      case (SNES_DIVERGED_TR_DELTA)
         s = "tr_delta"
      case (SNES_DIVERGED_FUNCTION_DOMAIN)
         s = "function_domain"
      case (SNES_DIVERGED_FUNCTION_COUNT)
         s = "function_count"
      case (SNES_DIVERGED_LINEAR_SOLVE)
         s = "linear_solver"
      case (SNES_DIVERGED_FNORM_NAN)
         s = "function_norm_NaN"
      case (SNES_DIVERGED_MAX_IT)
         s = "max_iterations"
      case (SNES_DIVERGED_LINE_SEARCH)
         s = "line_search"
      case (SNES_DIVERGED_INNER)
         s = "inner_solve"
      case (SNES_DIVERGED_LOCAL_MIN)
         s = "local_min"
      case (SNES_DIVERGED_DTOL)
         s = "divergence tolerance"
      case default
         s = "unknown"
      end select
    end function SNES_reason_str

    function KSP_reason_str(ksp_reason) result (s)
      KSPConvergedReason :: ksp_reason
      character(len = reason_str_len) :: s
      select case (ksp_reason)
         case (KSP_DIVERGED_NULL)
            s = "null"
         case (KSP_DIVERGED_ITS)
            s = "max_iterations"
         case (KSP_DIVERGED_DTOL)
            s = "dtol"
         case (KSP_DIVERGED_BREAKDOWN)
            s = "breakdown"
         case(KSP_DIVERGED_BREAKDOWN_BICG)
            s = "BiCG_breakdown"
         case (KSP_DIVERGED_NONSYMMETRIC)
            s = "non-symmetric"
         case (KSP_DIVERGED_INDEFINITE_PC)
            s = "indefinite PC"
         case (KSP_DIVERGED_NANORINF)
            s = "NaN_or_Inf"
         case (KSP_DIVERGED_INDEFINITE_MAT)
            s = "indefinite matrix"
         case (KSP_DIVERGED_PC_FAILED)
            s = "PC failed"
         case default
            s = "unknown"
      end select
    end function KSP_reason_str

  end subroutine timestepper_log_step_status

!------------------------------------------------------------------------

  subroutine timestepper_run(self)
    !! Runs the timestepper until finished.

    use iso_fortran_env, only: int32
    use utils_module, only : date_time_str, clock_elapsed_time

    class(timestepper_type), intent(in out) :: self
    ! Locals:
    PetscInt :: since_output
    PetscErrorCode :: err
    integer(int32) :: start_clock
    PetscReal :: elapsed_time

    call self%ode%logfile%write(LOG_LEVEL_INFO, 'timestepper', 'start', &
         real_keys = ['time'], real_values = [self%steps%current%time], &
         str_key = 'wall_time            ', str_value = date_time_str())
    call self%ode%logfile%flush()
    call system_clock(start_clock)

    err = 0
    self%steps%taken = 0
    self%output_index = 0
    since_output = 0

    call self%initial_function_calls(err)

    if (err == 0) then

       if (self%output_initial) then
          call self%ode%output(self%output_index, self%steps%current%time)
          self%output_index = self%output_index + 1
       end if

       do while (.not. self%steps%finished)

          if (associated(self%before_step_output)) then
             call self%before_step_output()
          end if

          call self%step()

          since_output = since_output + 1

          if (associated(self%after_step_output)) then
             call self%after_step_output()
          end if

          if ((self%steps%checkpoints%hit) .or. &
               (since_output == self%output_frequency)) then
             call self%ode%output(self%output_index, self%steps%current%time)
             self%output_index = self%output_index + 1
             if (self%steps%checkpoints%hit) then
                call self%steps%checkpoints%update()
             end if
             if (since_output == self%output_frequency) then
                since_output = 0
             end if
          end if

       end do

       if (self%output_final .and. (since_output > 0)) then
          call self%ode%output(self%output_index, self%steps%current%time)
       end if

    end if

    elapsed_time = clock_elapsed_time(start_clock)
    call self%ode%logfile%write_blank()
    call self%ode%logfile%write(LOG_LEVEL_INFO, 'timestepper', 'end', &
         real_keys = ['elapsed_seconds'], real_values = [elapsed_time], &
         str_key = 'wall_time', str_value = date_time_str())
    call self%ode%logfile%flush()

  end subroutine timestepper_run

!------------------------------------------------------------------------

  subroutine timestepper_input_summary(self)
    !! Writes summary of important inputs to the ODE logfile.

    class(timestepper_type), intent(in out) :: self
    ! Locals:
    character(max_ksp_type_str_len) :: ksp_type_str
    KSP :: ksp
    KSPType :: ksp_type
    character(max_pc_type_str_len) :: pc_type_str
    PC :: pc
    PCType :: pc_type
    PetscErrorCode :: ierr

    call self%ode%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'time.step.method', &
         str_value = self%method%name)

    call SNESGetKSP(self%solver, ksp, ierr); CHKERRQ(ierr)
    call KSPGetType(ksp, ksp_type, ierr); CHKERRQ(ierr)
    ksp_type_str = ksp_type_str_from_type(ksp_type)
    call self%ode%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'time.step.solver.linear.type', &
         str_value = ksp_type_str)

    call KSPGetPC(ksp, pc, ierr); CHKERRQ(ierr)
    call PCGetType(pc, pc_type, ierr); CHKERRQ(ierr)
    pc_type_str = pc_type_str_from_type(pc_type)
    call self%ode%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'time.step.solver.linear.preconditioner.type', &
         str_value = pc_type_str)

    call self%ode%logfile%write_blank()

  contains

    function ksp_type_str_from_type(ksp_type) result(str)
      character(max_ksp_type_str_len) :: str
      KSPType, intent(in) :: ksp_type
      select case (ksp_type)
      case (KSPGMRES)
         str = 'gmres'
      case (KSPLGMRES)
         str = 'lgmres'
      case (KSPBCGS)
         str = 'bcgs'
      case (KSPBCGSL)
         str = 'bcgsl'
      case default
         str = 'unknown'
      end select
    end function ksp_type_str_from_type

    function pc_type_str_from_type(pc_type) result(str)
      character(max_pc_type_str_len) :: str
      PCType, intent(in) :: pc_type
      select case (pc_type)
      case (PCBJACOBI)
         str = 'bjacobi'
      case (PCASM)
         str = 'asm'
      case (PCGASM)
         str = 'gasm'
      case (PCILU)
         str = 'ilu'
      case (PCSOR)
         str = 'sor'
      case (PCSPAI)
         str = 'spai'
      case default
         str = 'unknown'
      end select
    end function pc_type_str_from_type

  end subroutine timestepper_input_summary

!------------------------------------------------------------------------

end module timestepper_module
