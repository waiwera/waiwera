module timestepper_module
  !! Time stepping methods for solving \(\frac{d}{dt} L(t,y) = R(t,y)\), where y is a parallel vector.

  use kinds_module
  use mpi_module
  use ode_module
  use logfile_module

  implicit none

  private

#include <petsc/finclude/petsc.h90>

  ! Timestepping methods
     PetscInt, parameter, public :: TS_BEULER = 0, TS_BDF2 = 1, TS_DIRECTSS = 2

  ! Timestepping adaptor methods
     PetscInt, parameter, public :: TS_ADAPT_CHANGE = 0, TS_ADAPT_ITERATION = 1
     PetscInt, parameter, public :: timestepper_max_monitor_name_length = 12

  ! Timestep status
     PetscInt, parameter, public :: TIMESTEP_OK = 0, TIMESTEP_NOT_CONVERGED = 1, &
     TIMESTEP_TOO_SMALL = 2, TIMESTEP_TOO_BIG = 3, TIMESTEP_ABORTED = 4, &
     TIMESTEP_FINAL = 5

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
     PetscReal, public :: next_stepsize !! Step size to be used for next step
     PetscReal, public :: stop_time !! Maximum allowable time
     PetscReal, public :: nonlinear_solver_relative_tol !! Relative tolerance for non-linear solver
     PetscReal, public :: nonlinear_solver_abs_tol !! Absolute tolerance for non-linear solver
     PetscReal, public :: termination_tol = 1.e-6_dp !! Tolerance for detecting stop time reached
     PetscReal, allocatable, public :: sizes(:) !! Pre-specified time step sizes
     type(timestep_adaptor_type), public :: adaptor !! Time step size adaptor
     PetscBool, public :: stop_time_specified !! Whether stop time is specified or not
     PetscBool :: finished !! Whether simulation has yet finished
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
   contains
     private
     procedure, public :: init => timestepper_solver_context_init
     procedure, public :: destroy => timestepper_solver_context_destroy
  end type timestepper_solver_context_type

  type, public :: timestepper_type
     !! Timestepper class.
     private
     SNES :: solver
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
     PetscInt, public :: output_index !! Counter for determining when main output is needed
     PetscBool, public :: output_initial !! Whether to output initial conditions 
     PetscBool, public :: output_final !! Whether to output final results
   contains
     private
     procedure :: setup_jacobian => timestepper_setup_jacobian
     procedure :: setup_solver => timestepper_setup_solver
     procedure :: step => timestepper_step
     procedure :: log_step_status => timestepper_log_step_status
     procedure, public :: init => timestepper_init
     procedure, public :: destroy => timestepper_destroy
     procedure, public :: initial_function_calls => timestepper_initial_function_calls
     procedure, public :: run => timestepper_run
  end type timestepper_type

  interface

     subroutine method_residual(solver, y, residual, context, err)
       !! Residual routine to be minimised by nonlinear solver.
       import :: timestepper_solver_context_type
       SNES, intent(in) :: solver
       Vec, intent(in) :: y
       Vec, intent(out) :: residual
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

    !! Monitors relative change in solution. The parameter eps
    !! avoids division by zero if last%lhs is zero.

    type(timestepper_step_type), intent(in) :: current, last
    ! Locals:
    PetscReal, parameter :: eps = 1.e-3_dp
    NormType, parameter :: nt = NORM_2
    Vec :: diff
    PetscReal :: norm_diff, norm_current
    PetscErrorCode :: ierr

    call VecDuplicate(current%lhs, diff, ierr); CHKERRQ(ierr)
    call VecCopy(current%lhs, diff, ierr); CHKERRQ(ierr)
    call VecAXPY(diff, -1.0_dp, last%lhs, ierr); CHKERRQ(ierr)
    call VecNorm(diff, nt, norm_diff, ierr); CHKERRQ(ierr)
    call VecNorm(current%lhs, nt, norm_current, ierr); CHKERRQ(ierr)
    call VecDestroy(diff, ierr); CHKERRQ(ierr)

    relative_change_monitor = norm_diff / (norm_current + eps)

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
    Vec, intent(out) :: residual
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: t, dt
    PetscErrorCode :: ierr

    err = 0
    t = context%steps%current%time
    dt = context%steps%current%stepsize

    call context%ode%lhs(t, y, context%steps%current%lhs, err)
    if (err == 0) then
       call VecCopy(context%steps%current%lhs, residual, ierr); CHKERRQ(ierr)
       call VecAXPY(residual, -1.0_dp, context%steps%last%lhs, ierr); CHKERRQ(ierr)
       call context%ode%rhs(t, y, context%steps%current%rhs, err)
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
    Vec, intent(out) :: residual
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(timestepper_step_type), pointer :: last2
    PetscReal :: t, dt, dtlast
    PetscReal :: r, r1
    PetscErrorCode :: ierr

    if (context%steps%taken == 0) then

       ! Startup- use backwards Euler
       call backwards_Euler_residual(solver, y, residual, context, err)

    else

       err = 0
       t  = context%steps%current%time
       dt = context%steps%current%stepsize
       dtlast = context%steps%last%stepsize
       r = dt / dtlast
       r1 = r + 1._dp

       last2 => context%steps%pstore(3)%p

       call context%ode%lhs(t, y, context%steps%current%lhs, err)
       if (err == 0) then
          call VecCopy(context%steps%current%lhs, residual, ierr); CHKERRQ(ierr)
          call VecScale(residual, 1._dp + 2._dp * r, ierr); CHKERRQ(ierr)
          call VecAXPY(residual, -r1 * r1, context%steps%last%lhs, ierr)
          CHKERRQ(ierr)
          call VecAXPY(residual, r * r, last2%lhs, ierr); CHKERRQ(ierr)
          call context%ode%rhs(t, y, context%steps%current%rhs, err)
          if (err == 0) then
             call VecAXPY(residual, -dt * r1, context%steps%current%rhs, ierr)
             CHKERRQ(ierr)
          end if
       end if

    end if

  end subroutine BDF2_residual

!------------------------------------------------------------------------

  subroutine direct_ss_residual(solver, y, residual, context, err)
    !! Residual for direct solution of steady state equations R(y) = 0.
    !! Here we evaluate R() at the steps final time.

    SNES, intent(in) :: solver
    Vec, intent(in) :: y
    Vec, intent(out) :: residual
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr

    err = 0
    call context%ode%rhs(context%steps%stop_time, y, &
         context%steps%current%rhs, err)
    if (err == 0) then
       call VecCopy(context%steps%current%rhs, residual, ierr)
       CHKERRQ(ierr)
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

    SNES, intent(in) :: solver
    Vec, intent(in) :: y
    Vec, intent(out) :: residual
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr, ferr

    err = 0; ferr = 0
    call context%ode%pre_eval(context%steps%current%time, y, ferr)
    if (ferr == 0) then
       call context%residual(solver, y, residual, context, ferr)
    end if

    if (ferr > 0) then
       call SNESSetFunctionDomainError(solver, ierr); CHKERRQ(ierr)
    end if

  end subroutine SNES_residual

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
         ['tries           ', 'iters           '], &
         [self%num_tries, self%num_iterations], &
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
! Timestepper_steps procedures
!------------------------------------------------------------------------

  subroutine timestepper_steps_init(self, num_stored, &
       time, solution, initial_stepsize, stop_time_specified, &
       stop_time, max_num_steps, max_stepsize, &
       adapt_on, adapt_method, adapt_min, adapt_max, &
       adapt_reduction, adapt_amplification, step_sizes, &
       nonlinear_solver_relative_tol, nonlinear_solver_abs_tol, &
       max_num_tries)

    !! Sets up array of timesteps and pointers to them. This array stores the
    !! current step and one or more previous steps. The number of stored
    !! steps depends on the timestepping method (e.g. 2 for single-step methods,
    !! > 2 for multistep methods). 

    class(timestepper_steps_type), intent(in out) :: self
    PetscInt, intent(in) :: num_stored
    PetscReal, intent(in) :: time, initial_stepsize
    PetscBool, intent(in) :: stop_time_specified
    Vec, intent(in) :: solution
    PetscReal, intent(in) :: stop_time
    PetscInt, intent(in) :: max_num_steps
    PetscReal, intent(in) :: max_stepsize
    PetscBool, intent(in) :: adapt_on
    PetscInt, intent(in) :: adapt_method
    PetscReal, intent(in) :: adapt_min, adapt_max
    PetscReal, intent(in) :: adapt_reduction, adapt_amplification
    PetscReal, intent(in), allocatable, optional :: step_sizes(:)
    PetscReal, intent(in) :: nonlinear_solver_relative_tol, &
         nonlinear_solver_abs_tol
    PetscInt, intent(in) :: max_num_tries
    ! Locals:
    PetscInt :: i
    PetscErrorCode :: ierr

    self%taken = 0
    self%finished = PETSC_FALSE
    self%num_stored = num_stored
    allocate(self%store(num_stored), self%pstore(num_stored))

    do i = 1, num_stored
       call self%store(i)%init(solution)
       call self%set_pstore(self%store, i, i)
    end do

    call self%set_aliases()

    self%current%time = time
    self%next_stepsize = initial_stepsize
    call VecCopy(solution, self%current%solution, ierr); CHKERRQ(ierr)

    self%stop_time_specified = stop_time_specified
    self%stop_time = stop_time
    self%max_num = max_num_steps
    self%max_num_tries = max_num_tries

    self%adaptor%on = adapt_on
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

    self%nonlinear_solver_relative_tol = nonlinear_solver_relative_tol
    self%nonlinear_solver_abs_tol = nonlinear_solver_abs_tol

    if ((present(step_sizes) .and. (allocated(step_sizes)))) then
       ! Fixed time step sizes override adaptor:
       self%sizes = step_sizes
       self%adaptor%on = PETSC_FALSE
       self%next_stepsize = step_sizes(1)
    end if

  end subroutine timestepper_steps_init

!------------------------------------------------------------------------

  subroutine timestepper_initial_function_calls(self, err)
    !! Performs initial LHS function call at start of run (and pre-evaluation
    !! procedure if needed).

    class(timestepper_type), intent(in out) :: self
    PetscErrorCode, intent(out) :: err

    err = 0

    call self%ode%pre_solve(self%steps%current%time, &
         self%steps%current%solution, err)

    if (err == 0) then
       call self%ode%lhs(self%steps%current%time, &
            self%steps%current%solution, self%steps%current%lhs, err)
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

  subroutine timestepper_steps_check_finished(self)
    !! Checks if any termination criteria have been met, and reduces
    !! timestep if needed.

    class(timestepper_steps_type), intent(in out) :: self
    
    self%finished = PETSC_FALSE

    if (self%stop_time_specified .and. &
         (self%current%time > self%stop_time - self%termination_tol)) then
       self%current%stepsize = self%stop_time - self%last%time
       self%current%time = self%stop_time
       self%finished = PETSC_TRUE
    end if

    if (self%taken + 1 >= self%max_num) then
       self%finished = PETSC_TRUE
    end if

  end subroutine timestepper_steps_check_finished

!------------------------------------------------------------------------

  subroutine timestepper_steps_set_current_status(self, converged_reason)
    !! Sets status of current step.

    class(timestepper_steps_type), intent(in out) :: self
    SNESConvergedReason, intent(in) :: converged_reason
    ! Locals:
    PetscReal :: eta

    if (self%current%num_tries >= self%max_num_tries) then
       self%current%status = TIMESTEP_ABORTED
       self%finished = PETSC_TRUE
    else if (converged_reason >= 0 ) then
       if ((self%finished) .and. (self%current%status /= &
            TIMESTEP_ABORTED)) then
          self%current%status = TIMESTEP_FINAL
       else
          eta = self%adaptor%monitor(self%current, self%last)
          if (self%adaptor%on) then
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
    else
       self%current%status = TIMESTEP_NOT_CONVERGED
    end if

  end subroutine timestepper_steps_set_current_status

!------------------------------------------------------------------------

  subroutine timestepper_steps_get_next_fixed_stepsize(self)
    !! Gets next timestep size if using fixed time step sizes.

    class(timestepper_steps_type), intent(in out) :: self
    ! Locals:
    PetscInt :: step_index, next_index

    step_index = self%taken + 1
    next_index = step_index + 1

    if ((allocated(self%sizes)) .and. (next_index <= size(self%sizes))) then
       self%next_stepsize = self%sizes(next_index)
    else
       self%next_stepsize = self%current%stepsize
    end if

  end subroutine timestepper_steps_get_next_fixed_stepsize

!------------------------------------------------------------------------

  subroutine timestepper_steps_set_next_stepsize(self, accepted)

    class(timestepper_steps_type), intent(in out) :: self
    PetscBool, intent(out) :: accepted
    
    if (self%adaptor%on) then
       call self%adapt(accepted)
    else
       if (self%current%status == TIMESTEP_OK) then
          accepted = PETSC_TRUE
          call self%get_next_fixed_stepsize()
       else
          accepted = PETSC_FALSE
          self%adaptor%on = PETSC_TRUE
          self%next_stepsize = self%adaptor%reduce(self%current%stepsize)
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
! Timestepper procedures
!------------------------------------------------------------------------

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

    nullify(self%ode)
    nullify(self%steps)
    nullify(self%residual)

  end subroutine timestepper_solver_context_destroy

!------------------------------------------------------------------------
! Timestepper procedures
!------------------------------------------------------------------------

  subroutine timestepper_setup_solver(self, max_iterations)
    !! Sets up SNES nonlinear solver for the timestepper.

    class(timestepper_type), intent(in out) :: self
    PetscInt, intent(in) :: max_iterations
    ! Locals:
    PetscErrorCode :: ierr
    KSP :: ksp
    PC :: pc
    SNESLineSearch :: linesearch
    ! This tolerance needs to be set very small so it doesn't override
    ! time step reduction when primary variables go out of bounds:
    PetscReal, parameter :: stol = 1.e-99_dp

    call self%context%init(self%ode, self%steps, self%method%residual)

    call SNESCreate(mpi%comm, self%solver, ierr); CHKERRQ(ierr)
    call SNESSetApplicationContext(self%solver, self%context, ierr); CHKERRQ(ierr)
    call SNESSetFunction(self%solver, self%residual, &
         SNES_residual, self%context, ierr); CHKERRQ(ierr)
    call SNESSetJacobian(self%solver, self%jacobian, self%jacobian, &
         PETSC_NULL_FUNCTION, PETSC_NULL_OBJECT, ierr); CHKERRQ(ierr)
    call SNESSetDM(self%solver, self%ode%mesh%dm, ierr); CHKERRQ(ierr)

    ! Set nonlinear and linear solver options from command line options:
    call SNESSetFromOptions(self%solver, ierr); CHKERRQ(ierr)
    call SNESGetKSP(self%solver, ksp, ierr); CHKERRQ(ierr)
    call KSPSetFromOptions(ksp, ierr); CHKERRQ(ierr)
    call KSPGetPC(ksp, pc, ierr); CHKERRQ(ierr)
    call PCSetFromOptions(pc, ierr); CHKERRQ(ierr)

    call SNESSetTolerances(self%solver, PETSC_DEFAULT_REAL, &
         PETSC_DEFAULT_REAL, stol, max_iterations, &
         PETSC_DEFAULT_INTEGER, ierr); CHKERRQ(ierr)
    call SNESSetConvergenceTest(self%solver, SNES_convergence, self%context, &
         PETSC_NULL_FUNCTION, ierr); CHKERRQ(ierr)
    call SNESMonitorSet(self%solver, SNES_monitor, self%context, &
         PETSC_NULL_FUNCTION, ierr); CHKERRQ(ierr)

    ! Set function to be called at start of each solver iteration:
    call SNESSetUpdate(self%solver, SNES_pre_iteration_update, ierr)
    CHKERRQ(ierr)

    ! Nonlinear solver line search:
    call SNESGetLineSearch(self%solver, linesearch, ierr); CHKERRQ(ierr)
    call SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC, ierr)
    CHKERRQ(ierr)
    call SNESLineSearchSetPostCheck(linesearch, SNES_linesearch_post_check, &
         self%context, ierr); CHKERRQ(ierr)

  end subroutine timestepper_setup_solver

!------------------------------------------------------------------------

  subroutine SNES_monitor(solver, num_iterations, fnorm, context, ierr)
    !! SNES monitor routine for output at each nonlinear iteration.

    SNES, intent(in) :: solver
    PetscInt, intent(in) :: num_iterations
    PetscReal, intent(in) :: fnorm
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode :: ierr

    if ((num_iterations > 0)) then
       call context%ode%logfile%write(LOG_LEVEL_INFO, 'nonlinear_solver', &
            'iteration', &
            ['count           '], [num_iterations], &
            ['max_residual    '], [context%steps%current%max_residual])
    end if
    call context%ode%logfile%flush()

  end subroutine SNES_monitor

!------------------------------------------------------------------------

  subroutine SNES_convergence(solver, num_iterations, xnorm, pnorm, &
       fnorm, reason, context, ierr)
    !! Tests for convergence of nonlinear solver.

    SNES, intent(in) :: solver
    PetscInt, intent(in) :: num_iterations
    PetscReal, intent(in) :: xnorm, pnorm, fnorm
    SNESConvergedReason, intent(out) :: reason
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode :: ierr
    ! Locals:
    Vec :: unscaled_residual, residual
    PetscReal, pointer, contiguous :: unscaled_residual_array(:), &
         residual_array(:), lhs_array(:)
    PetscInt :: low, hi, i
    PetscReal :: scale

    call SNESGetFunction(solver, unscaled_residual, PETSC_NULL_OBJECT, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    call DMGetGlobalVector(context%ode%mesh%dm, residual, ierr)
    CHKERRQ(ierr)

    call VecGetArrayF90(unscaled_residual, unscaled_residual_array, ierr)
    CHKERRQ(ierr)
    call VecGetArrayF90(residual, residual_array, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(context%steps%last%lhs, lhs_array, ierr)
    CHKERRQ(ierr)

    call VecGetOwnershipRange(unscaled_residual, low, hi, ierr); CHKERRQ(ierr)
    do i = 1, hi - low
       scale = max(abs(lhs_array(i)), context%steps%nonlinear_solver_abs_tol)
       residual_array(i) = unscaled_residual_array(i) / scale
    end do

    call VecRestoreArrayF90(unscaled_residual, unscaled_residual_array, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(residual, residual_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(context%steps%current%lhs, lhs_array, ierr)
    CHKERRQ(ierr)

    call VecNorm(residual, NORM_INFINITY, &
         context%steps%current%max_residual, ierr); CHKERRQ(ierr)

    call DMRestoreGlobalVector(context%ode%mesh%dm, residual, ierr)
    CHKERRQ(ierr)

    if (context%steps%current%max_residual < &
         context%steps%nonlinear_solver_relative_tol) then
       reason = SNES_CONVERGED_FNORM_ABS
    else
       call SNESConvergedDefault(solver, num_iterations, xnorm, pnorm, &
            fnorm, reason, PETSC_NULL_OBJECT, ierr); CHKERRQ(ierr)
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

    call PetscOptionsSetValue(PETSC_NULL_OBJECT, "-mat_fd_type", "ds", ierr)
    CHKERRQ(ierr)
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
    use fson_value_m, only : TYPE_NULL

    class(timestepper_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    class(ode_type), intent(in), target :: ode
    ! Locals:
    PetscInt :: max_num_steps
    PetscReal, parameter :: default_stop_time = 1.0_dp
    PetscReal, parameter :: default_initial_stepsize = 0.1_dp
    PetscReal :: initial_stepsize, max_stepsize, stop_time
    PetscReal, parameter :: default_max_stepsize = 0.0_dp
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
         default_adapt_method_str = "change"
    PetscInt :: adapt_method
    PetscReal, parameter :: default_adapt_min = 0.01_dp, &
         default_adapt_max = 0.1_dp
    PetscReal :: adapt_min, adapt_max
    PetscReal, parameter :: default_adapt_reduction = 0.5_dp, &
         default_adapt_amplification = 2.0_dp
    PetscReal :: adapt_reduction, adapt_amplification
    PetscReal, allocatable :: step_sizes(:)
    PetscInt, parameter :: default_nonlinear_solver_max_iterations = 8
    PetscReal, parameter :: default_nonlinear_solver_relative_tol = 1.e-6_dp
    PetscReal, parameter :: default_nonlinear_solver_abs_tol = 1._dp
    PetscInt :: nonlinear_solver_max_iterations
    PetscReal :: nonlinear_solver_relative_tol, nonlinear_solver_abs_tol
    PetscInt :: max_num_tries
    PetscInt, parameter :: default_max_num_tries = 10
    PetscInt, parameter :: default_output_frequency = 1
    PetscBool, parameter :: default_output_initial = PETSC_TRUE
    PetscBool, parameter :: default_output_final = PETSC_TRUE
    PetscBool :: stop_time_specified
    PetscErrorCode :: ierr

    self%ode => ode
    call VecDuplicate(self%ode%solution, self%residual, ierr); CHKERRQ(ierr)
    call self%setup_jacobian()

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
       
    call fson_get_mpi(json, "time.step.method", &
         default_method_str, method_str, self%ode%logfile)
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

    call self%method%init(method)

    call fson_get_mpi(json, "time.step.initial", &
         default_initial_stepsize, initial_stepsize, self%ode%logfile)
    call fson_get_mpi(json, "time.step.maximum.size", &
         default_max_stepsize, max_stepsize, self%ode%logfile)
    call fson_get_mpi(json, "time.step.maximum.number", &
         default_max_num_steps, max_num_steps, self%ode%logfile)
    call fson_get_mpi(json, "time.step.maximum.tries", &
         default_max_num_tries, max_num_tries, self%ode%logfile)

    call fson_get_mpi(json, "time.step.adapt.on", &
         default_adapt_on, adapt_on, self%ode%logfile)

    call fson_get_mpi(json, "time.step.adapt.method", &
         default_adapt_method_str, adapt_method_str, self%ode%logfile)
    select case (str_to_lower(adapt_method_str))
    case ("change")
       adapt_method = TS_ADAPT_CHANGE
    case ("iteration")
       adapt_method = TS_ADAPT_ITERATION
    case default
       adapt_method = TS_ADAPT_CHANGE
    end select

    call fson_get_mpi(json, "time.step.adapt.min", &
         default_adapt_min, adapt_min, self%ode%logfile)
    call fson_get_mpi(json, "time.step.adapt.max", &
         default_adapt_max, adapt_max, self%ode%logfile)

    call fson_get_mpi(json, "time.step.adapt.reduction", &
         default_adapt_reduction, adapt_reduction, self%ode%logfile)
    call fson_get_mpi(json, "time.step.adapt.amplification", &
         default_adapt_amplification, adapt_amplification, self%ode%logfile)

    if (fson_has_mpi(json, "time.step.sizes")) then
       call fson_get_mpi(json, "time.step.sizes", val = step_sizes)
    else
       call self%ode%logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
            str_key = "time.step.sizes", str_value = "adaptive")
    end if

    call fson_get_mpi(json, &
         "time.step.solver.nonlinear.maximum.iterations", &
         default_nonlinear_solver_max_iterations, &
         nonlinear_solver_max_iterations, self%ode%logfile)
    call fson_get_mpi(json, &
         "time.step.solver.nonlinear.tolerance.relative", &
         default_nonlinear_solver_relative_tol, &
         nonlinear_solver_relative_tol, self%ode%logfile)
    call fson_get_mpi(json, &
         "time.step.solver.nonlinear.tolerance.absolute", &
         default_nonlinear_solver_abs_tol, nonlinear_solver_abs_tol, &
         self%ode%logfile)

    call self%steps%init(self%method%num_stored_steps, &
         self%ode%time, self%ode%solution, initial_stepsize, &
         stop_time_specified, stop_time, max_num_steps, max_stepsize, &
         adapt_on, adapt_method, adapt_min, adapt_max, &
         adapt_reduction, adapt_amplification, step_sizes, &
         nonlinear_solver_relative_tol, nonlinear_solver_abs_tol, &
         max_num_tries)

    call self%setup_solver(nonlinear_solver_max_iterations)

    call fson_get_mpi(json, "output.frequency", &
         default_output_frequency, self%output_frequency, self%ode%logfile)
    call fson_get_mpi(json, "output.initial", &
         default_output_initial, self%output_initial, self%ode%logfile)
    call fson_get_mpi(json, "output.final", &
         default_output_final, self%output_final, self%ode%logfile)

    call self%ode%logfile%write_blank()

  end subroutine timestepper_init

!------------------------------------------------------------------------

  subroutine timestepper_destroy(self)
    !! Destroys a timestepper.

    class(timestepper_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call self%context%destroy()

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

       call self%steps%check_finished()

       call SNESSolve(self%solver, PETSC_NULL_OBJECT, self%steps%current%solution, &
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

    if (self%steps%current%status == TIMESTEP_ABORTED) then
       call self%ode%logfile%write(LOG_LEVEL_WARN, 'timestep', 'aborted', &
            int_keys = ['num_tries'], &
            int_values = [self%steps%current%num_tries])
    end if

    if (.not. accepted) then
       if (converged_reason < 0) then
          call log_SNES_divergence()
       end if
       call self%ode%logfile%write(LOG_LEVEL_WARN, 'timestep', 'reduction', &
            real_keys = ['new_size        '], &
            real_values = [self%steps%next_stepsize])
    end if

  contains

    subroutine log_SNES_divergence()
      character(len = 80) :: reason_str
      select case (converged_reason)
         case (SNES_DIVERGED_FUNCTION_DOMAIN)
            reason_str = "function_domain"
         case (SNES_DIVERGED_FUNCTION_COUNT)
            reason_str = "function_count"
         case (SNES_DIVERGED_LINEAR_SOLVE)
            reason_str = "linear_solver"
            call log_linear_solver_divergence()
         case (SNES_DIVERGED_FNORM_NAN)
            reason_str = "function_norm_NaN"
         case (SNES_DIVERGED_MAX_IT)
            reason_str = "max_iterations"
         case (SNES_DIVERGED_LINE_SEARCH)
            reason_str = "line_search"
         case (SNES_DIVERGED_INNER)
            reason_str = "inner_solve"
         case (SNES_DIVERGED_LOCAL_MIN)
            reason_str = "local_min"
         case default
            reason_str = "unknown"
      end select
      call self%ode%logfile%write(LOG_LEVEL_WARN, 'nonlinear_solver', &
           'failed', str_key = 'reason', &
           str_value = trim(reason_str))
    end subroutine log_SNES_divergence

    subroutine log_linear_solver_divergence()
      ! Locals:
      KSP :: linear_solver
      KSPConvergedReason :: ksp_reason
      character(len = 80) :: ksp_reason_str
      PetscErrorCode :: ierr
      call SNESGetKSP(self%solver, linear_solver, ierr); CHKERRQ(ierr)
      call KSPGetConvergedReason(linear_solver, ksp_reason, ierr); CHKERRQ(ierr)
      select case (ksp_reason)
         case (KSP_DIVERGED_DTOL)
            ksp_reason_str = "dtol"
         case (KSP_DIVERGED_BREAKDOWN)
            ksp_reason_str = "breakdown"
         case (KSP_DIVERGED_ITS)
            ksp_reason_str = "max_iterations"
         case (KSP_DIVERGED_NANORINF)
            ksp_reason_str = "NaN_or_Inf"
         case(KSP_DIVERGED_BREAKDOWN_BICG)
            ksp_reason_str = "BiCG_breakdown"
         case default
            ksp_reason_str = "unknown"
      end select
      call self%ode%logfile%write(LOG_LEVEL_WARN, &
           'linear_solver', 'failed', &
           str_key = 'reason', &
           str_value = ksp_reason_str)
    end subroutine log_linear_solver_divergence

  end subroutine timestepper_log_step_status

!------------------------------------------------------------------------

  subroutine timestepper_run(self)
    !! Runs the timestepper until finished.

    use utils_module, only : date_time_str

    class(timestepper_type), intent(in out) :: self
    ! Locals:
    PetscInt :: since_output
    PetscErrorCode :: err, ierr
    PetscLogDouble :: start_wall_time, end_wall_time, elapsed_time

    call self%ode%logfile%write(LOG_LEVEL_INFO, 'timestepper', 'start', &
         str_key = 'time            ', str_value = date_time_str())
    call self%ode%logfile%flush()
    call PetscTime(start_wall_time, ierr); CHKERRQ(ierr)

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

          if (since_output == self%output_frequency) then
             call self%ode%output(self%output_index, self%steps%current%time)
             self%output_index = self%output_index + 1
             since_output = 0
          end if

       end do

       if (self%output_final .and. (since_output > 0)) then
          call self%ode%output(self%output_index, self%steps%current%time)
       end if

    end if

    call PetscTime(end_wall_time, ierr); CHKERRQ(ierr)
    elapsed_time = end_wall_time - start_wall_time

    call self%ode%logfile%write_blank()
    call self%ode%logfile%write(LOG_LEVEL_INFO, 'timestepper', 'end', &
         real_keys = ['elapsed_seconds'], real_values = [elapsed_time], &
         str_key = 'time', str_value = date_time_str())
    call self%ode%logfile%flush()

  end subroutine timestepper_run

!------------------------------------------------------------------------

end module timestepper_module
