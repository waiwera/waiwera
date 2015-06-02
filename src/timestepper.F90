module timestepper_module
  !! Time stepping methods for solving \(\frac{d}{dt} L(t,y) = R(t,y)\), where y is a parallel vector.

  use kinds_module
  use mpi_module
  use ode_module

  implicit none

  private

#include <petsc-finclude/petsc.h90>

  ! Timestepping methods
     PetscInt, parameter, public :: TS_BEULER = 0, TS_BDF2 = 1, TS_DIRECTSS = 2

  ! Timestepping adaptor methods
     PetscInt, parameter, public :: TS_ADAPT_CHANGE = 0, TS_ADAPT_ITERATION = 1
     PetscInt, parameter, public :: timestepper_max_monitor_name_length = 12

  ! Timestep status
     PetscInt, parameter, public :: TIMESTEP_OK = 0, TIMESTEP_NOT_CONVERGED = 1, &
     TIMESTEP_TOO_SMALL = 2, TIMESTEP_TOO_BIG = 3

  type timestepper_step_type
     !! Results of a time step.
     private
     PetscReal, public :: time, stepsize
     Vec, public :: solution, lhs, rhs
     PetscInt, public :: num_tries, num_iterations, status = TIMESTEP_OK
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
     PetscBool, public :: on
     character(timestepper_max_monitor_name_length), public :: name
     procedure(monitor_function), pointer, nopass, public :: monitor => relative_change_monitor
     PetscReal, public :: monitor_min, monitor_max
     PetscReal, public :: reduction, amplification
     PetscReal, public :: max_stepsize
   contains
     private
     procedure :: increase => timestep_adaptor_increase
     procedure :: reduce => timestep_adaptor_reduce
  end type timestep_adaptor_type

  type timestepper_steps_type
     !! Storage for current and immediate past steps in timestepper.
     private
     type(timestepper_step_type), allocatable :: store(:)
     type(ptimestepper_step_type), pointer :: pstore(:)
     type(timestepper_step_type), pointer, public :: current, last
     PetscBool :: finished
     PetscInt, public :: num_stored, taken, max_num
     PetscReal, public :: next_stepsize, start_time, final_time
     PetscReal, public :: termination_tol = 1.e-6_dp
     type(timestep_adaptor_type), public :: adaptor
   contains
     private
     procedure :: set_aliases => timestepper_steps_set_aliases
     procedure :: set_pstore => timestepper_steps_set_pstore
     procedure :: update => timestepper_steps_update
     procedure :: rotate => timestepper_steps_rotate
     procedure, public :: init => timestepper_steps_init
     procedure, public :: destroy => timestepper_steps_destroy
     procedure, public :: check_finished => timestepper_steps_check_finished
     procedure, public :: set_current_status => timestepper_steps_set_current_status
     procedure, public :: adapt => timestepper_steps_adapt
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
     class(ode_type), pointer, public :: ode
     type(timestepper_steps_type), public :: steps
     type(timestepper_method_type), public :: method
     procedure(step_output_routine), pointer, public :: &
          step_output => step_output_default
   contains
     private
     procedure :: setup_solver => timestepper_setup_solver
     procedure :: step => timestepper_step
     procedure, public :: init => timestepper_init
     procedure, public :: destroy => timestepper_destroy
     procedure, public :: initial_function_calls => timestepper_initial_function_calls
     procedure, public :: run => timestepper_run
  end type timestepper_type

  interface

     subroutine method_residual(solver, y, residual, context, ierr)
       !! Residual routine to be minimised by nonlinear solver.
       import :: timestepper_solver_context_type
       SNES, intent(in) :: solver
       Vec, intent(in) :: y
       Vec, intent(out) :: residual
       type(timestepper_solver_context_type), intent(in out) :: context
       PetscErrorCode, intent(out) :: ierr
     end subroutine method_residual

     PetscReal function monitor_function(current, last)
       !! Function for monitoring timestep acceptability.
       import :: timestepper_step_type
       type(timestepper_step_type), intent(in) :: current, last
     end function monitor_function

     subroutine step_output_routine(self)
       !! Routine for producing output at each time step.
       import :: timestepper_type
       class(timestepper_type), intent(in out) :: self
     end subroutine step_output_routine

  end interface

  ! Subroutines to be available outside this module:
  public :: iteration_monitor, relative_change_monitor
  public :: step_output_routine, step_output_default

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

  subroutine step_output_default(self)
    !! Default routine for printing diagnostic information at each
    !! time step.

    class(timestepper_type), intent(in out) :: self

    if (mpi%rank == mpi%output_rank) then
       write(*, '(a, i4)'), 'step:', self%steps%taken
       call self%steps%current%print()
    end if

  end subroutine step_output_default

!------------------------------------------------------------------------
! Residual routines
!------------------------------------------------------------------------

  subroutine backwards_Euler_residual(solver, y, residual, context, ierr)
    !! Residual for backwards Euler method.
    !! residual = L(1) - L(0)  - dt * R(1)

    SNES, intent(in) :: solver
    Vec, intent(in) :: y
    Vec, intent(out) :: residual
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode, intent(out) :: ierr
    ! Locals:
    PetscReal :: t, dt

    t = context%steps%current%time
    dt = context%steps%current%stepsize

    call context%ode%lhs(t, y, context%steps%current%lhs)
    call VecCopy(context%steps%current%lhs, residual, ierr); CHKERRQ(ierr)
    call VecAXPY(residual, -1.0_dp, context%steps%last%lhs, ierr); CHKERRQ(ierr)
    call context%ode%rhs(t, y, context%steps%current%rhs)
    call VecAXPY(residual, -dt, context%steps%current%rhs, ierr); CHKERRQ(ierr)

  end subroutine backwards_Euler_residual

!------------------------------------------------------------------------

  subroutine BDF2_residual(solver, y, residual, context, ierr)
    !! Residual for variable-stepsize BDF2 method.
    !! residual = (1 + 2r) * L(1) - (r+1)^2 * L(0) + r^2 * L(-1) - dt * (r+1) * R(1)
    !! where r = dt / (last dt)

    SNES, intent(in) :: solver
    Vec, intent(in) :: y
    Vec, intent(out) :: residual
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode, intent(out) :: ierr
    ! Locals:
    type(timestepper_step_type), pointer :: last2
    PetscReal :: t, dt, dtlast
    PetscReal :: r, r1

    if (context%steps%taken == 0) then

       ! Startup- use backwards Euler
       call backwards_Euler_residual(solver, y, residual, context, ierr)

    else

       t  = context%steps%current%time
       dt = context%steps%current%stepsize
       dtlast = context%steps%last%stepsize
       r = dt / dtlast
       r1 = r + 1._dp

       last2 => context%steps%pstore(3)%p

       call context%ode%lhs(t, y, context%steps%current%lhs)
       call VecCopy(context%steps%current%lhs, residual, ierr); CHKERRQ(ierr)
       call VecScale(residual, 1._dp + 2._dp * r, ierr); CHKERRQ(ierr)
       call VecAXPY(residual, -r1 * r1, context%steps%last%lhs, ierr); CHKERRQ(ierr)
       call VecAXPY(residual, r * r, last2%lhs, ierr); CHKERRQ(ierr)
       call context%ode%rhs(t, y, context%steps%current%rhs)
       call VecAXPY(residual, -dt * r1, context%steps%current%rhs, ierr)
       CHKERRQ(ierr)

    end if

  end subroutine BDF2_residual

!------------------------------------------------------------------------

  subroutine direct_ss_residual(solver, y, residual, context, ierr)
    !! Residual for direct solution of steady state equations R(y) = 0.
    !! Here we evaluate R() at the steps final time.

    SNES, intent(in) :: solver
    Vec, intent(in) :: y
    Vec, intent(out) :: residual
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode, intent(out) :: ierr

    call context%ode%rhs(context%steps%final_time, y, residual)

  end subroutine direct_ss_residual

!------------------------------------------------------------------------

  subroutine SNES_residual(solver, y, residual, &
       context, ierr)
    !! Residual routine to be minimized by SNES solver. This calls the
    !! pre-evaluation routine first (if needed) before calling the
    !! timestepper method residual routine.

    SNES, intent(in) :: solver
    Vec, intent(in) :: y
    Vec, intent(out) :: residual
    type(timestepper_solver_context_type), intent(in out) :: context
    PetscErrorCode, intent(out) :: ierr

    call context%ode%pre_eval(context%steps%current%time, y)
    call context%residual(solver, y, residual, context, ierr)

  end subroutine SNES_residual

!------------------------------------------------------------------------
! Timestepper_step procedures
!------------------------------------------------------------------------

  subroutine timestepper_step_init(self, template_vec)
    !! Initializes a timestep.

    class(timestepper_step_type), intent(inout) :: self
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

    class(timestepper_step_type), intent(inout) :: self
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
       case default
          timestepper_step_status_str = 'unknown'
    end select

  end function timestepper_step_status_str

!------------------------------------------------------------------------

  subroutine timestepper_step_print(self)
    !! Prints a timestep.

    class(timestepper_step_type), intent(in) :: self

    write(*, '(a, e12.4, a, e12.4, a, i2, a, a)') &
         'time: ', self%time, &
         ' stepsize: ',self%stepsize, &
         ' iters: ', self%num_iterations, &
         ' status: ', self%status_str()

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
       initial_time, solution, initial_stepsize, &
       final_time, max_num_steps, max_stepsize, &
       adapt_on, adapt_method, adapt_min, adapt_max, &
       adapt_reduction, adapt_amplification)

    !! Sets up array of timesteps and pointers to them. This array stores the
    !! current step and one or more previous steps. The number of stored
    !! steps depends on the timestepping method (e.g. 2 for single-step methods,
    !! > 2 for multistep methods). 

    class(timestepper_steps_type), intent(in out) :: self
    PetscInt, intent(in) :: num_stored
    PetscReal, intent(in) :: initial_time, initial_stepsize
    Vec, intent(in) :: solution
    PetscReal, intent(in) :: final_time
    PetscInt, intent(in) :: max_num_steps
    PetscReal, intent(in) :: max_stepsize
    PetscBool, intent(in) :: adapt_on
    PetscInt, intent(in) :: adapt_method
    PetscReal, intent(in) :: adapt_min, adapt_max
    PetscReal, intent(in) :: adapt_reduction, adapt_amplification
    ! Locals:
    PetscInt :: i
    PetscErrorCode :: ierr

    self%taken = 0
    self%finished = .false.
    self%num_stored = num_stored
    allocate(self%store(num_stored), self%pstore(num_stored))

    do i = 1, num_stored
       call self%store(i)%init(solution)
       call self%set_pstore(self%store, i, i)
    end do

    call self%set_aliases()

    self%start_time = initial_time
    self%current%time = self%start_time
    self%next_stepsize = initial_stepsize
    call VecCopy(solution, self%current%solution, ierr); CHKERRQ(ierr)

    self%final_time = final_time
    self%max_num = max_num_steps

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

  end subroutine timestepper_steps_init

!------------------------------------------------------------------------

  subroutine timestepper_initial_function_calls(self)
    !! Performs initial LHS function call at start of run (and pre-evaluation
    !! procedure if needed).

    class(timestepper_type), intent(in out) :: self

    call self%ode%pre_eval(self%steps%current%time, self%steps%current%solution)

    call self%ode%lhs(self%steps%current%time, &
         self%steps%current%solution, self%steps%current%lhs)

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

  subroutine timestepper_steps_update(self)
    !! Updates pointers to store array, when timestep is advanced.

    class(timestepper_steps_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call self%rotate()
    call self%set_aliases()

    if (self%num_stored > 1) then
       ! Last solution is initial guess for current solution:
       call VecCopy(self%last%solution, self%current%solution, ierr); CHKERRQ(ierr)
    end if

    self%current%num_tries = 0

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
    
    self%finished = .false.

    if (self%current%time > self%final_time - self%termination_tol) then
       self%current%stepsize = self%final_time - self%last%time
       self%current%time = self%final_time
       self%finished = .true.
    end if

    if (self%taken >= self%max_num) then
       self%finished = .true.
    end if

  end subroutine timestepper_steps_check_finished

!------------------------------------------------------------------------

  subroutine timestepper_steps_set_current_status(self, converged)
    !! Sets status of current step.

    class(timestepper_steps_type), intent(in out) :: self
    PetscBool, intent(in) :: converged
    ! Locals:
    PetscReal :: eta

    if (converged) then
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
          ! Not quite sure what the desired behaviour is for 
          ! non-adaptive stepping- this may not be it
          self%current%status = TIMESTEP_OK
       end if
    else
       self%current%status = TIMESTEP_NOT_CONVERGED
    end if

  end subroutine timestepper_steps_set_current_status

!------------------------------------------------------------------------

  subroutine timestepper_steps_adapt(self, accepted)
    !! Checks current time stepsize for acceptability and adapts if 
    !! necessary.

    class(timestepper_steps_type), intent(in out) :: self
    PetscBool, intent(out) :: accepted

    select case (self%current%status)
    case (TIMESTEP_TOO_SMALL)
       accepted = .true.
       self%next_stepsize = self%adaptor%increase(self%current%stepsize)
    case (TIMESTEP_TOO_BIG, TIMESTEP_NOT_CONVERGED)
       accepted = .false.
       self%next_stepsize = self%adaptor%reduce(self%current%stepsize)
    case default
       accepted = .true.
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

  subroutine timestepper_setup_solver(self)
    !! Sets up SNES nonlinear solver for the timestepper.

    class(timestepper_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    KSP :: ksp

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

  end subroutine timestepper_setup_solver

!------------------------------------------------------------------------

  subroutine timestepper_init(self, json, ode)

    !! Initializes a timestepper.

    use utils_module, only : str_to_lower
    use fson
    use fson_mpi_module

    class(timestepper_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    class(ode_type), intent(in), target :: ode
    ! Locals:
    PetscInt :: max_num_steps
    PetscReal, parameter :: default_start_time = 0.0_dp, &
         default_stop_time = 1.0_dp
    PetscReal, parameter :: default_initial_stepsize = 0.1_dp
    PetscReal :: initial_stepsize, max_stepsize
    PetscReal :: start_time, stop_time
    PetscReal, parameter :: default_max_stepsize = 0.0_dp
    PetscInt :: method
    PetscInt, parameter :: max_method_str_len = 12
    character(max_method_str_len) :: method_str
    character(max_method_str_len), parameter :: default_method_str = "beuler"
    PetscInt, parameter :: default_max_num_steps = 100
    PetscBool, parameter :: default_adapt_on = .false.
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
    PetscErrorCode :: ierr

    self%ode => ode

    call VecDuplicate(self%ode%solution, self%residual, ierr); CHKERRQ(ierr)
    call DMSetMatType(self%ode%mesh%dm, MATAIJ, ierr); CHKERRQ(ierr)
    call DMCreateMatrix(self%ode%mesh%dm, self%jacobian, ierr); CHKERRQ(ierr)
    call MatSetFromOptions(self%jacobian, ierr); CHKERRQ(ierr)
    call MatSetUp(self%jacobian, ierr); CHKERRQ(ierr)

    call fson_get_mpi(json, "time.start", default_start_time, start_time)
    call fson_get_mpi(json, "time.stop", default_stop_time, stop_time)

    call fson_get_mpi(json, "time.step.method", &
         default_method_str, method_str)
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
         default_initial_stepsize, initial_stepsize)
    call fson_get_mpi(json, "time.step.maximum.size", &
         default_max_stepsize, max_stepsize)
    call fson_get_mpi(json, "time.step.maximum.number", &
         default_max_num_steps, max_num_steps)

    call fson_get_mpi(json, "time.step.adapt.on", &
         default_adapt_on, adapt_on)

    call fson_get_mpi(json, "time.step.adapt.method", &
         default_adapt_method_str, adapt_method_str)
    select case (str_to_lower(adapt_method_str))
    case ("change")
       adapt_method = TS_ADAPT_CHANGE
    case ("iteration")
       adapt_method = TS_ADAPT_ITERATION
    case default
       adapt_method = TS_ADAPT_CHANGE
    end select

    call fson_get_mpi(json, "time.step.adapt.min", &
         default_adapt_min, adapt_min)
    call fson_get_mpi(json, "time.step.adapt.max", &
         default_adapt_max, adapt_max)

    call fson_get_mpi(json, "time.step.adapt.reduction", &
         default_adapt_reduction, adapt_reduction)
    call fson_get_mpi(json, "time.step.adapt.amplification", &
         default_adapt_amplification, adapt_amplification)

    call self%steps%init(self%method%num_stored_steps, &
         start_time, self%ode%solution, initial_stepsize, &
         stop_time, max_num_steps, max_stepsize, &
         adapt_on, adapt_method, adapt_min, adapt_max, &
         adapt_reduction, adapt_amplification)
       
    call self%setup_solver()

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
    PetscBool :: converged, accepted
    SNESConvergedReason :: converged_reason
    
    call self%steps%update()
    accepted = .false.

    do while (.not. accepted)

       self%steps%current%stepsize = self%steps%next_stepsize
       self%steps%current%time = self%steps%last%time + self%steps%current%stepsize
       call self%steps%check_finished()

       call SNESSolve(self%solver, PETSC_NULL_OBJECT, self%steps%current%solution, &
            ierr); CHKERRQ(ierr)
       call SNESGetIterationNumber(self%solver, self%steps%current%num_iterations, &
            ierr); CHKERRQ(ierr)
       call SNESGetConvergedReason(self%solver, converged_reason, ierr); CHKERRQ(ierr)
       converged = (converged_reason >= 0)

       self%steps%current%num_tries = self%steps%current%num_tries + 1
       call self%steps%set_current_status(converged)
       call self%steps%adapt(accepted)

    end do

    self%steps%taken = self%steps%taken + 1
    call VecCopy(self%steps%current%solution, self%ode%solution, ierr)
    CHKERRQ(ierr)

    ! This may possibly be needed- to make sure current LHS corresponds exactly
    ! to current solution:
    ! call self%ode%lhs(self%steps%current%time, self%steps%current%solution, &
    !      self%steps%current%lhs)

  end subroutine timestepper_step

!------------------------------------------------------------------------

  subroutine timestepper_run(self)
    !! Runs the timestepper until finished.

    class(timestepper_type), intent(in out) :: self

    self%steps%taken = 0
    call self%initial_function_calls()

    do while (.not. self%steps%finished)
       call self%step()
       if (associated(self%step_output)) then
          call self%step_output()
       end if
    end do

  end subroutine timestepper_run

!------------------------------------------------------------------------

end module timestepper_module
