module timestepper_test

  ! Tests for timestepper module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use ode_module
  use timestepper_module
  use fson
  use fson_mpi_module

  implicit none

  private

  type, extends(ode_type) :: test_ode_type
     private
     PetscInt, public :: dim = 8, dof = 1, stencil = 0
     PetscReal, allocatable, public :: initial_values(:)
     PetscReal, public :: start_time
     PetscReal, public :: maxdiff
     Vec, public :: exact_solution, diff
     class(unit_test_type), pointer, public :: test
     character(:), allocatable, public :: case_name
   contains
     private
     procedure, public :: init => init_test_ode
     procedure, public :: destroy => destroy_test_ode
     procedure, public :: lhs => lhs_test_ode
     procedure, public :: rhs => rhs_test_ode
     procedure, public :: pre_eval => pre_eval_test_ode
     procedure, public :: output => output_test_ode
     procedure, public :: exact => exact_test_ode
     procedure, public :: run_cases => run_cases_test_ode
     procedure, public :: read_start_time => read_start_time_test_ode
     procedure, public :: set_initial_conditions => set_initial_conditions_test_ode
  end type test_ode_type

  type, extends(test_ode_type) :: linear_ode_type
     private
     PetscReal :: k = -0.5_dp
   contains
     private
     procedure, public :: init => init_linear
     procedure, public :: rhs => rhs_const
     procedure, public :: exact => exact_linear
     procedure, public :: aux_rhs => aux_rhs_const
  end type linear_ode_type

  type, extends(test_ode_type) :: exponential_ode_type
     private
     PetscReal :: k = -5._dp
   contains
     private
     procedure, public :: init => init_exponential
     procedure, public :: rhs => rhs_linear
     procedure, public :: exact => exact_exponential
     procedure, public :: aux_rhs => aux_rhs_linear
  end type exponential_ode_type

  type, extends(test_ode_type) :: logistic_ode_type
     private
     PetscReal, allocatable, public :: c(:)
   contains
     private
     procedure, public :: rhs => rhs_logistic
     procedure, public :: exact => exact_logistic
  end type logistic_ode_type

  type, extends(exponential_ode_type) :: nontrivial_lhs_ode_type
   contains
     private
     procedure, public :: lhs => lhs_nontrivial_lhs
     procedure, public :: exact => exact_nontrivial_lhs
     procedure, public :: aux_lhs => aux_lhs_nontrivial_lhs
  end type nontrivial_lhs_ode_type

  type, extends(test_ode_type) :: nonlinear_lhs_ode_type
     private
     PetscReal :: k = -1._dp
   contains
     private
     procedure, public :: lhs => lhs_nonlinear_lhs
     procedure, public :: rhs => rhs_nonlinear_lhs
     procedure, public :: exact => exact_nonlinear_lhs
  end type nonlinear_lhs_ode_type

  type, extends(test_ode_type) :: heat1d_ode_type
     PetscReal, public :: L, a
   contains
     private
     procedure :: fn => fn_heat1d
     procedure, public :: init => init_heat1d
     procedure, public :: rhs => rhs_heat1d
     procedure, public :: exact => exact_heat1d
  end type heat1d_ode_type

  type, extends(heat1d_ode_type) :: heat1d_nonlinear_ode_type
   contains
     private
     procedure, public :: lhs => lhs_heat1d_nonlinear
     procedure, public :: rhs => rhs_heat1d_nonlinear
  end type heat1d_nonlinear_ode_type

  type, extends(nontrivial_lhs_ode_type) :: pre_eval_ode_type
     private
     Vec :: secondary
   contains
     private
     procedure, public :: init => init_pre_eval
     procedure, public :: destroy => destroy_pre_eval
     procedure, public :: pre_solve => calculate_secondary_pre_eval
     procedure, public :: pre_eval => calculate_secondary_pre_eval
     procedure, public :: lhs => lhs_pre_eval
  end type pre_eval_ode_type

  type, extends(test_ode_type) :: ss_ode_type
   contains
     private
     procedure, public :: rhs => rhs_ss
     procedure, public :: exact => exact_ss
     procedure, public :: set_initial_conditions => set_initial_conditions_ss
  end type ss_ode_type

  PetscInt, parameter :: max_json_len = 256
  PetscReal, parameter :: time_tolerance = 1.e-6_dp

  public :: setup, teardown
  public :: test_timestepper_read, &
       test_timestepper_linear, test_timestepper_exponential, &
       test_timestepper_logistic, test_timestepper_nontrivial_lhs, &
       test_timestepper_nonlinear_lhs, test_timestepper_heat1d, &
       test_timestepper_heat1d_nonlinear, test_timestepper_pre_eval, &
       test_timestepper_steady, test_checkpoints

contains

!------------------------------------------------------------------------

  subroutine setup()

    use profiling_module, only: init_profiling

    ! Locals:
    PetscErrorCode :: ierr

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call init_profiling()

  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown()

    PetscErrorCode :: ierr

    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------
! Utility routines
!------------------------------------------------------------------------

  subroutine VecSetArray(v, arr)

    ! Sets local part of vector v with values from array.

    Vec, intent(in out) :: v
    PetscReal, intent(in) :: arr(:)
    ! Locals:
    PetscInt :: hi, low
    PetscErrorCode :: ierr
    PetscReal, pointer :: va(:)

    call VecGetOwnershipRange(v, low, hi, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(v, va, ierr); CHKERRQ(ierr)
    va(1: hi-low) = arr(low+1: hi)
    call VecRestoreArrayF90(v, va, ierr); CHKERRQ(ierr)
    
  end subroutine VecSetArray

!------------------------------------------------------------------------
! Default test ode functions
!------------------------------------------------------------------------

  subroutine init_test_ode(self, initial_array, auxiliary, test, err)
    ! Default ode initialization.
    class(test_ode_type), intent(in out) :: self
    PetscReal, intent(in), optional :: initial_array(:)
    PetscBool, intent(in), optional :: auxiliary
    class(unit_test_type), target, intent(in out) :: test
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr

    if (present(auxiliary)) then
       self%auxiliary = auxiliary
    else
       self%auxiliary = PETSC_FALSE
    end if

    err = 0
    self%dim = 8
    self%dof = 1
    self%stencil = 0
    self%test => test

    call DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, self%dim, self%dof, &
         self%stencil, PETSC_NULL_INTEGER, self%mesh%interior_dm, ierr); CHKERRQ(ierr)
    call DMSetUp(self%mesh%interior_dm, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(self%mesh%interior_dm, self%solution, ierr); CHKERRQ(ierr)
    call VecDuplicate(self%solution, self%exact_solution, ierr); CHKERRQ(ierr)
    call VecDuplicate(self%solution, self%diff, ierr); CHKERRQ(ierr)
    if (self%auxiliary) then
       call VecDuplicate(self%solution, self%aux_solution, ierr); CHKERRQ(ierr)
    end if
    if (present(initial_array)) then
       self%initial_values = initial_array
    end if

    call self%setup_jacobian()
    call self%setup_auxiliary()

    call self%logfile%init('', echo = PETSC_FALSE)

  end subroutine init_test_ode

  subroutine destroy_test_ode(self)
    ! Default ode destroy method.
    class(test_ode_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    call VecDestroy(self%solution, ierr); CHKERRQ(ierr)
    call MatDestroy(self%jacobian, ierr); CHKERRQ(ierr)
    call VecDestroy(self%exact_solution, ierr); CHKERRQ(ierr)
    if (self%auxiliary) then
       call VecDestroy(self%aux_solution, ierr); CHKERRQ(ierr)
       call MatDestroy(self%A_aux, ierr); CHKERRQ(ierr)
       call VecDestroy(self%b_aux, ierr); CHKERRQ(ierr)
    end if
    call VecDestroy(self%diff, ierr); CHKERRQ(ierr)
    call DMDestroy(self%mesh%interior_dm, ierr); CHKERRQ(ierr)
    call self%logfile%destroy()
  end subroutine destroy_test_ode

  subroutine lhs_test_ode(self, t, interval, y, lhs, err)
    ! Default identity LHS function.
    class(test_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Vec, intent(in) :: y
    Vec, intent(in out) :: lhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    call VecCopy(y, lhs, ierr); CHKERRQ(ierr)
  end subroutine lhs_test_ode

  subroutine rhs_test_ode(self, t, interval, y, rhs, err)
    ! Default zero RHS function.
    class(test_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Vec, intent(in) :: y
    Vec, intent(in out) :: rhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    call VecSet(rhs, 0._dp, ierr); CHKERRQ(ierr)
  end subroutine rhs_test_ode

  subroutine pre_eval_test_ode(self, t, y, perturbed_columns, err)
    ! Default do-nothing pre-evaluation routine.
    class(test_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    PetscInt, intent(in), optional :: perturbed_columns(:)
    PetscErrorCode, intent(out) :: err
    continue
  end subroutine pre_eval_test_ode

  subroutine exact_test_ode(self, t, v)
    ! Default zero exact solution.
    class(test_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in out) :: v
    ! Locals:
    PetscErrorCode :: ierr
    call VecSet(v, 0._dp, ierr); CHKERRQ(ierr)
  end subroutine exact_test_ode

  subroutine output_test_ode(self, time_index, time)
    ! Output from test ode- compute difference betweeen solution and
    ! exact solution.
    class(test_ode_type), intent(in out) :: self
    PetscInt, intent(in) :: time_index
    PetscReal, intent(in) :: time
    ! Locals:
    PetscReal, pointer :: y(:), yex(:), yaux(:)
    PetscErrorCode :: ierr
    PetscInt :: local_size, low, hi
    character(len = 50) :: msg

    call self%exact(self%time, self%exact_solution)
    call VecGetLocalSize(self%exact_solution, local_size, ierr)
    CHKERRQ(ierr)
    call VecGetOwnershipRange(self%exact_solution, low, hi, ierr)
    CHKERRQ(ierr)
    call VecGetArrayReadF90(self%solution, y, ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(self%exact_solution, yex, ierr); CHKERRQ(ierr)
    write(msg, '(a, i0)') " solution at time index ", time_index
    call self%test%assert(yex, y, self%case_name // trim(msg))
    if (self%auxiliary) then
       call VecGetArrayReadF90(self%aux_solution, yaux, ierr); CHKERRQ(ierr)
       call self%test%assert(yex, yaux, self%case_name // " aux" // trim(msg))
       call VecRestoreArrayReadF90(self%aux_solution, yaux, ierr); CHKERRQ(ierr)
    end if
    call VecRestoreArrayReadF90(self%exact_solution, yex, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%solution, y, ierr)
    CHKERRQ(ierr)

  end subroutine output_test_ode

!------------------------------------------------------------------------
! Specific test ode functions
!------------------------------------------------------------------------

  subroutine init_linear(self, initial_array, auxiliary, test, err)
    ! Linear ode initialization.
    class(linear_ode_type), intent(in out) :: self
    PetscReal, intent(in), optional :: initial_array(:)
    PetscBool, intent(in), optional :: auxiliary
    class(unit_test_type), target, intent(in out) :: test
    PetscErrorCode, intent(out) :: err
    call self%test_ode_type%init(initial_array, auxiliary = PETSC_TRUE, &
         test = test, err = err)
  end subroutine init_linear

  subroutine rhs_const(self, t, interval, y, rhs, err)
    ! rhs(t, y) = k
    class(linear_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Vec, intent(in) :: y
    Vec, intent(in out) :: rhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call VecSet(rhs, self%k, ierr); CHKERRQ(ierr)
  end subroutine rhs_const

  subroutine exact_linear(self, t, v)
    ! Linear solution
    class(linear_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in out) :: v
    call VecSetArray(v, self%initial_values + &
         self%k * (t - self%start_time))
  end subroutine exact_linear

  subroutine aux_rhs_const(self, t, interval, Ar, br, err)
    class(linear_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Mat, intent(in out) :: Ar
    Vec, intent(in out) :: br
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call MatZeroEntries(Ar, ierr); CHKERRQ(ierr)
    call VecSet(br, self%k, ierr); CHKERRQ(ierr)
  end subroutine aux_rhs_const

!------------------------------------------------------------------------

  subroutine init_exponential(self, initial_array, auxiliary, test, err)
    ! Exponential ode initialization.
    class(exponential_ode_type), intent(in out) :: self
    PetscReal, intent(in), optional :: initial_array(:)
    PetscBool, intent(in), optional :: auxiliary
    class(unit_test_type), target, intent(in out) :: test
    PetscErrorCode, intent(out) :: err
    call self%test_ode_type%init(initial_array, auxiliary = PETSC_TRUE, &
         test = test, err = err)
  end subroutine init_exponential

!------------------------------------------------------------------------

  subroutine rhs_linear(self, t, interval, y, rhs, err)
    ! rhs(t, y) = k * y
    class(exponential_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Vec, intent(in) :: y
    Vec, intent(in out) :: rhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call VecCopy(y, rhs, ierr); CHKERRQ(ierr)
    call VecScale(rhs, self%k, ierr); CHKERRQ(ierr)
  end subroutine rhs_linear

!------------------------------------------------------------------------

  subroutine exact_exponential(self, t, v)
    ! Exponential solution
    class(exponential_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in out) :: v
    call VecSetArray(v, self%initial_values * &
         exp(self%k * (t - self%start_time)))
  end subroutine exact_exponential

!------------------------------------------------------------------------

  subroutine aux_rhs_linear(self, t, interval, Ar, br, err)
    ! Auxiliary rhs(t, y) = k * y
    class(exponential_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Mat, intent(in out) :: Ar
    Vec, intent(in out) :: br
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call MatZeroEntries(Ar, ierr); CHKERRQ(ierr)
    call MatShift(Ar, self%k, ierr); CHKERRQ(ierr)
    call VecSet(br, 0._dp, ierr); CHKERRQ(ierr)
  end subroutine aux_rhs_linear

!------------------------------------------------------------------------

  subroutine rhs_logistic(self, t, interval, y, rhs, err)
    ! rhs(t, y) = (3 - 2 * y) * y
    class(logistic_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Vec, intent(in) :: y
    Vec, intent(in out) :: rhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call VecSet(rhs, 3._dp, ierr); CHKERRQ(ierr)
    call VecAXPY(rhs, -2._dp, y, ierr); CHKERRQ(ierr)
    call VecPointwiseMult(rhs, rhs, y, ierr); CHKERRQ(ierr)
  end subroutine rhs_logistic

  subroutine exact_logistic(self, t, v)
    ! Logistic solution 3 * c * exp(3 * t) / (1 + 2 * c * exp(3 * t))
    class(logistic_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in out) :: v
    ! Locals:
    PetscReal :: ce3t
    PetscErrorCode :: ierr
    PetscInt :: low, hi, i, ig
    PetscReal, pointer :: va(:)
    call VecGetOwnershipRange(v, low, hi, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(v, va, ierr); CHKERRQ(ierr)
    do i = 1, hi-low
       ig = low + i
       ce3t = self%c(ig) * exp(3._dp * (t - self%start_time))
       va(i) = 3._dp * ce3t / (1._dp + 2._dp * ce3t)
    end do
    call VecRestoreArrayF90(v, va, ierr); CHKERRQ(ierr)
  end subroutine exact_logistic

!------------------------------------------------------------------------

  subroutine lhs_nontrivial_lhs(self, t, interval, y, lhs, err)
    ! lhs(t, y) = t * y
    class(nontrivial_lhs_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Vec, intent(in) :: y
    Vec, intent(in out) :: lhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call VecCopy(y, lhs, ierr); CHKERRQ(ierr)
    call VecScale(lhs, t, ierr); CHKERRQ(ierr)
  end subroutine lhs_nontrivial_lhs

  subroutine exact_nontrivial_lhs(self, t, v)
    ! Solution y = y0 * t ** (k-1)
    class(nontrivial_lhs_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in out) :: v
    call VecSetArray(v, self%initial_values * t ** (self%k - 1._dp))
  end subroutine exact_nontrivial_lhs

  subroutine aux_lhs_nontrivial_lhs(self, t, interval, Al, err)
    ! Al = t * I
    class(nontrivial_lhs_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Vec, intent(in out) :: Al
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call VecSet(Al, t, ierr); CHKERRQ(ierr)
  end subroutine aux_lhs_nontrivial_lhs

!------------------------------------------------------------------------
  
  subroutine lhs_nonlinear_lhs(self, t, interval, y, lhs, err)
    ! lhs(t, y) = y * (y - 2)
    class(nonlinear_lhs_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Vec, intent(in) :: y
    Vec, intent(in out) :: lhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call VecPointwiseMult(lhs, y, y, ierr); CHKERRQ(ierr)
    call VecAXPY(lhs, -2._dp, y, ierr); CHKERRQ(ierr)
  end subroutine lhs_nonlinear_lhs

  subroutine rhs_nonlinear_lhs(self, t, interval, y, rhs, err)
    ! rhs(t, y) = k * (y - 1)
    class(nonlinear_lhs_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Vec, intent(in) :: y
    Vec, intent(in out) :: rhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call VecSet(rhs, -1._dp, ierr); CHKERRQ(ierr)
    call VecAXPY(rhs, 1._dp, y, ierr); CHKERRQ(ierr)
    call VecScale(rhs, self%k, ierr); CHKERRQ(ierr)
  end subroutine rhs_nonlinear_lhs

  subroutine exact_nonlinear_lhs(self, t, v)
    ! Linear solution y = y0 + k * t / 2
    class(nonlinear_lhs_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in out) :: v
    call VecSetArray(v, self%initial_values + &
         0.5_dp * self%k * t)
  end subroutine exact_nonlinear_lhs

!------------------------------------------------------------------------

  PetscReal function fn_heat1d(self, x, t)
    ! Exact heat equation solution at given point and time.
    use utils_module, only: pi
    class(heat1d_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: x, t
    ! Locals:
    fn_heat1d = sin(pi * x) * exp(-pi * pi * t)
  end function fn_heat1d

  subroutine exact_heat1d(self, t, v)
    ! Exact heat equation solution vector on mesh.
    class(heat1d_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in out) :: v
    ! Locals:
    DM :: cdm
    Vec :: cv
    PetscReal, pointer :: va(:), coords(:)
    PetscInt :: i1, im, i
    PetscReal :: x
    PetscErrorCode :: ierr

    call DMGetCoordinateDM(self%mesh%interior_dm, cdm, ierr); CHKERRQ(ierr)
    call DMGetCoordinates(self%mesh%interior_dm, cv, ierr); CHKERRQ(ierr)
    call DMDAVecGetArrayF90(cdm, cv, coords, ierr); CHKERRQ(ierr)
    call DMDAVecGetArrayF90(self%mesh%interior_dm, v, va, ierr); CHKERRQ(ierr)
    call DMDAGetCorners(self%mesh%interior_dm, i1, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, im, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    do i = i1, i1+im-1
       x = coords(i)
       va(i) = self%fn(x, t)
    end do
    call DMDAVecRestoreArrayF90(self%mesh%interior_dm, v, va, ierr)
    CHKERRQ(ierr)
    call DMDAVecRestoreArrayF90(cdm, cv, coords, ierr); CHKERRQ(ierr)

  end subroutine exact_heat1d

  subroutine init_heat1d(self, initial_array, auxiliary, test, err)
    ! Initialization for heat equation.
    class(heat1d_ode_type), intent(in out) :: self
    PetscReal, intent(in), optional :: initial_array(:)
    PetscBool, intent(in), optional :: auxiliary
    class(unit_test_type), target, intent(in out) :: test
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: dx
    PetscErrorCode :: ierr

    if (present(auxiliary)) then
       self%auxiliary = auxiliary
    else
       self%auxiliary = PETSC_FALSE
    end if

    err = 0
    self%L = 1._dp
    self%dim = 21
    self%dof = 1
    self%stencil = 1
    self%test => test

    dx = self%L / (self%dim - 1._dp)
    self%a = 1._dp / (dx*dx)

    call DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_GHOSTED, self%dim-2, &
         self%dof, self%stencil, PETSC_NULL_INTEGER, self%mesh%interior_dm, &
         ierr); CHKERRQ(ierr)
    call DMSetUp(self%mesh%interior_dm, ierr); CHKERRQ(ierr)
    call DMDASetUniformCoordinates(self%mesh%interior_dm, dx, self%L - dx, &
         0._dp, 0._dp, 0._dp, 0._dp, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(self%mesh%interior_dm, self%solution, ierr)
    CHKERRQ(ierr)
    call VecDuplicate(self%solution, self%exact_solution, ierr)
    CHKERRQ(ierr)
    call VecDuplicate(self%solution, self%diff, ierr); CHKERRQ(ierr)

    call self%setup_jacobian()
    call self%setup_auxiliary()

    call self%logfile%init('', echo = PETSC_FALSE)

  end subroutine init_heat1d

  subroutine rhs_heat1d(self, t, interval, y, rhs, err)
    ! rhs = d2y/dx2
    class(heat1d_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Vec, intent(in) :: y
    Vec, intent(in out) :: rhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    PetscReal, pointer :: ya(:), rhsa(:)
    PetscInt :: i1, im, i, i1g, img, i2g
    Vec :: ylocal

    err = 0
    call DMGetLocalVector(self%mesh%interior_dm, ylocal, ierr); CHKERRQ(ierr)
    call DMGlobalToLocalBegin(self%mesh%interior_dm, y, INSERT_VALUES, ylocal, &
         ierr); CHKERRQ(ierr)
    call DMGlobalToLocalEnd(self%mesh%interior_dm, y, INSERT_VALUES, ylocal, &
         ierr); CHKERRQ(ierr)

    call DMDAGetCorners(self%mesh%interior_dm, i1, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, im, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    call DMDAVecGetArrayF90(self%mesh%interior_dm, ylocal, ya, ierr)
    CHKERRQ(ierr)
    call DMDAVecGetArrayF90(self%mesh%interior_dm, rhs, rhsa, ierr)
    CHKERRQ(ierr)

    ! BCs:
    call DMDAGetGhostCorners(self%mesh%interior_dm, i1g, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, img, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    i2g = i1g+img-1
    if (i1g == 0) ya(i1g) = self%fn(0._dp, t)
    if (i2g == self%dim-2) ya(i2g) = self%fn(self%L, t)

    do i = i1, i1+im-1
       rhsa(i) = self%a * (ya(i-1) - 2._dp * ya(i) + ya(i+1))
    end do

    call DMDAVecRestoreArrayF90(self%mesh%interior_dm, ylocal, ya, ierr)
    CHKERRQ(ierr)
    call DMRestoreLocalVector(self%mesh%interior_dm, ylocal, ierr); CHKERRQ(ierr)

    call DMDAVecRestoreArrayF90(self%mesh%interior_dm, rhs, rhsa, ierr)
    CHKERRQ(ierr)

  end subroutine rhs_heat1d

!------------------------------------------------------------------------

  subroutine lhs_heat1d_nonlinear(self, t, interval, y, lhs, err)
    ! lhs = y * y
    class(heat1d_nonlinear_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Vec, intent(in) :: y
    Vec, intent(in out) :: lhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    PetscReal, pointer :: ya(:), lhsa(:)
    PetscInt :: i1, im, i
    Vec :: ylocal

    err = 0
    call DMGetLocalVector(self%mesh%interior_dm, ylocal, ierr); CHKERRQ(ierr)
    call DMGlobalToLocalBegin(self%mesh%interior_dm, y, INSERT_VALUES, ylocal, &
         ierr); CHKERRQ(ierr)
    call DMGlobalToLocalEnd(self%mesh%interior_dm, y, INSERT_VALUES, ylocal, &
         ierr); CHKERRQ(ierr)

    call DMDAGetCorners(self%mesh%interior_dm, i1, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, im, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    call DMDAVecGetArrayF90(self%mesh%interior_dm, ylocal, ya, ierr)
    CHKERRQ(ierr)
    call DMDAVecGetArrayF90(self%mesh%interior_dm, lhs, lhsa, ierr)
    CHKERRQ(ierr)

    do i = i1, i1+im-1
       lhsa(i) = ya(i) * ya(i)
    end do

    call DMDAVecRestoreArrayF90(self%mesh%interior_dm, ylocal, ya, ierr)
    CHKERRQ(ierr)
    call DMRestoreLocalVector(self%mesh%interior_dm, ylocal, ierr); CHKERRQ(ierr)

    call DMDAVecRestoreArrayF90(self%mesh%interior_dm, lhs, lhsa, ierr)
    CHKERRQ(ierr)

  end subroutine lhs_heat1d_nonlinear
  
  subroutine rhs_heat1d_nonlinear(self, t, interval, y, rhs, err)
    ! rhs = 2 * y * d2y/dx2
    class(heat1d_nonlinear_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Vec, intent(in) :: y
    Vec, intent(in out) :: rhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    PetscReal, pointer :: ya(:), rhsa(:)
    PetscInt :: i1, im, i, i1g, img, i2g
    Vec :: ylocal

    err = 0
    call DMGetLocalVector(self%mesh%interior_dm, ylocal, ierr); CHKERRQ(ierr)
    call DMGlobalToLocalBegin(self%mesh%interior_dm, y, INSERT_VALUES, &
         ylocal, ierr); CHKERRQ(ierr)
    call DMGlobalToLocalEnd(self%mesh%interior_dm, y, INSERT_VALUES, ylocal, &
         ierr); CHKERRQ(ierr)

    call DMDAGetCorners(self%mesh%interior_dm, i1, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, im, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    call DMDAVecGetArrayF90(self%mesh%interior_dm, ylocal, ya, ierr)
    CHKERRQ(ierr)
    call DMDAVecGetArrayF90(self%mesh%interior_dm, rhs, rhsa, ierr)
    CHKERRQ(ierr)

    ! BCs:
    call DMDAGetGhostCorners(self%mesh%interior_dm, i1g, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, img, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    i2g = i1g+img-1
    if (i1g == 0) ya(i1g) = self%fn(0._dp, t)
    if (i2g == self%dim-2) ya(i2g) = self%fn(self%L, t)

    do i = i1, i1+im-1
       rhsa(i) = 2._dp * ya(i) * self%a * (ya(i-1) - 2._dp * ya(i) + ya(i+1))
    end do

    call DMDAVecRestoreArrayF90(self%mesh%interior_dm, ylocal, ya, ierr)
    CHKERRQ(ierr)
    call DMRestoreLocalVector(self%mesh%interior_dm, ylocal, ierr); CHKERRQ(ierr)

    call DMDAVecRestoreArrayF90(self%mesh%interior_dm, rhs, rhsa, ierr)
    CHKERRQ(ierr)

  end subroutine rhs_heat1d_nonlinear

!------------------------------------------------------------------------
  
  subroutine init_pre_eval(self, initial_array, auxiliary, test, err)
    class(pre_eval_ode_type), intent(in out) :: self
    PetscReal, intent(in), optional :: initial_array(:)
    PetscBool, intent(in), optional :: auxiliary
    class(unit_test_type), target, intent(in out) :: test
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call self%test_ode_type%init(initial_array, auxiliary, test, err)
    call VecDuplicate(self%solution, self%secondary, ierr); CHKERRQ(ierr)
    self%test => test
  end subroutine init_pre_eval

  subroutine destroy_pre_eval(self)
    class(pre_eval_ode_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    call VecDestroy(self%secondary, ierr); CHKERRQ(ierr)
    call self%test_ode_type%destroy()
  end subroutine destroy_pre_eval

  subroutine calculate_secondary_pre_eval(self, t, y, perturbed_columns, err)
    ! Calculates secondary vector = t * y
    class(pre_eval_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    PetscInt, intent(in), optional :: perturbed_columns(:)
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call VecCopy(y, self%secondary, ierr); CHKERRQ(ierr)
    call VecScale(self%secondary, t, ierr); CHKERRQ(ierr)
  end subroutine calculate_secondary_pre_eval

  subroutine lhs_pre_eval(self, t, interval, y, lhs, err)
    ! lhs(t, y) = secondary
    class(pre_eval_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Vec, intent(in) :: y
    Vec, intent(in out) :: lhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call VecCopy(self%secondary, lhs, ierr); CHKERRQ(ierr)
  end subroutine lhs_pre_eval

!------------------------------------------------------------------------

  subroutine set_initial_conditions_ss(self)
    class(ss_ode_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    self%time = self%start_time
    call VecSet(self%solution, 0.5_dp, ierr); CHKERRQ(ierr)
  end subroutine set_initial_conditions_ss

  subroutine rhs_ss(self, t, interval, y, rhs, err)
    ! rhs(t, y) = 1 - y * y
    class(ss_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    Vec, intent(in) :: y
    Vec, intent(in out) :: rhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    Vec :: y2
    PetscErrorCode :: ierr
    err = 0
    call VecDuplicate(y, y2, ierr); CHKERRQ(ierr)
    call VecPointwiseMult(y2, y, y, ierr); CHKERRQ(ierr)
    call VecSet(rhs, 1._dp, ierr); CHKERRQ(ierr)
    call VecAXPY(rhs, -1._dp, y2, ierr); CHKERRQ(ierr)
    call VecDestroy(y2, ierr); CHKERRQ(ierr)
  end subroutine rhs_ss

  subroutine exact_ss(self, t, v)
    ! Constant y = 1
    class(ss_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in out) :: v
    ! Locals:
    PetscErrorCode :: ierr
    call VecSet(v, 1._dp, ierr); CHKERRQ(ierr)
  end subroutine exact_ss

!------------------------------------------------------------------------

  subroutine read_start_time_test_ode(self, json)
    ! Reads start_time property from JSON.

    class(test_ode_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscReal, parameter :: default_start_time = 0.0_dp

    call fson_get_mpi(json, "time.start", default_start_time, &
         self%start_time)

  end subroutine read_start_time_test_ode

!------------------------------------------------------------------------

  subroutine set_initial_conditions_test_ode(self)
    ! Sets initial conditions for ODE.

    class(test_ode_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    self%time = self%start_time
    call self%exact(self%time, self%solution)
    if (self%auxiliary) then
       call VecCopy(self%solution, self%aux_solution, ierr); CHKERRQ(ierr)
    end if

  end subroutine set_initial_conditions_test_ode

!------------------------------------------------------------------------

  subroutine run_cases_test_ode(self, json_str, tol)

    ! Run tests on ODE, comparing results with exact solutions.

    class(test_ode_type), intent(in out) :: self
    character(len = max_json_len), intent(in) :: json_str(:)
    PetscReal, intent(in) :: tol(:)
    ! Locals:
    PetscInt :: num_cases, i
    type(timestepper_type) :: ts
    type(fson_value), pointer :: json
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    num_cases = size(json_str)

    do i = 1, num_cases

       json => fson_parse_mpi(str = trim(json_str(i)))

       call self%read_start_time(json)
       call self%set_initial_conditions()
       self%maxdiff = 0._dp

       call ts%init(json, self)
       ts%before_step_output => NULL()
       ts%after_step_output => timestepper_step_output
       call SNESMonitorSet(ts%solver, PETSC_NULL_FUNCTION, 0, &
            PETSC_NULL_FUNCTION, ierr); CHKERRQ(ierr)

       self%test%tolerance = real(tol(i))
       self%case_name = trim(ts%method%name)

       call ts%run()

       if (rank == 0) then
          call self%test%assert(ts%steps%stop_time, &
               ts%steps%current%time, &
               trim(ts%method%name) // ' stop time', &
               time_tolerance)
       end if

       call ts%destroy()
       call fson_destroy_mpi(json)

    end do

  end subroutine run_cases_test_ode

!------------------------------------------------------------------------

  subroutine timestepper_step_output(self)
    ! Output at each time.

    class(timestepper_type), intent(in out) :: self

    call self%ode%output(self%steps%taken, self%steps%current%time)

  end subroutine timestepper_step_output

!------------------------------------------------------------------------
! Test routines
!------------------------------------------------------------------------

  subroutine test_timestepper_read(test)
    ! Tests assigning timestepper parameters from JSON string.

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    character(len = max_json_len) :: json_str
    type(fson_value), pointer :: json
    type(timestepper_type) :: ts
    type(test_ode_type) :: test_ode
    PetscReal, allocatable :: initial(:)
    PetscErrorCode :: err
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    allocate(initial(test_ode%dim))
    initial = [-4._dp, -3.0_dp, -2.0_dp, -1.0_dp, &
         0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]

    json_str = '{"time": {"stop": 1.0, ' // &
         '"step": {"size": 0.01, "maximum": {"number": 200, ' // &
         '"size": 3.e6}, "method": "beuler", ' // &
         '"adapt": {"on": true, "method": "change", "minimum": 0.01, ' // &
         '"maximum": 0.2, "reduction": 0.6, "amplification": 1.9}}}}'

    json => fson_parse_mpi(str = trim(json_str))
    call test_ode%init(initial, test = test, err = err)
    call ts%init(json, test_ode)

    if (rank == 0) then
       call test%assert(1.0_dp, ts%steps%stop_time, &
            "Timestepper stop time", time_tolerance)
       call test%assert(0.01_dp, ts%steps%next_stepsize, &
            "Timestepper initial stepsize", time_tolerance)
       call test%assert(200, ts%steps%max_num, "Timestepper max. num steps")
       call test%assert("Backward Euler", ts%method%name, "Timestepper method")
       call test%assert(.not. ts%steps%fixed, "Timestepper steps fixed")
       call test%assert("change", trim(ts%steps%adaptor%name), "Timestepper adapt method")
       call test%assert(0.01_dp, ts%steps%adaptor%monitor_min, "Timestepper monitor min")
       call test%assert(0.2_dp, ts%steps%adaptor%monitor_max, "Timestepper monitor max")
       call test%assert(0.6_dp, ts%steps%adaptor%reduction, "Timestepper monitor reduction")
       call test%assert(1.9_dp, ts%steps%adaptor%amplification, "Timestepper monitor amplification")
       call test%assert(3.e6_dp, ts%steps%adaptor%max_stepsize, "Timestepper max stepsize")
    end if

    call ts%destroy()
    call test_ode%destroy()
    call fson_destroy_mpi(json)

  end subroutine test_timestepper_read

!------------------------------------------------------------------------

  subroutine test_timestepper_linear(test)

    ! Linear function

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(linear_ode_type), target :: linear
    PetscInt,  parameter :: num_cases = 3
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [1.5e-6_dp, 1.5e-6_dp, 1.5e-6_dp]
    PetscReal, allocatable :: initial(:)
    PetscErrorCode :: err

    allocate(initial(linear%dim))
    initial = [-4._dp, -3.0_dp, -2.0_dp, -1.0_dp, &
         0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    call linear%init(initial, test = test, err = err)

    json_str(1) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"size": 0.1, "maximum": {"number": 20, "size": null}, ' // &
         '"method": "beuler", "adapt": {"on": false}}}}'

    json_str(2) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"size": 0.1, "maximum": {"number": null, "size": null}, ' // &
         '"method": "bdf2", "adapt": {"on": false}}}}'

    json_str(3) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"maximum": {"number": 10, "size": null}, ' // &
         '"method": "beuler", ' // &
         '"size": [0.1, 0.1, 0.2, 0.2, 0.3]}}}'

    call linear%run_cases(json_str, tol)

    call linear%destroy()
    deallocate(initial)

  end subroutine test_timestepper_linear

!------------------------------------------------------------------------

  subroutine test_timestepper_exponential(test)

    ! Exponential function

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(exponential_ode_type), target :: exponential
    PetscInt,  parameter :: num_cases = 3
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [0.12_dp, 0.04_dp, 0.18_dp]
    PetscReal, allocatable :: initial(:)
    PetscErrorCode :: err

    allocate(initial(exponential%dim))
    initial = [-4._dp, -3.0_dp, -2.0_dp, -1.0_dp, &
         0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    call exponential%init(initial, test = test, err = err)

    json_str(1) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"size": 0.01, "maximum": {"number": 200}, ' // &
         '"method": "beuler", ' // &
         '"adapt": {"on": true, "method": "change", ' // &
         '"minimum": 0.01, "maximum": 0.2}}}}'

    json_str(2) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"size": 0.05, "maximum": {"number": 200}, ' // &
         '"method": "bdf2", ' // &
         '"adapt": {"on": true, "method": "change", ' // &
         '"minimum": 0.01, "maximum": 0.2}}}}'

    json_str(3) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"maximum": {"number": 200}, ' // &
         '"method": "beuler", "adapt": {"on": false}, ' // &
         '"size": [0.005, 0.007, 0.01, 0.012, 0.014, 0.015]}}}'

    call exponential%run_cases(json_str, tol)

    call exponential%destroy()
    deallocate(initial)

  end subroutine test_timestepper_exponential

!------------------------------------------------------------------------

  subroutine test_timestepper_logistic(test)

    ! Logistic equation

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(logistic_ode_type), target :: logistic
    PetscInt,  parameter :: num_cases = 2
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [0.05_dp, 0.006_dp]
    PetscReal, allocatable :: initial(:)
    PetscErrorCode :: err

    allocate(logistic%c(logistic%dim), initial(logistic%dim))
    logistic%c = [0.0_dp, 0.5_dp, 1._dp, 1.5_dp, 2._dp, &
         2.5_dp, 3._dp, 3.5_dp]
    initial = 3._dp * logistic%c / (1._dp + 2._dp * logistic%c)
    call logistic%init(initial, test = test, err = err)

    json_str(1) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"size": 0.1, "maximum": {"number": 100}, ' // &
         '"method": "beuler", ' // &
         '"adapt": {"on": true, "method": "change", ' // &
         '"minimum": 0.01, "maximum": 0.2}}}}'

    json_str(2) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"size": 0.1, "maximum": {"number": 100}, ' // &
         '"method": "bdf2", ' // &
         '"adapt": {"on": true, "method": "change", ' // &
         '"minimum": 0.01, "maximum": 0.2}}}}'

    call logistic%run_cases(json_str, tol)

    call logistic%destroy()
    deallocate(initial)

  end subroutine test_timestepper_logistic

!------------------------------------------------------------------------

  subroutine test_timestepper_nontrivial_lhs(test)

    ! Nontrivial LHS

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(nontrivial_lhs_ode_type), target :: nontrivial_lhs
    PetscInt,  parameter :: num_cases = 2
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [0.1_dp, 0.02_dp]
    PetscReal, allocatable :: initial(:)
    PetscErrorCode :: err

    nontrivial_lhs%k = -1._dp
    allocate(initial(nontrivial_lhs%dim))
    initial = [-4._dp, -3.0_dp, -2.0_dp, -1.0_dp, &
         0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    call nontrivial_lhs%init(initial, test = test, err = err)

    json_str(1) = '{"time": {"start": 1.0, "stop": 10.0, ' // &
         '"step": {"size": 0.01, "maximum": {"number": 100}, ' // &
         '"method": "beuler", ' // &
         '"adapt": {"on": true, "method": "change", ' // &
         '"minimum": 0.03, "maximum": 0.1}}}}'

    json_str(2) = '{"time": {"start": 1.0, "stop": 10.0, ' // &
         '"step": {"size": 0.1, "maximum": {"number": 100}, ' // &
         '"method": "bdf2", ' // &
         '"adapt": {"on": true, "method": "change", ' // &
         '"minimum": 0.05, "maximum": 0.1}}}}'

    call nontrivial_lhs%run_cases(json_str, tol)

    call nontrivial_lhs%destroy()
    deallocate(initial)

  end subroutine test_timestepper_nontrivial_lhs

!------------------------------------------------------------------------

  subroutine test_timestepper_nonlinear_lhs(test)

    ! Nonlinear LHS

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(nonlinear_lhs_ode_type), target :: nonlinear_lhs
    PetscInt,  parameter :: num_cases = 2
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [0.15_dp, 0.05_dp]
    PetscReal, allocatable :: initial(:)
    PetscErrorCode :: err

    allocate(initial(nonlinear_lhs%dim))
    initial = [-4._dp, -3.0_dp, -2.0_dp, -1.0_dp, &
         0.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    call nonlinear_lhs%init(initial, test = test, err = err)

    json_str(1) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"size": 0.01, "maximum": {"number": 200}, ' // &
         '"method": "beuler", ' // &
         '"adapt": {"on": true, "method": "change", ' // &
         '"minimum": 0.03, "maximum": 0.1}}}}'

    json_str(2) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"size": 0.04, "maximum": {"number": 200}, ' // &
         '"method": "bdf2", ' // &
         '"adapt": {"on": true, "method": "change", ' // &
         '"minimum": 0.05, "maximum": 0.1}}}}'

    call nonlinear_lhs%run_cases(json_str, tol)

    call nonlinear_lhs%destroy()
    deallocate(initial)

  end subroutine test_timestepper_nonlinear_lhs

!------------------------------------------------------------------------

  subroutine test_timestepper_heat1d(test)

    ! 1-D heat equation PDE
    ! Solved in its usual form: dc/dt = d2c/dx2

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(heat1d_ode_type), target :: heat1d
    PetscInt,  parameter :: num_cases = 2
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [0.10_dp, 0.05_dp]
    PetscErrorCode :: err

    call heat1d%init(test = test, err = err)

    json_str(1) = '{"time": {"start": 0.0, "stop": 0.2, ' // &
         '"step": {"size": 0.01, "maximum": {"number": 20}, ' // &
         '"method": "beuler", ' // &
         '"adapt": {"on": true, "method": "change", ' // &
         '"minimum": 0.01, "maximum": 0.1}}}}'

    json_str(2) = '{"time": {"start": 0.0, "stop": 0.2, ' // &
         '"step": {"size": 0.01, "maximum": {"number": 20}, ' // &
         '"method": "bdf2", ' // &
         '"adapt": {"on": true, "method": "change", ' // &
         ' "minimum": 0.01, "maximum": 0.2}}}}'

    call heat1d%run_cases(json_str, tol)

    call heat1d%destroy()

  end subroutine test_timestepper_heat1d

!------------------------------------------------------------------------

  subroutine test_timestepper_heat1d_nonlinear(test)

    ! 1-D heat equation PDE (nonlinear form)
    ! d/dt(y2) = 2 y d2y/dx2

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(heat1d_nonlinear_ode_type), target :: heat1d
    PetscInt,  parameter :: num_cases = 2
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [0.10_dp, 0.05_dp]
    PetscErrorCode :: err

    call heat1d%init(test = test, err = err)

    json_str(1) = '{"time": {"start": 0.0, "stop": 0.2, ' // &
         '"step": {"size": [0.005], "maximum": {"number": 100}, ' // &
         '"method": "beuler", ' // &
         '"adapt": {"on": true, "method": "change", ' // &
         '"minimum": 0.01, "maximum": 0.1}}}}'

    json_str(2) = '{"time": {"start": 0.0, "stop": 0.2, ' // &
         '"step": {"size": 0.01, "maximum": {"number": 100}, ' // &
         '"method": "bdf2", ' // &
         '"adapt": {"on": true, "method": "change", ' // &
         '"minimum": 0.01, "maximum": 0.2}}}}'

    call heat1d%run_cases(json_str, tol)

    call heat1d%destroy()

  end subroutine test_timestepper_heat1d_nonlinear

!------------------------------------------------------------------------
  
  subroutine test_timestepper_pre_eval(test)

    ! Pre-evaluation procedure problem
    ! (This is just a re-casting of test_timestepper_nontrivial_lhs,
    ! with the LHS function essentially computed in the pre-evaluation
    ! routine, and just copied in the actual LHS function routine.)

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(pre_eval_ode_type), target :: pre_eval
    PetscInt,  parameter :: num_cases = 2
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [0.1_dp, 0.02_dp]
    PetscReal, allocatable :: initial(:)
    PetscErrorCode :: err

    pre_eval%k = -1._dp
    allocate(initial(pre_eval%dim))
    initial = [-4._dp, -3.0_dp, -2.0_dp, -1.0_dp, &
         0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    call pre_eval%init(initial, test = test, err = err)

    json_str(1) = '{"time": {"start": 1.0, "stop": 10.0, ' // &
         '"step": {"size": 0.01, "maximum": {"number": 100}, ' // &
         '"method": "beuler", ' // &
         '"adapt": {"on": true, "method": "change", ' // &
         '"minimum": 0.03, "maximum": 0.1}}}}'

    json_str(2) = '{"time": {"start": 1.0, "stop": 10.0, ' // &
         '"step": {"size": 0.1, "maximum": {"number": 100}, ' // &
         '"method": "bdf2", ' // &
         '"adapt": {"on": true, "method": "change", ' // &
         '"minimum": 0.05, "maximum": 0.1}}}}'

    call pre_eval%run_cases(json_str, tol)

    call pre_eval%destroy()
    deallocate(initial)

  end subroutine test_timestepper_pre_eval

!------------------------------------------------------------------------

  subroutine test_timestepper_steady(test)

    ! Steady state solution of dy/dt = 1 - y2

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(ss_ode_type), target :: ss
    PetscInt,  parameter :: num_cases = 1
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [1.e-6_dp]
    PetscReal, allocatable :: initial(:)
    PetscErrorCode :: err

    allocate(initial(ss%dim))
    initial = 0._dp
    call ss%init(initial, test = test, err = err)

    json_str(1) = '{"time": {"step": {"method": "directss"}}, ' // &
         '"output": {"initial": false}}'

    call ss%run_cases(json_str, tol)

    call ss%destroy()
    deallocate(initial)

  end subroutine test_timestepper_steady

!------------------------------------------------------------------------

  subroutine test_checkpoints(test)

    ! Timestepper checkpoints

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(timestepper_checkpoints_type) :: checkpoints
    PetscReal, allocatable :: times(:)
    PetscInt :: repeat
    PetscReal :: start_time
    PetscReal, parameter :: tolerance = 0.01_dp
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       times = [1._dp, 3._dp, 4._dp]
       repeat = 2
       start_time = 0._dp
       call checkpoints%init(times, repeat, tolerance, start_time)

       call test%assert(.not. checkpoints%done, 'initial done')
       call test%assert(times(1), checkpoints%next_time, 'initial next_time')

       call checkpoints%check(0.4_dp, 0.4_dp)
       call test%assert(.not. checkpoints%hit, 't = 0.4 hit')

       call checkpoints%check(1.0_dp, 0.6_dp)
       call test%assert(checkpoints%hit, 't = 1.0 hit')
       call checkpoints%update()

       call checkpoints%check(2.5_dp, 0.8_dp)
       call test%assert(.not. checkpoints%hit, 't = 2.5 hit')

       call checkpoints%check(3.2_dp, 1.0_dp)
       call test%assert(checkpoints%hit, 't = 3.2 hit')
       call checkpoints%update()

       call checkpoints%check(4.5_dp, 1.3_dp)
       call test%assert(checkpoints%hit, 't = 4.5 hit')
       call checkpoints%update()
       call test%assert(5._dp, checkpoints%next_time, 't = 4.5 next_time')
       call test%assert(2, checkpoints%repeat_index, 't = 4.5 repeat_index')

       call checkpoints%check(5.1_dp, 1.2_dp)
       call test%assert(checkpoints%hit, 't = 5.1 hit')
       call checkpoints%update()

       call checkpoints%update()
       call test%assert(3, checkpoints%index, 'last index')
       call test%assert(8._dp, checkpoints%next_time, 'last next_time')
       call checkpoints%check(7.9_dp, 1.5_dp)
       call test%assert(.not. checkpoints%hit, 't = 7.9 hit')
       call checkpoints%check(7.99_dp, 1.5_dp)
       call test%assert(checkpoints%hit, 't = 7.99 hit')

       call checkpoints%update()
       call test%assert(checkpoints%done, 'last done')

       call checkpoints%destroy()

       ! Test indefinite repeating:
       times = [1._dp]
       repeat = -1
       call checkpoints%init(times, repeat, tolerance, start_time)
       call test%assert(.not. checkpoints%done, 'indefinite repeat done 1')
       call checkpoints%update()
       call test%assert(.not. checkpoints%done, 'indefinite repeat done 2')
       call checkpoints%update()
       call test%assert(.not. checkpoints%done, 'indefinite repeat done 3')
       call checkpoints%update()
       call checkpoints%check(4.1_dp, 0.5_dp)
       call test%assert(checkpoints%hit, 'indefinite repeat t = 4.1 hit')
       call checkpoints%destroy()

       ! Test start time after first checkpoint time:
       times = [1._dp, 2._dp, 4._dp, 8._dp]
       start_time = 3._dp
       repeat = -1
       call checkpoints%init(times, repeat, tolerance, start_time)
       call test%assert(.not. checkpoints%done, 'start time done 1')
       call test%assert(4._dp, checkpoints%next_time, 'start time first next time')
       call checkpoints%update()
       call test%assert(.not. checkpoints%done, 'start time done 2')
       call test%assert(8._dp, checkpoints%next_time, 'start time t = 4 next time')
       call checkpoints%destroy()

       ! Test start time after last checkpoint time:
       start_time = 10._dp
       repeat = 1
       call checkpoints%init(times, repeat, tolerance, start_time)
       call test%assert(checkpoints%done, 'start time 2 done 1')
       call checkpoints%destroy()

       ! Test start time after last checkpoint time with repeat:
       times = [1._dp]
       repeat = 10
       start_time = 5._dp
       call checkpoints%init(times, repeat, tolerance, start_time)
       call test%assert(.not. checkpoints%done, 'start time 3 done 1')
       call test%assert(5._dp, checkpoints%next_time, 'start time 3 first next time')
       call checkpoints%destroy()

       deallocate(times)

    end if

  end subroutine test_checkpoints

!------------------------------------------------------------------------

end module timestepper_test
