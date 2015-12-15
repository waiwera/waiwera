module timestepper_test

  ! Tests for timestepper module

  use kinds_module
  use mpi_module
  use fruit
  use ode_module
  use timestepper_module
  use fson
  use fson_mpi_module

  implicit none

  private

#include <petsc/finclude/petsc.h90>

  type, extends(ode_type) :: test_ode_type
     private
     PetscInt, public :: dim = 8, dof = 1, stencil = 0
     PetscReal, allocatable, public :: initial_values(:)
     PetscReal, public :: start_time
     PetscReal, public :: maxdiff
     Vec, public :: exact_solution, diff
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
     procedure, public :: rhs => rhs_const
     procedure, public :: exact => exact_linear
  end type linear_ode_type

  type, extends(test_ode_type) :: exponential_ode_type
     private
     PetscReal :: k = -5._dp
   contains
     private
     procedure, public :: rhs => rhs_linear
     procedure, public :: exact => exact_exponential
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

  PetscInt, parameter :: max_json_len = 256
  PetscReal, parameter :: time_tolerance = 1.e-6_dp

public :: test_timestepper_linear, test_timestepper_exponential, &
     test_timestepper_logistic, test_timestepper_nontrivial_lhs, &
     test_timestepper_nonlinear_lhs, test_timestepper_heat1d, &
     test_timestepper_heat1d_nonlinear, test_timestepper_pre_eval, &
     test_timestepper_read

contains

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

  subroutine init_test_ode(self, initial_array, err)
    ! Default ode initialization.
    class(test_ode_type), intent(in out) :: self
    PetscReal, intent(in), optional :: initial_array(:)
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr

    err = 0
    self%dim = 8
    self%dof = 1
    self%stencil = 0

    call DMDACreate1d(mpi%comm, DM_BOUNDARY_NONE, self%dim, self%dof, &
         self%stencil, PETSC_NULL_INTEGER, self%mesh%dm, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(self%mesh%dm, self%solution, ierr); CHKERRQ(ierr)
    call VecDuplicate(self%solution, self%exact_solution, ierr); CHKERRQ(ierr)
    call VecDuplicate(self%solution, self%diff, ierr); CHKERRQ(ierr)
    if (present(initial_array)) then
       self%initial_values = initial_array
    end if

  end subroutine init_test_ode

  subroutine destroy_test_ode(self)
    ! Default ode destroy method.
    class(test_ode_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    call VecDestroy(self%solution, ierr); CHKERRQ(ierr)
    call VecDestroy(self%exact_solution, ierr); CHKERRQ(ierr)
    call VecDestroy(self%diff, ierr); CHKERRQ(ierr)
    call DMDestroy(self%mesh%dm, ierr); CHKERRQ(ierr)
  end subroutine destroy_test_ode

  subroutine lhs_test_ode(self, t, y, lhs, err)
    ! Default identity LHS function.
    class(test_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    Vec, intent(out) :: lhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    call VecCopy(y, lhs, ierr); CHKERRQ(ierr)
  end subroutine lhs_test_ode

  subroutine rhs_test_ode(self, t, y, rhs, err)
    ! Default zero RHS function.
    class(test_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    Vec, intent(out) :: rhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    call VecSet(rhs, 0._dp, ierr); CHKERRQ(ierr)
  end subroutine rhs_test_ode

  subroutine pre_eval_test_ode(self, t, y, err)
    ! Default do-nothing pre-evaluation routine.
    class(test_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    PetscErrorCode, intent(out) :: err
    continue
  end subroutine pre_eval_test_ode

  subroutine exact_test_ode(self, t, v)
    ! Default zero exact solution.
    class(test_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(out) :: v
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
    PetscReal, pointer :: y(:), yex(:), diffa(:)
    PetscErrorCode :: ierr
    PetscReal :: normdiff
    PetscInt :: i, local_size, low, hi

    call self%exact(self%time, self%exact_solution)
    call VecGetLocalSize(self%exact_solution, local_size, ierr)
    CHKERRQ(ierr)
    call VecGetOwnershipRange(self%exact_solution, low, hi, ierr)
    CHKERRQ(ierr)
    call VecGetArrayReadF90(self%solution, y, ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(self%exact_solution, yex, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(self%diff, diffa, ierr); CHKERRQ(ierr)
    do i = 1, local_size
       diffa(i) = relerr(yex(i), y(i))
    end do
    call VecRestoreArrayF90(self%diff, diffa, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%exact_solution, yex, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%solution, y, ierr)
    CHKERRQ(ierr)
    call VecNorm(self%diff, NORM_INFINITY, normdiff, ierr)
    CHKERRQ(ierr)
    self%maxdiff = max(self%maxdiff, normdiff)

  contains

    PetscReal function relerr(exact, val)
      ! Calculates relative error of val with respect to exact, unless
      ! exact is near zero, in which case the absolute difference is returned.
      PetscReal, intent(in) :: exact, val
      ! Locals:
      PetscReal, parameter :: tol = 1.e-6_dp
      if (abs(exact) > tol) then
         relerr = abs((val - exact) / exact)
      else
         relerr = abs(val - exact)
      end if
    end function relerr

  end subroutine output_test_ode

!------------------------------------------------------------------------
! Specific test ode functions
!------------------------------------------------------------------------

  subroutine rhs_const(self, t, y, rhs, err)
    ! rhs(t, y) = k
    class(linear_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    Vec, intent(out) :: rhs
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
    Vec, intent(out) :: v
    call VecSetArray(v, self%initial_values + &
         self%k * (t - self%start_time))
  end subroutine exact_linear

!------------------------------------------------------------------------
  
  subroutine rhs_linear(self, t, y, rhs, err)
    ! rhs(t, y) = k * y
    class(exponential_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    Vec, intent(out) :: rhs
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
    Vec, intent(out) :: v
    call VecSetArray(v, self%initial_values * &
         exp(self%k * (t - self%start_time)))
  end subroutine exact_exponential

!------------------------------------------------------------------------

  subroutine rhs_logistic(self, t, y, rhs, err)
    ! rhs(t, y) = (3 - 2 * y) * y
    class(logistic_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    Vec, intent(out) :: rhs
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
    Vec, intent(out) :: v
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

  subroutine lhs_nontrivial_lhs(self, t, y, lhs, err)
    ! lhs(t, y) = t * y
    class(nontrivial_lhs_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    Vec, intent(out) :: lhs
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
    Vec, intent(out) :: v
    call VecSetArray(v, self%initial_values * t ** (self%k - 1._dp))
  end subroutine exact_nontrivial_lhs

!------------------------------------------------------------------------
  
  subroutine lhs_nonlinear_lhs(self, t, y, lhs, err)
    ! lhs(t, y) = y * (y - 2)
    class(nonlinear_lhs_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    Vec, intent(out) :: lhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call VecPointwiseMult(lhs, y, y, ierr); CHKERRQ(ierr)
    call VecAXPY(lhs, -2._dp, y, ierr); CHKERRQ(ierr)
  end subroutine lhs_nonlinear_lhs

  subroutine rhs_nonlinear_lhs(self, t, y, rhs, err)
    ! rhs(t, y) = k * (y - 1)
    class(nonlinear_lhs_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    Vec, intent(out) :: rhs
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
    Vec, intent(out) :: v
    call VecSetArray(v, self%initial_values + &
         0.5_dp * self%k * t)
  end subroutine exact_nonlinear_lhs

!------------------------------------------------------------------------

  PetscReal function fn_heat1d(self, x, t)
    ! Exact heat equation solution at given point and time.
    class(heat1d_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: x, t
    ! Locals:
    PetscReal, parameter :: pi = 4._dp * atan(1._dp)
    fn_heat1d = sin(pi * x) * exp(-pi * pi * t)
  end function fn_heat1d

  subroutine exact_heat1d(self, t, v)
    ! Exact heat equation solution vector on mesh.
    class(heat1d_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(out) :: v
    ! Locals:
    DM :: cdm
    Vec :: cv
    PetscReal, pointer :: va(:), coords(:)
    PetscInt :: i1, im, i
    PetscReal :: x
    PetscErrorCode :: ierr

    call DMGetCoordinateDM(self%mesh%dm, cdm, ierr); CHKERRQ(ierr)
    call DMGetCoordinates(self%mesh%dm, cv, ierr); CHKERRQ(ierr)
    call DMDAVecGetArrayF90(cdm, cv, coords, ierr); CHKERRQ(ierr)
    call DMDAVecGetArrayF90(self%mesh%dm, v, va, ierr); CHKERRQ(ierr)
    call DMDAGetCorners(self%mesh%dm, i1, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, im, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    do i = i1, i1+im-1
       x = coords(i)
       va(i) = self%fn(x, t)
    end do
    call DMDAVecRestoreArrayF90(self%mesh%dm, v, va, ierr)
    CHKERRQ(ierr)
    call DMDAVecRestoreArrayF90(cdm, cv, coords, ierr); CHKERRQ(ierr)

  end subroutine exact_heat1d

  subroutine init_heat1d(self, initial_array, err)
    ! Initialization for heat equation.
    class(heat1d_ode_type), intent(in out) :: self
    PetscReal, intent(in), optional :: initial_array(:)
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal :: dx
    PetscErrorCode :: ierr

    err = 0
    self%L = 1._dp
    self%dim = 21
    self%dof = 1
    self%stencil = 1

    dx = self%L / (self%dim - 1._dp)
    self%a = 1._dp / (dx*dx)

    call DMDACreate1d(mpi%comm, DM_BOUNDARY_GHOSTED, self%dim-2, &
         self%dof, self%stencil, PETSC_NULL_INTEGER, self%mesh%dm, &
         ierr); CHKERRQ(ierr)
    call DMDASetUniformCoordinates(self%mesh%dm, dx, self%L - dx, &
         0._dp, 0._dp, 0._dp, 0._dp, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(self%mesh%dm, self%solution, ierr)
    CHKERRQ(ierr)
    call VecDuplicate(self%solution, self%exact_solution, ierr)
    CHKERRQ(ierr)
    call VecDuplicate(self%solution, self%diff, ierr); CHKERRQ(ierr)

  end subroutine init_heat1d

  subroutine rhs_heat1d(self, t, y, rhs, err)
    ! rhs = d2y/dx2
    class(heat1d_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    Vec, intent(out) :: rhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    PetscReal, pointer :: ya(:), rhsa(:)
    PetscInt :: i1, im, i, i1g, img, i2g
    Vec :: ylocal

    err = 0
    call DMGetLocalVector(self%mesh%dm, ylocal, ierr); CHKERRQ(ierr)
    call DMGlobalToLocalBegin(self%mesh%dm, y, INSERT_VALUES, ylocal, &
         ierr); CHKERRQ(ierr)
    call DMGlobalToLocalEnd(self%mesh%dm, y, INSERT_VALUES, ylocal, &
         ierr); CHKERRQ(ierr)

    call DMDAGetCorners(self%mesh%dm, i1, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, im, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    call DMDAVecGetArrayF90(self%mesh%dm, ylocal, ya, ierr)
    CHKERRQ(ierr)
    call DMDAVecGetArrayF90(self%mesh%dm, rhs, rhsa, ierr)
    CHKERRQ(ierr)

    ! BCs:
    call DMDAGetGhostCorners(self%mesh%dm, i1g, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, img, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    i2g = i1g+img-1
    if (i1g == 0) ya(i1g) = self%fn(0._dp, t)
    if (i2g == self%dim-2) ya(i2g) = self%fn(self%L, t)

    do i = i1, i1+im-1
       rhsa(i) = self%a * (ya(i-1) - 2._dp * ya(i) + ya(i+1))
    end do

    call DMDAVecRestoreArrayF90(self%mesh%dm, ylocal, ya, ierr)
    CHKERRQ(ierr)
    call DMRestoreLocalVector(self%mesh%dm, ylocal, ierr); CHKERRQ(ierr)

    call DMDAVecRestoreArrayF90(self%mesh%dm, rhs, rhsa, ierr)
    CHKERRQ(ierr)

  end subroutine rhs_heat1d

!------------------------------------------------------------------------

  subroutine lhs_heat1d_nonlinear(self, t, y, lhs, err)
    ! lhs = y * y
    class(heat1d_nonlinear_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    Vec, intent(out) :: lhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    PetscReal, pointer :: ya(:), lhsa(:)
    PetscInt :: i1, im, i
    Vec :: ylocal

    err = 0
    call DMGetLocalVector(self%mesh%dm, ylocal, ierr); CHKERRQ(ierr)
    call DMGlobalToLocalBegin(self%mesh%dm, y, INSERT_VALUES, ylocal, &
         ierr); CHKERRQ(ierr)
    call DMGlobalToLocalEnd(self%mesh%dm, y, INSERT_VALUES, ylocal, &
         ierr); CHKERRQ(ierr)

    call DMDAGetCorners(self%mesh%dm, i1, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, im, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    call DMDAVecGetArrayF90(self%mesh%dm, ylocal, ya, ierr)
    CHKERRQ(ierr)
    call DMDAVecGetArrayF90(self%mesh%dm, lhs, lhsa, ierr)
    CHKERRQ(ierr)

    do i = i1, i1+im-1
       lhsa(i) = ya(i) * ya(i)
    end do

    call DMDAVecRestoreArrayF90(self%mesh%dm, ylocal, ya, ierr)
    CHKERRQ(ierr)
    call DMRestoreLocalVector(self%mesh%dm, ylocal, ierr); CHKERRQ(ierr)

    call DMDAVecRestoreArrayF90(self%mesh%dm, lhs, lhsa, ierr)
    CHKERRQ(ierr)

  end subroutine lhs_heat1d_nonlinear
  
  subroutine rhs_heat1d_nonlinear(self, t, y, rhs, err)
    ! rhs = 2 * y * d2y/dx2
    class(heat1d_nonlinear_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    Vec, intent(out) :: rhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    PetscReal, pointer :: ya(:), rhsa(:)
    PetscInt :: i1, im, i, i1g, img, i2g
    Vec :: ylocal

    err = 0
    call DMGetLocalVector(self%mesh%dm, ylocal, ierr); CHKERRQ(ierr)
    call DMGlobalToLocalBegin(self%mesh%dm, y, INSERT_VALUES, &
         ylocal, ierr); CHKERRQ(ierr)
    call DMGlobalToLocalEnd(self%mesh%dm, y, INSERT_VALUES, ylocal, &
         ierr); CHKERRQ(ierr)

    call DMDAGetCorners(self%mesh%dm, i1, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, im, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    call DMDAVecGetArrayF90(self%mesh%dm, ylocal, ya, ierr)
    CHKERRQ(ierr)
    call DMDAVecGetArrayF90(self%mesh%dm, rhs, rhsa, ierr)
    CHKERRQ(ierr)

    ! BCs:
    call DMDAGetGhostCorners(self%mesh%dm, i1g, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, img, PETSC_NULL_INTEGER, &
         PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
    i2g = i1g+img-1
    if (i1g == 0) ya(i1g) = self%fn(0._dp, t)
    if (i2g == self%dim-2) ya(i2g) = self%fn(self%L, t)

    do i = i1, i1+im-1
       rhsa(i) = 2._dp * ya(i) * self%a * (ya(i-1) - 2._dp * ya(i) + ya(i+1))
    end do

    call DMDAVecRestoreArrayF90(self%mesh%dm, ylocal, ya, ierr)
    CHKERRQ(ierr)
    call DMRestoreLocalVector(self%mesh%dm, ylocal, ierr); CHKERRQ(ierr)

    call DMDAVecRestoreArrayF90(self%mesh%dm, rhs, rhsa, ierr)
    CHKERRQ(ierr)

  end subroutine rhs_heat1d_nonlinear

!------------------------------------------------------------------------
  
  subroutine init_pre_eval(self, initial_array, err)
    class(pre_eval_ode_type), intent(in out) :: self
    PetscReal, intent(in), optional :: initial_array(:)
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call self%test_ode_type%init(initial_array, err)
    call VecDuplicate(self%solution, self%secondary, ierr); CHKERRQ(ierr)
  end subroutine init_pre_eval

  subroutine destroy_pre_eval(self)
    class(pre_eval_ode_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    call VecDestroy(self%secondary, ierr); CHKERRQ(ierr)
    call self%test_ode_type%destroy()
  end subroutine destroy_pre_eval

  subroutine calculate_secondary_pre_eval(self, t, y, err)
    ! Calculates secondary vector = t * y
    class(pre_eval_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call VecCopy(y, self%secondary, ierr); CHKERRQ(ierr)
    call VecScale(self%secondary, t, ierr); CHKERRQ(ierr)
  end subroutine calculate_secondary_pre_eval

  subroutine lhs_pre_eval(self, t, y, lhs, err)
    ! lhs(t, y) = secondary
    class(pre_eval_ode_type), intent(in out) :: self
    PetscReal, intent(in) :: t
    Vec, intent(in) :: y
    Vec, intent(out) :: lhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscErrorCode :: ierr
    err = 0
    call VecCopy(self%secondary, lhs, ierr); CHKERRQ(ierr)
  end subroutine lhs_pre_eval

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

    self%time = self%start_time
    call self%exact(self%time, self%solution)

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

    num_cases = size(json_str)

    do i = 1, num_cases

       json => fson_parse_mpi(str = trim(json_str(i)))

       call self%read_start_time(json)
       call self%set_initial_conditions()
       self%maxdiff = 0._dp

       call ts%init(json, self)
       ts%before_step_output => NULL()
       ts%after_step_output => timestepper_step_output
 
       call ts%run()

       if (mpi%rank == mpi%output_rank) then
          call assert_equals(ts%steps%stop_time, &
               ts%steps%current%time, time_tolerance, &
               trim(ts%method%name) // ' stop time')
          call assert_equals(0._dp, self%maxdiff, tol(i), &
               trim(ts%method%name) // ' max. relative error')
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

  subroutine test_timestepper_read
    ! Tests assigning timestepper parameters from JSON string.

    ! Locals:
    character(len = max_json_len) :: json_str
    type(fson_value), pointer :: json
    type(timestepper_type) :: ts
    type(test_ode_type) :: test_ode
    PetscReal, allocatable :: initial(:)
    PetscErrorCode :: err

    initial = [-4._dp, -3.0_dp, -2.0_dp, -1.0_dp, &
         0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]

    json_str = '{"time": {"stop": 1.0, ' // &
         '"step": {"initial": 0.01, "maximum": {"number": 200, ' // &
         '"size": 3.e6}, "method": "beuler", ' // &
         '"adapt": {"on": true, "method": "change", "min": 0.01, ' // &
         '"max": 0.2, "reduction": 0.6, "amplification": 1.9}}}}'

    json => fson_parse_mpi(str = trim(json_str))
    call test_ode%init(initial, err)
    call ts%init(json, test_ode)

    if (mpi%rank == mpi%output_rank) then
       call assert_equals(1.0_dp, ts%steps%stop_time, time_tolerance, &
            "Timestepper stop time")
       call assert_equals(0.01_dp, ts%steps%next_stepsize, &
            time_tolerance, "Timestepper initial stepsize")
       call assert_equals(200, ts%steps%max_num, "Timestepper max. num steps")
       call assert_equals("Backward Euler", ts%method%name, "Timestepper method")
       call assert_equals(.true., ts%steps%adaptor%on, "Timestepper adapt on")
       call assert_equals("change", trim(ts%steps%adaptor%name), "Timestepper adapt method")
       call assert_equals(0.01_dp, ts%steps%adaptor%monitor_min, "Timestepper monitor min")
       call assert_equals(0.2_dp, ts%steps%adaptor%monitor_max, "Timestepper monitor max")
       call assert_equals(0.6_dp, ts%steps%adaptor%reduction, "Timestepper monitor reduction")
       call assert_equals(1.9_dp, ts%steps%adaptor%amplification, "Timestepper monitor amplification")
       call assert_equals(3.e6_dp, ts%steps%adaptor%max_stepsize, "Timestepper max stepsize")
    end if

    call ts%destroy()
    call test_ode%destroy()
    call fson_destroy_mpi(json)

  end subroutine test_timestepper_read

!------------------------------------------------------------------------

  subroutine test_timestepper_linear

    ! Linear function

    type(linear_ode_type), target :: linear
    PetscInt,  parameter :: num_cases = 3
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [1.e-6_dp, 1.e-6_dp, 1.e-6_dp]
    PetscReal, allocatable :: initial(:)
    PetscErrorCode :: err

    initial = [-4._dp, -3.0_dp, -2.0_dp, -1.0_dp, &
         0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    call linear%init(initial, err)

    json_str(1) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"initial": 0.1, "maximum": {"number": 20}, ' // &
         '"method": "beuler", "adapt": {"on": false}}}}'

    json_str(2) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"initial": 0.1, "maximum": {"number": 20}, ' // &
         '"method": "bdf2", "adapt": {"on": false}}}}'

    json_str(3) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"maximum": {"number": 10}, ' // &
         '"method": "beuler", ' // &
         '"sizes": [0.1, 0.1, 0.2, 0.2, 0.3]}}}'

    call linear%run_cases(json_str, tol)

    call linear%destroy()
    deallocate(initial)

  end subroutine test_timestepper_linear

!------------------------------------------------------------------------

  subroutine test_timestepper_exponential

    ! Exponential function

    type(exponential_ode_type), target :: exponential
    PetscInt,  parameter :: num_cases = 2
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [0.15_dp, 0.05_dp]
    PetscReal, allocatable :: initial(:)
    PetscErrorCode :: err

    initial = [-4._dp, -3.0_dp, -2.0_dp, -1.0_dp, &
         0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    call exponential%init(initial, err)

    json_str(1) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"initial": 0.01, "maximum": {"number": 200}, ' // &
         '"method": "beuler", ' // &
         '"adapt": {"on": true, "min": 0.01, "max": 0.2}}}}'

    json_str(2) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"initial": 0.05, "maximum": {"number": 200}, ' // &
         '"method": "bdf2", ' // &
         '"adapt": {"on": true, "min": 0.01, "max": 0.2}}}}'

    call exponential%run_cases(json_str, tol)

    call exponential%destroy()
    deallocate(initial)

  end subroutine test_timestepper_exponential

!------------------------------------------------------------------------

  subroutine test_timestepper_logistic

    ! Logistic equation

    type(logistic_ode_type), target :: logistic
    PetscInt,  parameter :: num_cases = 2
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [0.05_dp, 0.006_dp]
    PetscReal, allocatable :: initial(:)
    PetscErrorCode :: err

    logistic%c = [0.0_dp, 0.5_dp, 1._dp, 1.5_dp, 2._dp, &
         2.5_dp, 3._dp, 3.5_dp]
    initial = 3._dp * logistic%c / (1._dp + 2._dp * logistic%c)
    call logistic%init(initial, err)

    json_str(1) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"initial": 0.1, "maximum": {"number": 100}, ' // &
         '"method": "beuler", ' // &
         '"adapt": {"on": true, "min": 0.01, "max": 0.2}}}}'

    json_str(2) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"initial": 0.1, "maximum": {"number": 100}, ' // &
         '"method": "bdf2", ' // &
         '"adapt": {"on": true, "min": 0.01, "max": 0.2}}}}'

    call logistic%run_cases(json_str, tol)

    call logistic%destroy()
    deallocate(initial)

  end subroutine test_timestepper_logistic

!------------------------------------------------------------------------

  subroutine test_timestepper_nontrivial_lhs

    ! Nontrivial LHS

    type(nontrivial_lhs_ode_type), target :: nontrivial_lhs
    PetscInt,  parameter :: num_cases = 2
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [0.1_dp, 0.02_dp]
    PetscReal, allocatable :: initial(:)
    PetscErrorCode :: err

    nontrivial_lhs%k = -1._dp
    initial = [-4._dp, -3.0_dp, -2.0_dp, -1.0_dp, &
         0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    call nontrivial_lhs%init(initial, err)

    json_str(1) = '{"time": {"start": 1.0, "stop": 10.0, ' // &
         '"step": {"initial": 0.01, "maximum": {"number": 100}, ' // &
         '"method": "beuler", ' // &
         '"adapt": {"on": true, "min": 0.03, "max": 0.1}}}}'

    json_str(2) = '{"time": {"start": 1.0, "stop": 10.0, ' // &
         '"step": {"initial": 1.0, "maximum": {"number": 100}, ' // &
         '"method": "bdf2", ' // &
         '"adapt": {"on": true, "min": 0.05, "max": 0.1}}}}'

    call nontrivial_lhs%run_cases(json_str, tol)

    call nontrivial_lhs%destroy()
    deallocate(initial)

  end subroutine test_timestepper_nontrivial_lhs

!------------------------------------------------------------------------

  subroutine test_timestepper_nonlinear_lhs

    ! Nonlinear LHS

    type(nonlinear_lhs_ode_type), target :: nonlinear_lhs
    PetscInt,  parameter :: num_cases = 2
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [0.15_dp, 0.05_dp]
    PetscReal, allocatable :: initial(:)
    PetscErrorCode :: err

    initial = [-4._dp, -3.0_dp, -2.0_dp, -1.0_dp, &
         0.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    call nonlinear_lhs%init(initial, err)

    json_str(1) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"initial": 0.01, "maximum": {"number": 200}, ' // &
         '"method": "beuler", ' // &
         '"adapt": {"on": true, "min": 0.03, "max": 0.1}}}}'

    json_str(2) = '{"time": {"start": 0.0, "stop": 1.0, ' // &
         '"step": {"initial": 0.04, "maximum": {"number": 200}, ' // &
         '"method": "bdf2", ' // &
         '"adapt": {"on": true, "min": 0.05, "max": 0.1}}}}'

    call nonlinear_lhs%run_cases(json_str, tol)

    call nonlinear_lhs%destroy()
    deallocate(initial)

  end subroutine test_timestepper_nonlinear_lhs

!------------------------------------------------------------------------

  subroutine test_timestepper_heat1d

    ! 1-D heat equation PDE
    ! Solved in its usual form: dc/dt = d2c/dx2

    type(heat1d_ode_type), target :: heat1d
    PetscInt,  parameter :: num_cases = 2
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [0.10_dp, 0.05_dp]
    PetscErrorCode :: err

    call heat1d%init(err = err)

    json_str(1) = '{"time": {"start": 0.0, "stop": 0.2, ' // &
         '"step": {"initial": 0.01, "maximum": {"number": 20}, ' // &
         '"method": "beuler", ' // &
         '"adapt": {"on": true, "min": 0.01, "max": 0.1}}}}'

    json_str(2) = '{"time": {"start": 0.0, "stop": 0.2, ' // &
         '"step": {"initial": 0.01, "maximum": {"number": 20}, ' // &
         '"method": "bdf2", ' // &
         '"adapt": {"on": true, "min": 0.01, "max": 0.2}}}}'

    call heat1d%run_cases(json_str, tol)

    call heat1d%destroy()

  end subroutine test_timestepper_heat1d

!------------------------------------------------------------------------

  subroutine test_timestepper_heat1d_nonlinear

    ! 1-D heat equation PDE (nonlinear form)
    ! d/dt(y2) = 2 y d2y/dx2

    type(heat1d_nonlinear_ode_type), target :: heat1d
    PetscInt,  parameter :: num_cases = 2
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [0.10_dp, 0.05_dp]
    PetscErrorCode :: err

    call heat1d%init(err = err)

    json_str(1) = '{"time": {"start": 0.0, "stop": 0.2, ' // &
         '"step": {"initial": 0.01, "maximum": {"number": 40}, ' // &
         '"method": "beuler", ' // &
         '"adapt": {"on": true, "min": 0.01, "max": 0.1}}}}'

    json_str(2) = '{"time": {"start": 0.0, "stop": 0.2, ' // &
         '"step": {"initial": 0.01, "maximum": {"number": 40}, ' // &
         '"method": "bdf2", ' // &
         '"adapt": {"on": true, "min": 0.01, "max": 0.2}}}}'

    call heat1d%run_cases(json_str, tol)

    call heat1d%destroy()

  end subroutine test_timestepper_heat1d_nonlinear

!------------------------------------------------------------------------
  
  subroutine test_timestepper_pre_eval

    ! Pre-evaluation procedure problem
    ! (This is just a re-casting of test_timestepper_nontrivial_lhs,
    ! with the LHS function essentially computed in the pre-evaluation
    ! routine, and just copied in the actual LHS function routine.)

    type(pre_eval_ode_type), target :: pre_eval
    PetscInt,  parameter :: num_cases = 2
    character(len = max_json_len) :: json_str(num_cases)
    PetscReal, parameter :: tol(num_cases) = [0.1_dp, 0.02_dp]
    PetscReal, allocatable :: initial(:)
    PetscErrorCode :: err

    pre_eval%k = -1._dp
    initial = [-4._dp, -3.0_dp, -2.0_dp, -1.0_dp, &
         0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    call pre_eval%init(initial, err)

    json_str(1) = '{"time": {"start": 1.0, "stop": 10.0, ' // &
         '"step": {"initial": 0.01, "maximum": {"number": 100}, ' // &
         '"method": "beuler", ' // &
         '"adapt": {"on": true, "min": 0.03, "max": 0.1}}}}'

    json_str(2) = '{"time": {"start": 1.0, "stop": 10.0, ' // &
         '"step": {"initial": 1.0, "maximum": {"number": 100}, ' // &
         '"method": "bdf2", ' // &
         '"adapt": {"on": true, "min": 0.05, "max": 0.1}}}}'

    call pre_eval%run_cases(json_str, tol)

    call pre_eval%destroy()
    deallocate(initial)

  end subroutine test_timestepper_pre_eval

!------------------------------------------------------------------------

end module timestepper_test
