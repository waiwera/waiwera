module timestepping_test

  ! Tests for timestepping module

  use kinds_module
  use mpi_module
  use fruit
  use timestepping_module

  implicit none

  private

#include <petsc-finclude/petscsys.h>
#include <petsc-finclude/petscdef.h>
#include <petsc-finclude/petscvec.h>
#include <petsc-finclude/petscvec.h90>
#include <petsc-finclude/petscsnes.h>
#include <petsc-finclude/petscdm.h>
#include <petsc-finclude/petscdmda.h>
#include <petsc-finclude/petscdmda.h90>

  PetscReal :: solution_tolerance, maxdiff
  Vec :: exact, diff

  interface
     subroutine exact_function(t, v)
       PetscReal, intent(in) :: t
       Vec, intent(out) :: v
     end subroutine exact_function
  end interface

procedure(lhs_function), pointer :: lhs_func
procedure(rhs_function), pointer :: rhs_func
procedure(exact_function), pointer :: exact_func

public :: test_linear, test_exponential, test_logistic, &
     test_nontrivial_lhs, test_nonlinear_lhs, test_heat1d

contains

!------------------------------------------------------------------------
! Utility routines
!------------------------------------------------------------------------

  subroutine VecSetArray(v, arr)

    ! Sets local part of vector v with values from array.

    Vec, intent(in out) :: v
    PetscReal, intent(in) :: arr(:)
    ! Locals:
    PetscInt :: hi, low, i
    PetscErrorCode :: ierr
    PetscReal, pointer :: va(:)

    call VecGetOwnershipRange(v, low, hi, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(v, va, ierr); CHKERRQ(ierr)
    va(1: hi-low) = arr(low+1: hi)
    call VecRestoreArrayF90(v, va, ierr); CHKERRQ(ierr)
    
  end subroutine VecSetArray

!------------------------------------------------------------------------

  subroutine odetest(dm, methods, t0, t1, dt, initial, max_steps, tol, &
       adaptive, eta_min, eta_max)

    ! Test ODE solver with specified dimension, methods, time parameters, initial
    ! conditions, solution tolerance, LHS / RHS functions and subroutine for output.

    DM, intent(in) :: dm
    PetscInt, intent(in) :: methods(:)
    PetscReal, intent(in) :: t0, t1, dt(:)
    Vec :: initial
    PetscInt, intent(in) :: max_steps(:)
    PetscReal, intent(in) :: tol(:)
    PetscBool, intent(in) :: adaptive
    PetscReal, intent(in) :: eta_min(:), eta_max(:)
    ! Locals:
    type(timestepper_type) :: ts
    PetscErrorCode :: ierr
    PetscReal, parameter :: time_tolerance = 1.e-6_dp
    PetscInt :: i
    PetscReal, parameter :: max_stepsize = 0._dp

    call VecDuplicate(initial, exact, ierr); CHKERRQ(ierr)
    call VecDuplicate(initial, diff, ierr); CHKERRQ(ierr)

    do i = 1, size(methods)
       solution_tolerance = tol(i)
       maxdiff = 0._dp
       call ts%init(methods(i), dm, lhs_func, rhs_func, t0, initial, dt(i), &
            t1, max_steps(i), max_stepsize)
       ts%step_output => step_output_compare
       ts%steps%adaptor%on = adaptive
       ts%steps%adaptor%monitor_min = eta_min(i)
       ts%steps%adaptor%monitor_max = eta_max(i)
       call ts%run()
       if (mpi%rank == mpi%output_rank) then
          call assert_equals(t1, ts%steps%current%time, time_tolerance, &
               ts%method%name // 'final time')
          call assert_equals(0._dp, maxdiff, solution_tolerance, &
               ts%method%name // 'max. relative error')
       end if
       call ts%destroy()
    end do

    call VecDestroy(exact, ierr); CHKERRQ(ierr)
    call VecDestroy(diff, ierr); CHKERRQ(ierr)

  end subroutine odetest

!------------------------------------------------------------------------

  subroutine step_output_compare(self)

    ! Compare result with exact solution vector at each time

    class(timestepper_type), intent(in out) :: self
    ! Locals:
    PetscReal, pointer :: y(:), yex(:), diffa(:)
    PetscErrorCode :: ierr
    PetscReal :: t, normdiff
    character(8) :: s
    PetscInt :: i, local_size, low, hi

    t = self%steps%current%time
    call VecGetArrayF90(self%steps%current%solution, y, ierr); CHKERRQ(ierr)
    call exact_func(t, exact)
    call VecGetArrayF90(exact, yex, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(diff, diffa, ierr); CHKERRQ(ierr)
    call VecGetLocalSize(exact, local_size, ierr); CHKERRQ(ierr)
    call VecGetOwnershipRange(exact, low, hi, ierr); CHKERRQ(ierr)
    do i = 1, local_size
       diffa(i) = relerr(yex(i), y(i))
    end do
    call VecRestoreArrayF90(diff, diffa, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(exact, yex, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(self%steps%current%solution, y, ierr); CHKERRQ(ierr)
    call VecNorm(diff, NORM_INFINITY, normdiff, ierr); CHKERRQ(ierr)
    maxdiff = max(maxdiff, normdiff)

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

  end subroutine step_output_compare

!------------------------------------------------------------------------
! Test routines
!------------------------------------------------------------------------

  subroutine test_linear

    ! Linear function

    PetscInt,  parameter :: dim = 8
    PetscReal, parameter :: k = -0.5_dp
    PetscReal, parameter :: t0 = 0._dp, t1 = 1._dp
    PetscReal, parameter :: initial(dim) = [-4._dp, -3.0_dp, -2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    PetscInt,  parameter :: num_methods = 2, methods(num_methods) = [TS_BEULER, TS_BDF2]
    PetscInt,  parameter :: max_steps(num_methods) = [20, 20]
    PetscReal, parameter :: tol(num_methods) = [1.e-6_dp, 1.e-6_dp]
    PetscReal, parameter :: dt(num_methods) = [0.1_dp, 0.1_dp]
    PetscBool, parameter :: adaptive = .false.
    PetscReal, parameter :: eta_min(num_methods) = [0.01_dp, 0.01_dp]
    PetscReal, parameter :: eta_max(num_methods) = [0.1_dp, 0.1_dp]
    PetscInt,  parameter :: dof = 1, stencil = 0
    PetscErrorCode :: ierr
    DM :: dm
    Vec :: initialv

    lhs_func => lhs_identity
    rhs_func => rhs_const
    exact_func => exact_linear

    call DMDACreate1d(mpi%comm, DM_BOUNDARY_NONE, dim, dof, stencil, &
         PETSC_NULL_INTEGER, dm, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(dm, initialv, ierr); CHKERRQ(ierr)
    call VecSetArray(initialv, initial)
    
    call odetest(dm, methods, t0, t1, dt, initialv, max_steps, &
         tol, adaptive, eta_min, eta_max)

    call VecDestroy(initialv, ierr); CHKERRQ(ierr)
    call DMDestroy(dm, ierr); CHKERRQ(ierr)

  contains

    subroutine rhs_const(t, y, rhs)
      ! rhs(t, y) = k
      PetscReal, intent(in) :: t
      Vec, intent(in) :: y
      Vec, intent(out) :: rhs
      ! Locals:
      PetscErrorCode :: ierr
      call VecSet(rhs, k, ierr); CHKERRQ(ierr)
    end subroutine rhs_const

    subroutine exact_linear(t, v)
      ! Linear solution
      PetscReal, intent(in) :: t
      Vec, intent(out) :: v
      call VecSetArray(v, initial + k * (t - t0))
    end subroutine exact_linear

  end subroutine test_linear

!------------------------------------------------------------------------

  subroutine test_exponential

    ! Exponential function

    PetscInt,  parameter :: dim = 8
    PetscReal, parameter :: k = -5._dp
    PetscReal, parameter :: t0 = 0._dp, t1 = 1._dp
    PetscReal, parameter :: initial(dim) = [-4._dp, -3.0_dp, -2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    PetscInt,  parameter :: num_methods = 2, methods(num_methods) = [TS_BEULER, TS_BDF2]
    PetscInt,  parameter :: max_steps(num_methods) = [200, 200]
    PetscReal, parameter :: tol(num_methods) = [0.15_dp, 0.05_dp]
    PetscReal, parameter :: dt(num_methods) = [0.01_dp, 0.05_dp]
    PetscBool, parameter :: adaptive = .true.
    PetscReal, parameter :: eta_min(num_methods) = [0.01_dp, 0.01_dp]
    PetscReal, parameter :: eta_max(num_methods) = [0.2_dp, 0.2_dp]
    PetscInt,  parameter :: dof = 1, stencil = 0
    PetscErrorCode :: ierr
    DM :: dm
    Vec :: initialv

    lhs_func => lhs_identity
    rhs_func => rhs_linear
    exact_func => exact_exp

    call DMDACreate1d(mpi%comm, DM_BOUNDARY_NONE, dim, dof, stencil, &
         PETSC_NULL_INTEGER, dm, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(dm, initialv, ierr); CHKERRQ(ierr)
    call VecSetArray(initialv, initial)
    
    call odetest(dm, methods, t0, t1, dt, initialv, max_steps, &
         tol, adaptive, eta_min, eta_max)

    call VecDestroy(initialv, ierr); CHKERRQ(ierr)
    call DMDestroy(dm, ierr); CHKERRQ(ierr)

  contains

    subroutine rhs_linear(t, y, rhs)
      ! rhs(t, y) = k * y
      PetscReal, intent(in) :: t
      Vec, intent(in) :: y
      Vec, intent(out) :: rhs
      ! Locals:
      PetscErrorCode :: ierr
      call VecCopy(y, rhs, ierr); CHKERRQ(ierr)
      call VecScale(rhs, k, ierr); CHKERRQ(ierr)
    end subroutine rhs_linear

    subroutine exact_exp(t, v)
      ! Exponential solution
      PetscReal, intent(in) :: t
      Vec, intent(out) :: v
      call VecSetArray(v, initial * exp(k * (t - t0)))
    end subroutine exact_exp

  end subroutine test_exponential

!------------------------------------------------------------------------

  subroutine test_logistic

    ! Logistic equation

    PetscInt,  parameter :: dim = 8
    PetscReal, parameter :: t0 = 0._dp, t1 = 1._dp
    PetscReal, parameter :: c(dim) = [0.0_dp, 0.5_dp, 1._dp, 1.5_dp, 2._dp, 2.5_dp, 3._dp, 3.5_dp]
    PetscReal, parameter :: initial(dim) = 3._dp * c / (1._dp + 2._dp * c)
    PetscInt,  parameter :: num_methods = 2, methods(num_methods) = [TS_BEULER, TS_BDF2]
    PetscInt,  parameter :: max_steps(num_methods) = [100, 100]
    PetscReal, parameter :: tol(num_methods) = [0.05_dp, 0.006_dp]
    PetscReal, parameter :: dt(num_methods) = [0.1_dp, 0.1_dp]
    PetscBool, parameter :: adaptive = .true.
    PetscReal, parameter :: eta_min(num_methods) = [0.01_dp, 0.01_dp]
    PetscReal, parameter :: eta_max(num_methods) = [0.2_dp, 0.2_dp]
    PetscInt,  parameter :: dof = 1, stencil = 0
    PetscErrorCode :: ierr
    DM :: dm
    Vec :: initialv

    lhs_func => lhs_identity
    rhs_func => rhs_logistic
    exact_func => exact_logistic

    call DMDACreate1d(mpi%comm, DM_BOUNDARY_NONE, dim, dof, stencil, &
         PETSC_NULL_INTEGER, dm, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(dm, initialv, ierr); CHKERRQ(ierr)
    call VecSetArray(initialv, initial)
    
    call odetest(dm, methods, t0, t1, dt, initialv, max_steps, &
         tol, adaptive, eta_min, eta_max)

    call VecDestroy(initialv, ierr); CHKERRQ(ierr)
    call DMDestroy(dm, ierr); CHKERRQ(ierr)

  contains

    subroutine rhs_logistic(t, y, rhs)
      ! rhs(t, y) = (3 - 2 * y) * y
      PetscReal, intent(in) :: t
      Vec, intent(in) :: y
      Vec, intent(out) :: rhs
      ! Locals:
      PetscErrorCode :: ierr
      call VecSet(rhs, 3._dp, ierr); CHKERRQ(ierr)
      call VecAXPY(rhs, -2._dp, y, ierr); CHKERRQ(ierr)
      call VecPointwiseMult(rhs, rhs, y, ierr); CHKERRQ(ierr)
    end subroutine rhs_logistic

    subroutine exact_logistic(t, v)
      ! Logistic solution 3 * c * exp(3 * t) / (1 + 2 * c * exp(3 * t))
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
         ce3t = c(ig) * exp(3._dp * (t - t0))
         va(i) = 3._dp * ce3t / (1._dp + 2._dp * ce3t)
      end do
      call VecRestoreArrayF90(v, va, ierr); CHKERRQ(ierr)
    end subroutine exact_logistic

  end subroutine test_logistic

!------------------------------------------------------------------------

  subroutine test_nontrivial_lhs

    ! Nontrivial LHS

    PetscInt,  parameter :: dim = 8
    PetscReal, parameter :: t0 = 1._dp, t1 = 10._dp
    PetscReal, parameter :: k = -1._dp
    PetscReal, parameter :: initial(dim) = [-4._dp, -3.0_dp, -2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    PetscInt,  parameter :: num_methods = 2, methods(num_methods) = [TS_BEULER, TS_BDF2]
    PetscInt,  parameter :: max_steps(num_methods) = [100, 100]
    PetscReal, parameter :: tol(num_methods) = [0.1_dp, 0.02_dp]
    PetscReal, parameter :: dt(num_methods) = [0.01_dp, 1._dp]
    PetscBool, parameter :: adaptive = .true.
    PetscReal, parameter :: eta_min(num_methods) = [0.03_dp, 0.05_dp]
    PetscReal, parameter :: eta_max(num_methods) = [0.1_dp, 0.1_dp]
    PetscInt,  parameter :: dof = 1, stencil = 0
    PetscErrorCode :: ierr
    DM :: dm
    Vec :: initialv

    lhs_func => lhs_fn
    rhs_func => rhs_fn
    exact_func => exact_fn

    call DMDACreate1d(mpi%comm, DM_BOUNDARY_NONE, dim, dof, stencil, &
         PETSC_NULL_INTEGER, dm, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(dm, initialv, ierr); CHKERRQ(ierr)
    call VecSetArray(initialv, initial)
    
    call odetest(dm, methods, t0, t1, dt, initialv, max_steps, &
         tol, adaptive, eta_min, eta_max)

    call VecDestroy(initialv, ierr); CHKERRQ(ierr)
    call DMDestroy(dm, ierr); CHKERRQ(ierr)

  contains

    subroutine lhs_fn(t, y, rhs)
      ! lhs(t, y) = t * y
      PetscReal, intent(in) :: t
      Vec, intent(in) :: y
      Vec, intent(out) :: rhs
      ! Locals:
      PetscErrorCode :: ierr
      call VecCopy(y, rhs, ierr); CHKERRQ(ierr)
      call VecScale(rhs, t, ierr); CHKERRQ(ierr)
    end subroutine lhs_fn

    subroutine rhs_fn(t, y, rhs)
      ! rhs(t, y) = k * y
      PetscReal, intent(in) :: t
      Vec, intent(in) :: y
      Vec, intent(out) :: rhs
      ! Locals:
      PetscErrorCode :: ierr
      call VecCopy(y, rhs, ierr); CHKERRQ(ierr)
      call VecScale(rhs, k, ierr); CHKERRQ(ierr)
    end subroutine rhs_fn

    subroutine exact_fn(t, v)
      ! Solution y = y0 * t ** (k-1)
      PetscReal, intent(in) :: t
      Vec, intent(out) :: v
      call VecSetArray(v, initial * t ** (k - 1._dp))
    end subroutine exact_fn

  end subroutine test_nontrivial_lhs

!------------------------------------------------------------------------

  subroutine test_nonlinear_lhs

    ! Nonlinear LHS

    PetscInt,  parameter :: dim = 8
    PetscReal, parameter :: t0 = 0._dp, t1 = 1._dp
    PetscReal, parameter :: k = -1._dp
    PetscReal, parameter :: initial(dim) = [-4._dp, -3._dp, -2.0_dp, -1.0_dp, 0.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    PetscInt,  parameter :: num_methods = 2, methods(num_methods) = [TS_BEULER, TS_BDF2]
    PetscInt,  parameter :: max_steps(num_methods) = [200, 200]
    PetscReal, parameter :: tol(num_methods) = [0.15_dp, 0.05_dp]
    PetscReal, parameter :: dt(num_methods) = [0.01_dp, 0.04_dp]
    PetscBool, parameter :: adaptive = .true.
    PetscReal, parameter :: eta_min(num_methods) = [0.03_dp, 0.05_dp]
    PetscReal, parameter :: eta_max(num_methods) = [0.1_dp, 0.1_dp]
    PetscInt,  parameter :: dof = 1, stencil = 0
    PetscErrorCode :: ierr
    DM :: dm
    Vec :: initialv

    lhs_func => lhs_fn
    rhs_func => rhs_fn
    exact_func => exact_fn

    call DMDACreate1d(mpi%comm, DM_BOUNDARY_NONE, dim, dof, stencil, &
         PETSC_NULL_INTEGER, dm, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(dm, initialv, ierr); CHKERRQ(ierr)
    call VecSetArray(initialv, initial)
    
    call odetest(dm, methods, t0, t1, dt, initialv, max_steps, &
         tol, adaptive, eta_min, eta_max)

    call VecDestroy(initialv, ierr); CHKERRQ(ierr)
    call DMDestroy(dm, ierr); CHKERRQ(ierr)

  contains

    subroutine lhs_fn(t, y, rhs)
      ! lhs(t, y) = y * (y - 2)
      PetscReal, intent(in) :: t
      Vec, intent(in) :: y
      Vec, intent(out) :: rhs
      ! Locals:
      PetscErrorCode :: ierr
      call VecPointwiseMult(rhs, y, y, ierr); CHKERRQ(ierr)
      call VecAXPY(rhs, -2._dp, y, ierr); CHKERRQ(ierr)
    end subroutine lhs_fn

    subroutine rhs_fn(t, y, rhs)
      ! rhs(t, y) = k * (y - 1)
      PetscReal, intent(in) :: t
      Vec, intent(in) :: y
      Vec, intent(out) :: rhs
      ! Locals:
      PetscErrorCode :: ierr
      call VecSet(rhs, -1._dp, ierr); CHKERRQ(ierr)
      call VecAXPY(rhs, 1._dp, y, ierr); CHKERRQ(ierr)
      call VecScale(rhs, k, ierr); CHKERRQ(ierr)
    end subroutine rhs_fn

    subroutine exact_fn(t, v)
      ! Linear solution y = y0 + k * t / 2
      PetscReal, intent(in) :: t
      Vec, intent(out) :: v
      call VecSetArray(v, initial + 0.5_dp * k * t)
    end subroutine exact_fn

  end subroutine test_nonlinear_lhs

!------------------------------------------------------------------------

  subroutine test_heat1d

    ! 1-D heat equation PDE
    ! Solved in its usual form: dc/dt = d2c/dx2
    ! and also in a nonlinear form: d/dt(y2) = 2 y d2y/dx2

    PetscInt,  parameter :: dim = 21
    PetscReal, parameter :: t0 = 0._dp, t1 = 0.2_dp
    PetscInt,  parameter :: num_methods = 2, methods(num_methods) = [TS_BEULER, TS_BDF2]
    PetscInt             :: max_steps(num_methods) = [20, 20]
    PetscReal, parameter :: tol(num_methods) = [0.1_dp, 0.05_dp]
    PetscReal, parameter :: dt(num_methods) = [0.01_dp, 0.01_dp]
    PetscBool, parameter :: adaptive = .true.
    PetscReal, parameter :: eta_min(num_methods) = [0.01_dp, 0.01_dp]
    PetscReal, parameter :: eta_max(num_methods) = [0.1_dp, 0.2_dp]
    PetscInt,  parameter :: dof = 1, stencil = 1
    PetscReal :: initial(dim)
    PetscErrorCode :: ierr
    DM :: dm
    Vec :: initialv
    PetscReal, parameter :: L = 1._dp
    PetscReal, parameter :: dx = L / (dim - 1._dp), a = 1._dp / (dx*dx)

    ! Usual form
    lhs_func => lhs_identity
    rhs_func => rhs_fn
    exact_func => exact_fn

    call DMDACreate1d(mpi%comm, DM_BOUNDARY_GHOSTED, dim-2, dof, stencil, &
         PETSC_NULL_INTEGER, dm, ierr); CHKERRQ(ierr)
    call DMDASetUniformCoordinates(dm, dx, L - dx, 0._dp, 0._dp, &
         0._dp, 0._dp, ierr); CHKERRQ(ierr)
    call DMCreateGlobalVector(dm, initialv, ierr); CHKERRQ(ierr)
    call exact_fn(0._dp, initialv)

    call odetest(dm, methods, t0, t1, dt, initialv, max_steps, &
         tol, adaptive, eta_min, eta_max)

    ! Nonlinear form
    lhs_func => lhs_fn_nonlinear
    rhs_func => rhs_fn_nonlinear
    max_steps = [40, 40]

    call odetest(dm, methods, t0, t1, dt, initialv, max_steps, &
         tol, adaptive, eta_min, eta_max)

    call VecDestroy(initialv, ierr); CHKERRQ(ierr)
    call DMDestroy(dm, ierr); CHKERRQ(ierr)

  contains

    real(dp) function fn(x, t)
      PetscReal, intent(in) :: x, t
      ! Locals:
      PetscReal, parameter :: pi = 4._dp * atan(1._dp)
      fn = sin(pi * x) * exp(-pi * pi * t)
    end function fn

    subroutine rhs_fn(t, y, rhs)

      PetscReal, intent(in) :: t
      Vec, intent(in) :: y
      Vec, intent(out) :: rhs
      ! Locals:
      PetscErrorCode :: ierr
      PetscReal, pointer :: ya(:), rhsa(:)
      PetscInt :: i1, im, i, i1g, img, i2g
      Vec :: ylocal, rlocal

      call DMGetLocalVector(dm, ylocal, ierr); CHKERRQ(ierr)
      call DMGetLocalVector(dm, rlocal, ierr); CHKERRQ(ierr)
      call DMGlobalToLocalBegin(dm, y, INSERT_VALUES, ylocal, ierr); CHKERRQ(ierr)
      call DMGlobalToLocalEnd(dm, y, INSERT_VALUES, ylocal, ierr); CHKERRQ(ierr)

      call DMDAGetCorners(dm, i1, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
           im, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
      call DMDAVecGetArrayF90(dm, ylocal, ya, ierr); CHKERRQ(ierr)
      call DMDAVecGetArrayF90(dm, rlocal, rhsa, ierr); CHKERRQ(ierr)

      ! BCs:
      call DMDAGetGhostCorners(dm, i1g, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
           img, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
      i2g = i1g+img-1
      if (i1g == 0) ya(i1g) = fn(0._dp, t)
      if (i2g == dim-2) ya(i2g) = fn(L, t)

      do i = i1, i1+im-1
         rhsa(i) = a * (ya(i-1) - 2._dp * ya(i) + ya(i+1))
      end do

      call DMDAVecRestoreArrayF90(dm, ylocal, ya, ierr); CHKERRQ(ierr)
      call DMRestoreLocalVector(dm, ylocal, ierr); CHKERRQ(ierr)

      call DMDAVecRestoreArrayF90(dm, rlocal, rhsa, ierr); CHKERRQ(ierr)
      call DMLocalToGlobalBegin(dm, rlocal, INSERT_VALUES, rhs, ierr); CHKERRQ(ierr)
      call DMLocalToGlobalEnd(dm, rlocal, INSERT_VALUES, rhs, ierr); CHKERRQ(ierr)
      call DMRestoreLocalVector(dm, rlocal, ierr); CHKERRQ(ierr)

    end subroutine rhs_fn

    subroutine lhs_fn_nonlinear(t, y, lhs)

      PetscReal, intent(in) :: t
      Vec, intent(in) :: y
      Vec, intent(out) :: lhs
      ! Locals:
      PetscErrorCode :: ierr
      PetscReal, pointer :: ya(:), lhsa(:)
      PetscInt :: i1, im, i
      Vec :: ylocal, llocal

      call DMGetLocalVector(dm, ylocal, ierr); CHKERRQ(ierr)
      call DMGetLocalVector(dm, llocal, ierr); CHKERRQ(ierr)
      call DMGlobalToLocalBegin(dm, y, INSERT_VALUES, ylocal, ierr); CHKERRQ(ierr)
      call DMGlobalToLocalEnd(dm, y, INSERT_VALUES, ylocal, ierr); CHKERRQ(ierr)

      call DMDAGetCorners(dm, i1, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
           im, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
      call DMDAVecGetArrayF90(dm, ylocal, ya, ierr); CHKERRQ(ierr)
      call DMDAVecGetArrayF90(dm, llocal, lhsa, ierr); CHKERRQ(ierr)

      do i = i1, i1+im-1
         lhsa(i) = ya(i) * ya(i)
      end do

      call DMDAVecRestoreArrayF90(dm, ylocal, ya, ierr); CHKERRQ(ierr)
      call DMRestoreLocalVector(dm, ylocal, ierr); CHKERRQ(ierr)

      call DMDAVecRestoreArrayF90(dm, llocal, lhsa, ierr); CHKERRQ(ierr)
      call DMLocalToGlobalBegin(dm, llocal, INSERT_VALUES, lhs, ierr); CHKERRQ(ierr)
      call DMLocalToGlobalEnd(dm, llocal, INSERT_VALUES, lhs, ierr); CHKERRQ(ierr)
      call DMRestoreLocalVector(dm, llocal, ierr); CHKERRQ(ierr)

    end subroutine lhs_fn_nonlinear

    subroutine rhs_fn_nonlinear(t, y, rhs)

      PetscReal, intent(in) :: t
      Vec, intent(in) :: y
      Vec, intent(out) :: rhs
      ! Locals:
      PetscErrorCode :: ierr
      PetscReal, pointer :: ya(:), rhsa(:)
      PetscInt :: i1, im, i, i1g, img, i2g
      Vec :: ylocal, rlocal

      call DMGetLocalVector(dm, ylocal, ierr); CHKERRQ(ierr)
      call DMGetLocalVector(dm, rlocal, ierr); CHKERRQ(ierr)
      call DMGlobalToLocalBegin(dm, y, INSERT_VALUES, ylocal, ierr); CHKERRQ(ierr)
      call DMGlobalToLocalEnd(dm, y, INSERT_VALUES, ylocal, ierr); CHKERRQ(ierr)

      call DMDAGetCorners(dm, i1, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
           im, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
      call DMDAVecGetArrayF90(dm, ylocal, ya, ierr); CHKERRQ(ierr)
      call DMDAVecGetArrayF90(dm, rlocal, rhsa, ierr); CHKERRQ(ierr)

      ! BCs:
      call DMDAGetGhostCorners(dm, i1g, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
           img, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
      i2g = i1g+img-1
      if (i1g == 0) ya(i1g) = fn(0._dp, t)
      if (i2g == dim-2) ya(i2g) = fn(L, t)

      do i = i1, i1+im-1
         rhsa(i) = 2._dp * ya(i) * a * (ya(i-1) - 2._dp * ya(i) + ya(i+1))
      end do

      call DMDAVecRestoreArrayF90(dm, ylocal, ya, ierr); CHKERRQ(ierr)
      call DMRestoreLocalVector(dm, ylocal, ierr); CHKERRQ(ierr)

      call DMDAVecRestoreArrayF90(dm, rlocal, rhsa, ierr); CHKERRQ(ierr)
      call DMLocalToGlobalBegin(dm, rlocal, INSERT_VALUES, rhs, ierr); CHKERRQ(ierr)
      call DMLocalToGlobalEnd(dm, rlocal, INSERT_VALUES, rhs, ierr); CHKERRQ(ierr)
      call DMRestoreLocalVector(dm, rlocal, ierr); CHKERRQ(ierr)

    end subroutine rhs_fn_nonlinear

    subroutine exact_fn(t, v)

      PetscReal, intent(in) :: t
      Vec, intent(out) :: v
      ! Locals:
      DM :: cdm
      Vec :: cv, vlocal
      PetscReal, pointer :: va(:), coords(:)
      PetscInt :: i1, im, i
      PetscReal :: x

      call DMGetCoordinateDM(dm, cdm, ierr); CHKERRQ(ierr)
      call DMGetCoordinatesLocal(dm, cv, ierr); CHKERRQ(ierr)
      call DMDAVecGetArrayF90(cdm, cv, coords, ierr); CHKERRQ(ierr)
      call DMGetLocalVector(dm, vlocal, ierr); CHKERRQ(ierr)
      call DMDAVecGetArrayF90(dm, vlocal, va, ierr); CHKERRQ(ierr)
      call DMDAGetCorners(dm, i1, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
           im, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
      do i = i1, i1+im-1
         x = coords(i)
         va(i) = fn(x, t)
      end do
      call DMDAVecRestoreArrayF90(dm, vlocal, va, ierr); CHKERRQ(ierr)
      call DMDAVecRestoreArrayF90(cdm, cv, coords, ierr); CHKERRQ(ierr)
      call DMLocalToGlobalBegin(dm, vlocal, INSERT_VALUES, v, ierr); CHKERRQ(ierr)
      call DMLocalToGlobalEnd(dm, vlocal, INSERT_VALUES, v, ierr); CHKERRQ(ierr)

    end subroutine exact_fn

  end subroutine test_heat1d

!------------------------------------------------------------------------

end module timestepping_test
