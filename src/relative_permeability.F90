module relative_permeability_module
  !! Relative permeability functions.

  use kinds_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>

  PetscInt, parameter, public :: max_relative_permeability_name_length = 12

!------------------------------------------------------------------------

  type, public, abstract :: relative_permeability_type
     !! Abstract type for relative permeability objects.
     private
     character(max_relative_permeability_name_length), &
          public :: name !! Name of relative permeability curves
   contains
     private
     procedure(relative_permeability_function), public, deferred :: values
  end type relative_permeability_type

!------------------------------------------------------------------------

  abstract interface

     function relative_permeability_function(self, sl) result(rp)
       !! Relative permeability function, returning array of values for both phases.
       import :: relative_permeability_type
       class(relative_permeability_type), intent(in) :: self
       PetscReal, intent(in) :: sl  ! Liquid saturation
       PetscReal, dimension(2) :: rp ! Relative permeabilities for liquid and vapour
     end function relative_permeability_function

  end interface

!------------------------------------------------------------------------

  type, extends(relative_permeability_type), &
       public :: relative_permeability_linear_type
     !! Linear relative permeability functions.
     private
     PetscReal, public :: liquid_limits(2)
     PetscReal, public :: vapour_limits(2)
   contains
     procedure, public :: init => relative_permeability_linear_init
     procedure, public :: values => relative_permeability_linear_values
  end type relative_permeability_linear_type

!------------------------------------------------------------------------

  type, extends(relative_permeability_type), &
       public :: relative_permeability_corey_type
     !! Corey's relative permeability curves.
     private
     PetscReal, public :: slr, ssr
   contains
     procedure, public :: init => relative_permeability_corey_init
     procedure, public :: values => relative_permeability_corey_values
     procedure :: sstar => relative_permeability_corey_sstar
  end type relative_permeability_corey_type

!------------------------------------------------------------------------

  type, extends(relative_permeability_corey_type), &
       public :: relative_permeability_grant_type
     !! Grant's relative permeability curves.
   contains
     procedure, public :: init => relative_permeability_grant_init
     procedure, public :: values => relative_permeability_grant_values
  end type relative_permeability_grant_type

!------------------------------------------------------------------------

contains

!------------------------------------------------------------------------
! Linear functions
!------------------------------------------------------------------------

  subroutine relative_permeability_linear_init(self, liquid_limits, &
       vapour_limits)
    !! Initialize linear relative permeability function.

    class(relative_permeability_linear_type), intent(in out) :: self
    PetscReal, intent(in) :: liquid_limits(2), vapour_limits(2)

    self%name = "Linear"
    self%liquid_limits = liquid_limits
    self%vapour_limits = vapour_limits

  end subroutine relative_permeability_linear_init

!------------------------------------------------------------------------

  function relative_permeability_linear_values(self, sl) result(rp)
    !! Evaluate linear relative permeability function.

    class(relative_permeability_linear_type), intent(in) :: self
    PetscReal, intent(in) :: sl !! Liquid saturation
    PetscReal, dimension(2) :: rp !! Relative permeabilities
    ! Locals:
    PetscReal :: sv

    sv = 1._dp - sl
    rp(1) = linear_interpolate(sl, self%liquid_limits)
    rp(2) = linear_interpolate(sv, self%vapour_limits)

  contains

    PetscReal function linear_interpolate(x, x_values) result(y)
      !! Interpolates linearly between 0 and 1 within specified
      !! limits.

      PetscReal, intent(in) :: x, x_values(2)

      if (x < x_values(1)) then
         y = 0._dp
      else if (x > x_values(2)) then
         y = 1._dp
      else
         y = (x - x_values(1)) / (x_values(2) - x_values(1))
      end if

    end function linear_interpolate

  end function relative_permeability_linear_values

!------------------------------------------------------------------------
! Corey's curves
!------------------------------------------------------------------------

  subroutine relative_permeability_corey_init(self, slr, ssr)
    !! Initialize Corey's relative permeability function.

    class(relative_permeability_corey_type), intent(in out) :: self
    PetscReal, intent(in) :: slr, ssr

    self%name = "Corey"
    self%slr = slr
    self%ssr = ssr

  end subroutine relative_permeability_corey_init

!------------------------------------------------------------------------

  PetscReal function relative_permeability_corey_sstar(self, sl) &
       result(sstar)
    !! Return sstar value used to calculate relative permeabilities.

    class(relative_permeability_corey_type), intent(in) :: self
    PetscReal, intent(in) :: sl !! Liquid saturation

    sstar = (sl - self%slr) / (1._dp - self%slr - self%ssr)

  end function relative_permeability_corey_sstar

!------------------------------------------------------------------------

  function relative_permeability_corey_values(self, sl) result(rp)
    !! Evaluate Corey's relative permeability function.

    class(relative_permeability_corey_type), intent(in) :: self
    PetscReal, intent(in) :: sl !! Liquid saturation
    PetscReal, dimension(2) :: rp !! Relative permeabilities
    ! Locals:
    PetscReal :: sv, sstar, sstar2

    sv = 1._dp - sl
    if (sv < self%ssr) then
       rp = [1._dp, 0._dp]
    else if (sv > 1._dp - self%slr) then
       rp = [0._dp, 1._dp]
    else
       sstar = self%sstar(sl)
       sstar2 = sstar * sstar
       rp(1) = sstar2 * sstar2
       rp(2) = (1._dp - 2._dp * sstar + sstar2) * (1._dp - sstar2)
    end if

  end function relative_permeability_corey_values

!------------------------------------------------------------------------
! Grant's curves
!------------------------------------------------------------------------

  subroutine relative_permeability_grant_init(self, slr, ssr)
    !! Initialize Grant's relative permeability function.

    class(relative_permeability_grant_type), intent(in out) :: self
    PetscReal, intent(in) :: slr, ssr

    self%name = "Grant"
    self%slr = slr
    self%ssr = ssr

  end subroutine relative_permeability_grant_init

!------------------------------------------------------------------------

  function relative_permeability_grant_values(self, sl) result(rp)
    !! Evaluate Grant's relative permeability function.

    class(relative_permeability_grant_type), intent(in) :: self
    PetscReal, intent(in) :: sl !! Liquid saturation
    PetscReal, dimension(2) :: rp !! Relative permeabilities
    ! Locals:
    PetscReal :: sv, sstar, sstar2

    sv = 1._dp - sl
    if (sv < self%ssr) then
       rp = [1._dp, 0._dp]
    else if (sv > 1._dp - self%slr) then
       rp = [0._dp, 1._dp]
    else
       sstar = self%sstar(sl)
       sstar2 = sstar * sstar
       rp(1) = sstar2 * sstar2
       rp(2) = 1._dp - rp(1)
    end if

  end function relative_permeability_grant_values

!------------------------------------------------------------------------

end module relative_permeability_module
