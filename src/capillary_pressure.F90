!   Copyright 2017 University of Auckland.

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

module capillary_pressure_module
  !! Capillary pressure functions.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use fson

  implicit none
  private

  PetscInt, parameter, public :: max_capillary_pressure_name_length = 12

!------------------------------------------------------------------------

  type, public, abstract :: capillary_pressure_type
     !! Abstract type for capillary pressure objects.
     private
     character(max_capillary_pressure_name_length), &
          public :: name !! Name of capillary pressure function
   contains
     private
     procedure(capillary_pressure_init_routine), public, deferred :: init
     procedure(capillary_pressure_function), public, deferred :: value
  end type capillary_pressure_type

!------------------------------------------------------------------------

  abstract interface

     subroutine capillary_pressure_init_routine(self, json, logfile)
       !! Initializes capillary pressure object from JSON data.
       use logfile_module
       import :: capillary_pressure_type, fson_value
       class(capillary_pressure_type), intent(in out) :: self
       type(fson_value), pointer, intent(in) :: json
       type(logfile_type), intent(in out), optional :: logfile
     end subroutine capillary_pressure_init_routine

     PetscReal function capillary_pressure_function(self, sl)
       !! Capillary pressure function.
       import :: capillary_pressure_type
       class(capillary_pressure_type), intent(in) :: self
       PetscReal, intent(in) :: sl  ! Liquid saturation
     end function capillary_pressure_function

  end interface

!------------------------------------------------------------------------

  type, extends(capillary_pressure_type), public :: &
       capillary_pressure_zero_type
     !! Zero capillary pressure function.
   contains
     procedure, public :: init => capillary_pressure_zero_init
     procedure, public :: value => capillary_pressure_zero_value
  end type capillary_pressure_zero_type

!------------------------------------------------------------------------

!   type, extends(relative_permeability_type), &
!        public :: relative_permeability_linear_type
!      !! Linear relative permeability functions.
!      private
!      PetscReal, public :: liquid_limits(2)
!      PetscReal, public :: vapour_limits(2)
!    contains
!      procedure, public :: init => relative_permeability_linear_init
!      procedure, public :: values => relative_permeability_linear_values
!   end type relative_permeability_linear_type

! !------------------------------------------------------------------------

!   type, extends(relative_permeability_type), &
!        public :: relative_permeability_pickens_type
!      !! Pickens relative permeability curves.
!      private
!      PetscReal, public :: power
!    contains
!      procedure, public :: init => relative_permeability_pickens_init
!      procedure, public :: values => relative_permeability_pickens_values
!   end type relative_permeability_pickens_type

! !------------------------------------------------------------------------

!   type, extends(relative_permeability_type), &
!        public :: relative_permeability_corey_type
!      !! Corey's relative permeability curves.
!      private
!      PetscReal, public :: slr, ssr
!    contains
!      procedure, public :: init => relative_permeability_corey_init
!      procedure, public :: values => relative_permeability_corey_values
!      procedure :: sstar => relative_permeability_corey_sstar
!   end type relative_permeability_corey_type

! !------------------------------------------------------------------------

!   type, extends(relative_permeability_corey_type), &
!        public :: relative_permeability_grant_type
!      !! Grant's relative permeability curves.
!    contains
!      procedure, public :: init => relative_permeability_grant_init
!      procedure, public :: values => relative_permeability_grant_values
!   end type relative_permeability_grant_type

!------------------------------------------------------------------------

  public :: setup_capillary_pressures

contains

!------------------------------------------------------------------------
! Zero function
!------------------------------------------------------------------------

  subroutine capillary_pressure_zero_init(self, json, logfile)
    !! Initialize zero capillary pressure function.

    use fson_mpi_module
    use logfile_module

    class(capillary_pressure_zero_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    type(logfile_type), intent(in out), optional :: logfile

    self%name = "Zero"

  end subroutine capillary_pressure_zero_init

!------------------------------------------------------------------------

  PetscReal function capillary_pressure_zero_value(self, sl) result(cp)
    !! Evaluate zero capillary pressure function.

    class(capillary_pressure_zero_type), intent(in) :: self
    PetscReal, intent(in) :: sl !! Liquid saturation

    cp = 0._dp

  end function capillary_pressure_zero_value

!------------------------------------------------------------------------
! Linear functions
!------------------------------------------------------------------------

!   subroutine relative_permeability_linear_init(self, json, logfile)
!     !! Initialize linear relative permeability function.

!     use fson_mpi_module
!     use logfile_module

!     class(relative_permeability_linear_type), intent(in out) :: self
!     type(fson_value), pointer, intent(in) :: json
!     type(logfile_type), intent(in out), optional :: logfile
!     ! Locals:
!     PetscReal, allocatable :: liquid_limits(:), vapour_limits(:)
!     PetscReal, parameter :: default_liquid_limits(2) = [0._dp, 1._dp]
!     PetscReal, parameter :: default_vapour_limits(2) = [0._dp, 1._dp]

!     self%name = "Linear"

!     call fson_get_mpi(json, "liquid", default_liquid_limits, &
!          liquid_limits, logfile, "rock.relative_permeability.liquid")
!     call fson_get_mpi(json, "vapour", default_vapour_limits, &
!          vapour_limits, logfile, "rock.relative_permeability.vapour")

!     self%liquid_limits = liquid_limits
!     self%vapour_limits = vapour_limits

!     deallocate(liquid_limits, vapour_limits)

!   end subroutine relative_permeability_linear_init

! !------------------------------------------------------------------------

!   function relative_permeability_linear_values(self, sl) result(rp)
!     !! Evaluate linear relative permeability function.

!     class(relative_permeability_linear_type), intent(in) :: self
!     PetscReal, intent(in) :: sl !! Liquid saturation
!     PetscReal, dimension(2) :: rp !! Relative permeabilities
!     ! Locals:
!     PetscReal :: sv

!     sv = 1._dp - sl
!     rp(1) = linear_interpolate(sl, self%liquid_limits)
!     rp(2) = linear_interpolate(sv, self%vapour_limits)

!   contains

!     PetscReal function linear_interpolate(x, x_values) result(y)
!       !! Interpolates linearly between 0 and 1 within specified
!       !! limits.

!       PetscReal, intent(in) :: x, x_values(2)

!       if (x < x_values(1)) then
!          y = 0._dp
!       else if (x > x_values(2)) then
!          y = 1._dp
!       else
!          y = (x - x_values(1)) / (x_values(2) - x_values(1))
!       end if

!     end function linear_interpolate

!   end function relative_permeability_linear_values

! !------------------------------------------------------------------------
! ! Pickens curves
! !------------------------------------------------------------------------

!   subroutine relative_permeability_pickens_init(self, json, logfile)
!     !! Initialize Pickens relative permeability function.

!     use fson_mpi_module
!     use logfile_module

!     class(relative_permeability_pickens_type), intent(in out) :: self
!     type(fson_value), pointer, intent(in) :: json
!     type(logfile_type), intent(in out), optional :: logfile
!     ! Locals:
!     PetscReal, parameter :: default_power = 1._dp

!     self%name = "Pickens"

!     call fson_get_mpi(json, "power", default_power, self%power, &
!          logfile, "rock.relative_permeability.power")

!   end subroutine relative_permeability_pickens_init

! !------------------------------------------------------------------------

!   function relative_permeability_pickens_values(self, sl) result(rp)
!     !! Evaluate Pickens relative permeability function.

!     class(relative_permeability_pickens_type), intent(in) :: self
!     PetscReal, intent(in) :: sl !! Liquid saturation
!     PetscReal, dimension(2) :: rp !! Relative permeabilities

!     rp(1) = sl ** self%power
!     rp(2) = 1._dp

!   end function relative_permeability_pickens_values

! !------------------------------------------------------------------------
! ! Corey's curves
! !------------------------------------------------------------------------

!   subroutine relative_permeability_corey_init(self, json, logfile)
!     !! Initialize Corey's relative permeability function.

!     use fson_mpi_module
!     use logfile_module

!     class(relative_permeability_corey_type), intent(in out) :: self
!     type(fson_value), pointer, intent(in) :: json
!     type(logfile_type), intent(in out), optional :: logfile
!     ! Locals:
!     PetscReal, parameter :: default_slr = 0.3_dp, default_ssr = 0.6_dp

!     self%name = "Corey"

!     call fson_get_mpi(json, "slr", default_slr, self%slr, logfile, &
!          "rock.relative_permeability.slr")
!     call fson_get_mpi(json, "ssr", default_ssr, self%ssr, logfile, &
!          "rock.relative_permeability.ssr")

!   end subroutine relative_permeability_corey_init

! !------------------------------------------------------------------------

!   PetscReal function relative_permeability_corey_sstar(self, sl) &
!        result(sstar)
!     !! Return sstar value used to calculate relative permeabilities.

!     class(relative_permeability_corey_type), intent(in) :: self
!     PetscReal, intent(in) :: sl !! Liquid saturation

!     sstar = (sl - self%slr) / (1._dp - self%slr - self%ssr)

!   end function relative_permeability_corey_sstar

! !------------------------------------------------------------------------

!   function relative_permeability_corey_values(self, sl) result(rp)
!     !! Evaluate Corey's relative permeability function.

!     class(relative_permeability_corey_type), intent(in) :: self
!     PetscReal, intent(in) :: sl !! Liquid saturation
!     PetscReal, dimension(2) :: rp !! Relative permeabilities
!     ! Locals:
!     PetscReal :: sv, sstar, sstar2

!     sv = 1._dp - sl
!     if (sv < self%ssr) then
!        rp = [1._dp, 0._dp]
!     else if (sv > 1._dp - self%slr) then
!        rp = [0._dp, 1._dp]
!     else
!        sstar = self%sstar(sl)
!        sstar2 = sstar * sstar
!        rp(1) = sstar2 * sstar2
!        rp(2) = (1._dp - 2._dp * sstar + sstar2) * (1._dp - sstar2)
!     end if

!   end function relative_permeability_corey_values

! !------------------------------------------------------------------------
! ! Grant's curves
! !------------------------------------------------------------------------

!   subroutine relative_permeability_grant_init(self, json, logfile)
!     !! Initialize Grant's relative permeability function.

!     use fson_mpi_module
!     use logfile_module

!     class(relative_permeability_grant_type), intent(in out) :: self
!     type(fson_value), pointer, intent(in) :: json
!     type(logfile_type), intent(in out), optional :: logfile
!     ! Locals:
!     PetscReal, parameter :: default_slr = 0.3_dp, default_ssr = 0.6_dp

!     self%name = "Grant"

!     call fson_get_mpi(json, "slr", default_slr, self%slr, logfile, &
!          "rock.relative_permeability.slr")
!     call fson_get_mpi(json, "ssr", default_ssr, self%ssr, logfile, &
!          "rock.relative_permeability.ssr")

!   end subroutine relative_permeability_grant_init

! !------------------------------------------------------------------------

!   function relative_permeability_grant_values(self, sl) result(rp)
!     !! Evaluate Grant's relative permeability function.

!     class(relative_permeability_grant_type), intent(in) :: self
!     PetscReal, intent(in) :: sl !! Liquid saturation
!     PetscReal, dimension(2) :: rp !! Relative permeabilities
!     ! Locals:
!     PetscReal :: sv, sstar, sstar2

!     sv = 1._dp - sl
!     if (sv < self%ssr) then
!        rp = [1._dp, 0._dp]
!     else if (sv > 1._dp - self%slr) then
!        rp = [0._dp, 1._dp]
!     else
!        sstar = self%sstar(sl)
!        sstar2 = sstar * sstar
!        rp(1) = sstar2 * sstar2
!        rp(2) = 1._dp - rp(1)
!     end if

!   end function relative_permeability_grant_values

!------------------------------------------------------------------------
! Setup procedures
!------------------------------------------------------------------------

  subroutine setup_capillary_pressure(json, cp, logfile)
    !! Sets up single capillary pressure object from JSON object
    !! for rock data in input.

    use fson_mpi_module
    use utils_module, only: str_to_lower
    use logfile_module

    type(fson_value), pointer, intent(in) :: json
    class(capillary_pressure_type), allocatable, intent(out) :: cp
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    character(max_capillary_pressure_name_length), parameter :: &
         default_capillary_type = "zero"
    character(max_capillary_pressure_name_length) :: capillary_type

    call fson_get_mpi(json, "type", default_capillary_type, &
         capillary_type, logfile, "rock.capillary_pressure.type")

    select case (str_to_lower(capillary_type))
    case ("zero")
       allocate(capillary_pressure_zero_type :: cp)
    case default
       allocate(capillary_pressure_zero_type :: cp)
    end select
    call cp%init(json, logfile)

  end subroutine setup_capillary_pressure

!------------------------------------------------------------------------

  subroutine setup_capillary_pressures(json, cp, logfile)
    !! Sets up capillary pressure objects from rock data in JSON input.
    !! Eventually cp will be an array of objects.

    use fson_mpi_module
    use logfile_module

    type(fson_value), pointer, intent(in) :: json
    class(capillary_pressure_type), allocatable, intent(out) :: cp
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    type(fson_value), pointer :: capillary
    PetscBool :: default_present

    if (fson_has_mpi(json, "rock.capillary_pressure")) then
       call fson_get_mpi(json, "rock.capillary_pressure", capillary)
       default_present = PETSC_TRUE
    else
       capillary => fson_parse(str = "{}")
       default_present = PETSC_FALSE
    end if

    call setup_capillary_pressure(capillary, cp, logfile)

    if (.not. (default_present)) then
       call fson_destroy(capillary)
    end if

  end subroutine setup_capillary_pressures

!------------------------------------------------------------------------

end module capillary_pressure_module
