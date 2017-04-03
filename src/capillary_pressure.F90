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
  use interpolation_module
  use fson

  implicit none
  private

  PetscInt, parameter, public :: max_capillary_pressure_name_length = 16

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
     procedure, public :: destroy => capillary_pressure_destroy
  end type capillary_pressure_type

!------------------------------------------------------------------------

  type, public, extends(capillary_pressure_type) :: &
       capillary_pressure_zero_type
     !! Zero capillary pressure function.
   contains
     procedure, public :: init => capillary_pressure_zero_init
     procedure, public :: value => capillary_pressure_zero_value
  end type capillary_pressure_zero_type

!------------------------------------------------------------------------

  type, public, extends(capillary_pressure_type) :: &
       capillary_pressure_linear_type
     !! Linear capillary pressure function.
     private
     PetscReal, public :: saturation_limits(2)
     PetscReal, public :: pressure
   contains
     procedure, public :: init => capillary_pressure_linear_init
     procedure, public :: value => capillary_pressure_linear_value
  end type capillary_pressure_linear_type

!------------------------------------------------------------------------

  type, public, extends(capillary_pressure_type) :: &
       capillary_pressure_van_genuchten_type
     !! Van Genuchten capillary pressure function.
     private
     PetscReal, public :: P0, Pmax
     PetscReal, public :: lambda
     PetscReal, public :: slr, sls
     PetscBool, public :: apply_Pmax
   contains
     procedure, public :: init => capillary_pressure_van_genuchten_init
     procedure, public :: value => capillary_pressure_van_genuchten_value
  end type capillary_pressure_van_genuchten_type

!------------------------------------------------------------------------

  type, public, extends(capillary_pressure_type) :: &
       capillary_pressure_table_type
     !! Table capillary pressure function.
     private
     type(interpolation_table_type), public :: pressure
   contains
     procedure, public :: init => capillary_pressure_table_init
     procedure, public :: value => capillary_pressure_table_value
     procedure, public :: destroy => capillary_pressure_table_destroy
  end type capillary_pressure_table_type

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

     PetscReal function capillary_pressure_function(self, sl, t)
       !! Capillary pressure function.
       import :: capillary_pressure_type
       class(capillary_pressure_type), intent(in out) :: self
       PetscReal, intent(in) :: sl  !! Liquid saturation
       PetscReal, intent(in) :: t   !! Temperature
     end function capillary_pressure_function

  end interface

  public :: setup_capillary_pressures

contains

!------------------------------------------------------------------------
!  capillary_pressure_type
!------------------------------------------------------------------------

  subroutine capillary_pressure_destroy(self)
    !! Destroys capillary pressure object. Dummy method, to be
    !! overridden (as needed) by derived types.

    class(capillary_pressure_type), intent(in out) :: self

    continue

  end subroutine capillary_pressure_destroy

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

    self%name = "zero"

  end subroutine capillary_pressure_zero_init

!------------------------------------------------------------------------

  PetscReal function capillary_pressure_zero_value(self, sl, t) result(cp)
    !! Evaluate zero capillary pressure function.

    class(capillary_pressure_zero_type), intent(in out) :: self
    PetscReal, intent(in) :: sl !! Liquid saturation
    PetscReal, intent(in) :: t  !! Temperature

    cp = 0._dp

  end function capillary_pressure_zero_value

!------------------------------------------------------------------------
! Linear function
!------------------------------------------------------------------------

  subroutine capillary_pressure_linear_init(self, json, logfile)
    !! Initialize linear capillary pressure function.

    use fson_mpi_module
    use logfile_module

    class(capillary_pressure_linear_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscReal, allocatable :: saturation_limits(:)
    PetscReal, parameter :: default_saturation_limits(2) = [0._dp, 1._dp]
    PetscReal, parameter :: default_pressure = 0.125e5_dp

    self%name = "linear"

    call fson_get_mpi(json, "saturation_limits", default_saturation_limits, &
         saturation_limits, logfile, "rock.capillary_pressure.saturation_limits")
    call fson_get_mpi(json, "pressure", default_pressure, &
         self%pressure, logfile, "rock.capillary_pressure.pressure")
    self%pressure = abs(self%pressure)

    self%saturation_limits = saturation_limits

    deallocate(saturation_limits)

  end subroutine capillary_pressure_linear_init

!------------------------------------------------------------------------

  PetscReal function capillary_pressure_linear_value(self, sl, t) result(cp)
    !! Evaluate linear capillary pressure function.

    class(capillary_pressure_linear_type), intent(in out) :: self
    PetscReal, intent(in) :: sl !! Liquid saturation
    PetscReal, intent(in) :: t  !! Temperature

    cp = ramp_interpolate(sl, self%saturation_limits, &
         [-self%pressure, 0._dp])

  end function capillary_pressure_linear_value

!------------------------------------------------------------------------
! Van Genuchten function
!------------------------------------------------------------------------

  subroutine capillary_pressure_van_genuchten_init(self, json, logfile)
    !! Initialize van Genuchten capillary pressure function.

    use fson_mpi_module
    use logfile_module

    class(capillary_pressure_van_genuchten_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscReal, parameter :: default_P0 = 0.125e5_dp
    PetscReal, parameter :: default_lambda = 0.45_dp
    PetscReal, parameter :: default_slr = 1.e-3_dp
    PetscReal, parameter :: default_sls = 1._dp

    self%name = "van Genuchten"

    call fson_get_mpi(json, "P0", default_P0, &
         self%P0, logfile, "rock.capillary_pressure.P0")
    self%P0 = abs(self%P0)
    call fson_get_mpi(json, "lambda", default_lambda, &
         self%lambda, logfile, "rock.capillary_pressure.lambda")
    call fson_get_mpi(json, "slr", default_slr, &
         self%slr, logfile, "rock.capillary_pressure.slr")
    call fson_get_mpi(json, "sls", default_sls, &
         self%sls, logfile, "rock.capillary_pressure.sls")

    if (fson_has_mpi(json, "Pmax")) then
       call fson_get_mpi(json, "Pmax", val = self%Pmax)
       self%Pmax = abs(self%Pmax)
       self%apply_Pmax = PETSC_TRUE
    else
       self%apply_Pmax = PETSC_FALSE
    end if
    
  end subroutine capillary_pressure_van_genuchten_init

!------------------------------------------------------------------------

  PetscReal function capillary_pressure_van_genuchten_value(self, sl, t) &
       result(cp)
    !! Evaluate van Genuchten capillary pressure function.

    class(capillary_pressure_van_genuchten_type), intent(in out) :: self
    PetscReal, intent(in) :: sl !! Liquid saturation
    PetscReal, intent(in) :: t  !! Temperature

    associate(sstar => (sl - self%slr) / (self%sls - self%slr))
      if (sstar < 0._dp) then
         cp = -self%P0
      else if (sstar < 1._dp) then
         cp = -self%P0 * (sstar ** (-1._dp / self%lambda) - 1._dp) ** &
              (1._dp - self%lambda)
      else
         cp = 0._dp
      end if
      cp = min(0._dp, cp)
      if (self%apply_Pmax) then
         cp = max(-self%Pmax, cp)
      end if
    end associate

  end function capillary_pressure_van_genuchten_value

!------------------------------------------------------------------------
! Table function
!------------------------------------------------------------------------

  subroutine capillary_pressure_table_init(self, json, logfile)
    !! Initialize table capillary pressure function.

    use fson_mpi_module
    use logfile_module

    class(capillary_pressure_table_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscReal, allocatable :: pressure_array(:,:)
    PetscReal, parameter :: default_pressure_array(2, 2) = reshape( &
         [0._dp, 1._dp, 0._dp, 0._dp], [2, 2])

    self%name = "table"

    call fson_get_mpi(json, "pressure", default_pressure_array, &
         pressure_array, logfile, "rock.capillary_pressure.pressure")
    call self%pressure%init(pressure_array)
    deallocate(pressure_array)

  end subroutine capillary_pressure_table_init

!------------------------------------------------------------------------

  PetscReal function capillary_pressure_table_value(self, sl, t) result(cp)
    !! Evaluate table capillary pressure function.

    class(capillary_pressure_table_type), intent(in out) :: self
    PetscReal, intent(in) :: sl !! Liquid saturation
    PetscReal, intent(in) :: t  !! Temperature

    cp = self%pressure%interpolate(sl)

  end function capillary_pressure_table_value

!------------------------------------------------------------------------

  subroutine capillary_pressure_table_destroy(self)
    !! Destroy table capillary pressure function.

    class(capillary_pressure_table_type), intent(in out) :: self

    call self%pressure%destroy()

  end subroutine capillary_pressure_table_destroy

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
    case ("linear")
       allocate(capillary_pressure_linear_type :: cp)
    case ("van_genuchten", "van genuchten")
       allocate(capillary_pressure_van_genuchten_type :: cp)
    case ("table")
       allocate(capillary_pressure_table_type :: cp)
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
    use fson_value_m, only: TYPE_OBJECT, TYPE_NULL
    use logfile_module

    type(fson_value), pointer, intent(in) :: json
    class(capillary_pressure_type), allocatable, intent(out) :: cp
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    type(fson_value), pointer :: capillary
    PetscInt :: cp_type
    PetscBool :: default_present

    if (fson_has_mpi(json, "rock.capillary_pressure")) then
       cp_type = fson_type_mpi(json, "rock.capillary_pressure")
       if (cp_type == TYPE_OBJECT) then
          call fson_get_mpi(json, "rock.capillary_pressure", capillary)
          default_present = PETSC_TRUE
       else if (cp_type == TYPE_NULL) then
          capillary => fson_parse(str = "{}")
          default_present = PETSC_FALSE
       end if
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
