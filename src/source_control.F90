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

module source_control_module
  !! Module for source controls- for controlling source parameters (e.g. flow rate, enthalpy) over time.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use eos_module
  use thermodynamics_module
  use list_module
  use source_module
  use interpolation_module

  implicit none
  private

  PetscReal, parameter, public :: default_source_control_separator_pressure = 0.55e6_dp
  PetscInt, parameter, public :: max_limiter_type_length = 5
  character(max_limiter_type_length), parameter, public :: &
       default_source_control_limiter_type_str = "total"
  PetscInt, parameter, public :: SRC_CONTROL_LIMITER_TYPE_TOTAL = 1, &
       SRC_CONTROL_LIMITER_TYPE_WATER = 2, SRC_CONTROL_LIMITER_TYPE_STEAM = 3
  PetscReal, parameter, public :: default_source_control_limiter_limit = 1._dp
  PetscReal, parameter, public :: default_deliverability_productivity = 1.e-11_dp
  PetscReal, parameter, public :: default_deliverability_reference_pressure = 1.e5_dp
  PetscReal, parameter, public :: default_recharge_coefficient = 1.e-2_dp
  PetscInt, parameter, public :: SRC_DIRECTION_PRODUCTION = 1, &
       SRC_DIRECTION_INJECTION = 2, SRC_DIRECTION_BOTH = 3
  PetscInt, parameter, public :: default_source_direction = SRC_DIRECTION_BOTH
  PetscInt, parameter, public :: SRC_PRESSURE_TABLE_COORD_TIME = 1, &
       SRC_PRESSURE_TABLE_COORD_ENTHALPY = 2
  PetscInt, parameter, public :: default_source_pressure_table_coordinate = &
       SRC_PRESSURE_TABLE_COORD_TIME

  type, public, abstract :: source_control_type
     !! Abstract type for source control, controlling source
     !! parameters over time, for one or more sources.
     private
   contains
     procedure(source_control_destroy_procedure), public, deferred :: destroy
     procedure(source_control_update_procedure), public, deferred :: update
  end type source_control_type

  type, abstract, extends(source_control_type) :: source_control_table_type
     !! Controls a source parameter (e.g. rate or enthalpy) via a
     !! table of values vs. time.
     private
     PetscInt, allocatable, public :: source_indices(:)
     type(interpolation_table_type), public :: table !! Table of values vs. time
   contains
     private
     procedure, public :: init => source_control_table_init
     procedure, public :: destroy => source_control_table_destroy
  end type source_control_table_type

  type, public, extends(source_control_table_type) :: source_control_rate_table_type
     !! Controls source rate via a table of values vs. time.
   contains
     private
     procedure, public :: update => source_control_rate_table_update
  end type source_control_rate_table_type

  type, public, extends(source_control_table_type) :: source_control_enthalpy_table_type
     !! Controls source injection enthalpy via a table of values vs. time.
   contains
     private
     procedure, public :: update => source_control_enthalpy_table_update
  end type source_control_enthalpy_table_type

  type, public, extends(source_control_table_type) :: source_control_rate_factor_type
     !! Multiplies source rate by a factor from a table of values vs. time.
   contains
     private
     procedure, public :: update => source_control_rate_factor_update
  end type source_control_rate_factor_type

  type, public, abstract, extends(source_control_type) :: source_control_pressure_reference_type
     !! Controls a source by comparing fluid pressure with a reference
     !! pressure (e.g. deliverability or recharge).
     private
     PetscInt, public :: local_source_index !! Index of source
     type(interpolation_table_type), public :: reference_pressure !! Reference pressure vs. time
   contains
     procedure, public :: set_reference_pressure_initial => &
          source_control_set_reference_pressure_initial
  end type source_control_pressure_reference_type

  type, public, extends(source_control_pressure_reference_type) :: source_control_deliverability_type
     !! Controls a source on deliverability.
     private
     PetscInt, public :: pressure_table_coordinate !! Coordinate variable of pressure table
     type(interpolation_table_type), public :: productivity !! Productivity index vs. time
     PetscReal, public :: threshold !! Pressure threshold below which deliverability is switched on (< 0 for always on)
     PetscReal, public :: threshold_productivity !! Productivity index computed from flow rate and used when pressure drops below threshold
   contains
     procedure, public :: init => source_control_deliverability_init
     procedure, public :: destroy => source_control_deliverability_destroy
     procedure, public :: update => source_control_deliverability_update
     procedure, public :: calculate_PI_from_rate => &
          source_control_deliverability_calculate_PI_from_rate
  end type source_control_deliverability_type

  type, public, extends(source_control_pressure_reference_type) :: source_control_recharge_type
     !! Controls a source simulating recharge through a model boundary.
     private
     type(interpolation_table_type), public :: coefficient !! Recharge coefficient vs. time
   contains
     procedure, public :: init => source_control_recharge_init
     procedure, public :: destroy => source_control_recharge_destroy
     procedure, public :: update => source_control_recharge_update
  end type source_control_recharge_type

  type, public, extends(source_control_type) :: source_control_separator_type
     !! Takes flow from a single source and outputs water and steam flow.
     private
     PetscInt :: local_source_index
     class(thermodynamics_type), pointer :: thermo
     PetscReal, public :: separator_pressure !! Separator pressure
     PetscReal :: water_enthalpy !! Enthalpy of water at separator pressure
     PetscReal :: steam_enthalpy !! Enthalpy of steam at separator pressure
     PetscReal, public :: steam_fraction  !! Output steam fraction
     PetscReal, public :: water_flow_rate !! Output separated water flow rate
     PetscReal, public :: steam_flow_rate !! Output separated steam flow rate
   contains
     procedure, public :: init => source_control_separator_init
     procedure, public :: destroy => source_control_separator_destroy
     procedure, public :: update => source_control_separator_update
     procedure :: update_steam_fraction => source_control_separator_update_steam_fraction
  end type source_control_separator_type

  type, public, extends(source_control_type) :: source_control_limiter_type
     !! Limits total flow, or separated steam or water flow (output
     !! from a separator) through a source.
     private
     PetscInt, allocatable, public :: local_source_indices(:) !! Local source indices
     PetscInt, public :: input_local_source_index !! Local index of source to be used as input
     class(source_control_type), pointer :: input_source_control !! Source control to be used as input
     PetscInt, public :: type !! Type of limiter (water, steam or total)
     PetscReal, public :: limit !! Flow limit
   contains
     private
     procedure :: rate_scale => source_control_limiter_rate_scale
     generic, public :: init => init_source, init_control
     procedure :: init_source => source_control_limiter_init_source
     procedure :: init_control => source_control_limiter_init_control
     procedure, public :: destroy => source_control_limiter_destroy
     procedure, public :: update => source_control_limiter_update
     procedure :: get_rate => source_control_limiter_get_rate
  end type source_control_limiter_type

  type, public, extends(source_control_type) :: source_control_direction_type
     !! Allows source to flow only in a specified direction
     !! (production or injection).
     private
     PetscInt, allocatable, public :: local_source_indices(:) !! Local source indices
     PetscInt, public :: direction !! Flow direction
   contains
     private
     procedure, public :: init => source_control_direction_init
     procedure, public :: destroy => source_control_direction_destroy
     procedure, public :: update => source_control_direction_update
  end type source_control_direction_type

  abstract interface

     subroutine source_control_init_procedure(self)
       !! Initialises a source control object.
       import :: source_control_type
       class(source_control_type), intent(in out) :: self
     end subroutine source_control_init_procedure

     subroutine source_control_destroy_procedure(self)
       !! Destroys a source control object.
       import :: source_control_type
       class(source_control_type), intent(in out) :: self
     end subroutine source_control_destroy_procedure

     subroutine source_control_update_procedure(self, t, interval, &
          source_data, source_section, source_range_start, &
          local_fluid_data, local_fluid_section, eos)
       use petscis
       !! Updates sources at the specified time.
       import :: source_control_type, eos_type
       class(source_control_type), intent(in out) :: self
       PetscReal, intent(in) :: t, interval(2)
       PetscReal, pointer, contiguous, intent(in) :: source_data(:)
       PetscSection, intent(in) :: source_section
       PetscInt, intent(in) :: source_range_start
       PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
       PetscSection, intent(in) :: local_fluid_section
       class(eos_type), intent(in) :: eos
     end subroutine source_control_update_procedure

  end interface

contains

!------------------------------------------------------------------------
! Table source control:
!------------------------------------------------------------------------

  subroutine source_control_table_init(self, data, interpolation_type, &
       averaging_type, source_indices, err)
    !! Initialises source_control_table object.

    class(source_control_table_type), intent(in out) :: self
    PetscReal, intent(in) :: data(:,:)
    PetscInt, intent(in) :: interpolation_type
    PetscInt, intent(in) :: averaging_type
    PetscInt, intent(in) :: source_indices(:)
    PetscErrorCode, intent(out) :: err

    call self%table%init(data, interpolation_type, averaging_type, err)
    if (err == 0) then
       self%source_indices = source_indices
    end if

  end subroutine source_control_table_init

!------------------------------------------------------------------------

  subroutine source_control_table_destroy(self)
    !! Destroys source_control_table_type object.

    class(source_control_table_type), intent(in out) :: self

    call self%table%destroy()
    if (allocated(self%source_indices)) deallocate(self%source_indices)

  end subroutine source_control_table_destroy

!------------------------------------------------------------------------
! Rate table source control:
!------------------------------------------------------------------------

  subroutine source_control_rate_table_update(self, t, interval, &
       source_data, source_section, source_range_start, &
       local_fluid_data, local_fluid_section, eos)
    !! Update flow rate for source_control_rate_table_type.

    use dm_utils_module, only: global_section_offset

    class(source_control_rate_table_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    ! Locals:
    PetscReal :: rate
    type(source_type) :: source
    PetscInt :: i, s, source_offset

    call source%init(eos)
    rate = self%table%average(interval, 1)

    do i = 1, size(self%source_indices)
       s = self%source_indices(i)
       source_offset = global_section_offset(source_section, s, &
            source_range_start)
       call source%assign(source_data, source_offset)
       source%rate = rate
    end do

    call source%destroy()

  end subroutine source_control_rate_table_update

!------------------------------------------------------------------------
! Enthalpy table source control:
!------------------------------------------------------------------------

  subroutine source_control_enthalpy_table_update(self, t, interval, &
       source_data, source_section, source_range_start, &
       local_fluid_data, local_fluid_section, eos)
    !! Update injection enthalpy for source_control_enthalpy_table_type.

    use dm_utils_module, only: global_section_offset

    class(source_control_enthalpy_table_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    ! Locals:
    PetscReal :: enthalpy
    type(source_type) :: source
    PetscInt :: i, s, source_offset

    call source%init(eos)
    enthalpy = self%table%average(interval, 1)

    do i = 1, size(self%source_indices)
       s = self%source_indices(i)
       source_offset = global_section_offset(source_section, s, &
            source_range_start)
       call source%assign(source_data, source_offset)
       source%injection_enthalpy = enthalpy
    end do

    call source%destroy()

  end subroutine source_control_enthalpy_table_update

!------------------------------------------------------------------------
! Rate factor source control:
!------------------------------------------------------------------------

  subroutine source_control_rate_factor_update(self, t, interval, &
       source_data, source_section, source_range_start, &
       local_fluid_data, local_fluid_section, eos)
    !! Update flow rate for source_control_rate_factor_type.

    use dm_utils_module, only: global_section_offset

    class(source_control_rate_factor_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    ! Locals:
    PetscReal :: factor
    type(source_type) :: source
    PetscInt :: i, s, source_offset

    call source%init(eos)
    factor = self%table%average(interval, 1)

    do i = 1, size(self%source_indices)
       s = self%source_indices(i)
       source_offset = global_section_offset(source_section, s, &
            source_range_start)
       call source%assign(source_data, source_offset)
       source%rate = source%rate * factor
    end do

    call source%destroy()

  end subroutine source_control_rate_factor_update

!------------------------------------------------------------------------
! Pressure reference source control:
!------------------------------------------------------------------------

  subroutine source_control_set_reference_pressure_initial(self, &
       source_data, source_section, source_range_start, &
       global_fluid_data, global_fluid_section, fluid_range_start, eos)
    !! Sets reference pressure for pressure reference control to be
    !! the initial fluid pressure in the source cell.

    use dm_utils_module, only: global_section_offset

    class(source_control_pressure_reference_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: global_fluid_data(:)
    PetscSection, intent(in) :: global_fluid_section
    PetscInt, intent(in) :: fluid_range_start
    class(eos_type), intent(in) :: eos
    ! Locals:
    type(source_type) :: source
    PetscInt :: c, fluid_offset, source_offset

    call source%init(eos)
    source_offset = global_section_offset(source_section, &
         self%local_source_index, source_range_start)
    call source%assign(source_data, source_offset)

    c = nint(source%local_cell_index)
    fluid_offset = global_section_offset(global_fluid_section, c, &
         fluid_range_start)
    call source%fluid%assign(global_fluid_data, fluid_offset)

    self%reference_pressure%val(1, 1) = source%fluid%pressure

    call source%destroy()

  end subroutine source_control_set_reference_pressure_initial

!------------------------------------------------------------------------
! Deliverability source control:
!------------------------------------------------------------------------

  subroutine source_control_deliverability_init(self, productivity_data, &
       interpolation_type, averaging_type, reference_pressure_data, &
       pressure_table_coordinate, threshold, local_source_index, err)
    !! Initialises source_control_deliverability object. Error flag err
    !! returns 1 if there are problems with the productivity index
    !! array, or 2 if there are problems with the reference pressure
    !! array.

    class(source_control_deliverability_type), intent(in out) :: self
    PetscReal, intent(in) :: productivity_data(:,:)
    PetscInt, intent(in) :: interpolation_type, averaging_type
    PetscReal, intent(in) :: reference_pressure_data(:,:)
    PetscInt, intent(in) :: pressure_table_coordinate
    PetscReal, intent(in) :: threshold
    PetscInt, intent(in) :: local_source_index
    PetscErrorCode, intent(out) :: err

    call self%productivity%init(productivity_data, &
         interpolation_type, averaging_type, err)
    if (err == 0) then
       call self%reference_pressure%init(reference_pressure_data, &
            interpolation_type, averaging_type, err)
       if (err == 0) then
          self%pressure_table_coordinate = pressure_table_coordinate
          self%threshold = threshold
          self%local_source_index = local_source_index
       else
          err = 2
       end if
    end if

  end subroutine source_control_deliverability_init

!------------------------------------------------------------------------

  subroutine source_control_deliverability_destroy(self)
    !! Destroys source_control_deliverability object.

    class(source_control_deliverability_type), intent(in out) :: self

    call self%productivity%destroy()
    call self%reference_pressure%destroy()

  end subroutine source_control_deliverability_destroy

!------------------------------------------------------------------------

  subroutine source_control_deliverability_update(self, t, interval, &
       source_data, source_section, source_range_start, &
       local_fluid_data, local_fluid_section, eos)
    !! Update flow rate for source_control_deliverability_type.

    use dm_utils_module, only: global_section_offset

    class(source_control_deliverability_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    ! Locals:
    type(source_type) :: source
    PetscInt :: source_offset
    PetscReal :: productivity

    call source%init(eos)

    source_offset = global_section_offset(source_section, &
         self%local_source_index, source_range_start)
    call source%assign(source_data, source_offset)
    call source%assign_fluid(local_fluid_data, local_fluid_section)

    if (self%threshold <= 0._dp) then
       productivity = self%productivity%average(interval, 1)
       call update_flow_rate(source, productivity)
    else
       if (source%fluid%pressure < self%threshold) then
          call update_flow_rate(source, self%threshold_productivity)
       else
          call self%calculate_PI_from_rate(t, source%rate, &
               source_data, source_section, source_range_start, &
               local_fluid_data, local_fluid_section, -1, &
               eos, self%threshold_productivity)
       end if
    end if

    call source%destroy()

  contains

    subroutine update_flow_rate(source, productivity)
      !! Updates flow rate using deliverability relation and the
      !! specified productivity index.
      type(source_type), intent(in out)  :: source
      PetscReal, intent(in) :: productivity
      ! Locals:
      PetscInt :: p, phases
      PetscReal :: h, reference_pressure, pressure_difference
      PetscReal, allocatable :: phase_mobilities(:), phase_flow_fractions(:)

      allocate(phase_mobilities(source%fluid%num_phases))
      allocate(phase_flow_fractions(source%fluid%num_phases))
      phase_mobilities = source%fluid%phase_mobilities()
      phase_flow_fractions = phase_mobilities / sum(phase_mobilities)

      select case (self%pressure_table_coordinate)
      case (SRC_PRESSURE_TABLE_COORD_TIME)
         reference_pressure = self%reference_pressure%average(interval, 1)
      case (SRC_PRESSURE_TABLE_COORD_ENTHALPY)
         h = source%fluid%specific_enthalpy(phase_flow_fractions)
         reference_pressure = self%reference_pressure%interpolate(h, 1)
      end select

      pressure_difference = source%fluid%pressure - reference_pressure
      source%rate = 0._dp

      phases = nint(source%fluid%phase_composition)
      do p = 1, source%fluid%num_phases
         if (btest(phases, p - 1)) then
            source%rate = source%rate - productivity * &
                 phase_mobilities(p) * pressure_difference
         end if
      end do

      deallocate(phase_mobilities, phase_flow_fractions)

    end subroutine update_flow_rate

  end subroutine source_control_deliverability_update

!------------------------------------------------------------------------

  subroutine source_control_deliverability_calculate_PI_from_rate(&
       self, time, rate, source_data, source_section, source_range_start, &
       fluid_data, fluid_section, fluid_range_start, eos, &
       productivity)
    !! Calculates productivity index for deliverability control, from
    !! specified initial flow rate.

    use dm_utils_module, only: global_section_offset, section_offset

    class(source_control_deliverability_type), intent(in out) :: self
    PetscReal, intent(in) :: time
    PetscReal, intent(in) :: rate
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start !! Specify -1 for local data rather than global
    class(eos_type), intent(in) :: eos
    PetscReal, intent(out) :: productivity
    ! Locals:
    type(source_type) :: source
    PetscInt :: c, fluid_offset, source_offset
    PetscReal, allocatable :: phase_mobilities(:)
    PetscReal :: reference_pressure, pressure_difference, factor
    PetscReal, parameter :: tol = 1.e-9_dp

    call source%init(eos)
    source_offset = global_section_offset(source_section, self%local_source_index, &
         source_range_start)
    call source%assign(source_data, source_offset)

    c = nint(source%local_cell_index)
    if (fluid_range_start >= 0) then ! global
       fluid_offset = global_section_offset(fluid_section, c, &
            fluid_range_start)
    else ! local
       fluid_offset = section_offset(fluid_section, c)
    end if

    call source%fluid%assign(fluid_data, fluid_offset)
    allocate(phase_mobilities(source%fluid%num_phases))
    phase_mobilities = source%fluid%phase_mobilities()

    reference_pressure = self%reference_pressure%interpolate(time, 1)
    pressure_difference = source%fluid%pressure - reference_pressure
    factor = sum(phase_mobilities) * pressure_difference

    if (abs(factor) > tol) then
       productivity = abs(rate) / factor
    end if

    deallocate(phase_mobilities)
    call source%destroy()

  end subroutine source_control_deliverability_calculate_PI_from_rate

!------------------------------------------------------------------------
! Recharge source control:
!------------------------------------------------------------------------

  subroutine source_control_recharge_init(self, recharge_data, &
       interpolation_type, averaging_type, reference_pressure_data, &
       local_source_index, err)
    !! Initialises source_control_recharge object. Error flag err
    !! returns 1 if there are problems with the recharge coefficient
    !! array, or 2 if there are problems with the reference pressure
    !! array.

    class(source_control_recharge_type), intent(in out) :: self
    PetscReal, intent(in) :: recharge_data(:,:)
    PetscInt, intent(in) :: interpolation_type, averaging_type
    PetscReal, intent(in) :: reference_pressure_data(:,:)
    PetscInt, intent(in) :: local_source_index
    PetscErrorCode, intent(out) :: err

    call self%coefficient%init(recharge_data, &
         interpolation_type, averaging_type, err)
    if (err == 0) then
       call self%reference_pressure%init(reference_pressure_data, &
            interpolation_type, averaging_type, err)
       if (err > 0) then
          err = 2
       else
          self%local_source_index = local_source_index
       end if
    end if

  end subroutine source_control_recharge_init

!------------------------------------------------------------------------

  subroutine source_control_recharge_destroy(self)
    !! Destroys source_control_recharge object.

    class(source_control_recharge_type), intent(in out) :: self

    call self%coefficient%destroy()
    call self%reference_pressure%destroy()

  end subroutine source_control_recharge_destroy

!------------------------------------------------------------------------

  subroutine source_control_recharge_update(self, t, interval, &
       source_data, source_section, source_range_start, &
       local_fluid_data, local_fluid_section, eos)
    !! Update flow rate for source_control_recharge_type.

    use dm_utils_module, only: global_section_offset

    class(source_control_recharge_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    ! Locals:
    type(source_type) :: source
    PetscInt :: source_offset
    PetscReal :: reference_pressure, pressure_difference
    PetscReal :: recharge_coefficient

    call source%init(eos)

    source_offset = global_section_offset(source_section, &
         self%local_source_index, source_range_start)
    call source%assign(source_data, source_offset)
    call source%assign_fluid(local_fluid_data, local_fluid_section)

    reference_pressure = self%reference_pressure%average(interval, 1)
    pressure_difference = source%fluid%pressure - reference_pressure
    recharge_coefficient = self%coefficient%average(interval, 1)
    source%rate = -recharge_coefficient * pressure_difference

    call source%destroy()

  end subroutine source_control_recharge_update

!------------------------------------------------------------------------
! Separator source control:
!------------------------------------------------------------------------

  subroutine source_control_separator_init(self, local_source_index, &
       thermo, separator_pressure)
    !! Initialises source_control_separator object.

    class(source_control_separator_type), intent(in out) :: self
    PetscInt, intent(in) :: local_source_index
    class(thermodynamics_type), target, intent(in) :: thermo
    PetscReal, intent(in) :: separator_pressure
    ! Locals:
    PetscReal :: saturation_temperature
    PetscReal :: params(2), water_props(2), steam_props(2)
    PetscErrorCode :: err

    self%local_source_index = local_source_index
    self%thermo => thermo
    self%separator_pressure = separator_pressure

    call self%thermo%saturation%temperature(self%separator_pressure, &
         saturation_temperature, err)
    params = [self%separator_pressure, saturation_temperature]
    call self%thermo%water%properties(params, water_props, err)
    call self%thermo%steam%properties(params, steam_props, err)

    associate(water_density => water_props(1), &
         water_internal_energy => water_props(2), &
         steam_density => steam_props(1), &
         steam_internal_energy => steam_props(2))
      self%water_enthalpy = water_internal_energy + &
           self%separator_pressure / water_density
      self%steam_enthalpy = steam_internal_energy + &
           self%separator_pressure / steam_density
    end associate

  end subroutine source_control_separator_init

!------------------------------------------------------------------------

  subroutine source_control_separator_destroy(self)
    !! Destroys source_control_separator object.

    class(source_control_separator_type), intent(in out) :: self

    self%thermo => null()

  end subroutine source_control_separator_destroy

!------------------------------------------------------------------------

  subroutine source_control_separator_update_steam_fraction(self, &
       enthalpy)
    !! Updates steam fraction for given enthalpy, based on the
    !! separator pressure.

    class(source_control_separator_type), intent(in out) :: self
    PetscReal, intent(in) :: enthalpy

    if (enthalpy <= self%water_enthalpy) then
       self%steam_fraction = 0._dp
    else if (enthalpy <= self%steam_enthalpy) then
       self%steam_fraction = (enthalpy - self%water_enthalpy) / &
            (self%steam_enthalpy - self%water_enthalpy)
    else
       self%steam_fraction = 1._dp
    end if

  end subroutine source_control_separator_update_steam_fraction

!------------------------------------------------------------------------

  subroutine source_control_separator_update(self, t, interval, &
       source_data, source_section, source_range_start, &
       local_fluid_data, local_fluid_section, eos)
    !! Update separated water and steam flow rates for
    !! source_control_separator_type.

    use dm_utils_module, only: global_section_offset

    class(source_control_separator_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    ! Locals:
    type(source_type) :: source
    PetscInt :: source_offset
    PetscReal, allocatable :: phase_flow_fractions(:)
    PetscReal :: enthalpy

    call source%init(eos)
    source_offset = global_section_offset(source_section, &
         self%local_source_index, source_range_start)
    call source%assign(source_data, source_offset)

    if ((source%rate < 0._dp) .and. (source%production_component < &
         eos%num_primary_variables)) then

       call source%assign_fluid(local_fluid_data, local_fluid_section)
       allocate(phase_flow_fractions(source%fluid%num_phases))
       phase_flow_fractions = source%fluid%phase_flow_fractions()
       enthalpy = source%fluid%specific_enthalpy(phase_flow_fractions)
       deallocate(phase_flow_fractions)

       call self%update_steam_fraction(enthalpy)
       self%water_flow_rate = (1._dp - self%steam_fraction) * source%rate
       self%steam_flow_rate = self%steam_fraction * source%rate

    else

       self%water_flow_rate = 0._dp
       self%steam_flow_rate = 0._dp

    end if

    call source%destroy()

  end subroutine source_control_separator_update

!------------------------------------------------------------------------
! Limiter source control:
!------------------------------------------------------------------------

  subroutine source_control_limiter_init_source(self, limiter_type, &
       input_local_source_index, limit, local_source_indices)
    !! Initialises source_control_limiter object, with input taken
    !! from a source.

    class(source_control_limiter_type), intent(in out) :: self
    PetscInt, intent(in) :: limiter_type
    PetscInt, intent(in) :: input_local_source_index
    PetscReal, intent(in) :: limit
    PetscInt, intent(in) :: local_source_indices(:)

    self%type = limiter_type
    self%input_local_source_index = input_local_source_index
    self%limit = limit
    self%local_source_indices = local_source_indices

  end subroutine source_control_limiter_init_source

!------------------------------------------------------------------------

  subroutine source_control_limiter_init_control(self, &
       limiter_type, input_source_control, limit, local_source_indices)
    !! Initialises source_control_limiter object, with input taken
    !! from another source control.

    class(source_control_limiter_type), intent(in out) :: self
    PetscInt, intent(in) :: limiter_type
    class(source_control_type), target, intent(in) :: input_source_control
    PetscReal, intent(in) :: limit
    PetscInt, intent(in) :: local_source_indices(:)

    self%type = limiter_type
    self%input_source_control => input_source_control
    self%limit = limit
    self%local_source_indices = local_source_indices

  end subroutine source_control_limiter_init_control

!------------------------------------------------------------------------

  subroutine source_control_limiter_destroy(self)
    !! Destroys source_control_limiter_type object.

    class(source_control_limiter_type), intent(in out) :: self

    if (allocated(self%local_source_indices)) &
         deallocate(self%local_source_indices)
    self%input_source_control => null()

  end subroutine source_control_limiter_destroy

!------------------------------------------------------------------------

  PetscReal function source_control_limiter_get_rate(self, &
       source_data, source_section, source_range_start, eos) result (rate)
    !! Gets rate to limit from input source or source control,
    !! depending on limiter type.

    use dm_utils_module, only: global_section_offset

    class(source_control_limiter_type), intent(in) :: self
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    class(eos_type), intent(in) :: eos
    ! Locals:
    type(source_type) :: source
    PetscInt :: source_offset

    select case (self%type)
    case (SRC_CONTROL_LIMITER_TYPE_TOTAL)
       call source%init(eos)
       source_offset = global_section_offset(source_section, &
            self%input_local_source_index, source_range_start)
       call source%assign(source_data, source_offset)
       rate = source%rate
       call source%destroy()
    case (SRC_CONTROL_LIMITER_TYPE_WATER)
       select type (separator => self%input_source_control)
       type is (source_control_separator_type)
          rate = separator%water_flow_rate
       end select
    case (SRC_CONTROL_LIMITER_TYPE_STEAM)
       select type (separator => self%input_source_control)
       type is (source_control_separator_type)
          rate = separator%steam_flow_rate
       end select
    end select

  end function source_control_limiter_get_rate

!------------------------------------------------------------------------

  PetscReal function source_control_limiter_rate_scale(self, &
       source_data, source_section, source_range_start, eos) &
       result(scale)
    !! Returns factor by which flow should be scaled to avoid
    !! exceededing limit.

    class(source_control_limiter_type), intent(in) :: self
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    class(eos_type), intent(in) :: eos
    ! Locals:
    PetscReal :: rate, abs_rate
    PetscReal, parameter :: small = 1.e-6_dp

    rate = self%get_rate(source_data, source_section, &
         source_range_start, eos)
    abs_rate = abs(rate)

    if ((abs_rate > self%limit) .and. (abs_rate > small)) then
       scale = self%limit / abs_rate
    else
       scale = 1._dp
    end if

  end function source_control_limiter_rate_scale

!------------------------------------------------------------------------

  subroutine source_control_limiter_update(self, t, interval, &
       source_data, source_section, source_range_start, &
       local_fluid_data, local_fluid_section, eos)
    !! Update flow rate for source_control_limiter_type.

    use dm_utils_module, only: global_section_offset

    class(source_control_limiter_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    ! Locals:
    type(source_type) :: source
    PetscReal :: scale
    PetscInt :: i, s, source_offset

    call source%init(eos)
    scale = self%rate_scale(source_data, source_section, &
         source_range_start, eos)

    do i = 1, size(self%local_source_indices)

       s = self%local_source_indices(i)
       source_offset = global_section_offset(source_section, s, source_range_start)
       call source%assign(source_data, source_offset)
       source%rate = source%rate * scale

    end do

    call source%destroy()

  end subroutine source_control_limiter_update

!------------------------------------------------------------------------
! Direction source control:
!------------------------------------------------------------------------

  subroutine source_control_direction_init(self, direction, &
       local_source_indices)
    !! Initialises source_control_direction object.

    class(source_control_direction_type), intent(in out) :: self
    PetscInt, intent(in) :: direction
    PetscInt, intent(in) :: local_source_indices(:)

    self%direction = direction
    self%local_source_indices = local_source_indices

  end subroutine source_control_direction_init

!------------------------------------------------------------------------

  subroutine source_control_direction_destroy(self)
    !! Destroys source_control_direction_type object.

    class(source_control_direction_type), intent(in out) :: self

    if (allocated(self%local_source_indices)) &
         deallocate(self%local_source_indices)

  end subroutine source_control_direction_destroy

!------------------------------------------------------------------------

  subroutine source_control_direction_update(self, t, interval, &
       source_data, source_section, source_range_start, &
       local_fluid_data, local_fluid_section, eos)
    !! Update flow rate for source_control_direction_type.

    use dm_utils_module, only: global_section_offset

    class(source_control_direction_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    ! Locals:
    type(source_type) :: source
    PetscInt :: i, s, source_offset
    PetscBool :: flowing

    call source%init(eos)

    do i = 1, size(self%local_source_indices)

       s = self%local_source_indices(i)
       source_offset = global_section_offset(source_section, s, source_range_start)
       call source%assign(source_data, source_offset)

         select case (self%direction)
         case (SRC_DIRECTION_PRODUCTION)
            flowing = (source%rate < 0._dp)
         case (SRC_DIRECTION_INJECTION)
            flowing = (source%rate > 0._dp)
         case default
            flowing = PETSC_TRUE
         end select
         if (.not. flowing) source%rate = 0._dp

      end do

    call source%destroy()

  end subroutine source_control_direction_update

!------------------------------------------------------------------------

end module source_control_module
