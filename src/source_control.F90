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

  PetscInt, parameter, public :: max_limiter_type_length = 5
  character(max_limiter_type_length), parameter, public :: &
       default_source_control_limiter_type_str = "total"
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
     type(list_type), public :: sources
   contains
     procedure(source_control_destroy_procedure), public, deferred :: destroy
     procedure(source_control_update_procedure), public, deferred :: update
  end type source_control_type

  type, abstract, extends(source_control_type) :: source_control_table_type
     !! Controls a source parameter (e.g. rate or enthalpy) via a
     !! table of values vs. time.
     private
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

  type, public, extends(source_control_table_type) :: source_control_tracer_table_type
     !! Controls source tracer injection rate via a table of values
     !! vs. time. Each control applies to only one tracer.
     private
     PetscInt, public :: tracer_index !! Index of tracer
   contains
     private
     procedure, public :: update => source_control_tracer_table_update
  end type source_control_tracer_table_type

  type, public, abstract, extends(source_control_type) :: source_control_pressure_reference_type
     !! Controls a source by comparing fluid pressure with a reference
     !! pressure (e.g. deliverability or recharge).
     private
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

  type, abstract, public, extends(source_control_type) :: source_control_limiter_type
     !! Limits total flow, or separated steam or water flow (from a
     !! separator) through a source.
     private
     PetscReal, public :: limit !! Flow limit
   contains
     private
     procedure :: rate_scale => source_control_limiter_rate_scale
     procedure, public :: init => source_control_limiter_init
     procedure, public :: destroy => source_control_limiter_destroy
     procedure, public :: update => source_control_limiter_update
     procedure(source_control_limiter_get_rate_procedure), public, deferred :: get_rate
  end type source_control_limiter_type

  type, public, extends(source_control_limiter_type) :: source_control_total_limiter_type
     !! Limiter based on total flow in source.
   contains
     private
     procedure, public :: get_rate => source_control_total_limiter_get_rate
  end type source_control_total_limiter_type

  type, public, extends(source_control_limiter_type) :: source_control_water_limiter_type
     !! Limiter based on separated water flow in source.
   contains
     private
     procedure, public :: get_rate => source_control_water_limiter_get_rate
  end type source_control_water_limiter_type

  type, public, extends(source_control_limiter_type) :: source_control_steam_limiter_type
     !! Limiter based on separated steam flow in source.
   contains
     private
     procedure, public :: get_rate => source_control_steam_limiter_get_rate
  end type source_control_steam_limiter_type

  type, public, extends(source_control_type) :: source_control_direction_type
     !! Allows source to flow only in a specified direction
     !! (production or injection).
     private
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
          local_fluid_data, local_fluid_section, eos, num_tracers)
       use petscis
       !! Updates sources at the specified time.
       import :: source_control_type, eos_type
       class(source_control_type), intent(in out) :: self
       PetscReal, intent(in) :: t, interval(2)
       PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
       PetscSection, intent(in) :: local_fluid_section
       class(eos_type), intent(in) :: eos
       PetscInt, intent(in) :: num_tracers
     end subroutine source_control_update_procedure

     PetscReal function source_control_limiter_get_rate_procedure(self, eos)
       use petscis
       import :: source_control_limiter_type, eos_type
       class(source_control_limiter_type), intent(in out) :: self
       class(eos_type), intent(in) :: eos
     end function source_control_limiter_get_rate_procedure

  end interface

contains

!------------------------------------------------------------------------
! Table source control:
!------------------------------------------------------------------------

  subroutine source_control_table_init(self, data, interpolation_type, &
       averaging_type, sources)
    !! Initialises source_control_table object.

    class(source_control_table_type), intent(in out) :: self
    PetscReal, intent(in) :: data(:,:)
    PetscInt, intent(in) :: interpolation_type
    PetscInt, intent(in) :: averaging_type
    type(list_type), intent(in) :: sources

    call self%table%init(data, interpolation_type, averaging_type)
    self%sources = sources

  end subroutine source_control_table_init

!------------------------------------------------------------------------

  subroutine source_control_table_destroy(self)
    !! Destroys source_control_table_type object.

    class(source_control_table_type), intent(in out) :: self

    call self%table%destroy()
    call self%sources%destroy()

  end subroutine source_control_table_destroy

!------------------------------------------------------------------------
! Rate table source control:
!------------------------------------------------------------------------

  subroutine source_control_rate_table_update(self, t, interval, &
       local_fluid_data, local_fluid_section, eos, num_tracers)
    !! Update flow rate for source_control_rate_table_type.

    class(source_control_rate_table_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    PetscInt, intent(in) :: num_tracers
    ! Locals:
    PetscReal :: rate

    rate = self%table%average(interval, 1)
    call self%sources%traverse(rate_iterator)

  contains

    subroutine rate_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         call source%set_rate(rate)
      end select

    end subroutine rate_iterator

  end subroutine source_control_rate_table_update

!------------------------------------------------------------------------
! Enthalpy table source control:
!------------------------------------------------------------------------

  subroutine source_control_enthalpy_table_update(self, t, interval, &
       local_fluid_data, local_fluid_section, eos, num_tracers)
    !! Update injection enthalpy for source_control_enthalpy_table_type.

    class(source_control_enthalpy_table_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    PetscInt, intent(in) :: num_tracers
    ! Locals:
    PetscReal :: enthalpy

    enthalpy = self%table%average(interval, 1)
    call self%sources%traverse(enthalpy_iterator)

  contains

    subroutine enthalpy_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         source%injection_enthalpy = enthalpy
      end select

    end subroutine enthalpy_iterator

  end subroutine source_control_enthalpy_table_update

!------------------------------------------------------------------------
! Rate factor source control:
!------------------------------------------------------------------------

  subroutine source_control_rate_factor_update(self, t, interval, &
       local_fluid_data, local_fluid_section, eos, num_tracers)
    !! Update flow rate for source_control_rate_factor_type.

    class(source_control_rate_factor_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    PetscInt, intent(in) :: num_tracers
    ! Locals:
    PetscReal :: factor

    factor = self%table%average(interval, 1)
    call self%sources%traverse(rate_factor_iterator)

  contains

    subroutine rate_factor_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         call source%set_rate(source%rate * factor)
      end select

    end subroutine rate_factor_iterator

  end subroutine source_control_rate_factor_update

!------------------------------------------------------------------------
! Tracer table source control:
!------------------------------------------------------------------------

  subroutine source_control_tracer_table_update(self, t, interval, &
       local_fluid_data, local_fluid_section, eos, num_tracers)
    !! Update tracer injection rate for source_control_tracer_table_type.

    class(source_control_tracer_table_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    PetscInt, intent(in) :: num_tracers
    ! Locals:
    PetscReal :: tracer_injection_rate

    tracer_injection_rate = self%table%average(interval, 1)
    call self%sources%traverse(tracer_iterator)

  contains

    subroutine tracer_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         source%tracer_injection_rate(self%tracer_index) = &
              tracer_injection_rate
      end select

    end subroutine tracer_iterator

  end subroutine source_control_tracer_table_update

!------------------------------------------------------------------------
! Pressure reference source control:
!------------------------------------------------------------------------

  subroutine source_control_set_reference_pressure_initial(self, &
       global_fluid_data, global_fluid_section, fluid_range_start, eos)
    !! Sets reference pressure for pressure reference control to be
    !! the initial fluid pressure in the source cell.

    use dm_utils_module, only: global_section_offset

    class(source_control_pressure_reference_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: global_fluid_data(:)
    PetscSection, intent(in) :: global_fluid_section
    PetscInt, intent(in) :: fluid_range_start
    class(eos_type), intent(in) :: eos

    call self%sources%traverse(ref_pressure_iterator)

  contains

    subroutine ref_pressure_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: c, fluid_offset

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         c = source%local_cell_index
         fluid_offset = global_section_offset(global_fluid_section, c, &
              fluid_range_start)
         call source%fluid%assign(global_fluid_data, fluid_offset)
         self%reference_pressure%val(1, 1) = source%fluid%pressure
      end select

    end subroutine ref_pressure_iterator

  end subroutine source_control_set_reference_pressure_initial

!------------------------------------------------------------------------
! Deliverability source control:
!------------------------------------------------------------------------

  subroutine source_control_deliverability_init(self, productivity_data, &
       interpolation_type, averaging_type, reference_pressure_data, &
       pressure_table_coordinate, threshold, sources)
    !! Initialises source_control_deliverability object.

    class(source_control_deliverability_type), intent(in out) :: self
    PetscReal, intent(in) :: productivity_data(:,:)
    PetscInt, intent(in) :: interpolation_type, averaging_type
    PetscReal, intent(in) :: reference_pressure_data(:,:)
    PetscInt, intent(in) :: pressure_table_coordinate
    PetscReal, intent(in) :: threshold
    type(list_type), intent(in) :: sources

    call self%productivity%init(productivity_data, &
         interpolation_type, averaging_type)
    call self%reference_pressure%init(reference_pressure_data, &
         interpolation_type, averaging_type)
    self%pressure_table_coordinate = pressure_table_coordinate
    self%threshold = threshold
    self%sources = sources

  end subroutine source_control_deliverability_init

!------------------------------------------------------------------------

  subroutine source_control_deliverability_destroy(self)
    !! Destroys source_control_deliverability object.

    class(source_control_deliverability_type), intent(in out) :: self

    call self%productivity%destroy()
    call self%reference_pressure%destroy()
    call self%sources%destroy()

  end subroutine source_control_deliverability_destroy

!------------------------------------------------------------------------

  subroutine source_control_deliverability_update(self, t, interval, &
       local_fluid_data, local_fluid_section, eos, num_tracers)
    !! Update flow rate for source_control_deliverability_type.

    class(source_control_deliverability_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    PetscInt, intent(in) :: num_tracers

    call self%sources%traverse(deliverability_iterator)

  contains

    subroutine deliverability_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscReal :: productivity, qd

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)

         call source%assign_fluid_local(local_fluid_data, local_fluid_section)
         
         if (self%threshold <= 0._dp) then
            productivity = self%productivity%average(interval, 1)
            call source%set_rate(flow_rate(source, productivity))
         else
            if (source%fluid%pressure < self%threshold) then
               qd = flow_rate(source, self%threshold_productivity)
               if (qd > source%rate) then
                  call source%set_rate(qd)
               else ! don't use qd, but update PI
                  call self%calculate_PI_from_rate(t, source%rate, &
                       local_fluid_data, local_fluid_section, -1, &
                       eos, self%threshold_productivity)
               end if
            else
               call self%calculate_PI_from_rate(t, source%rate, &
                    local_fluid_data, local_fluid_section, -1, &
                    eos, self%threshold_productivity)
            end if
         end if

      end select

    end subroutine deliverability_iterator

!........................................................................

    PetscReal function flow_rate(source, productivity)
      !! Computes flow rate using deliverability relation and the
      !! specified productivity index.

      type(source_type), intent(in)  :: source
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
      flow_rate = 0._dp

      phases = nint(source%fluid%phase_composition)
      do p = 1, source%fluid%num_phases
         if (btest(phases, p - 1)) then
            flow_rate = flow_rate - productivity * &
                 phase_mobilities(p) * pressure_difference
         end if
      end do

      deallocate(phase_mobilities, phase_flow_fractions)

    end function flow_rate

  end subroutine source_control_deliverability_update

!------------------------------------------------------------------------

  subroutine source_control_deliverability_calculate_PI_from_rate(&
       self, time, rate, fluid_data, fluid_section, fluid_range_start, &
       eos, productivity)
    !! Calculates productivity index for deliverability control, from
    !! specified initial flow rate.

    use dm_utils_module, only: global_section_offset, section_offset

    class(source_control_deliverability_type), intent(in out) :: self
    PetscReal, intent(in) :: time
    PetscReal, intent(in) :: rate
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start !! Specify -1 for local data rather than global
    class(eos_type), intent(in) :: eos
    PetscReal, intent(out) :: productivity

    call self%sources%traverse(PI_iterator)

  contains

    subroutine PI_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: c, fluid_offset
      PetscReal, allocatable :: phase_mobilities(:)
      PetscReal :: reference_pressure, pressure_difference, factor
      PetscReal, parameter :: tol = 1.e-9_dp

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)

         c = source%local_cell_index
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

      end select

    end subroutine PI_iterator

  end subroutine source_control_deliverability_calculate_PI_from_rate

!------------------------------------------------------------------------
! Recharge source control:
!------------------------------------------------------------------------

  subroutine source_control_recharge_init(self, recharge_data, &
       interpolation_type, averaging_type, reference_pressure_data, &
       sources)
    !! Initialises source_control_recharge object.

    class(source_control_recharge_type), intent(in out) :: self
    PetscReal, intent(in) :: recharge_data(:,:)
    PetscInt, intent(in) :: interpolation_type, averaging_type
    PetscReal, intent(in) :: reference_pressure_data(:,:)
    type(list_type), intent(in) :: sources

    call self%coefficient%init(recharge_data, &
         interpolation_type, averaging_type)
    call self%reference_pressure%init(reference_pressure_data, &
         interpolation_type, averaging_type)
    self%sources = sources

  end subroutine source_control_recharge_init

!------------------------------------------------------------------------

  subroutine source_control_recharge_destroy(self)
    !! Destroys source_control_recharge object.

    class(source_control_recharge_type), intent(in out) :: self

    call self%coefficient%destroy()
    call self%reference_pressure%destroy()
    call self%sources%destroy()

  end subroutine source_control_recharge_destroy

!------------------------------------------------------------------------

  subroutine source_control_recharge_update(self, t, interval, &
       local_fluid_data, local_fluid_section, eos, num_tracers)
    !! Update flow rate for source_control_recharge_type.

    class(source_control_recharge_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    PetscInt, intent(in) :: num_tracers

    call self%sources%traverse(recharge_iterator)

  contains

    subroutine recharge_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscReal :: reference_pressure, pressure_difference
      PetscReal :: recharge_coefficient

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         call source%assign_fluid_local(local_fluid_data, local_fluid_section)
         reference_pressure = self%reference_pressure%average(interval, 1)
         pressure_difference = source%fluid%pressure - reference_pressure
         recharge_coefficient = self%coefficient%average(interval, 1)
         call source%set_rate(-recharge_coefficient * pressure_difference)
      end select

    end subroutine recharge_iterator

  end subroutine source_control_recharge_update

!------------------------------------------------------------------------
! Limiter source controls:
!------------------------------------------------------------------------

  subroutine source_control_limiter_init(self, limit, sources)
    !! Initialises source_control_limiter object.

    class(source_control_limiter_type), intent(in out) :: self
    PetscReal, intent(in) :: limit
    type(list_type), intent(in) :: sources

    self%limit = limit
    self%sources = sources

  end subroutine source_control_limiter_init

!------------------------------------------------------------------------

  subroutine source_control_limiter_destroy(self)
    !! Destroys source_control_limiter_type object.

    class(source_control_limiter_type), intent(in out) :: self

    call self%sources%destroy()

  end subroutine source_control_limiter_destroy

!------------------------------------------------------------------------

  PetscReal function source_control_limiter_rate_scale(self, eos) &
       result(scale)
    !! Returns factor by which flow should be scaled to avoid
    !! exceededing limit.

    class(source_control_limiter_type), intent(in out) :: self
    class(eos_type), intent(in) :: eos
    ! Locals:
    PetscReal :: rate, abs_rate
    PetscReal, parameter :: small = 1.e-6_dp

    rate = self%get_rate(eos)
    abs_rate = abs(rate)

    if ((abs_rate > self%limit) .and. (abs_rate > small)) then
       scale = self%limit / abs_rate
    else
       scale = 1._dp
    end if

  end function source_control_limiter_rate_scale

!------------------------------------------------------------------------

  subroutine source_control_limiter_update(self, t, interval, &
       local_fluid_data, local_fluid_section, eos, num_tracers)
    !! Update flow rate for source_control_limiter_type.

    class(source_control_limiter_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    PetscInt, intent(in) :: num_tracers
    ! Locals:
    PetscReal :: scale

    scale = self%rate_scale(eos)
    call self%sources%traverse(limiter_iterator)

  contains

    subroutine limiter_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         call source%set_rate(source%rate * scale)
      end select

    end subroutine limiter_iterator

  end subroutine source_control_limiter_update

!------------------------------------------------------------------------

  PetscReal function source_control_total_limiter_get_rate(self, eos) &
       result (rate)
    !! Gets total rate to limit from input source.

    class(source_control_total_limiter_type), intent(in out) :: self
    class(eos_type), intent(in) :: eos

    call self%sources%traverse(total_limiter_rate_iterator)

  contains

    subroutine total_limiter_rate_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         rate = source%rate
      end select

    end subroutine total_limiter_rate_iterator

  end function source_control_total_limiter_get_rate

!------------------------------------------------------------------------

  PetscReal function source_control_water_limiter_get_rate(self, eos) &
       result (rate)
    !! Gets separated water rate to limit from input source.

    class(source_control_water_limiter_type), intent(in out) :: self
    class(eos_type), intent(in) :: eos

    call self%sources%traverse(water_limiter_rate_iterator)

  contains

    subroutine water_limiter_rate_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         rate = source%water_rate
      end select

    end subroutine water_limiter_rate_iterator

  end function source_control_water_limiter_get_rate

!------------------------------------------------------------------------

  PetscReal function source_control_steam_limiter_get_rate(self, eos) &
       result (rate)
    !! Gets separated steam rate to limit from input source.

    class(source_control_steam_limiter_type), intent(in out) :: self
    class(eos_type), intent(in) :: eos

    call self%sources%traverse(steam_limiter_rate_iterator)

  contains

    subroutine steam_limiter_rate_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         rate = source%steam_rate
      end select

    end subroutine steam_limiter_rate_iterator

  end function source_control_steam_limiter_get_rate

!------------------------------------------------------------------------
! Direction source control:
!------------------------------------------------------------------------

  subroutine source_control_direction_init(self, direction, sources)
    !! Initialises source_control_direction object.

    class(source_control_direction_type), intent(in out) :: self
    PetscInt, intent(in) :: direction
    type(list_type), intent(in) :: sources

    self%direction = direction
    self%sources = sources

  end subroutine source_control_direction_init

!------------------------------------------------------------------------

  subroutine source_control_direction_destroy(self)
    !! Destroys source_control_direction_type object.

    class(source_control_direction_type), intent(in out) :: self

    call self%sources%destroy()

  end subroutine source_control_direction_destroy

!------------------------------------------------------------------------

  subroutine source_control_direction_update(self, t, interval, &
       local_fluid_data, local_fluid_section, eos, num_tracers)
    !! Update flow rate for source_control_direction_type.

    class(source_control_direction_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    class(eos_type), intent(in) :: eos
    PetscInt, intent(in) :: num_tracers

    call self%sources%traverse(direction_iterator)

  contains

    subroutine direction_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscBool :: flowing

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         select case (self%direction)
         case (SRC_DIRECTION_PRODUCTION)
            flowing = (source%rate < 0._dp)
         case (SRC_DIRECTION_INJECTION)
            flowing = (source%rate > 0._dp)
         case default
            flowing = PETSC_TRUE
         end select
         if (.not. flowing) call source%set_rate(0._dp)
      end select

    end subroutine direction_iterator

  end subroutine source_control_direction_update

!------------------------------------------------------------------------

end module source_control_module
