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
  use eos_module, only: max_phase_name_length
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
     type(list_type), public :: sources
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

  type, public, abstract, extends(source_control_type) :: source_control_pressure_reference_type
     !! Controls a source by comparing fluid pressure with a reference
     !! pressure (e.g. deliverability or recharge).
     private
     type(list_type), public :: sources
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
     type(source_type), pointer :: source
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
     type(list_type), public :: sources
     PetscReal, pointer, public :: rate !! Pointer to flow rate variable
     PetscInt, public :: type !! Type of limiter (water, steam or total)
     PetscReal, public :: limit !! Flow limit
   contains
     private
     procedure :: rate_scale => source_control_limiter_rate_scale
     procedure, public :: init => source_control_limiter_init
     procedure, public :: destroy => source_control_limiter_destroy
     procedure, public :: update => source_control_limiter_update
  end type source_control_limiter_type

  type, public, extends(source_control_type) :: source_control_direction_type
     !! Allows source to flow only in a specified direction
     !! (production or injection).
     private
     type(list_type), public :: sources
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
          local_fluid_data, local_fluid_section)
       use petscis
       !! Updates sources at the specified time.
       import :: source_control_type
       class(source_control_type), intent(in out) :: self
       PetscReal, intent(in) :: t, interval(2)
       PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
       PetscSection, intent(in) :: local_fluid_section
     end subroutine source_control_update_procedure

  end interface

contains

!------------------------------------------------------------------------
! Table source control:
!------------------------------------------------------------------------

  subroutine source_control_table_init(self, data, interpolation_type, &
       averaging_type, sources, err)
    !! Initialises source_control_table object.

    class(source_control_table_type), intent(in out) :: self
    PetscReal, intent(in) :: data(:,:)
    PetscInt, intent(in) :: interpolation_type
    PetscInt, intent(in) :: averaging_type
    type(list_type), intent(in out) :: sources
    PetscErrorCode, intent(out) :: err

    call self%sources%init()
    call self%table%init(data, interpolation_type, averaging_type, err)
    if (err == 0) then
       call self%sources%add(sources)
    end if

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
       local_fluid_data, local_fluid_section)
    !! Update flow rate for source_control_rate_table_type.

    class(source_control_rate_table_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    ! Locals:
    PetscReal :: rate

    rate = self%table%average(interval, 1)
    call self%sources%traverse(source_control_rate_table_update_iterator)

  contains

    subroutine source_control_rate_table_update_iterator(node, stopped)
      !! Sets source rate at a list node.
      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      select type (source => node%data)
      type is (source_type)
         source%rate = rate
      end select
      stopped = PETSC_FALSE
    end subroutine source_control_rate_table_update_iterator

  end subroutine source_control_rate_table_update

!------------------------------------------------------------------------
! Enthalpy table source control:
!------------------------------------------------------------------------

  subroutine source_control_enthalpy_table_update(self, t, interval, &
       local_fluid_data, local_fluid_section)
    !! Update injection enthalpy for source_control_enthalpy_table_type.

    class(source_control_enthalpy_table_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    ! Locals:
    PetscReal :: enthalpy

    enthalpy = self%table%average(interval, 1)
    call self%sources%traverse(source_control_enthalpy_table_update_iterator)

  contains

    subroutine source_control_enthalpy_table_update_iterator(node, stopped)
      !! Sets source injection enthalpy at a list node.
      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      select type (source => node%data)
      type is (source_type)
         source%injection_enthalpy = enthalpy
      end select
      stopped = PETSC_FALSE
    end subroutine source_control_enthalpy_table_update_iterator

  end subroutine source_control_enthalpy_table_update

!------------------------------------------------------------------------
! Pressure reference source control:
!------------------------------------------------------------------------

  subroutine source_control_set_reference_pressure_initial(self, &
       global_fluid_data, global_fluid_section, fluid_range_start)
    !! Sets reference pressure for pressure reference control to be the
    !! initial fluid pressure.

    use dm_utils_module, only: global_section_offset

    class(source_control_pressure_reference_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: global_fluid_data(:)
    PetscSection, intent(in) :: global_fluid_section
    PetscInt, intent(in) :: fluid_range_start
    ! Locals:
    type(list_node_type), pointer :: node
    PetscInt :: c, fluid_offset
    PetscErrorCode :: ierr

    if (self%sources%count == 1) then
       node => self%sources%head
       select type (source => node%data)
       type is (source_type)

          c = source%cell_index
          call global_section_offset(global_fluid_section, c, &
               fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)

          call source%fluid%assign(global_fluid_data, fluid_offset)

          self%reference_pressure%val(1, 1) = source%fluid%pressure

       end select
    end if

  end subroutine source_control_set_reference_pressure_initial

!------------------------------------------------------------------------
! Deliverability source control:
!------------------------------------------------------------------------

  subroutine source_control_deliverability_init(self, productivity_data, &
       interpolation_type, averaging_type, reference_pressure_data, &
       pressure_table_coordinate, threshold, sources, err)
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
    type(list_type), intent(in out) :: sources
    PetscErrorCode, intent(out) :: err

    call self%sources%init()
    call self%sources%add(sources)
    call self%productivity%init(productivity_data, &
         interpolation_type, averaging_type, err)
    if (err == 0) then
       call self%reference_pressure%init(reference_pressure_data, &
            interpolation_type, averaging_type, err)
       if (err == 0) then
          self%pressure_table_coordinate = pressure_table_coordinate
          self%threshold = threshold
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
    call self%sources%destroy()

  end subroutine source_control_deliverability_destroy

!------------------------------------------------------------------------

  subroutine source_control_deliverability_update(self, t, interval, &
       local_fluid_data, local_fluid_section)
    !! Update flow rate for source_control_deliverability_type.

    class(source_control_deliverability_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section

    call self%sources%traverse(source_control_deliverability_update_iterator)

  contains

    subroutine source_control_deliverability_update_iterator(node, stopped)
      !! Applies deliverability control at a list node.
      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: p, phases
      PetscReal :: h, reference_pressure, pressure_difference, productivity
      PetscReal, allocatable :: phase_mobilities(:), phase_flow_fractions(:)

      select type (source => node%data)
      type is (source_type)

         call source%update_fluid(local_fluid_data, local_fluid_section)

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

         productivity = self%productivity%average(interval, 1)

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

      end select

      stopped = PETSC_FALSE

    end subroutine source_control_deliverability_update_iterator

  end subroutine source_control_deliverability_update

!------------------------------------------------------------------------

  subroutine source_control_deliverability_calculate_PI_from_rate(&
       self, time, rate, global_fluid_data, &
       global_fluid_section, fluid_range_start)
    !! Calculates productivity index for deliverability control, from
    !! specified initial flow rate. This only works if the control has
    !! exactly one source, otherwise the correct productivity index
    !! would not be well-defined. If the productivity index can't be
    !! calculated, it is left at its initial value.

    use dm_utils_module, only: global_section_offset

    class(source_control_deliverability_type), intent(in out) :: self
    PetscReal, intent(in) :: time
    PetscReal, intent(in) :: rate
    PetscReal, pointer, contiguous, intent(in) :: global_fluid_data(:)
    PetscSection, intent(in) :: global_fluid_section
    PetscInt, intent(in) :: fluid_range_start
    ! Locals:
    type(list_node_type), pointer :: node
    PetscInt :: c, fluid_offset
    PetscReal, allocatable :: phase_mobilities(:)
    PetscReal :: reference_pressure, pressure_difference, factor
    PetscErrorCode :: ierr
    PetscReal, parameter :: tol = 1.e-9_dp

    if (self%sources%count == 1) then
       node => self%sources%head
       select type (source => node%data)
       type is (source_type)

          c = source%cell_index
          call global_section_offset(global_fluid_section, c, &
               fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)

          call source%fluid%assign(global_fluid_data, fluid_offset)
          allocate(phase_mobilities(source%fluid%num_phases))
          phase_mobilities = source%fluid%phase_mobilities()

          reference_pressure = self%reference_pressure%interpolate(time, 1)
          pressure_difference = source%fluid%pressure - reference_pressure
          factor = sum(phase_mobilities) * pressure_difference

          if (abs(factor) > tol) then
             self%productivity%val(1, 1) = abs(rate) / factor
          end if

          deallocate(phase_mobilities)

       end select
    end if

  end subroutine source_control_deliverability_calculate_PI_from_rate

!------------------------------------------------------------------------
! Recharge source control:
!------------------------------------------------------------------------

  subroutine source_control_recharge_init(self, recharge_data, &
       interpolation_type, averaging_type, reference_pressure_data, &
       sources, err)
    !! Initialises source_control_recharge object. Error flag err
    !! returns 1 if there are problems with the recharge coefficient
    !! array, or 2 if there are problems with the reference pressure
    !! array.

    class(source_control_recharge_type), intent(in out) :: self
    PetscReal, intent(in) :: recharge_data(:,:)
    PetscInt, intent(in) :: interpolation_type, averaging_type
    PetscReal, intent(in) :: reference_pressure_data(:,:)
    type(list_type), intent(in out) :: sources
    PetscErrorCode, intent(out) :: err

    call self%sources%init()
    call self%sources%add(sources)
    call self%coefficient%init(recharge_data, &
         interpolation_type, averaging_type, err)
    if (err == 0) then
       call self%reference_pressure%init(reference_pressure_data, &
            interpolation_type, averaging_type, err)
       if (err > 0) err = 2
    end if

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
       local_fluid_data, local_fluid_section)
    !! Update flow rate for source_control_recharge_type.

    class(source_control_recharge_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section

    call self%sources%traverse(source_control_recharge_update_iterator)

  contains

    subroutine source_control_recharge_update_iterator(node, stopped)
      !! Applies recharge control at a list node.
      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscReal :: reference_pressure, pressure_difference, &
           recharge_coefficient

      select type (source => node%data)
      type is (source_type)

         call source%update_fluid(local_fluid_data, local_fluid_section)
         reference_pressure = self%reference_pressure%average(interval, 1)
         pressure_difference = source%fluid%pressure - reference_pressure
         recharge_coefficient = self%coefficient%average(interval, 1)
         source%rate = -recharge_coefficient * pressure_difference

      end select

      stopped = PETSC_FALSE

    end subroutine source_control_recharge_update_iterator

  end subroutine source_control_recharge_update

!------------------------------------------------------------------------
! Separator source control:
!------------------------------------------------------------------------

  subroutine source_control_separator_init(self, source, &
       thermo, separator_pressure)
    !! Initialises source_control_separator object.

    class(source_control_separator_type), intent(in out) :: self
    type(source_type), target, intent(in) :: source
    class(thermodynamics_type), target, intent(in) :: thermo
    PetscReal, intent(in) :: separator_pressure
    ! Locals:
    PetscReal :: saturation_temperature
    PetscReal :: params(2), water_props(2), steam_props(2)
    PetscErrorCode :: err

    self%source => source
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

    self%source => null()
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
       local_fluid_data, local_fluid_section)
    !! Update separated water and steam flow rates for
    !! source_control_separator_type.

    class(source_control_separator_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    ! Locals:
    PetscReal :: phase_flow_fractions(self%source%fluid%num_phases)
    PetscReal :: enthalpy

    associate(np => size(self%source%flow))

      if ((self%source%rate < 0._dp) .and. &
           (self%source%production_component < np)) then

         call self%source%update_fluid(local_fluid_data, local_fluid_section)
         phase_flow_fractions = self%source%fluid%phase_flow_fractions()
         enthalpy = self%source%fluid%specific_enthalpy(phase_flow_fractions)

         call self%update_steam_fraction(enthalpy)
         self%water_flow_rate = (1._dp - self%steam_fraction) * self%source%rate
         self%steam_flow_rate = self%steam_fraction * self%source%rate

      else

         self%water_flow_rate = 0._dp
         self%steam_flow_rate = 0._dp

      end if

    end associate

  end subroutine source_control_separator_update

!------------------------------------------------------------------------
! Limiter source control:
!------------------------------------------------------------------------

  subroutine source_control_limiter_init(self, limiter_type, rate, &
       limit, sources)
    !! Initialises source_control_limiter object.

    class(source_control_limiter_type), intent(in out) :: self
    PetscInt, intent(in) :: limiter_type
    PetscReal, target, intent(in) :: rate
    PetscReal, intent(in) :: limit
    type(list_type), intent(in out) :: sources

    call self%sources%init()
    call self%sources%add(sources)
    self%type = limiter_type
    self%rate => rate
    self%limit = limit

  end subroutine source_control_limiter_init

!------------------------------------------------------------------------

  subroutine source_control_limiter_destroy(self)
    !! Destroys source_control_limiter_type object.

    class(source_control_limiter_type), intent(in out) :: self

    self%rate => null()
    call self%sources%destroy()

  end subroutine source_control_limiter_destroy

!------------------------------------------------------------------------

  PetscReal function source_control_limiter_rate_scale(self) &
       result(scale)
    !! Returns factor by which flow should be scaled to avoid
    !! exceededing limit.

    class(source_control_limiter_type), intent(in) :: self
    ! Locals:
    PetscReal :: abs_rate
    PetscReal, parameter :: small = 1.e-6_dp

    abs_rate = abs(self%rate)

    if ((abs_rate > self%limit) .and. (abs_rate > small)) then
       scale = self%limit / abs_rate
    else
       scale = 1._dp
    end if

  end function source_control_limiter_rate_scale

!------------------------------------------------------------------------

  subroutine source_control_limiter_update(self, t, interval, &
       local_fluid_data, local_fluid_section)
    !! Update flow rate for source_control_limiter_type.

    class(source_control_limiter_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    PetscReal :: scale

    scale = self%rate_scale()
    call self%sources%traverse(source_control_limiter_update_iterator)

  contains

    subroutine source_control_limiter_update_iterator(node, stopped)
      !! Limits source rate at a list node.

      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped

      select type (source => node%data)
      type is (source_type)
         source%rate = source%rate * scale
      end select

      stopped = PETSC_FALSE

    end subroutine source_control_limiter_update_iterator

  end subroutine source_control_limiter_update

!------------------------------------------------------------------------
! Direction source control:
!------------------------------------------------------------------------

  subroutine source_control_direction_init(self, direction, sources)
    !! Initialises source_control_direction object.

    class(source_control_direction_type), intent(in out) :: self
    PetscInt, intent(in) :: direction
    type(list_type), intent(in out) :: sources

    call self%sources%init()
    call self%sources%add(sources)
    self%direction = direction

  end subroutine source_control_direction_init

!------------------------------------------------------------------------

  subroutine source_control_direction_destroy(self)
    !! Destroys source_control_direction_type object.

    class(source_control_direction_type), intent(in out) :: self

    call self%sources%destroy()

  end subroutine source_control_direction_destroy

!------------------------------------------------------------------------

  subroutine source_control_direction_update(self, t, interval, &
       local_fluid_data, local_fluid_section)
    !! Update flow rate for source_control_direction_type.

    class(source_control_direction_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section

    call self%sources%traverse(source_control_direction_update_iterator)

  contains

    subroutine source_control_direction_update_iterator(node, stopped)
      !! Controls source rate direction at a list node.

      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscBool :: flowing

      select type (source => node%data)
      type is (source_type)

         select case (self%direction)
         case (SRC_DIRECTION_PRODUCTION)
            flowing = (source%rate < 0._dp)
         case (SRC_DIRECTION_INJECTION)
            flowing = (source%rate > 0._dp)
         case default
            flowing = PETSC_TRUE
         end select
         if (.not. flowing) source%rate = 0._dp

      end select

      stopped = PETSC_FALSE

    end subroutine source_control_direction_update_iterator

  end subroutine source_control_direction_update

!------------------------------------------------------------------------

end module source_control_module
