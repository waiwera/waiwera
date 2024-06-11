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
  use control_module
  use eos_module
  use thermodynamics_module
  use list_module
  use source_module
  use interpolation_module

  implicit none
  private

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

  type, public, extends(table_object_control_type) :: rate_table_source_control_type
     !! Controls source rate via a table of values vs. time.
   contains
     private
     procedure, public :: iterator => rate_table_source_control_iterator
  end type rate_table_source_control_type

  type, public, extends(table_object_control_type) :: enthalpy_table_source_control_type
     !! Controls source injection enthalpy via a table of values vs. time.
   contains
     private
     procedure, public :: iterator => enthalpy_table_source_control_iterator
  end type enthalpy_table_source_control_type

  type, public, extends(table_object_control_type) :: rate_factor_table_source_control_type
     !! Multiplies source rate by a factor from a table of values vs. time.
   contains
     private
     procedure, public :: iterator => rate_factor_table_source_control_iterator
  end type rate_factor_table_source_control_type

  type, public, extends(table_object_control_type) :: tracer_table_source_control_type
     !! Controls source tracer injection rate via a table of values
     !! vs. time. Each control applies to only one tracer.
     private
     PetscInt, public :: tracer_index !! Index of tracer
   contains
     private
     procedure, public :: iterator => tracer_table_source_control_iterator
  end type tracer_table_source_control_type


  type, public, extends(object_control_type) :: pressure_reference_source_control_type
     private
     class(interpolation_table_type), allocatable, public :: reference_pressure !! Reference pressure vs. time
     PetscReal :: time !! Time
     PetscReal :: interval(2) !! Time interval
     PetscReal, pointer, contiguous :: local_fluid_data(:)
     PetscSection, pointer :: local_fluid_section
   contains
     procedure, public :: set_reference_pressure_initial => &
          pressure_reference_source_control_set_initial
     procedure, public :: iterator => pressure_reference_source_control_iterator
     procedure, public :: update => pressure_reference_source_control_update
     procedure, public :: destroy => pressure_reference_source_control_destroy
  end type pressure_reference_source_control_type

  type, public, extends(pressure_reference_source_control_type) :: &
       deliverability_source_control_type
     !! Controls a source on deliverability.
     private
     PetscInt, public :: pressure_table_coordinate !! Coordinate variable of pressure table
     class(interpolation_table_type), allocatable, public :: productivity !! Productivity index vs. time
     PetscReal, public :: threshold !! Pressure threshold below which deliverability is switched on (< 0 for always on)
     PetscReal, public :: threshold_productivity !! Productivity index computed from flow rate and used when pressure drops below threshold
   contains
     private
     procedure :: flow_rate => deliverability_source_control_flow_rate
     procedure, public :: calculate_PI_from_rate => &
          deliverability_source_control_calculate_PI_from_rate
     procedure, public :: init => deliverability_source_control_init
     procedure, public :: iterator => deliverability_source_control_iterator
     procedure, public :: destroy => deliverability_source_control_destroy
  end type

  type, public, extends(pressure_reference_source_control_type) :: &
       recharge_source_control_type
     !! Controls a source simulating recharge through a model boundary.
     private
     class(interpolation_table_type), allocatable, public :: coefficient !! Recharge coefficient vs. time
   contains
     private
     procedure, public :: init => recharge_source_control_init
     procedure, public :: iterator => recharge_source_control_iterator
     procedure, public :: destroy => recharge_source_control_destroy
  end type

  type, public, extends(integer_object_control_type) :: direction_source_control_type
     !! Allows source to flow only in a specified direction
     !! (production or injection).
   contains
     private
     procedure, public :: iterator => direction_source_control_iterator
  end type direction_source_control_type

contains

!------------------------------------------------------------------------
! Rate table source control:
!------------------------------------------------------------------------

  subroutine rate_table_source_control_iterator(self, node, stopped)
    !! Update source rate.
    class(rate_table_source_control_type), intent(in out) :: self
    type(list_node_type), pointer, intent(in out) :: node
    PetscBool, intent(out) :: stopped

    stopped = PETSC_FALSE
    select type(source => node%data)
    type is (source_type)
       associate(rate => self%value(1))
         call source%set_rate(rate)
       end associate
    end select

  end subroutine rate_table_source_control_iterator

!------------------------------------------------------------------------
! Enthalpy table source control:
!------------------------------------------------------------------------

  subroutine enthalpy_table_source_control_iterator(self, node, stopped)
    !! Update source enthalpy.

    class(enthalpy_table_source_control_type), intent(in out) :: self
    type(list_node_type), pointer, intent(in out) :: node
    PetscBool, intent(out) :: stopped

    stopped = PETSC_FALSE
    select type(source => node%data)
    type is (source_type)
       associate(enthalpy => self%value(1))
         call source%set_enthalpy(enthalpy)
       end associate
    end select

  end subroutine enthalpy_table_source_control_iterator

!------------------------------------------------------------------------
! Rate factor source control:
!------------------------------------------------------------------------

  subroutine rate_factor_table_source_control_iterator(self, node, stopped)
    !! Update source rate using rate factor.

    class(rate_factor_table_source_control_type), intent(in out) :: self
    type(list_node_type), pointer, intent(in out) :: node
    PetscBool, intent(out) :: stopped

    stopped = PETSC_FALSE
    select type(source => node%data)
    type is (source_type)
       associate(factor => self%value(1))
         call source%set_rate(source%rate * factor)
       end associate
    end select

  end subroutine rate_factor_table_source_control_iterator

!------------------------------------------------------------------------
! Tracer table source control:
!------------------------------------------------------------------------

  subroutine tracer_table_source_control_iterator(self, node, stopped)
    !! Update tracer injection rate.

    class(tracer_table_source_control_type), intent(in out) :: self
    type(list_node_type), pointer, intent(in out) :: node
    PetscBool, intent(out) :: stopped

    stopped = PETSC_FALSE
    select type(source => node%data)
    type is (source_type)
       associate(tracer_injection_rate => self%value(1))
         source%tracer_injection_rate(self%tracer_index) = &
              tracer_injection_rate
       end associate
    end select

  end subroutine tracer_table_source_control_iterator

!------------------------------------------------------------------------
! Pressure reference source control:
!------------------------------------------------------------------------

  subroutine pressure_reference_source_control_set_initial(self, &
       global_fluid_data, global_fluid_section, fluid_range_start)
    !! Sets reference pressure for pressure reference control to be
    !! the initial fluid pressure in the source cell.

    use dm_utils_module, only: global_section_offset

    class(pressure_reference_source_control_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: global_fluid_data(:)
    PetscSection, intent(in) :: global_fluid_section
    PetscInt, intent(in) :: fluid_range_start

    call self%objects%traverse(ref_pressure_iterator)

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

  end subroutine pressure_reference_source_control_set_initial

!------------------------------------------------------------------------

  subroutine pressure_reference_source_control_iterator(self, node, stopped)
    !! Update iterator for pressure reference source control. Derived
    !! types override this routine to update their sources.

    class(pressure_reference_source_control_type), intent(in out) :: self
    type(list_node_type), pointer, intent(in out) :: node
    PetscBool, intent(out) :: stopped

    stopped = PETSC_FALSE

  end subroutine pressure_reference_source_control_iterator

!------------------------------------------------------------------------

  subroutine pressure_reference_source_control_update(self, time, interval, &
       local_fluid_data, local_fluid_section)
    !! Updates source flow rates for pressure reference source controls.

    class(pressure_reference_source_control_type), intent(in out) :: self
    PetscReal, intent(in) :: time !! Time
    PetscReal, intent(in) :: interval(2) !! Time interval
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:) !! Fluid data array
    PetscSection, target, intent(in) :: local_fluid_section !! Fluid section

    self%time = time
    self%interval = interval
    self%local_fluid_data => local_fluid_data
    self%local_fluid_section => local_fluid_section
    call self%objects%traverse(iterator)

  contains

    subroutine iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      call self%iterator(node, stopped)

    end subroutine iterator

  end subroutine pressure_reference_source_control_update

!------------------------------------------------------------------------

  subroutine pressure_reference_source_control_destroy(self)
    !! Destroys pressure reference source control.

    class(pressure_reference_source_control_type), intent(in out) :: self

    call self%objects%destroy()
    call self%reference_pressure%destroy()
    deallocate(self%reference_pressure)
    self%local_fluid_data => null()
    self%local_fluid_section => null()

  end subroutine pressure_reference_source_control_destroy

!------------------------------------------------------------------------
! Deliverability source control  
!------------------------------------------------------------------------

  subroutine deliverability_source_control_init(self, objects, &
       productivity_data, interpolation_type, averaging_type, &
       reference_pressure_data, pressure_table_coordinate, threshold)
    !! Initialises deliverability source_control.

    class(deliverability_source_control_type), intent(in out) :: self
    type(list_type), intent(in) :: objects
    PetscReal, intent(in) :: productivity_data(:,:)
    PetscInt, intent(in) :: interpolation_type, averaging_type
    PetscReal, intent(in) :: reference_pressure_data(:,:)
    PetscInt, intent(in) :: pressure_table_coordinate
    PetscReal, intent(in) :: threshold

    self%objects = objects

    select case (interpolation_type)
    case (INTERP_STEP)
       allocate(interpolation_table_step_type :: self%reference_pressure)
       allocate(interpolation_table_step_type :: self%productivity)
    case default
       allocate(interpolation_table_type :: self%reference_pressure)
       allocate(interpolation_table_type :: self%productivity)
    end select

    call self%productivity%init(productivity_data, averaging_type)
    call self%reference_pressure%init(reference_pressure_data, &
         averaging_type)
    self%pressure_table_coordinate = pressure_table_coordinate
    self%threshold = threshold

  end subroutine deliverability_source_control_init

!------------------------------------------------------------------------

  PetscReal function deliverability_source_control_flow_rate(self, source, &
       productivity) result(flow_rate)
    !! Computes source flow rate using the deliverability relation and
    !! the specified productivity index.

    class(deliverability_source_control_type), intent(in out) :: self
    type(source_type), intent(in)  :: source !! Source
    PetscReal, intent(in) :: productivity !! Productivity index
    ! Locals:
    PetscInt :: p, phases
    PetscReal :: h, reference_pressure, pressure_difference
    PetscReal, allocatable :: phase_mobilities(:), phase_flow_fractions(:)
    PetscReal :: effective_productivity

    allocate(phase_mobilities(source%fluid%num_phases))
    allocate(phase_flow_fractions(source%fluid%num_phases))
    phase_mobilities = source%fluid%phase_mobilities()
    phase_flow_fractions = phase_mobilities / sum(phase_mobilities)

    select case (self%pressure_table_coordinate)
    case (SRC_PRESSURE_TABLE_COORD_TIME)
       reference_pressure = self%reference_pressure%average(self%interval, 1)
    case (SRC_PRESSURE_TABLE_COORD_ENTHALPY)
       h = source%fluid%specific_enthalpy(phase_flow_fractions)
       reference_pressure = self%reference_pressure%interpolate(h, 1)
    end select

    effective_productivity = productivity * source%fluid%permeability_factor
    pressure_difference = source%fluid%pressure - reference_pressure
    flow_rate = 0._dp

    phases = nint(source%fluid%phase_composition)
    do p = 1, source%fluid%num_phases
       if (btest(phases, p - 1)) then
          flow_rate = flow_rate - effective_productivity * &
               phase_mobilities(p) * pressure_difference
       end if
    end do

    deallocate(phase_mobilities, phase_flow_fractions)

  end function deliverability_source_control_flow_rate

!------------------------------------------------------------------------

  subroutine deliverability_source_control_calculate_PI_from_rate(&
       self, time, rate, fluid_data, fluid_section, fluid_range_start, &
       productivity)
    !! Calculates productivity index for deliverability control, from
    !! specified initial flow rate.

    use dm_utils_module, only: global_section_offset, section_offset

    class(deliverability_source_control_type), intent(in out) :: self
    PetscReal, intent(in) :: time
    PetscReal, intent(in) :: rate
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start !! Specify -1 for local data rather than global
    PetscReal, intent(out) :: productivity

    call self%objects%traverse(PI_iterator)

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
         factor = sum(phase_mobilities) * pressure_difference * &
              source%fluid%permeability_factor

         if (abs(factor) > tol) then
            productivity = abs(rate) / factor
         end if

         deallocate(phase_mobilities)

      end select

    end subroutine PI_iterator

  end subroutine deliverability_source_control_calculate_PI_from_rate

!------------------------------------------------------------------------

  subroutine deliverability_source_control_iterator(self, node, stopped)
    !! Updates source flow rates for deliverability source control.

    class(deliverability_source_control_type), intent(in out) :: self
    type(list_node_type), pointer, intent(in out) :: node
    PetscBool, intent(out) :: stopped

    ! Locals:
    PetscReal :: productivity, qd

    stopped = PETSC_FALSE
    select type(source => node%data)
    type is (source_type)

       call source%assign_fluid_local(self%local_fluid_data, &
            self%local_fluid_section)

       if (self%threshold <= 0._dp) then
          productivity = self%productivity%average(self%interval, 1)
          call source%set_rate(self%flow_rate(source, productivity))
       else
          if (source%fluid%pressure < self%threshold) then
             qd = self%flow_rate(source, self%threshold_productivity)
             if (qd > source%rate) then
                call source%set_rate(qd)
             end if
          else
             call self%calculate_PI_from_rate(self%time, source%rate, &
                  self%local_fluid_data, self%local_fluid_section, -1, &
                  self%threshold_productivity)
          end if
       end if

    end select

  end subroutine deliverability_source_control_iterator

!------------------------------------------------------------------------

  subroutine deliverability_source_control_destroy(self)
    !! Destroys deliverability source control.
    class(deliverability_source_control_type), intent(in out) :: self

    call self%pressure_reference_source_control_type%destroy()
    call self%productivity%destroy()
    deallocate(self%productivity)

  end subroutine deliverability_source_control_destroy
  
!------------------------------------------------------------------------
! Recharge source control:
!------------------------------------------------------------------------

  subroutine recharge_source_control_init(self, objects, recharge_data, &
       interpolation_type, averaging_type, reference_pressure_data)
    !! Initialises recharge source control.

    class(recharge_source_control_type), intent(in out) :: self
    type(list_type), intent(in) :: objects
    PetscReal, intent(in) :: recharge_data(:,:)
    PetscInt, intent(in) :: interpolation_type, averaging_type
    PetscReal, intent(in) :: reference_pressure_data(:,:)

    self%objects = objects

    select case (interpolation_type)
    case (INTERP_STEP)
       allocate(interpolation_table_step_type :: self%reference_pressure)
       allocate(interpolation_table_step_type :: self%coefficient)
    case default
       allocate(interpolation_table_type :: self%reference_pressure)
       allocate(interpolation_table_type :: self%coefficient)
    end select

    call self%coefficient%init(recharge_data, averaging_type)
    call self%reference_pressure%init(reference_pressure_data, &
         averaging_type)

  end subroutine recharge_source_control_init

!------------------------------------------------------------------------

  subroutine recharge_source_control_iterator(self, node, stopped)
    !! Updates source flow rates for recharge source control.

    class(recharge_source_control_type), intent(in out) :: self
    type(list_node_type), pointer, intent(in out) :: node
    PetscBool, intent(out) :: stopped
    ! Locals:
    PetscReal :: reference_pressure, pressure_difference
    PetscReal :: recharge_coefficient

    stopped = PETSC_FALSE

    select type(source => node%data)
    type is (source_type)

       call source%assign_fluid_local(self%local_fluid_data, self%local_fluid_section)
       reference_pressure = self%reference_pressure%average(self%interval, 1)
       pressure_difference = source%fluid%pressure - reference_pressure
       recharge_coefficient = self%coefficient%average(self%interval, 1)
       call source%set_rate(-recharge_coefficient * pressure_difference)

    end select

  end subroutine recharge_source_control_iterator

!------------------------------------------------------------------------

  subroutine recharge_source_control_destroy(self)
    !! Destroys recharge source control.

    class(recharge_source_control_type), intent(in out) :: self

    call self%pressure_reference_source_control_type%destroy()
    call self%coefficient%destroy()
    deallocate(self%coefficient)

  end subroutine recharge_source_control_destroy

!------------------------------------------------------------------------
! Direction source control:
!------------------------------------------------------------------------

  subroutine direction_source_control_iterator(self, node, stopped)

    class(direction_source_control_type), intent(in out) :: self
    type(list_node_type), pointer, intent(in out) :: node
    PetscBool, intent(out) :: stopped
    ! Locals:
    PetscBool :: flowing

    stopped = PETSC_FALSE
    select type(source => node%data)
    type is (source_type)
       associate(direction => self%value)
         select case (direction)
         case (SRC_DIRECTION_PRODUCTION)
            flowing = (source%rate < 0._dp)
         case (SRC_DIRECTION_INJECTION)
            flowing = (source%rate > 0._dp)
         case default
            flowing = PETSC_TRUE
         end select
         if (.not. flowing) call source%set_rate(0._dp)
       end associate
    end select

  end subroutine direction_source_control_iterator

!------------------------------------------------------------------------

end module source_control_module
