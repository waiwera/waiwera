!   Copyright 2021 University of Auckland.

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

module source_network_node_module
  !! Module for source network nodes.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use list_module
  use separator_module

  implicit none
  private

  PetscInt, parameter, public :: max_source_network_node_name_length = 32
  PetscInt, parameter, public :: num_source_network_node_variables = 6
  PetscInt, parameter, public :: max_source_network_variable_name_length = 24
  character(max_source_network_variable_name_length), parameter, public :: &
       source_network_variable_names(num_source_network_node_variables) = [ &
       "rate                ", "enthalpy            ", &
       "water_rate          ", "water_enthalpy      ", &
       "steam_rate          ", "steam_enthalpy      " ]

  type, public :: source_network_node_type
     !! Type for node in source network, e.g. source or source group.
     private
     character(max_source_network_node_name_length), public :: name !! Name of source network node
     PetscReal, pointer, public :: rate !! Flow rate
     PetscReal, pointer, public :: enthalpy !! Enthalpy of fluid
     PetscReal, pointer, public :: water_rate !! Separated water mass flow rate
     PetscReal, pointer, public :: water_enthalpy !! Separated water enthalpy
     PetscReal, pointer, public :: steam_rate !! Separated steam mass flow rate
     PetscReal, pointer, public :: steam_enthalpy !! Separated steam enthalpy
     type(separator_type), public :: separator !! Separator
     PetscBool, public :: heat !! Whether flow through node is heat-only (not mass)
     PetscInt, public :: link_index !! Index of node in group it outputs to, or reinjector it is output from (or -1)
   contains
     private
     procedure, public :: assign => source_network_node_assign
     procedure, public :: get_separated_flows => source_network_node_get_separated_flows
     procedure, public :: default_separated_flows => source_network_node_default_separated_flows
     procedure, public :: zero_separated_flows => source_network_node_zero_separated_flows
     procedure, public :: separate => source_network_node_separate
     procedure, public :: set_rate => source_network_node_set_rate
     procedure, public :: scale_rate => source_network_node_scale_rate
     procedure, public :: get_rate_by_type => source_network_node_get_rate_by_type
     procedure, public :: limit_rate => source_network_node_limit_rate
     procedure, public :: limit_inputs => source_network_node_limit_inputs
     procedure, public :: is_over => source_network_node_is_over
     procedure, public :: get_limit_scale => source_network_node_get_limit_scale
     procedure, public :: get_minimum_limit_scale => source_network_node_get_minimum_limit_scale
     procedure, public :: add_flows => source_network_node_add_flows
     procedure, public :: add_separated_flows => source_network_node_add_separated_flows
     procedure, public :: destroy => source_network_node_destroy
  end type source_network_node_type

  public :: source_network_node_list_node_data_destroy

contains

!------------------------------------------------------------------------

  subroutine source_network_node_list_node_data_destroy(node)
    !! Destroys source network node in each list node.

    type(list_node_type), pointer, intent(in out) :: node

    select type (source_network_node => node%data)
    class is (source_network_node_type)
       call source_network_node%destroy()
    end select

  end subroutine source_network_node_list_node_data_destroy

!------------------------------------------------------------------------
! Source network node
!------------------------------------------------------------------------

  subroutine source_network_node_assign(self, data, offset)
    !! Assigns pointers in network node object to elements in the data
    !! array, starting from the specified offset.

    class(source_network_node_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! source data array
    PetscInt, intent(in) :: offset  !! source array offset

    self%rate => data(offset)
    self%enthalpy => data(offset + 1)
    self%water_rate => data(offset + 2)
    self%water_enthalpy => data(offset + 3)
    self%steam_rate => data(offset + 4)
    self%steam_enthalpy => data(offset + 5)

    call self%separator%assign(data, offset + 6)

  end subroutine source_network_node_assign

!------------------------------------------------------------------------

    subroutine source_network_node_get_separated_flows(self)
      !! Gets separated water and steam flows.

      class(source_network_node_type), intent(in out) :: self

      if (self%rate < 0._dp) then
         if (self%separator%on) then
            call self%separate()
         else
            call self%default_separated_flows()
         end if
      else
         call self%zero_separated_flows()
      end if

    end subroutine source_network_node_get_separated_flows

!------------------------------------------------------------------------

    subroutine source_network_node_default_separated_flows(self)
      !! Gets default separated water and steam flows when the node
      !! does not have its own separator.

      class(source_network_node_type), intent(in out) :: self

      call self%zero_separated_flows()

    end subroutine source_network_node_default_separated_flows

!------------------------------------------------------------------------

    subroutine source_network_node_zero_separated_flows(self)
      !! Zeroes separated water and steam flows.

      class(source_network_node_type), intent(in out) :: self

      self%water_rate = 0._dp
      self%water_enthalpy = 0._dp
      self%steam_rate = 0._dp
      self%steam_enthalpy = 0._dp
      call self%separator%zero()

    end subroutine source_network_node_zero_separated_flows

!------------------------------------------------------------------------

    subroutine source_network_node_separate(self)
    !! Uses the node's separator to calculate separated water and
    !! steam flows.

    class(source_network_node_type), intent(in out) :: self

    call self%separator%separate(self%rate, self%enthalpy, &
         self%water_rate, self%water_enthalpy, &
         self%steam_rate, self%steam_enthalpy)

  end subroutine source_network_node_separate

!------------------------------------------------------------------------

  subroutine source_network_node_set_rate(self, rate, relax)
    !! Sets network node flow rate to specified value, with given
    !! relaxation coefficient.

    class(source_network_node_type), intent(in out) :: self
    PetscReal, intent(in) :: rate !! Flow rate
    PetscReal, intent(in), optional :: relax !! Relaxation coefficient

    if (present(relax)) then
       self%rate = (1._dp - relax) * self%rate + relax * rate
    else
       self%rate = rate
    end if
    call self%get_separated_flows()

  end subroutine source_network_node_set_rate

!------------------------------------------------------------------------

  subroutine source_network_node_scale_rate(self, scale)
    !! Scales network node flow rate by specified scale factor.

    class(source_network_node_type), intent(in out) :: self
    PetscReal, intent(in) :: scale !! Flow rate scale factor

    call self%set_rate(self%rate * scale)

  end subroutine source_network_node_scale_rate

!------------------------------------------------------------------------

  PetscReal function source_network_node_get_rate_by_type(self, flow_type) &
       result(rate)
    !! Gets rate corresponding to the specified flow type (total,
    !! water or steam).

    class(source_network_node_type), intent(in) :: self
    PetscInt, intent(in) :: flow_type

    select case (flow_type)
    case (SEPARATED_FLOW_TYPE_WATER)
       rate = self%water_rate
    case (SEPARATED_FLOW_TYPE_STEAM)
       rate = self%steam_rate
    case default
       rate = self%rate
    end select

  end function source_network_node_get_rate_by_type

!------------------------------------------------------------------------

  PetscBool function source_network_node_is_over(self, flow_type, limit) &
       result(over)
    !! Returns true if any of the specified flows are over their
    !! respective limit.

    class(source_network_node_type), intent(in) :: self
    PetscInt, intent(in) :: flow_type(:) !! Flow types
    PetscReal, intent(in) :: limit(:) !! Flow rate limits
    ! Locals:
    PetscInt :: i
    PetscReal :: rate

    over = PETSC_FALSE
    do i = 1, size(limit)
       rate = self%get_rate_by_type(flow_type(i))
       if (abs(rate) > limit(i)) then
          over = PETSC_TRUE
          exit
       end if
    end do

  end function source_network_node_is_over

!------------------------------------------------------------------------

  subroutine source_network_node_get_limit_scale(self, rate, limit, &
       over, scale)
    !! If rate is over limit, over returns true and scale is reduced
    !! (if necessary) to the value required to reduce the rate to the
    !! limit. Otherwise, over and scale are unchanged.

    class(source_network_node_type), intent(in out) :: self
    PetscReal, intent(in) :: rate !! Rate to limit
    PetscReal, intent(in) :: limit !! Rate limit
    PetscBool, intent(in out) :: over !! Whether rate is over limit
    PetscReal, intent(in out) :: scale !! Scale factor to reduce rate to limit
    ! Locals:
    PetscReal :: abs_rate
    PetscReal, parameter :: small = 1.e-6_dp

    abs_rate = abs(rate)
    if (abs_rate > limit) then
       over = PETSC_TRUE
       if (abs_rate > small) then
          scale = min(scale, limit / abs_rate)
       end if
    end if

  end subroutine source_network_node_get_limit_scale

!------------------------------------------------------------------------

  subroutine source_network_node_get_minimum_limit_scale(self, &
       flow_type, limit, over, scale)
    !! If any rates are over limit, over returns true and scale is the
    !! minimum value required to reduce all rates below limit.

    class(source_network_node_type), intent(in out) :: self
    PetscInt, intent(in) :: flow_type(:) !! Flow types
    PetscReal, intent(in) :: limit(:) !! Flow rate limits
    PetscBool, intent(out) :: over !! Whether any rate is over limit
    PetscReal, intent(out) :: scale !! Scale factor to reduce all rates to limits
    ! Locals:
    PetscInt :: i
    PetscReal :: rate

    over = PETSC_FALSE
    scale = 1._dp
    do i = 1, size(limit)
       rate = self%get_rate_by_type(flow_type(i))
       call self%get_limit_scale(rate, limit(i), over, scale)
    end do

  end subroutine source_network_node_get_minimum_limit_scale

!------------------------------------------------------------------------

  subroutine source_network_node_limit_rate(self, flow_type, limit)
    !! Limits network node flow rates (total, water or steam as
    !! specified by flow_type) to specified limits.

    class(source_network_node_type), intent(in out) :: self
    PetscInt, intent(in) :: flow_type(:) !! Flow types
    PetscReal, intent(in) :: limit(:) !! Flow rate limits
    ! Locals:
    PetscReal :: scale
    PetscBool :: over

    call self%get_minimum_limit_scale(flow_type, limit, over, scale)
    if (over) then
       call self%scale_rate(scale)
    end if

  end subroutine source_network_node_limit_rate

!------------------------------------------------------------------------

  subroutine source_network_node_limit_inputs(self, flow_type, limit)
    !! Limits network node input flow rates (total, water or steam as
    !! specified by flow_type) to specified limits.

    class(source_network_node_type), intent(in out) :: self
    PetscInt, intent(in) :: flow_type(:) !! Flow types
    PetscReal, intent(in) :: limit(:) !! Flow rate limits

    call self%limit_rate(flow_type, limit)

  end subroutine source_network_node_limit_inputs

!------------------------------------------------------------------------

  subroutine source_network_node_add_flows(self, mass_flow, energy_flow)
    !! Adds mass and energy flows from the node to the specified totals.

    class(source_network_node_type), intent(in out) :: self
    PetscReal, intent(in out) :: mass_flow !! Total mass flow rate
    PetscReal, intent(in out) :: energy_flow !! Total energy flow rate

    mass_flow = mass_flow + self%rate
    energy_flow = energy_flow + self%rate * self%enthalpy

  end subroutine source_network_node_add_flows

!------------------------------------------------------------------------

  subroutine source_network_node_add_separated_flows(self, water_mass_flow, &
       water_energy_flow, steam_mass_flow, steam_energy_flow)
    !! Adds mass and energy flows for separated water and steam from
    !! the node to the specified totals.

    class(source_network_node_type), intent(in out) :: self
    PetscReal, intent(in out) :: water_mass_flow !! Total water mass flow rate
    PetscReal, intent(in out) :: water_energy_flow !! Total water energy flow rate
    PetscReal, intent(in out) :: steam_mass_flow !! Total steam mass flow rate
    PetscReal, intent(in out) :: steam_energy_flow !! Total steam energy flow rate

    water_mass_flow = water_mass_flow + self%water_rate
    water_energy_flow = water_energy_flow + self%water_rate * self%water_enthalpy
    steam_mass_flow = steam_mass_flow + self%steam_rate
    steam_energy_flow = steam_energy_flow + self%steam_rate * self%steam_enthalpy

  end subroutine source_network_node_add_separated_flows

!------------------------------------------------------------------------

  subroutine source_network_node_destroy(self)
    !! Destroys source network node.

    class(source_network_node_type), intent(in out) :: self

    self%rate => null()
    self%enthalpy => null()
    self%water_rate => null()
    self%water_enthalpy => null()
    self%steam_rate => null()
    self%steam_enthalpy => null()

  end subroutine source_network_node_destroy

!------------------------------------------------------------------------

end module source_network_node_module
