!   Copyright 2022 University of Auckland.

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

module source_network_reinjector_module
  !! Module for source network reinjectors.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use list_module
  use source_network_node_module
  use source_network_group_module
  use source_module
  use separator_module, only: SEPARATED_FLOW_TYPE_WATER, &
       SEPARATED_FLOW_TYPE_STEAM
  use hdf5io_module, only: max_field_name_length

  implicit none
  private

  PetscInt, parameter, public :: num_source_network_reinjector_variables = &
       num_source_network_node_variables + 5
  PetscInt, parameter, public :: max_source_network_reinjector_variable_name_length = 24
  character(max_source_network_reinjector_variable_name_length), parameter, public :: &
       source_network_reinjector_variable_names(num_source_network_reinjector_variables) = [ &
       source_network_variable_names,  ["reinjector_index        ", &
       "overflow_water_rate     ", "overflow_water_enthalpy ", &
       "overflow_steam_rate     ", "overflow_steam_enthalpy " ]]
  PetscInt, parameter, public :: num_source_network_reinjector_constant_integer_variables = 1
  character(max_source_network_reinjector_variable_name_length), parameter, public :: &
       source_network_reinjector_constant_integer_variables( &
       num_source_network_reinjector_constant_integer_variables) = [ &
       "reinjector_index        "]
  character(max_field_name_length), parameter, public :: &
       required_output_source_network_reinjector_fields(0) = [&
       character(max_field_name_length)::]
  character(max_field_name_length), parameter, public :: &
       default_output_source_network_reinjector_fields(2) = [&
       "overflow_water_rate", "overflow_steam_rate"]

  type :: reinjector_output_type
     !! Type for reinjector outputs, distributing part of the
     !! reinjector flow to another network node.
     private
     class(source_network_reinjector_type), pointer :: reinjector !! Reinjector to output from
     class(source_network_node_type), pointer, public :: out !! Output network node
   contains
     private
     procedure, public :: node_limit => reinjector_output_node_limit
     procedure, public :: update => reinjector_output_update
     procedure, public :: destroy => reinjector_output_destroy
  end type reinjector_output_type

  type, public, extends(reinjector_output_type) :: specified_reinjector_output_type
     !! Type for reinjector outputs which have specified output
     !! rates/enthalpies, and output a particular flow type
     !! (e.g. water or steam).
     private
     PetscInt, public :: flow_type !! Type of output flow (water or steam)
     ! TODO: could make enthalpy time-dependent?
     PetscReal, public :: enthalpy !! Specified enthalpy for output (-1 to use input enthalpy)
   contains
     private
     procedure, public :: init => specified_reinjector_output_init
     procedure, public :: specified_enthalpies => specified_reinjector_output_specified_enthalpies
     procedure, public :: default_enthalpies => specified_reinjector_output_default_enthalpies
     procedure, public :: rates => specified_reinjector_output_rates
     procedure, public :: get_specified_flows => specified_reinjector_output_get_specified_flows
     procedure, public :: subtract_flow => specified_reinjector_output_subtract_flow
  end type specified_reinjector_output_type

  type, public, extends(specified_reinjector_output_type) :: rate_reinjector_output_type
     !! Type for reinjector output which has a specified output rate.
     private
     ! TODO: make this an interpolation table for time-dependent rate
     PetscReal, public :: rate !! Specified output rate
   contains
     private
     procedure, public :: rates => rate_reinjector_output_rates
  end type rate_reinjector_output_type

  type, public, extends(specified_reinjector_output_type) :: proportion_reinjector_output_type
     !! Type for reinjector output with an output rate which is a
     !! specified proportion of the input rate.
     private
     ! TODO: make this an interpolation table for time-dependent proportion
     PetscReal, public :: proportion !! Specified output rate
   contains
     private
     procedure, public :: rates => proportion_reinjector_output_rates
  end type proportion_reinjector_output_type

  type, public, extends(reinjector_output_type) :: overflow_reinjector_output_type
     !! Type for overflow reinjector output, which stores and can
     !! display mass flow balances left over after specified outputs
     !! have been reinjected, and optionally direct them to another
     !! network node.
     private
     PetscReal, pointer, public :: water_rate !! Separated water mass flow rate
     PetscReal, pointer, public :: water_enthalpy !! Separated water enthalpy
     PetscReal, pointer, public :: steam_rate !! Separated steam mass flow rate
     PetscReal, pointer, public :: steam_enthalpy !! Separated steam enthalpy
   contains
     private
     procedure, public :: init => overflow_reinjector_output_init
     procedure, public :: set_flows => overflow_reinjector_output_set_flows
  end type overflow_reinjector_output_type

  type, public, extends(source_network_node_type) :: source_network_reinjector_type
     !! Type for reinjectors, taking a single input and distributing
     !! it to one or more reinjection sources or other reinjectors.
     private
     class(source_network_node_type), pointer, public :: in !! Input node for reinjector
     type(list_type), public :: out !! List of outputs for reinjector
     type(overflow_reinjector_output_type), public :: overflow !! Output for overflow
     MPI_Comm :: comm !! MPI communicator for reinjector
     PetscMPIInt, public :: rank !! Rank of reinjector in its own communicator
     PetscMPIInt :: root_world_rank !! Rank in world communicator of reinjector root rank
     PetscInt, public :: local_reinjector_index !! Index of reinjector in local part of reinjector vector (-1 if not a root reinjector)
     PetscReal, pointer, public :: reinjector_index !! Index of reinjector in input
   contains
     private
     procedure, public :: init => source_network_reinjector_init
     procedure, public :: init_comm => source_network_reinjector_init_comm
     procedure, public :: assign => source_network_reinjector_assign
     procedure, public :: init_data => source_network_reinjector_init_data
     procedure, public :: balance_limit => source_network_reinjector_balance_limit
     procedure, public :: distribute => source_network_reinjector_distribute
     procedure, public :: destroy => source_network_reinjector_destroy
  end type source_network_reinjector_type

contains

!------------------------------------------------------------------------

  subroutine reinjector_output_list_node_data_destroy(node)
    ! Destroys reinjector output in a list node.

    use list_module, only: list_node_type

    type(list_node_type), pointer, intent(in out) :: node

    select type (output => node%data)
    class is (reinjector_output_type)
       call output%destroy()
    end select

  end subroutine reinjector_output_list_node_data_destroy

!------------------------------------------------------------------------
! Reinjector output type
!------------------------------------------------------------------------

  subroutine reinjector_output_node_limit(self, rate)
    !! Limits injection rate to what can be injected into the output
    !! node.

    class(reinjector_output_type), intent(in out) :: self
    PetscReal, intent(in out) :: rate

    if (associated(self%out)) then

       select type (n => self%out)
       class is (source_network_node_type)
          if (n%rate > 0._dp) then
             if (rate > 0._dp) then
                ! rates specified in both reinjector output and node:
                rate = min(rate, n%rate)
             else
                ! rate specified in node only:
                rate = n%rate
             end if
          else
             ! check if no rate specified in reinjector output or node:
             rate = max(rate, 0._dp)
          end if
       end select
    end if

  end subroutine reinjector_output_node_limit

!------------------------------------------------------------------------

  subroutine reinjector_output_update(self, water_rate, &
       water_enthalpy, steam_rate, steam_enthalpy)
    !! Assigns rates and enthalpies to the output node.

    class(reinjector_output_type), intent(in out) :: self
    PetscReal, target, intent(in) :: water_rate, water_enthalpy
    PetscReal, target, intent(in) :: steam_rate, steam_enthalpy
    ! Locals:
    PetscReal :: rate, enthalpy
    PetscReal, parameter :: small = 1.e-6_dp

    ! TODO: add recursive case of class is (source_network_reinjector_type)

    if (associated(self%out)) then

       rate = water_rate + steam_rate
       if (rate > small) then
          enthalpy = (water_rate * water_enthalpy + &
               steam_rate * steam_enthalpy) / rate
       else
          enthalpy = 0._dp
       end if

       select type (n => self%out)
       class is (source_type)
          call n%set_rate(rate)
          n%injection_enthalpy = enthalpy
       end select

       self%out%water_rate = water_rate
       self%out%water_enthalpy = water_enthalpy
       self%out%steam_rate = steam_rate
       self%out%steam_enthalpy = steam_enthalpy

    end if

  end subroutine reinjector_output_update

!------------------------------------------------------------------------

  subroutine reinjector_output_destroy(self)
    !! Destroys a reinjector output.

    class(reinjector_output_type), intent(in out) :: self

    self%reinjector => null()
    self%out => null()

  end subroutine reinjector_output_destroy

!------------------------------------------------------------------------
! Specified reinjector output type
!------------------------------------------------------------------------

  subroutine specified_reinjector_output_init(self, reinjector, &
       flow_type, enthalpy)
    !! Initialises specified reinjector output.

    class(specified_reinjector_output_type), intent(in out) :: self
    type(source_network_reinjector_type), target, intent(in) :: reinjector
    PetscInt, intent(in) :: flow_type
    PetscReal, intent(in) :: enthalpy

    self%reinjector => reinjector
    self%flow_type = flow_type
    self%enthalpy = enthalpy

  end subroutine specified_reinjector_output_init

!------------------------------------------------------------------------

  subroutine specified_reinjector_output_specified_enthalpies(self, &
       water_enthalpy, steam_enthalpy)
    !! Returns enthalpies, taken from specified enthalpy, depending on
    !! flow type.

    class(specified_reinjector_output_type), intent(in out) :: self
    PetscReal, intent(out) :: water_enthalpy, steam_enthalpy

    select case (self%flow_type)
    case (SEPARATED_FLOW_TYPE_WATER)
       water_enthalpy = self%enthalpy
       steam_enthalpy = 0._dp
    case (SEPARATED_FLOW_TYPE_STEAM)
       water_enthalpy = 0._dp
       steam_enthalpy = self%enthalpy
    end select

  end subroutine specified_reinjector_output_specified_enthalpies

!------------------------------------------------------------------------

  subroutine specified_reinjector_output_default_enthalpies(self, &
       water_enthalpy, steam_enthalpy)
    !! Returns default enthalpies, taken from input, depending on flow type.

    class(specified_reinjector_output_type), intent(in out) :: self
    PetscReal, intent(out) :: water_enthalpy, steam_enthalpy

    select case (self%flow_type)
    case (SEPARATED_FLOW_TYPE_WATER)
       water_enthalpy = self%reinjector%in%water_enthalpy
       steam_enthalpy = 0._dp
    case (SEPARATED_FLOW_TYPE_STEAM)
       water_enthalpy = 0._dp
       steam_enthalpy = self%reinjector%in%steam_enthalpy
    end select

  end subroutine specified_reinjector_output_default_enthalpies

!------------------------------------------------------------------------

  subroutine specified_reinjector_output_rates(self, water_rate, steam_rate)
    !! Gets rates for specified reinjector output. Derived types
    !! override this function.

    class(specified_reinjector_output_type), intent(in out) :: self
    PetscReal, intent(out) :: water_rate, steam_rate

    water_rate = 0._dp
    steam_rate = 0._dp

  end subroutine specified_reinjector_output_rates

!------------------------------------------------------------------------

  subroutine specified_reinjector_output_get_specified_flows(self, water_rate, &
       water_enthalpy, steam_rate, steam_enthalpy)
    !! Gets specified rates and enthalpies for specified reinjector output.

    class(specified_reinjector_output_type), intent(in out) :: self
    PetscReal, intent(out) :: water_rate, water_enthalpy
    PetscReal, intent(out) :: steam_rate, steam_enthalpy

    call self%rates(water_rate, steam_rate)

    if (self%enthalpy > 0._dp) then
       call self%specified_enthalpies(water_enthalpy, steam_enthalpy)
    else
       call self%default_enthalpies(water_enthalpy, steam_enthalpy)
    end if

  end subroutine specified_reinjector_output_get_specified_flows

!------------------------------------------------------------------------

  subroutine specified_reinjector_output_subtract_flow(self, water_balance, &
       steam_balance)
    !! Subtracts output phase flows from water and steam balances,
    !! according to the output flow type.

    class(specified_reinjector_output_type), intent(in) :: self
    PetscReal, intent(in out) :: water_balance !! Water mass flow balance
    PetscReal, intent(in out) :: steam_balance !! Steam mass flow balance

    if (associated(self%out)) then
       select type (n => self%out)
       class is (source_network_node_type)
          water_balance = water_balance - n%water_rate
          steam_balance = steam_balance - n%steam_rate
       end select
    end if

  end subroutine specified_reinjector_output_subtract_flow

!------------------------------------------------------------------------
! Rate reinjector output type
!------------------------------------------------------------------------

  subroutine rate_reinjector_output_rates(self, water_rate, steam_rate)
    !! Gets rates for rate reinjector output.

    class(rate_reinjector_output_type), intent(in out) :: self
    PetscReal, intent(out) :: water_rate, steam_rate

    select case (self%flow_type)
    case (SEPARATED_FLOW_TYPE_WATER)
       water_rate = self%rate
       steam_rate = 0._dp
    case (SEPARATED_FLOW_TYPE_STEAM)
       water_rate = 0._dp
       steam_rate = self%rate
    end select

  end subroutine rate_reinjector_output_rates

!------------------------------------------------------------------------
! Proportion reinjector output type
!------------------------------------------------------------------------

  subroutine proportion_reinjector_output_rates(self, water_rate, steam_rate)
    !! Gets rates for proportion reinjector output.

    class(proportion_reinjector_output_type), intent(in out) :: self
    PetscReal, intent(out) :: water_rate, steam_rate

    select case (self%flow_type)
    case (SEPARATED_FLOW_TYPE_WATER)
       water_rate = self%proportion * abs(self%reinjector%in%water_rate)
       steam_rate = 0._dp
    case (SEPARATED_FLOW_TYPE_STEAM)
       water_rate = 0._dp
       steam_rate = self%proportion * abs(self%reinjector%in%steam_rate)
    end select

  end subroutine proportion_reinjector_output_rates

!------------------------------------------------------------------------
! Overflow reinjector output type
!------------------------------------------------------------------------

  subroutine overflow_reinjector_output_init(self, reinjector)
    !! Initialises overflow reinjector output.

    class(overflow_reinjector_output_type), intent(in out) :: self
    type(source_network_reinjector_type), target, intent(in) :: reinjector

    self%reinjector => reinjector

  end subroutine overflow_reinjector_output_init

!------------------------------------------------------------------------

  subroutine overflow_reinjector_output_set_flows(self, water_rate, &
       water_enthalpy, steam_rate, steam_enthalpy)
    !! Sets overflow output flow rates and enthalpies.

    class(overflow_reinjector_output_type), intent(in out) :: self
    PetscReal, intent(in) :: water_rate !! Separated water mass flow rate
    PetscReal, intent(in) :: water_enthalpy !! Separated water enthalpy
    PetscReal, intent(in) :: steam_rate !! Separated steam mass flow rate
    PetscReal, intent(in) :: steam_enthalpy !! Separated steam enthalpy

    self%water_rate = water_rate
    self%water_enthalpy = water_enthalpy
    self%steam_rate = steam_rate
    self%steam_enthalpy = steam_enthalpy

  end subroutine overflow_reinjector_output_set_flows

!------------------------------------------------------------------------
! Reinjector type
!------------------------------------------------------------------------

  subroutine source_network_reinjector_init(self, name)
    !! Initialises a source network reinjector. Only variables stored in
    !! the object itself are initialised, not those stored in the
    !! reinjector data vector and accessed via pointers.

    class(source_network_reinjector_type), intent(in out) :: self
    character(*), intent(in) :: name !! Reinjector name

    self%name = name

    self%in => null()
    call self%out%init(owner = PETSC_TRUE)
    self%overflow%out => null()

  end subroutine source_network_reinjector_init

!------------------------------------------------------------------------

  subroutine source_network_reinjector_init_comm(self)
    !! Initialises MPI communicator, rank and root world rank for the
    !! reinjector. It is assumed that the output node list has already
    !! been populated and overflow (if any) has been assigned.

    use mpi_utils_module, only: mpi_comm_root_world_rank

    class(source_network_reinjector_type), intent(in out) :: self
    ! Locals:
    PetscInt :: colour
    PetscErrorCode :: ierr

    colour = MPI_UNDEFINED
    call self%out%traverse(reinjector_comm_iterator)
    call overflow_comm()

    call MPI_comm_split(PETSC_COMM_WORLD, colour, 0, self%comm, ierr)
    if (self%comm /= MPI_COMM_NULL) then
       call MPI_comm_rank(self%comm, self%rank, ierr)
    else
       self%rank = -1
    end if
    self%root_world_rank = mpi_comm_root_world_rank(self%comm)

  contains

    subroutine reinjector_comm_iterator(node, stopped)
      !! Stops and sets colour if there are any local reinjection
      !! sources, or reinjectors with local reinjection sources, in
      !! the reinjector output node list.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (output => node%data)
      class is (reinjector_output_type)
         select type (n => output%out)
         type is (source_type)
            colour = 1
            stopped = PETSC_TRUE
         class is (source_network_reinjector_type)
            if (n%rank >= 0) then
               colour = 1
               stopped = PETSC_TRUE
            end if
         end select
      end select

    end subroutine reinjector_comm_iterator

!........................................................................

    subroutine overflow_comm()
      !! Sets colour if overflow node is assigned to a local source or
      !! a reinjector with a local source as output or overflow.

      if (associated(self%overflow%out)) then
         select type (node => self%overflow%out)
         type is (source_type)
            colour = 1
         class is (source_network_reinjector_type)
            if (node%rank >= 0) then
               colour = 1
            end if
         end select
      end if

    end subroutine overflow_comm

  end subroutine source_network_reinjector_init_comm

!------------------------------------------------------------------------

  subroutine source_network_reinjector_assign(self, data, offset)
    !! Assigns pointers in source network reinjector object to elements in
    !! the data array, starting from the specified offset.

    class(source_network_reinjector_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! reinjector data array
    PetscInt, intent(in) :: offset  !! source array offset
    ! Locals:
    PetscInt :: reinjector_offset

    call self%source_network_node_type%assign(data, offset)

    reinjector_offset = offset + num_source_network_node_variables

    self%reinjector_index => data(reinjector_offset)

    self%overflow%water_rate => data(reinjector_offset + 1)
    self%overflow%water_enthalpy => data(reinjector_offset + 2)
    self%overflow%steam_rate => data(reinjector_offset + 3)
    self%overflow%steam_enthalpy => data(reinjector_offset + 4)

  end subroutine source_network_reinjector_assign

!------------------------------------------------------------------------

  subroutine source_network_reinjector_init_data(self, reinjector_index)
    !! Initialised source network reinjector variables accessed via
    !! pointers to the reinjector vector. The reinjector assign() method
    !! must be called first.

    class(source_network_reinjector_type), intent(in out) :: self
    PetscInt, intent(in) :: reinjector_index !! Index of reinjector in input

    self%reinjector_index = dble(reinjector_index)
    self%overflow%water_rate = 0._dp
    self%overflow%water_enthalpy = 0._dp
    self%overflow%steam_rate = 0._dp
    self%overflow%steam_enthalpy = 0._dp

  end subroutine source_network_reinjector_init_data

!------------------------------------------------------------------------

  subroutine source_network_reinjector_balance_limit(self, water_balance, &
       steam_balance, water_rate, steam_rate)
    !! Limits water and steam rates so they do not exceed the current
    !! remaining water and steam balances.

    class(source_network_reinjector_type), intent(in) :: self
    PetscReal, intent(in) :: water_balance !! Current reinjector water balance
    PetscReal, intent(in) :: steam_balance !! Current reinjector steam balance
    PetscReal, intent(in out) :: water_rate !! Water rate for next output
    PetscReal, intent(in out) :: steam_rate !! Steam rate for next output

    water_rate = min(water_rate, water_balance)
    steam_rate = min(steam_rate, steam_balance)

  end subroutine source_network_reinjector_balance_limit

!------------------------------------------------------------------------

  subroutine source_network_reinjector_distribute(self)
    !! Distributes reinjector input flow to outputs (and overflow if
    !! needed).

    class(source_network_reinjector_type), intent(in out) :: self
    ! Locals:
    PetscReal :: water_balance, steam_balance

    if (associated(self%in)) then
       water_balance = abs(self%in%water_rate)
       steam_balance = abs(self%in%steam_rate)
    else
       water_balance = 0._dp
       steam_balance = 0._dp
    end if

    call self%out%traverse(reinjector_distribute_iterator)

    call self%overflow%set_flows(water_balance, self%in%water_enthalpy, &
         steam_balance, self%in%steam_enthalpy)
    if (associated(self%overflow%out)) then
       call self%overflow%update(water_balance, self%in%water_enthalpy, &
            steam_balance, self%in%steam_enthalpy)
    end if

  contains

    subroutine reinjector_distribute_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscReal, target :: water_rate, steam_rate
      PetscReal :: water_enthalpy, steam_enthalpy
      PetscReal, pointer :: rate

      stopped = PETSC_FALSE
      select type (output => node%data)
      class is (specified_reinjector_output_type)
         call output%get_specified_flows(water_rate, water_enthalpy, &
              steam_rate, steam_enthalpy)
         select case (output%flow_type)
         case (SEPARATED_FLOW_TYPE_WATER)
            rate => water_rate
         case (SEPARATED_FLOW_TYPE_STEAM)
            rate => steam_rate
         end select
         call output%node_limit(rate)
         call self%balance_limit(water_balance, steam_balance, &
              water_rate, steam_rate)
         call output%update(water_rate, water_enthalpy, steam_rate, &
              steam_enthalpy)
         call output%subtract_flow(water_balance, steam_balance)
      end select

      ! TODO use MPI_send()/recv() to pass output flow rate to self%in%comm root?

    end subroutine reinjector_distribute_iterator

  end subroutine source_network_reinjector_distribute

!------------------------------------------------------------------------

  subroutine source_network_reinjector_destroy(self)
    !! Destroys a source network reinjector.

    class(source_network_reinjector_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    if (associated(self%in)) self%in => null()
    call self%out%destroy(reinjector_output_list_node_data_destroy)
    call MPI_comm_free(self%comm, ierr)

    call self%source_network_node_type%destroy()

  end subroutine source_network_reinjector_destroy

!------------------------------------------------------------------------

end module source_network_reinjector_module
