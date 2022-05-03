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

  implicit none
  private

  type, public :: reinjector_output_type
     !! Type for reinjector outputs, distributing part of the
     !! reinjector flow to another network node.
     private
     PetscInt, public :: flow_type !! Type of output flow (water or steam)
     ! TODO: could make enthalpy time-dependent?
     PetscReal, public :: enthalpy !! Specified enthalpy for output (-1 to use input enthalpy)
     class(source_network_node_type), pointer :: in !! Input network node
     class(source_network_node_type), pointer, public :: out !! Output network node
   contains
     private
     procedure, public :: init => reinjector_output_init
     procedure, public :: get_rate => reinjector_output_get_rate
     procedure, public :: update => reinjector_output_update
     procedure, public :: subtract_flow => reinjector_output_subtract_flow
     procedure, public :: destroy => reinjector_output_destroy
  end type reinjector_output_type

  type, public, extends(reinjector_output_type) :: rate_reinjector_output_type
     !! Type for reinjector output which has a specified output rate.
     private
     ! TODO: make this an interpolation table for time-dependent rate
     PetscReal, public :: rate !! Specified output rate
   contains
     procedure, public :: get_rate => rate_reinjector_output_get_rate
  end type rate_reinjector_output_type

  type, public, extends(reinjector_output_type) :: proportion_reinjector_output_type
     !! Type for reinjector output with an output rate which is a
     !! specified proportion of the input rate.
     private
     ! TODO: make this an interpolation table for time-dependent proportion
     PetscReal, public :: proportion !! Specified output rate
   contains
     procedure, public :: get_rate => proportion_reinjector_output_get_rate
  end type proportion_reinjector_output_type

  type, public, extends(source_network_node_type) :: source_network_reinjector_type
     !! Type for reinjectors, taking a single input and distributing
     !! it to one or more reinjection sources or other reinjectors.
     private
     class(source_network_node_type), pointer, public :: in !! Input node for reinjector
     type(list_type), public :: out !! List of outputs for reinjector
     type(reinjector_output_type), public :: overflow !! Output for overflow
     MPI_Comm :: comm !! MPI communicator for reinjector
     PetscMPIInt, public :: rank !! Rank of reinjector in its own communicator
     PetscMPIInt :: root_world_rank !! Rank in world communicator of reinjector root rank
   contains
     private
     procedure, public :: init => source_network_reinjector_init
     procedure, public :: init_comm => source_network_reinjector_init_comm
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

  subroutine reinjector_output_init(self, reinjector, flow_type, enthalpy)
    !! Initialises reinjector output.

    class(reinjector_output_type), intent(in out) :: self
    type(source_network_reinjector_type), target, intent(in) :: reinjector
    PetscInt, intent(in) :: flow_type
    PetscReal, intent(in) :: enthalpy

    self%in => reinjector
    self%flow_type = flow_type
    self%enthalpy = enthalpy

  end subroutine reinjector_output_init

!------------------------------------------------------------------------

  PetscReal function reinjector_output_get_rate(self) result(rate)
    !! Gets rate for reinjector output. Derived types override this
    !! function.

    class(reinjector_output_type), intent(in out) :: self

    rate = 0._dp

  end function reinjector_output_get_rate

!------------------------------------------------------------------------

  subroutine reinjector_output_update(self, rate)
    !! Assigns rate to the output node, and enthalpy from the input
    !! node.

    class(reinjector_output_type), intent(in out) :: self
    PetscReal, intent(in) :: rate !! Specified flow rate
    ! Locals:
    PetscReal :: effective_rate, effective_enthalpy

    ! TODO: add recursive case of class is (source_network_reinjector_type)

    select type (n => self%out)
    class is (source_type)

       if (n%rate > 0._dp) then
          if (rate > 0._dp) then
             effective_rate = min(rate, n%rate)
          else
             effective_rate = n%rate
          end if
       else
          if (rate > 0._dp) then
             effective_rate = rate
          else
             effective_rate = 0._dp
          end if
       end if
       call n%set_rate(effective_rate)

       if (self%enthalpy > 0._dp) then
          effective_enthalpy = self%enthalpy
       else
          select case (self%flow_type)
          case (SEPARATED_FLOW_TYPE_WATER)
             effective_enthalpy = self%in%water_enthalpy
          case (SEPARATED_FLOW_TYPE_STEAM)
             effective_enthalpy = self%in%steam_enthalpy
          end select
       end if
       n%injection_enthalpy = effective_enthalpy

       select case (self%flow_type)
       case (SEPARATED_FLOW_TYPE_WATER)
          n%water_rate = effective_rate
          n%steam_rate = 0._dp
          n%water_enthalpy = effective_enthalpy
          n%steam_enthalpy = 0._dp
       case (SEPARATED_FLOW_TYPE_STEAM)
          n%water_rate = 0._dp
          n%steam_rate = effective_rate
          n%water_enthalpy = 0._dp
          n%steam_enthalpy = effective_enthalpy
       end select

    end select

  end subroutine reinjector_output_update

!------------------------------------------------------------------------

  subroutine reinjector_output_subtract_flow(self, rate, water_rate, &
       steam_rate)
    !! Subtracts output flow from totals, according to the output flow
    !! type.

    class(reinjector_output_type), intent(in) :: self
    PetscReal, intent(in out) :: rate !! Total mass flow rate
    PetscReal, intent(in out) :: water_rate !! Water mass flow rate
    PetscReal, intent(in out) :: steam_rate !! Steam mass flow rate

    if (associated(self%out)) then
       select type (n => self%out)
       class is (source_network_node_type)
          rate = rate - n%rate
          water_rate = water_rate - n%water_rate
          steam_rate = steam_rate - n%steam_rate
       end select
    end if

  end subroutine reinjector_output_subtract_flow

!------------------------------------------------------------------------

  subroutine reinjector_output_destroy(self)
    !! Destroys a reinjector output.

    class(reinjector_output_type), intent(in out) :: self

    self%in => null()
    self%out => null()

  end subroutine reinjector_output_destroy

!------------------------------------------------------------------------
! Rate reinjector output type
!------------------------------------------------------------------------

  PetscReal function rate_reinjector_output_get_rate(self) result(rate)
    !! Gets rate for rate reinjector output.

    class(rate_reinjector_output_type), intent(in out) :: self

    rate = self%rate

  end function rate_reinjector_output_get_rate

!------------------------------------------------------------------------
! Proportion reinjector output type
!------------------------------------------------------------------------

  PetscReal function proportion_reinjector_output_get_rate(self) result(rate)
    !! Gets rate for proportion reinjector output.

    class(proportion_reinjector_output_type), intent(in out) :: self

    select case (self%flow_type)
    case (SEPARATED_FLOW_TYPE_WATER)
       rate = self%proportion * self%in%water_rate
    case (SEPARATED_FLOW_TYPE_STEAM)
       rate = self%proportion * self%in%steam_rate
    end select

  end function proportion_reinjector_output_get_rate

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

  subroutine source_network_reinjector_distribute(self)
    !! Distributes reinjector input flow to outputs (and overflow if
    !! needed).

    class(source_network_reinjector_type), intent(in out) :: self
    ! Locals:
    PetscReal :: rate_balance, water_balance, steam_balance

    if (associated(self%in)) then
       select type (n => self%in)
       ! Reverse flows from production side of network (sources and
       ! groups):
       type is (source_type)
          rate_balance = -n%rate
          water_balance = -n%water_rate
          steam_balance = -n%steam_rate
       type is (source_network_group_type)
          rate_balance = -n%rate
          water_balance = -n%water_rate
          steam_balance = -n%steam_rate
       type is (source_network_reinjector_type)
          rate_balance = n%rate
          water_balance = n%water_rate
          steam_balance = n%steam_rate
       end select
    else
       rate_balance = 0._dp
       water_balance = 0._dp
       steam_balance = 0._dp
    end if

    call self%out%traverse(reinjector_distribute_iterator)

    if (associated(self%overflow%out)) then
       select type (n => self%overflow%out)
       class is (source_network_node_type)
          call n%set_rate(rate_balance)
          n%water_rate = water_balance
          n%water_enthalpy = self%in%water_enthalpy
          n%steam_rate = steam_balance
          n%steam_enthalpy = self%in%steam_enthalpy
       end select
    end if

  contains

    subroutine reinjector_distribute_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscReal :: rate

      stopped = PETSC_FALSE
      select type (output => node%data)
      class is (reinjector_output_type)
         rate = output%get_rate()
         call output%update(rate)
         call output%subtract_flow(rate_balance, water_balance, steam_balance)
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
