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
       num_source_network_node_variables + 10
  PetscInt, parameter, public :: max_source_network_reinjector_variable_name_length = 24
  character(max_source_network_reinjector_variable_name_length), parameter, public :: &
       source_network_reinjector_variable_names(num_source_network_reinjector_variables) = [ &
       source_network_variable_names,  ["reinjector_index        ", &
       "output_rate             ", &
       "output_water_rate       ", "output_steam_rate       ", &
       "overflow_rate           ", "overflow_enthalpy       ", &
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
       default_output_source_network_reinjector_fields(4) = [&
       "output_water_rate  ", "output_steam_rate  ", &
       "overflow_water_rate", "overflow_steam_rate"]
  PetscReal, parameter, public :: default_reinjector_output_rate = -1._dp
  PetscReal, parameter, public :: default_reinjector_output_proportion = 0._dp
  PetscReal, parameter, public :: default_reinjector_output_enthalpy = -1._dp

  type, extends(source_network_node_type) :: reinjector_output_type
     !! Type for reinjector outputs, distributing part of the
     !! reinjector flow to another network node.
     private
     class(source_network_reinjector_type), pointer :: reinjector !! Reinjector to output from
     class(source_network_node_type), pointer, public :: out !! Output network node
   contains
     private
     procedure, public :: allocate_variables => reinjector_output_allocate_variables
     procedure, public :: deallocate_variables => reinjector_output_deallocate_variables
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
     PetscReal, public :: specified_enthalpy !! Specified enthalpy for output (-1 to use input enthalpy)
   contains
     private
     procedure, public :: init => specified_reinjector_output_init
     procedure, public :: specified_enthalpies => specified_reinjector_output_specified_enthalpies
     procedure, public :: default_enthalpies => specified_reinjector_output_default_enthalpies
     procedure, public :: rates => specified_reinjector_output_rates
     procedure, public :: enthalpies => specified_reinjector_output_enthalpies
     procedure, public :: node_limit => specified_reinjector_output_node_limit
     procedure, public :: destroy => specified_reinjector_destroy
  end type specified_reinjector_output_type

  type, public, extends(specified_reinjector_output_type) :: rate_reinjector_output_type
     !! Type for reinjector output which has a specified output rate.
     private
     PetscReal, public :: specified_rate !! Specified output rate
   contains
     private
     procedure, public :: rates => rate_reinjector_output_rates
  end type rate_reinjector_output_type

  type, public, extends(specified_reinjector_output_type) :: proportion_reinjector_output_type
     !! Type for reinjector output with an output rate which is a
     !! specified proportion of the input rate.
     private
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
     PetscMPIInt :: out_world_rank !! Rank in world communicator of output node
   contains
     private
     procedure, public :: init => overflow_reinjector_output_init
     procedure, public :: get_world_rank => overflow_reinjector_output_get_world_rank
     procedure, public :: allocate_variables => overflow_reinjector_output_allocate_variables
     procedure, public :: deallocate_variables => overflow_reinjector_output_deallocate_variables
     procedure, public :: assign => overflow_reinjector_output_assign
     procedure, public :: set_flows => overflow_reinjector_output_set_flows
     procedure, public :: node_limit => overflow_reinjector_output_node_limit
     procedure, public :: destroy => overflow_reinjector_output_destroy
  end type overflow_reinjector_output_type

  type, public, extends(source_network_node_type) :: source_network_reinjector_type
     !! Type for reinjectors, taking a single input and distributing
     !! it to one or more reinjection sources or other reinjectors.
     private
     class(source_network_node_type), pointer, public :: in !! Input node for reinjector
     type(list_type), public :: out !! List of outputs for reinjector
     type(overflow_reinjector_output_type), pointer, public :: overflow !! Output for overflow
     MPI_Comm :: comm !! MPI communicator for reinjector
     MPI_Comm :: in_comm !! MPI communicator for reinjector, including input node process
     PetscMPIInt, public :: rank !! Rank of reinjector in its own communicator
     PetscMPIInt, public :: in_comm_rank !! Rank of reinjector in self%in_comm
     PetscMPIInt, public :: in_comm_input_rank !! Rank of input node process in self%in_comm
     PetscMPIInt :: root_world_rank !! Rank in world communicator of reinjector root rank
     PetscInt, public :: local_reinjector_index !! Index of reinjector in local part of reinjector vector (-1 if not a root reinjector)
     PetscReal, pointer, public :: reinjector_index !! Index of reinjector in input
     PetscReal, pointer, public :: output_rate !! Sum of output rates
     PetscReal, pointer, public :: output_water_rate !! Sum of water output rates
     PetscReal, pointer, public :: output_steam_rate !! Sum of steam output rates
     PetscReal :: in_water_rate !! Water rate in input node
     PetscReal :: in_steam_rate !! Steam rate in input node
     PetscReal :: in_water_enthalpy !! Water enthalpy in input node
     PetscReal :: in_steam_enthalpy !! Steam enthalpy in input node
     PetscInt :: local_gather_count !! How many local outputs are included in gather operations
     PetscInt, allocatable :: gather_counts(:) !! Process counts for gather operations
     PetscInt :: gather_count !! Total count for gather operations (only computed on root rank)
     PetscInt, allocatable :: gather_displacements(:) !! Process displacements for gather operations
     PetscInt, allocatable :: gather_index(:) !! Sort index for gather operations
   contains
     private
     procedure, public :: init => source_network_reinjector_init
     procedure, public :: init_comm => source_network_reinjector_init_comm
     procedure, public :: init_in_comm => source_network_reinjector_init_in_comm
     procedure, public :: assign => source_network_reinjector_assign
     procedure, public :: init_data => source_network_reinjector_init_data
     procedure, public :: overflow_output => source_network_reinjector_overflow_output
     procedure, public :: capacity => source_network_reinjector_capacity
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

  subroutine node_limit_rate(node_rate, rate)
    !! Limits specified rate according to node_rate.

    PetscReal, intent(in) :: node_rate
    PetscReal, intent(in out) :: rate

    if (node_rate > -1._dp) then
       if (rate > -1._dp) then
          ! rates specified in both reinjector output and node:
          rate = min(rate, node_rate)
       else
          ! rate specified in node only:
          rate = node_rate
       end if
    end if

  end subroutine node_limit_rate

!------------------------------------------------------------------------
! Reinjector output type
!------------------------------------------------------------------------

  subroutine reinjector_output_allocate_variables(self)
    !! Allocates overflow rate, enthalpy etc. variables.

    class(reinjector_output_type), intent(in out) :: self

    allocate(self%rate, self%enthalpy)
    allocate(self%water_rate, self%water_enthalpy)
    allocate(self%steam_rate, self%steam_enthalpy)

  end subroutine reinjector_output_allocate_variables

!------------------------------------------------------------------------

  subroutine reinjector_output_deallocate_variables(self)
    !! Deallocates overflow rate, enthalpy etc. variables.

    class(reinjector_output_type), intent(in out) :: self

    deallocate(self%rate, self%enthalpy)
    deallocate(self%water_rate, self%water_enthalpy)
    deallocate(self%steam_rate, self%steam_enthalpy)

  end subroutine reinjector_output_deallocate_variables

!------------------------------------------------------------------------

  subroutine reinjector_output_node_limit(self, rate)
    !! Limits injection rate to what can be injected into the output
    !! node. Derived types override this routine.

    class(reinjector_output_type), intent(in out) :: self
    PetscReal, intent(in out) :: rate

    continue

  end subroutine reinjector_output_node_limit

!------------------------------------------------------------------------

  subroutine reinjector_output_update(self, water_rate, &
       water_enthalpy, steam_rate, steam_enthalpy, enthalpy_specified)
    !! Updates output rates and enthalpies, and assigns them to the
    !! output node if it is a source.

    class(reinjector_output_type), intent(in out) :: self
    PetscReal, target, intent(in) :: water_rate, water_enthalpy
    PetscReal, target, intent(in) :: steam_rate, steam_enthalpy
    PetscBool, intent(in) :: enthalpy_specified
    ! Locals:
    PetscReal :: rate, enthalpy
    PetscReal, parameter :: small = 1.e-6_dp

    if (associated(self%out)) then

       rate = water_rate + steam_rate
       if (rate > small) then
          enthalpy = (water_rate * water_enthalpy + &
               steam_rate * steam_enthalpy) / rate
       else
          enthalpy = 0._dp
       end if

       self%rate = rate
       self%enthalpy = enthalpy
       self%water_rate = water_rate
       self%water_enthalpy = water_enthalpy
       self%steam_rate = steam_rate
       self%steam_enthalpy = steam_enthalpy

       select type (n => self%out)
       class is (source_type)
          n%rate = rate
          n%water_rate = water_rate
          n%steam_rate = steam_rate
          if (.not. n%enthalpy_specified) then
             n%injection_enthalpy = enthalpy
             n%water_enthalpy = water_enthalpy
             n%steam_enthalpy = steam_enthalpy
          end if
       end select

    end if

  end subroutine reinjector_output_update

!------------------------------------------------------------------------

  subroutine reinjector_output_destroy(self)
    !! Destroys a reinjector output.

    class(reinjector_output_type), intent(in out) :: self

    call self%deallocate_variables()
    call self%source_network_node_type%destroy()

    self%reinjector => null()
    self%out => null()

  end subroutine reinjector_output_destroy

!------------------------------------------------------------------------
! Specified reinjector output type
!------------------------------------------------------------------------

  subroutine specified_reinjector_output_init(self, reinjector, &
       flow_type)
    !! Initialises specified reinjector output.

    class(specified_reinjector_output_type), intent(in out) :: self
    type(source_network_reinjector_type), target, intent(in) :: reinjector
    PetscInt, intent(in) :: flow_type

    self%reinjector => reinjector
    self%out => null()
    self%flow_type = flow_type
    call self%allocate_variables()

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
       water_enthalpy = self%specified_enthalpy
       steam_enthalpy = 0._dp
    case (SEPARATED_FLOW_TYPE_STEAM)
       water_enthalpy = 0._dp
       steam_enthalpy = self%specified_enthalpy
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
       water_enthalpy = self%reinjector%in_water_enthalpy
       steam_enthalpy = 0._dp
    case (SEPARATED_FLOW_TYPE_STEAM)
       water_enthalpy = 0._dp
       steam_enthalpy = self%reinjector%in_steam_enthalpy
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

  subroutine specified_reinjector_output_enthalpies(self, &
       water_enthalpy, steam_enthalpy, enthalpy_specified)
    !! Gets enthalpies for specified reinjector output, and
    !! whether these are specified or defaults.

    class(specified_reinjector_output_type), intent(in out) :: self
    PetscReal, intent(out) :: water_enthalpy, steam_enthalpy
    PetscBool, intent(out) :: enthalpy_specified

    if (self%specified_enthalpy > 0._dp) then
       call self%specified_enthalpies(water_enthalpy, steam_enthalpy)
       enthalpy_specified = PETSC_TRUE
    else
       call self%default_enthalpies(water_enthalpy, steam_enthalpy)
       enthalpy_specified = PETSC_FALSE
    end if

  end subroutine specified_reinjector_output_enthalpies

!------------------------------------------------------------------------

  subroutine specified_reinjector_output_node_limit(self, rate)
    !! Limits injection rate to what can be injected into the output
    !! node.

    class(specified_reinjector_output_type), intent(in out) :: self
    PetscReal, intent(in out) :: rate
    ! Locals:
    PetscReal :: node_rate

    if (associated(self%out)) then

       select type (n => self%out)
       class is (source_type)
          node_rate = n%specified_injection_rate()
          call node_limit_rate(node_rate, rate)
       class is (source_network_reinjector_type)
          select case (self%flow_type)
          case (SEPARATED_FLOW_TYPE_WATER)
             call node_limit_rate(n%water_rate, rate)
          case (SEPARATED_FLOW_TYPE_STEAM)
             call node_limit_rate(n%steam_rate, rate)
          end select
       end select

    end if

  end subroutine specified_reinjector_output_node_limit

!------------------------------------------------------------------------
! Rate reinjector output type
!------------------------------------------------------------------------

  subroutine rate_reinjector_output_rates(self, water_rate, steam_rate)
    !! Gets rates for rate reinjector output.

    class(rate_reinjector_output_type), intent(in out) :: self
    PetscReal, intent(out) :: water_rate, steam_rate

    select case (self%flow_type)
    case (SEPARATED_FLOW_TYPE_WATER)
       water_rate = self%specified_rate
       steam_rate = 0._dp
    case (SEPARATED_FLOW_TYPE_STEAM)
       water_rate = 0._dp
       steam_rate = self%specified_rate
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
       water_rate = self%proportion * abs(self%reinjector%in_water_rate)
       steam_rate = 0._dp
    case (SEPARATED_FLOW_TYPE_STEAM)
       water_rate = 0._dp
       steam_rate = self%proportion * abs(self%reinjector%in_steam_rate)
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
    self%out => null()

  end subroutine overflow_reinjector_output_init

!------------------------------------------------------------------------

  subroutine overflow_reinjector_output_get_world_rank(self)
    !! Gets rank in world communicator of output node, if it exists,
    !! otherwise -1. This is used to send flows from a reinjector's
    !! root rank to the overflow output rank.

    class(overflow_reinjector_output_type), intent(in out) :: self
    ! Locals:
    PetscMPIInt :: rank, overflow_world_rank, max_overflow_world_rank
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    if (associated(self%out)) then
       overflow_world_rank = rank
    else
       overflow_world_rank = -1
    end if
    call MPI_allreduce(overflow_world_rank, max_overflow_world_rank, 1, &
         MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr)
    self%out_world_rank = max_overflow_world_rank

  end subroutine overflow_reinjector_output_get_world_rank

!------------------------------------------------------------------------

  subroutine overflow_reinjector_output_allocate_variables(self)
    !! Allocates overflow rate, enthalpy etc. variables for ranks
    !! other than the reinjector root rank (on which these are
    !! pointers into the reinjection vector, for output purposes.)

    class(overflow_reinjector_output_type), intent(in out) :: self

    if (self%reinjector%rank > 0) then
       call self%reinjector_output_type%allocate_variables()
    end if

  end subroutine overflow_reinjector_output_allocate_variables

!------------------------------------------------------------------------

  subroutine overflow_reinjector_output_deallocate_variables(self)
    !! Deallocates overflow rate, enthalpy etc. variables for ranks
    !! other than the reinjector root rank.

    class(overflow_reinjector_output_type), intent(in out) :: self

    if (self%reinjector%rank > 0) then
       call self%reinjector_output_type%deallocate_variables()
    end if

  end subroutine overflow_reinjector_output_deallocate_variables

!------------------------------------------------------------------------

  subroutine overflow_reinjector_output_assign(self, data, offset)
    !! Assigns pointers in overflow reinjector output object to
    !! elements in the data array, starting from the specified offset.

    class(overflow_reinjector_output_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! reinjector data array
    PetscInt, intent(in) :: offset  !! source array offset

    self%water_rate => data(offset)
    self%water_enthalpy => data(offset + 1)
    self%steam_rate => data(offset + 2)
    self%steam_enthalpy => data(offset + 3)
    self%rate => data(offset)
    self%enthalpy => data(offset + 1)

  end subroutine overflow_reinjector_output_assign

!------------------------------------------------------------------------

  subroutine overflow_reinjector_output_set_flows(self, water_rate, &
       water_enthalpy, steam_rate, steam_enthalpy)
    !! Sets overflow output flow rates and enthalpies.

    class(overflow_reinjector_output_type), intent(in out) :: self
    PetscReal, intent(in) :: water_rate !! Separated water mass flow rate
    PetscReal, intent(in) :: water_enthalpy !! Separated water enthalpy
    PetscReal, intent(in) :: steam_rate !! Separated steam mass flow rate
    PetscReal, intent(in) :: steam_enthalpy !! Separated steam enthalpy

    self%rate = rate
    self%enthalpy = enthalpy
    self%water_rate = water_rate
    self%water_enthalpy = water_enthalpy
    self%steam_rate = steam_rate
    self%steam_enthalpy = steam_enthalpy

  end subroutine overflow_reinjector_output_set_flows

!------------------------------------------------------------------------

  subroutine overflow_reinjector_output_node_limit(self, rate)
    !! Limits injection rate from overflow to what can be injected
    !! into the output node. The overflow is limited only if
    !! reinjecting it directly to a source. If reinjecting it to
    !! another reinjector, it is permitted to exceed that reinjector's
    !! capacity, and the excess will be assigned to its own overflow.

    class(overflow_reinjector_output_type), intent(in out) :: self
    PetscReal, intent(in out) :: rate
    ! Locals:
    PetscReal :: node_rate

    if (associated(self%out)) then

       select type (n => self%out)
       class is (source_type)
          node_rate = n%specified_injection_rate()
          call node_limit_rate(node_rate, rate)
       end select

    end if

  end subroutine overflow_reinjector_output_node_limit

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
    allocate(self%overflow)
    call self%overflow%init(self)

  end subroutine source_network_reinjector_init

!------------------------------------------------------------------------

  subroutine source_network_reinjector_init_comm(self)
    !! Initialises MPI communicator, rank and root world rank for the
    !! reinjector. It is assumed that the output node list has already
    !! been populated and overflow (if any) has been assigned.

    use mpi_utils_module, only: mpi_comm_root_world_rank

    class(source_network_reinjector_type), intent(in out) :: self
    ! Locals:
    PetscInt :: colour, i
    PetscInt, allocatable :: local_index(:)
    PetscErrorCode :: ierr

    colour = MPI_UNDEFINED
    call self%out%traverse(reinjector_comm_iterator)

    call MPI_comm_split(PETSC_COMM_WORLD, colour, 0, self%comm, ierr)
    if (self%comm /= MPI_COMM_NULL) then
       call MPI_comm_rank(self%comm, self%rank, ierr)
    else
       self%rank = -1
    end if
    self%root_world_rank = mpi_comm_root_world_rank(self%comm)

    call get_gather_parameters()

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

    subroutine get_gather_parameters()
      !! Finds counts, displacements and order (to match the order of
      !! the self%out list) for data gathered over the reinjector
      !! communicator.

      use mpi_utils_module, only: get_mpi_int_gather_array
      use utils_module, only: array_cumulative_sum

      ! Locals:
      PetscMPIInt :: comm_size
      PetscInt, allocatable :: indices_all(:)
      PetscErrorCode :: ierr

      self%local_gather_count = 0
      self%gather_count = 0

      if (self%rank >= 0) then

         call self%out%traverse(local_gather_count_iterator)

         call MPI_comm_size(self%comm, comm_size, ierr)

         self%gather_counts = get_mpi_int_gather_array(self%comm)
         self%gather_displacements = get_mpi_int_gather_array(self%comm)
         call MPI_gather(self%local_gather_count, 1, MPI_INTEGER, &
              self%gather_counts, 1, MPI_INTEGER, 0, self%comm, ierr)
         if (self%rank == 0) then
            self%gather_displacements = [[0], &
                 array_cumulative_sum(self%gather_counts(1: comm_size - 1))]
            self%gather_count = sum(self%gather_counts)
         else
            self%gather_count = 1
         end if

         allocate(local_index(self%local_gather_count))
         local_index = -1
         i = 1
         call self%out%traverse(output_index_iterator)

         allocate(self%gather_index(self%gather_count), &
              indices_all(self%gather_count))
         call MPI_gatherv(local_index, self%local_gather_count, MPI_INTEGER, &
              indices_all, self%gather_counts, self%gather_displacements, &
              MPI_INTEGER, 0, self%comm, ierr)

         if (self%rank == 0) then
            self%gather_index = [(i, i = 0, self%gather_count - 1)]
            call PetscSortIntWithPermutation(self%gather_count, indices_all, &
                 self%gather_index, ierr); CHKERRQ(ierr)
            self%gather_index = self%gather_index + 1
         end if

         deallocate(indices_all)

      end if

    end subroutine get_gather_parameters

!........................................................................

    subroutine local_gather_count_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (output => node%data)
      class is (reinjector_output_type)
         if (output%out%link_index >= 0) then
            self%local_gather_count = self%local_gather_count + 1
         end if
      end select

    end subroutine local_gather_count_iterator

!........................................................................

    subroutine output_index_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (output => node%data)
      class is (reinjector_output_type)
         if (output%out%link_index >= 0) then
            local_index(i) = output%out%link_index
            i = i + 1
         end if
      end select

    end subroutine output_index_iterator

  end subroutine source_network_reinjector_init_comm

!------------------------------------------------------------------------

  subroutine source_network_reinjector_init_in_comm(self)
    !! Creates a communicator which includes all processes in the main
    !! self%comm communicator, plus the process containing the input
    !! node for the reinjector. This is used for broadcasting inflow
    !! parameters to the outputs. Also identifies the rank of the
    !! input node in this communicator.

    class(source_network_reinjector_type), intent(in out) :: self
    ! Locals:
    PetscInt :: colour
    PetscMPIInt :: in_comm_input_rank, max_in_comm_input_rank
    PetscErrorCode :: ierr

    if ((self%rank >= 0) .or. (associated(self%in))) then
       colour = 1
    else
       colour = MPI_UNDEFINED
    end if

    call MPI_comm_split(PETSC_COMM_WORLD, colour, 0, self%in_comm, ierr)

    if (self%in_comm /= MPI_COMM_NULL) then

       call MPI_comm_rank(self%in_comm, self%in_comm_rank, ierr)

       if (associated(self%in)) then
          in_comm_input_rank = self%in_comm_rank
       else
          in_comm_input_rank = -1
       end if
       call MPI_allreduce(in_comm_input_rank, max_in_comm_input_rank, 1, &
            MPI_INTEGER, MPI_MAX, self%in_comm, ierr)
       self%in_comm_input_rank = max_in_comm_input_rank

    else
       self%in_comm_rank = -1
       self%in_comm_input_rank = -1
    end if

  end subroutine source_network_reinjector_init_in_comm

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
    self%output_rate => data(reinjector_offset + 1)
    self%output_water_rate => data(reinjector_offset + 2)
    self%output_steam_rate => data(reinjector_offset + 3)
    call self%overflow%assign(data, reinjector_offset + 4)

  end subroutine source_network_reinjector_assign

!------------------------------------------------------------------------

  subroutine source_network_reinjector_init_data(self, reinjector_index)
    !! Initialised source network reinjector variables accessed via
    !! pointers to the reinjector vector. The reinjector assign() method
    !! must be called first.

    class(source_network_reinjector_type), intent(in out) :: self
    PetscInt, intent(in) :: reinjector_index !! Index of reinjector in input

    self%reinjector_index = dble(reinjector_index)

    self%output_rate = 0._dp
    self%output_water_rate = 0._dp
    self%output_steam_rate = 0._dp

    self%overflow%rate = 0._dp
    self%overflow%enthalpy = 0._dp
    self%overflow%water_rate = 0._dp
    self%overflow%water_enthalpy = 0._dp
    self%overflow%steam_rate = 0._dp
    self%overflow%steam_enthalpy = 0._dp

  end subroutine source_network_reinjector_init_data

!------------------------------------------------------------------------

  subroutine source_network_reinjector_overflow_output(self, water_balance, &
       water_enthalpy, steam_balance, steam_enthalpy)
    !! Updates overflow data in reinjection vector, and assigns to
    !! overflow output node if there is one.

    use mpi_utils_module, only: mpi_comm_send

    class(source_network_reinjector_type), intent(in out) :: self
    PetscReal, intent(in out) :: water_balance, water_enthalpy
    PetscReal, intent(in out) :: steam_balance, steam_enthalpy
    ! Locals:
    PetscBool, parameter :: output_enthalpy_specified = PETSC_FALSE

    if (self%rank == 0) then
       call self%overflow%set_flows(water_balance, water_enthalpy, &
            steam_balance, steam_enthalpy)
    end if

    if (self%overflow%out_world_rank >= 0) then

       call mpi_comm_send(PETSC_COMM_WORLD, water_balance, &
            self%root_world_rank, self%overflow%out_world_rank)
       call mpi_comm_send(PETSC_COMM_WORLD, water_enthalpy, &
            self%root_world_rank, self%overflow%out_world_rank)
       call mpi_comm_send(PETSC_COMM_WORLD, steam_balance, &
            self%root_world_rank, self%overflow%out_world_rank)
       call mpi_comm_send(PETSC_COMM_WORLD, steam_enthalpy, &
            self%root_world_rank, self%overflow%out_world_rank)

       call self%overflow%update(water_balance, water_enthalpy, &
            steam_balance, steam_enthalpy, output_enthalpy_specified)

    end if

  end subroutine source_network_reinjector_overflow_output

!------------------------------------------------------------------------

  subroutine source_network_reinjector_capacity(self)
    !! Calculates output capacity of a reinjector and assigns this
    !! capacity to its rate property. If any outputs are unrated
    !! (i.e. have no capacity specified, flagged by their flow rates
    !! set to -1) then the reinjector capacity is also set to -1.

    class(source_network_reinjector_type), intent(in out) :: self
    ! Locals:
    PetscReal :: local_water_capacity, local_steam_capacity
    PetscReal :: water_capacity, steam_capacity
    PetscErrorCode :: ierr

    local_water_capacity = 0._dp
    local_steam_capacity = 0._dp

    call self%out%traverse(output_capacity_iterator)

    water_capacity = total_capacity(local_water_capacity)
    if (self%rank == 0) self%water_rate = water_capacity
    steam_capacity = total_capacity(local_steam_capacity)
    if (self%rank == 0) self%steam_rate = steam_capacity

  contains

    PetscReal function total_capacity(local_capacity)
      !! Calculates total capacity (-1 if unrated).

      PetscReal, intent(in) :: local_capacity
      ! Locals:
      PetscBool :: unrated

      call MPI_allreduce(local_capacity < 0._dp, unrated, 1, &
           MPI_LOGICAL, MPI_LOR, self%comm, ierr)
      if (unrated) then
         total_capacity = -1._dp
      else
         call MPI_reduce(local_capacity, total_capacity, 1, &
              MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
      end if

    end function total_capacity

!........................................................................

    subroutine update_capacity(rate, capacity)
      !! Updates capacity according to rate. If rate = -1, capacity is
      !! also set to -1.

      PetscReal, intent(in) :: rate
      PetscReal, intent(in out) :: capacity

      if (rate > -1._dp) then
         if (capacity > -1._dp) capacity = capacity + rate
      else
         capacity = -1._dp
      end if

    end subroutine update_capacity

!........................................................................

    subroutine output_capacity_iterator(node, stopped)
      !! Updates water and steam capacities from reinjector output.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscReal :: node_rate

      stopped = PETSC_FALSE
      select type (output => node%data)
      class is (specified_reinjector_output_type)
         if (output%out%link_index >= 0) then
            select type (n => output%out)
            class is (source_type)
               node_rate = n%specified_injection_rate()
               select case (output%flow_type)
               case (SEPARATED_FLOW_TYPE_WATER)
                  call update_capacity(node_rate, local_water_capacity)
               case (SEPARATED_FLOW_TYPE_STEAM)
                  call update_capacity(node_rate, local_steam_capacity)
               end select
            class is (source_network_reinjector_type)
               select case (output%flow_type)
               case (SEPARATED_FLOW_TYPE_WATER)
                  call update_capacity(n%water_rate, local_water_capacity)
               case (SEPARATED_FLOW_TYPE_STEAM)
                  call update_capacity(n%steam_rate, local_steam_capacity)
               end select
            end select
         end if
      end select

    end subroutine output_capacity_iterator

  end subroutine source_network_reinjector_capacity

!------------------------------------------------------------------------

  subroutine source_network_reinjector_distribute(self)
    !! Distributes reinjector input flow to outputs (and overflow if
    !! needed).

    class(source_network_reinjector_type), intent(in out) :: self
    ! Locals:
    PetscReal :: water_balance
    PetscReal :: steam_balance
    PetscReal :: local_qw(self%local_gather_count), qw(self%gather_count)
    PetscReal :: local_qs(self%local_gather_count), qs(self%gather_count)
    PetscInt :: i, j
    PetscErrorCode :: ierr

    if (associated(self%in)) then
       self%in_water_rate = abs(self%in%water_rate)
       self%in_water_enthalpy = self%in%water_enthalpy
       self%in_steam_rate = abs(self%in%steam_rate)
       self%in_steam_enthalpy = self%in%steam_enthalpy
    end if
    if ((self%in_comm_rank >= 0) .and. (self%in_comm_input_rank >= 0)) then
       call MPI_bcast(self%in_water_rate, 1, MPI_DOUBLE_PRECISION, &
            self%in_comm_input_rank, self%in_comm, ierr)
       call MPI_bcast(self%in_water_enthalpy, 1, MPI_DOUBLE_PRECISION, &
            self%in_comm_input_rank, self%in_comm, ierr)
       call MPI_bcast(self%in_steam_rate, 1, MPI_DOUBLE_PRECISION, &
            self%in_comm_input_rank, self%in_comm, ierr)
       call MPI_bcast(self%in_steam_enthalpy, 1, MPI_DOUBLE_PRECISION, &
            self%in_comm_input_rank, self%in_comm, ierr)
    else
       self%in_water_rate = 0._dp
       self%in_water_enthalpy = 0._dp
       self%in_steam_rate = 0._dp
       self%in_steam_enthalpy = 0._dp
    end if

    i = 1
    call self%out%traverse(local_rates_iterator)

    if (self%rank >= 0) then

       call MPI_gatherv(local_qw, self%local_gather_count, &
            MPI_DOUBLE_PRECISION, qw, self%gather_counts, &
            self%gather_displacements, MPI_DOUBLE_PRECISION, 0, &
            self%comm, ierr)
       call MPI_gatherv(local_qs, self%local_gather_count, &
            MPI_DOUBLE_PRECISION, qs, self%gather_counts, &
            self%gather_displacements, MPI_DOUBLE_PRECISION, 0, &
            self%comm, ierr)

       if (self%rank == 0) then

          ! Set rates and enthalpies for output:
          self%water_rate = self%in_water_rate
          self%water_enthalpy = self%in_water_enthalpy
          self%steam_rate = self%in_steam_rate
          self%steam_enthalpy = self%in_steam_enthalpy

          self%output_water_rate = 0._dp
          self%output_steam_rate = 0._dp
          water_balance = self%in_water_rate
          steam_balance = self%in_steam_rate

          do i = 1, self%gather_count
             j = self%gather_index(i)
             call limit_rate(qw(j), water_balance, self%output_water_rate)
             call limit_rate(qs(j), steam_balance, self%output_steam_rate)
          end do
          self%output_rate = self%output_water_rate + self%output_steam_rate

       end if

       call MPI_scatterv(qw, self%gather_counts, &
            self%gather_displacements, MPI_DOUBLE_PRECISION, &
            local_qw, self%local_gather_count, &
            MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
       call MPI_scatterv(qs, self%gather_counts, &
            self%gather_displacements, MPI_DOUBLE_PRECISION, &
            local_qs, self%local_gather_count, &
            MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
    end if

    i = 1
    call self%out%traverse(local_update_iterator)

    call self%overflow_output(water_balance, self%in_water_enthalpy, &
         steam_balance, self%in_steam_enthalpy)

  contains

    subroutine local_rates_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (output => node%data)
      class is (specified_reinjector_output_type)
         if (output%out%link_index >= 0) then
            call output%rates(local_qw(i), local_qs(i))
            select case (output%flow_type)
            case (SEPARATED_FLOW_TYPE_WATER)
               call output%node_limit(local_qw(i))
            case (SEPARATED_FLOW_TYPE_STEAM)
               call output%node_limit(local_qs(i))
            end select
            i = i + 1
         end if
      end select

    end subroutine local_rates_iterator

!........................................................................

    subroutine limit_rate(rate, balance, total)
      !! Limits rate to remaining balance, and updates balance and
      !! total.

      PetscReal, intent(in out) :: rate, balance, total

      if (rate < 0._dp) then
         ! No limit on flow rate - set to remaining balance:
         rate = balance
      end if

      ! Rate cannot exceed remaining balance:
      rate = min(rate, balance)

      ! Update balance:
      balance = max(balance - rate, 0._dp)
      total = total + rate

    end subroutine limit_rate

!........................................................................

    subroutine local_update_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscReal :: water_enthalpy, steam_enthalpy
      PetscBool :: output_enthalpy_specified

      stopped = PETSC_FALSE
      select type (output => node%data)
      class is (specified_reinjector_output_type)
         if (output%out%link_index >= 0) then
            call output%enthalpies(water_enthalpy, steam_enthalpy, &
                 output_enthalpy_specified)
            call output%update(local_qw(i), water_enthalpy, local_qs(i), &
                 steam_enthalpy, output_enthalpy_specified)
            i = i + 1
         end if
      end select

    end subroutine local_update_iterator

  end subroutine source_network_reinjector_distribute

!------------------------------------------------------------------------

  subroutine source_network_reinjector_destroy(self)
    !! Destroys a source network reinjector.

    class(source_network_reinjector_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    if (associated(self%in)) self%in => null()
    call self%out%destroy(reinjector_output_list_node_data_destroy)
    call self%overflow%destroy()
    deallocate(self%overflow)

    call MPI_comm_free(self%comm, ierr)
    call MPI_comm_free(self%in_comm, ierr)

    if (allocated(self%gather_index)) deallocate(self%gather_index)
    if (allocated(self%gather_counts)) deallocate(self%gather_counts)
    if (allocated(self%gather_displacements)) deallocate(self%gather_displacements)

    call self%source_network_node_type%destroy()

  end subroutine source_network_reinjector_destroy

!------------------------------------------------------------------------

end module source_network_reinjector_module
