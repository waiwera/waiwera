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

module source_network_group_module
  !! Module for source groups.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use source_network_module
  use source_module, only: source_type
  use separator_module, only: num_separator_variables, separator_variable_names
  use list_module, only: list_type, list_node_type
  use hdf5io_module, only: max_field_name_length
  use thermodynamics_module

  implicit none
  private

  PetscInt, parameter, public :: num_source_network_group_variables = &
       num_source_network_node_variables + num_separator_variables + 1
  PetscInt, parameter, public :: max_source_network_group_variable_name_length = 24
  character(max_source_network_group_variable_name_length), parameter, public :: &
       source_network_group_variable_names(num_source_network_group_variables) = [ &
       source_network_variable_names, separator_variable_names, [ &
       "group_index         "]]
  PetscInt, parameter, public :: num_source_network_group_constant_integer_variables = 1
  character(max_source_network_group_variable_name_length), parameter, public :: &
       source_network_group_constant_integer_variables( &
       num_source_network_group_constant_integer_variables) = [ &
       "group_index       "]
  character(max_field_name_length), parameter, public :: &
       required_output_source_network_group_fields(0) = [&
       character(max_field_name_length)::]
  character(max_field_name_length), parameter, public :: &
       default_output_source_network_group_fields(2) = [&
       "rate              ", "enthalpy          "]

  type, public, extends(source_network_node_type) :: source_network_group_type
     !! Type for group of source network nodes, e.g. multi-feed well
     !! or group of makeup wells.
     private
     type(list_type), public :: in !! List of input nodes in group
     class(source_network_node_type), pointer, public :: out !! Output node for group (if any)
     MPI_Comm :: comm !! MPI communicator for group
     PetscMPIInt, public :: rank !! Rank of group in its own communicator
     PetscMPIInt :: root_world_rank !! Rank in world communicator of group root rank
     PetscInt, public :: local_group_index !! Index of group in local part of source group vector (-1 if not a root group)
     PetscReal, pointer, public :: group_index !! Index of source group in input
   contains
     private
     procedure, public :: init => source_network_group_init
     procedure, public :: init_comm => source_network_group_init_comm
     procedure, public :: assign => source_network_group_assign
     procedure, public :: init_data => source_network_group_init_data
     procedure, public :: set_rate => source_network_group_set_rate
     procedure, public :: limit_rate => source_network_group_limit_rate
     procedure, public :: sum => source_network_group_sum
     procedure, public :: sum_out => source_network_group_sum_out
     procedure, public :: default_separated_flows => source_network_group_default_separated_flows
     procedure, public :: get_separated_flows => source_network_group_get_separated_flows
     procedure, public :: add_separated_flows => source_network_group_add_separated_flows
     procedure, public :: add_flows => source_network_group_add_flows
     procedure, public :: destroy => source_network_group_destroy
  end type source_network_group_type

  type, public, extends(source_network_group_type) :: &
       uniform_scaling_source_network_group_type
     !! Type for source network group in which inputs are scaled
     !! uniformly by a constant factor.
   contains
     procedure, public :: scale_rate => uniform_scaling_source_network_group_scale_rate
  end type uniform_scaling_source_network_group_type

  type, public, extends(source_network_group_type) :: &
       progressive_scaling_source_network_group_type
     !! Type for source network group in which inputs are scaled
     !! progressively, starting from the last input in the group and
     !! proceeding back to the first.
     private
     PetscInt :: local_gather_count !! How many local inputs are included in gather operations
     PetscInt, allocatable :: gather_counts(:) !! Process counts for gather operations
     PetscInt :: gather_count !! Total count for gather operations (only computed on root rank)
     PetscInt, allocatable :: gather_displacements(:) !! Process displacements for gather operations
     PetscInt, allocatable :: gather_order(:) !! Sort order for gather operations
   contains
     procedure, public :: init_comm => progressive_scaling_source_network_group_init_comm
     procedure, public :: scale_rate => progressive_scaling_source_network_group_scale_rate
     procedure, public :: destroy => progressive_scaling_source_network_group_destroy
  end type progressive_scaling_source_network_group_type

contains

!------------------------------------------------------------------------
!  Source network group
!------------------------------------------------------------------------

  subroutine source_network_group_init(self, name)
    !! Initialises a source network group. Only variables stored in
    !! the object itself are initialised, not those stored in the
    !! group data vector and accessed via pointers.

    class(source_network_group_type), intent(in out) :: self
    character(*), intent(in) :: name !! Group name

    self%name = name

    call self%in%init(owner = PETSC_FALSE)
    self%out => null()
    self%out_input_index = -1

  end subroutine source_network_group_init

!------------------------------------------------------------------------

  subroutine source_network_group_init_comm(self)
    !! Initialises MPI communicator, rank and root world rank for the
    !! group. It is assumed that the input nodes list has already been
    !! populated.

    class(source_network_group_type), intent(in out) :: self
    ! Locals:
    PetscInt :: colour
    PetscErrorCode :: ierr

    colour = MPI_UNDEFINED
    call self%in%traverse(group_comm_iterator)
    call MPI_comm_split(PETSC_COMM_WORLD, colour, 0, self%comm, ierr)

    call get_root_world_rank()

  contains

    subroutine group_comm_iterator(node, stopped)
      !! Stops and sets colour if there are any local sources, or
      !! groups with local sources, in the group node list.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (n => node%data)
      type is (source_type)
         colour = 1
         stopped = PETSC_TRUE
      class is (source_network_group_type)
         if (n%rank >= 0) then
            colour = 1
            stopped = PETSC_TRUE
         end if
      end select

    end subroutine group_comm_iterator

!........................................................................

    subroutine get_root_world_rank()
      !! Finds root world rank, i.e. rank in world communicator of
      !! root rank of the group.

      ! Locals:
      PetscMPIInt :: rank, root_world_rank

      call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
      root_world_rank = -1
      if (self%comm /= MPI_COMM_NULL) then
         call MPI_comm_rank(self%comm, self%rank, ierr)
         if (self%rank == 0) then
            root_world_rank = rank
         end if
      else
         self%rank = -1
      end if

      call MPI_allreduce(root_world_rank, self%root_world_rank, 1, &
           MPI_INTEGER, MPI_MAX, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)

    end subroutine get_root_world_rank

  end subroutine source_network_group_init_comm

!------------------------------------------------------------------------

  subroutine source_network_group_assign(self, data, offset)
    !! Assigns pointers in source network group object to elements in
    !! the data array, starting from the specified offset.

    class(source_network_group_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! source data array
    PetscInt, intent(in) :: offset  !! source array offset
    ! Locals:
    PetscInt :: group_offset

    call self%source_network_node_type%assign(data, offset)

    group_offset = offset + num_source_network_node_variables + &
         num_separator_variables

    self%group_index => data(group_offset)

  end subroutine source_network_group_assign

!------------------------------------------------------------------------

  subroutine source_network_group_init_data(self, group_index, separator_pressure, &
       thermo)
    !! Initialised source network group variables accessed via
    !! pointers to the source group vector. The group assign() method
    !! must be called first.

    class(source_network_group_type), intent(in out) :: self
    PetscInt, intent(in) :: group_index !! Index of group in input
    PetscReal, intent(in) :: separator_pressure(:) !! Separator pressures ([-1] for no separator)
    class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics

    self%group_index = dble(group_index)
    call self%separator%init(separator_pressure, thermo)

  end subroutine source_network_group_init_data

!------------------------------------------------------------------------

  subroutine source_network_group_set_rate(self, rate)
    !! Sets source network group flow rate to specified value, on the
    !! root process of the group.

    class(source_network_group_type), intent(in out) :: self
    PetscReal, intent(in) :: rate !! Flow rate

    if (self%rank == 0) then
       self%rate = rate
    end if
    call self%get_separated_flows()

  end subroutine source_network_group_set_rate

!------------------------------------------------------------------------

  subroutine source_network_group_limit_rate(self, flow_type, limit)
    !! Limits source network group flow rates (total, water or steam
    !! as specified by flow_type) to specified limits.

    class(source_network_group_type), intent(in out) :: self
    PetscInt, intent(in) :: flow_type(:) !! Flow types
    PetscReal, intent(in) :: limit(:) !! Flow rate limits
    ! Locals:
    PetscReal :: scale
    PetscBool :: over
    PetscErrorCode :: ierr

    if (self%rank == 0) then
       call self%get_minimum_limit_scale(flow_type, limit, over, scale)
    end if
    call MPI_bcast(over, 1, MPI_LOGICAL, self%root_world_rank, &
         PETSC_COMM_WORLD, ierr)
    if (over) then
       call self%scale_rate(scale)
       call self%sum_out()
    end if

  end subroutine source_network_group_limit_rate

!------------------------------------------------------------------------

  subroutine source_network_group_sum(self)
    !! Computes total flow rate, enthalpy etc. in a source group. The
    !! results are stored only on the root rank of the group
    !! communicator.

    use dm_utils_module, only: global_section_offset

    class(source_network_group_type), intent(in out) :: self
    ! Locals:
    PetscReal :: local_q, local_qh
    PetscReal :: total_q, total_qh
    PetscErrorCode :: ierr
    PetscReal, parameter :: rate_tol = 1.e-9_dp

    local_q = 0._dp
    local_qh = 0._dp

    call self%in%traverse(group_sum_iterator)

    call MPI_reduce(local_q, total_q, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
    call MPI_reduce(local_qh, total_qh, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)

    if (self%rank == 0) then
       if (abs(total_q) > rate_tol) then
          self%enthalpy = total_qh / total_q
       else
          self%enthalpy = 0._dp
       end if
    end if
    call self%set_rate(total_q)

  contains

    subroutine group_sum_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (s => node%data)
      class is (source_network_node_type)
         call s%add_flows(local_q, local_qh)
      end select

    end subroutine group_sum_iterator

  end subroutine source_network_group_sum

!------------------------------------------------------------------------

  recursive subroutine source_network_group_sum_out(self)
    !! Re-computes downstream flow rate, enthalpy etc. in output
    !! network group (if there is one).

    class(source_network_group_type), intent(in out) :: self

    if (associated(self%out)) then
       select type (group => self%out)
       class is (source_network_group_type)
         call group%sum()
         call group%sum_out()
       end select
    end if

  end subroutine source_network_group_sum_out

!------------------------------------------------------------------------

  subroutine source_network_group_default_separated_flows(self)
      !! Gets default separated water and steam flows when the group
      !! does not have its own separator. These are calculated by
      !! summing the separated flows from the inputs.

    class(source_network_group_type), intent(in out) :: self

    ! Locals:
    PetscReal :: local_water_q, local_water_qh
    PetscReal :: total_water_q, total_water_qh
    PetscReal :: local_steam_q, local_steam_qh
    PetscReal :: total_steam_q, total_steam_qh
    PetscErrorCode :: ierr
    PetscReal, parameter :: rate_tol = 1.e-9_dp

    local_water_q = 0._dp
    local_water_qh = 0._dp
    local_steam_q = 0._dp
    local_steam_qh = 0._dp

    call self%in%traverse(group_sum_separated_iterator)

    call MPI_reduce(local_water_q, total_water_q, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
    call MPI_reduce(local_water_qh, total_water_qh, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
    call MPI_reduce(local_steam_q, total_steam_q, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)
    call MPI_reduce(local_steam_qh, total_steam_qh, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, 0, self%comm, ierr)

    if (self%rank == 0) then
       self%water_rate = total_water_q
       self%steam_rate = total_steam_q
       if (abs(total_water_q) > rate_tol) then
          self%water_enthalpy = total_water_qh / total_water_q
       else
          self%water_enthalpy = 0._dp
       end if
       if (abs(total_steam_q) > rate_tol) then
          self%steam_enthalpy = total_steam_qh / total_steam_q
       else
          self%steam_enthalpy = 0._dp
       end if
    end if

  contains

    subroutine group_sum_separated_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (s => node%data)
      class is (source_network_node_type)
         call s%add_separated_flows(local_water_q, local_water_qh, &
              local_steam_q, local_steam_qh)
      end select

    end subroutine group_sum_separated_iterator

  end subroutine source_network_group_default_separated_flows

!------------------------------------------------------------------------

  subroutine source_network_group_get_separated_flows(self)
      !! Gets separated water and steam flows on the root process of
      !! the group.

    class(source_network_group_type), intent(in out) :: self
    ! Locals:
    PetscBool :: use_default
    PetscErrorCode :: ierr

    use_default = PETSC_FALSE

    if (self%rank == 0) then
       if (self%rate < 0._dp) then
          if (self%separator%on) then
             call self%separate()
          else
             use_default = PETSC_TRUE
          end if
       else
          call self%zero_separated_flows()
       end if
    end if

    call MPI_bcast(use_default, 1, MPI_LOGICAL, 0, self%comm, ierr)
    if (use_default) then
       call self%default_separated_flows()
    end if

  end subroutine source_network_group_get_separated_flows

!------------------------------------------------------------------------

  subroutine source_network_group_add_flows(self, mass_flow, energy_flow)
    !! Adds mass and energy flows from the group to the specified
    !! totals (only on the root process of the group).

    class(source_network_group_type), intent(in out) :: self
    PetscReal, intent(in out) :: mass_flow !! Total mass flow rate
    PetscReal, intent(in out) :: energy_flow !! Total energy flow rate

    if (self%rank == 0) then
       call self%source_network_node_type%add_flows(mass_flow, energy_flow)
    end if

  end subroutine source_network_group_add_flows

!------------------------------------------------------------------------

  subroutine source_network_group_add_separated_flows(self, water_mass_flow, &
       water_energy_flow, steam_mass_flow, steam_energy_flow)
    !! Adds mass and energy flows for separated water and steam from
    !! the group to the specified totals, on the group's root process
    !! only.

    class(source_network_group_type), intent(in out) :: self
    PetscReal, intent(in out) :: water_mass_flow !! Total water mass flow rate
    PetscReal, intent(in out) :: water_energy_flow !! Total water energy flow rate
    PetscReal, intent(in out) :: steam_mass_flow !! Total steam mass flow rate
    PetscReal, intent(in out) :: steam_energy_flow !! Total steam energy flow rate

    if (self%rank == 0) then
       call self%source_network_node_type%add_separated_flows(water_mass_flow, &
       water_energy_flow, steam_mass_flow, steam_energy_flow)
    end if

  end subroutine source_network_group_add_separated_flows

!------------------------------------------------------------------------

  subroutine source_network_group_destroy(self)
    !! Destroys a source network group.

    class(source_network_group_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call self%in%destroy()
    if (associated(self%out)) self%out => null()
    call MPI_comm_free(self%comm, ierr)

    call self%source_network_node_type%destroy()

  end subroutine source_network_group_destroy
    
!------------------------------------------------------------------------
! Uniform scaling source network group
!------------------------------------------------------------------------

  recursive subroutine uniform_scaling_source_network_group_scale_rate(self, scale)
    !! Scales source network group flow rate by specified scale
    !! factor, by scaling flows in group input nodes uniformly, all by
    !! the same scale factor.

    class(uniform_scaling_source_network_group_type), intent(in out) :: self
    PetscReal, intent(in) :: scale !! Flow rate scale factor
    ! Locals:
    PetscErrorCode :: ierr

    call MPI_bcast(scale, 1, MPI_DOUBLE_PRECISION, 0, self%comm, ierr)
    call self%in%traverse(group_scale_iterator)
    call self%sum()

  contains

    recursive subroutine group_scale_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (s => node%data)
      class is (source_network_node_type)
         call s%scale_rate(scale)
      end select

    end subroutine group_scale_iterator

  end subroutine uniform_scaling_source_network_group_scale_rate

!------------------------------------------------------------------------
! Progressive scaling source network group
!------------------------------------------------------------------------

  subroutine progressive_scaling_source_network_group_init_comm(self)
    !! Initialises MPI communicator, rank and root world rank for the
    !! group, as well as parameters related to gathering data from the
    !! relevant input nodes.

    class(progressive_scaling_source_network_group_type), intent(in out) :: self
    ! Locals:
    PetscInt :: i
    PetscInt, allocatable :: local_index(:)

    call self%source_network_group_type%init_comm()
    call get_gather_parameters()

  contains

    subroutine get_gather_parameters()
      !! Finds counts, displacements and sort order (to match the
      !! order of the self%in list) for data gathered over the group
      !! communicator.

      use mpi_utils_module, only: get_mpi_int_gather_array
      use utils_module, only: array_cumulative_sum

      ! Locals:
      PetscMPIInt :: comm_size
      PetscInt, allocatable :: indices_all(:)
      PetscErrorCode :: ierr

      self%local_gather_count = 0

      if (self%rank >= 0) then

         call self%in%traverse(local_gather_count_iterator)

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
         call self%in%traverse(input_index_iterator)

         allocate(self%gather_order(self%gather_count), &
              indices_all(self%gather_count))
         call MPI_gatherv(local_index, self%local_gather_count, MPI_INTEGER, &
              indices_all, self%gather_counts, self%gather_displacements, &
              MPI_INTEGER, 0, self%comm, ierr)

         if (self%rank == 0) then
            self%gather_order = [(i, i = 0, self%gather_count - 1)]
            call PetscSortIntWithPermutation(self%gather_count, indices_all, &
                 self%gather_order, ierr); CHKERRQ(ierr)
            self%gather_order = self%gather_order + 1
         end if

         deallocate(indices_all)

      end if

    end subroutine get_gather_parameters

!........................................................................

    subroutine local_gather_count_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (network_node => node%data)
      class is (source_network_node_type)
         if (network_node%out_input_index >= 0) then
            self%local_gather_count = self%local_gather_count + 1
         end if
      end select

    end subroutine local_gather_count_iterator

!........................................................................

    subroutine input_index_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (network_node => node%data)
      class is (source_network_node_type)
         if (network_node%out_input_index >= 0) then
            local_index(i) = network_node%out_input_index
            i = i + 1
         end if
      end select

    end subroutine input_index_iterator

  end subroutine progressive_scaling_source_network_group_init_comm

!------------------------------------------------------------------------

  recursive subroutine progressive_scaling_source_network_group_scale_rate(self, &
       scale)
    !! Scales source network group flow rate by specified scale
    !! factor, by scaling flows in group input nodes progressively.

    use utils_module, only: array_progressive_scale

    class(progressive_scaling_source_network_group_type), intent(in out) :: self
    PetscReal, intent(in) :: scale !! Flow rate scale factor
    ! Locals:
    PetscInt :: i
    PetscReal :: local_q(self%local_gather_count), q(self%gather_count)
    PetscReal :: local_s(self%local_gather_count), s(self%gather_count)
    PetscReal :: target_rate
    PetscErrorCode :: ierr

    i = 1
    call self%in%traverse(get_local_rate_iterator)
    call MPI_gatherv(local_q, self%local_gather_count, &
         MPI_DOUBLE_PRECISION, q, self%gather_counts, &
         self%gather_displacements, MPI_INTEGER, 0, self%comm, ierr)

    if (self%rank == 0) then
       target_rate = scale * self%rate
       s = array_progressive_scale(q, target_rate, PETSC_TRUE)
    end if

    ! TODO:
    ! scatterv s to local_s
    ! call self%in%traverse() to scale inputs

    call self%sum()

  contains

    subroutine get_local_rate_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (network_node => node%data)
      class is (source_network_node_type)
         if (network_node%out_input_index >= 0) then
            local_q(i) = network_node%rate
            i = i + 1
         end if
      end select

    end subroutine get_local_rate_iterator

  end subroutine progressive_scaling_source_network_group_scale_rate

!------------------------------------------------------------------------

  subroutine progressive_scaling_source_network_group_destroy(self)
    !! Destroys a progressive scaling source network group.

    class(progressive_scaling_source_network_group_type), intent(in out) :: self

    if (allocated(self%gather_order)) deallocate(self%gather_order)
    if (allocated(self%gather_counts)) deallocate(self%gather_counts)
    if (allocated(self%gather_displacements)) deallocate(self%gather_displacements)

    call self%source_network_group_type%destroy()

  end subroutine progressive_scaling_source_network_group_destroy

!------------------------------------------------------------------------

end module source_network_group_module
