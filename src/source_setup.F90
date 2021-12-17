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

module source_setup_module
  !! Module for setting up sources and sinks.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use fson
  use fson_value_m, only: TYPE_INTEGER, TYPE_REAL, TYPE_ARRAY, &
       TYPE_STRING, TYPE_OBJECT, TYPE_LOGICAL
  use fson_mpi_module
  use list_module
  use logfile_module
  use eos_module
  use thermodynamics_module, only: thermodynamics_type
  use source_module
  use source_control_module
  use separator_module
  use source_group_module
  use source_network_module

  implicit none
  private

  public :: setup_sources

contains

!------------------------------------------------------------------------

  subroutine setup_sources(json, dm, ao, eos, tracer_names, thermo, start_time, &
       fluid_vector, fluid_range_start, source_vector, source_range_start, &
       sources, num_sources, source_controls, source_index, &
       separated_sources, source_groups, logfile, err)
    !! Sets up sinks / sources, source controls and source groups.

    use dm_utils_module
    use dictionary_module

    type(fson_value), pointer, intent(in) :: json !! JSON file object
    DM, intent(in) :: dm !! Mesh DM
    AO, intent(in) :: ao !! Application ordering for natural to global cell indexing
    class(eos_type), intent(in) :: eos !! Equation of state
    character(*), intent(in) :: tracer_names(:) !! Tracer names
    class(thermodynamics_type), intent(in out) :: thermo !! Thermodynamics formulation
    PetscReal, intent(in) :: start_time
    Vec, intent(in) :: fluid_vector !! Fluid vector
    PetscInt, intent(in) :: fluid_range_start !! Range start for global fluid vector
    Vec, intent(out) :: source_vector !! Source vector
    PetscInt, intent(out) :: source_range_start !! Range start for global source vector
    type(list_type), intent(in out) :: sources !! List of local source objects
    PetscInt, intent(out) :: num_sources !! Total number of sources created on all processes
    type(list_type), intent(in out) :: source_controls !! List of source controls
    IS, intent(in out) :: source_index !! IS defining natural-to-global source ordering
    type(list_type), intent(out) :: separated_sources !! List of sources with separators
    type(list_type), intent(in out) :: source_groups !! List of source groups
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscMPIInt :: rank
    DM :: dm_source
    PetscInt :: num_local_sources, source_spec_index, local_source_index
    type(fson_value), pointer :: sources_json, source_json, group_spec
    PetscInt :: num_source_specs, num_tracers
    PetscReal, pointer, contiguous :: fluid_data(:), source_data(:)
    PetscSection :: fluid_section, source_section
    type(list_type) :: group_specs
    type(dictionary_type) :: source_dict
    PetscInt, allocatable :: indices(:)
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    num_sources = 0
    num_tracers = size(tracer_names)
    err = 0

    call label_source_zones(json, num_local_sources, logfile, err)
    if (err == 0) then

       call create_path_dm(num_local_sources, dm_source)
       call setup_source_dm_data_layout(dm_source)
       call DMCreateGlobalVector(dm_source, source_vector, ierr); CHKERRQ(ierr)
       call PetscObjectSetName(source_vector, "source", ierr); CHKERRQ(ierr)
       call global_vec_range_start(source_vector, source_range_start)

       call global_vec_section(fluid_vector, fluid_section)
       call VecGetArrayReadF90(fluid_vector, fluid_data, ierr); CHKERRQ(ierr)
       call global_vec_section(source_vector, source_section)
       call VecGetArrayF90(source_vector, source_data, ierr); CHKERRQ(ierr)

       call sources%init(owner = PETSC_TRUE)
       call source_dict%init(owner = PETSC_FALSE)
       call separated_sources%init(owner = PETSC_FALSE)
       call source_controls%init(owner = PETSC_TRUE)
       call source_groups%init(owner = PETSC_TRUE)
       call group_specs%init(owner = PETSC_TRUE)

       if (fson_has_mpi(json, "source")) then
          call fson_get_mpi(json, "source", sources_json)
          num_source_specs = fson_value_count_mpi(sources_json, ".")
          source_json => fson_value_children_mpi(sources_json)
          local_source_index = 0
          do source_spec_index = 0, num_source_specs - 1
             if (fson_has_mpi(source_json, "source")) then ! source group
                allocate(group_spec)
                if (rank == 0) then
                   group_spec = source_json
                end if
                call group_specs%append(group_spec)
             else
                call setup_source(source_spec_index, local_source_index, &
                     source_json, ao, tracer_names, num_sources, &
                     separated_sources, thermo, sources, source_dict, err)
             end if
             if (err > 0) exit
             source_json => fson_value_next_mpi(source_json)
          end do
       else
          if (present(logfile)) then
             call logfile%write(LOG_LEVEL_INFO, "input", "no_sources")
          end if
       end if

       if (err == 0) then
          allocate(indices(num_local_sources))
          call sources%traverse(source_indices_iterator)
          call setup_source_index(indices, source_index)
          deallocate(indices)
          call group_specs%traverse(setup_group_iterator)
       end if

       call VecRestoreArrayF90(source_vector, source_data, ierr); CHKERRQ(ierr)
       call VecRestoreArrayReadF90(fluid_vector, fluid_data, ierr)
       CHKERRQ(ierr)
       call group_specs%destroy()
       call source_dict%destroy()

    end if

  contains

!........................................................................

    subroutine label_source_zones(json, num_local_sources, logfile, err)
      !! Set DM source label for cells defined on zones. Also returns
      !! total number of local sources.

      use zone_label_module
      use mpi_utils_module, only: mpi_broadcast_error_flag

      type(fson_value), pointer, intent(in) :: json
      PetscInt, intent(out) :: num_local_sources
      type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
      PetscErrorCode, intent(out) :: err !! Error code
      ! Locals:
      type(fson_value), pointer :: sources_json, source_json
      PetscInt :: num_source_specs, source_index, i, iz, ghost, num_cells
      PetscInt :: zones_type, num_zone_cells
      character(max_zone_name_length), allocatable :: zones(:)
      character(:), allocatable :: label_name
      PetscBool :: has_label
      IS :: cell_IS
      PetscInt, pointer, contiguous :: zone_cells(:), cells(:)
      PetscInt, allocatable :: labelled_cells(:)
      DMLabel :: ghost_label, source_label
      PetscErrorCode :: ierr

      num_local_sources = 0
      call DMGetLabel(dm, "ghost", ghost_label, ierr)
      call DMGetLabel(dm, source_label_name, source_label, ierr)

      if (fson_has_mpi(json, "source")) then

         call fson_get_mpi(json, "source", sources_json)
         num_source_specs = fson_value_count_mpi(sources_json, ".")
         source_json => fson_value_children_mpi(sources_json)
         do source_index = 0, num_source_specs -1
            if (fson_has_mpi(source_json, "zones")) then
               zones_type = fson_type_mpi(source_json, "zones")
               select case (zones_type)
               case (TYPE_STRING)
                  allocate(zones(1))
                  call fson_get_mpi(source_json, "zones", val = zones(1))
               case (TYPE_ARRAY)
                  call fson_get_mpi(source_json, "zones", &
                       string_length = max_zone_name_length, val = zones)
               end select
               CHKERRQ(ierr)
               associate(num_zones => size(zones))
                 do iz = 1, num_zones
                    label_name = zone_label_name(zones(iz))
                    call DMHasLabel(dm, label_name, has_label, ierr); CHKERRQ(ierr)
                    if (has_label) then
                       call DMGetStratumSize(dm, label_name, 1, num_zone_cells, &
                            ierr); CHKERRQ(ierr)
                       if (num_zone_cells > 0) then
                          call DMGetStratumIS(dm, label_name, 1, cell_IS, &
                               ierr); CHKERRQ(ierr)
                          call ISGetIndicesF90(cell_IS, zone_cells, ierr); CHKERRQ(ierr)
                          do i = 1, num_zone_cells
                             associate(c => zone_cells(i))
                               call DMLabelGetValue(ghost_label, c, ghost, ierr)
                               if (ghost < 0) then
                                  call DMSetLabelValue(dm, source_label_name, &
                                       c, source_index, ierr); CHKERRQ(ierr)
                               end if
                             end associate
                          end do
                          call ISRestoreIndicesF90(cell_IS, zone_cells, ierr); CHKERRQ(ierr)
                          call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)
                       end if
                    else
                       err = 1
                       if (present(logfile)) then
                          call logfile%write(LOG_LEVEL_ERR, "input", "unrecognised zone", &
                               str_key = "name", str_value = zones(iz))
                       end if
                       exit
                    end if
                 end do
               end associate
               deallocate(zones)
            end if

            call mpi_broadcast_error_flag(err)
            if (err > 0) exit

            ! Unlabel any ghost cells:
            call DMGetStratumSize(dm, source_label_name, source_index, &
                 num_cells, ierr); CHKERRQ(ierr)
            if (num_cells > 0) then
               call DMGetStratumIS(dm, source_label_name, source_index, cell_IS, &
                    ierr); CHKERRQ(ierr)
               call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
               labelled_cells = cells ! make copy of indices to check
               call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
               do i = 1, num_cells
                  associate(c => labelled_cells(i))
                    call DMLabelGetValue(ghost_label, c, ghost, ierr)
                    if (ghost < 0) then
                       num_local_sources = num_local_sources + 1
                    else
                       call DMLabelClearValue(source_label, c, source_index, &
                            ierr); CHKERRQ(ierr)
                    end if
                  end associate
               end do
               deallocate(labelled_cells)
            end if

            source_json => fson_value_next_mpi(source_json)

         end do

      end if

    end subroutine label_source_zones

!........................................................................

    subroutine setup_source_dm_data_layout(dm_source)
      !! Sets up data layout on source DM.

      DM, intent(in out) :: dm_source
      ! Locals:
      type(source_type) :: source
      PetscInt :: i, j
      PetscInt, allocatable :: num_field_components(:), field_dim(:)
      PetscInt, parameter :: max_field_name_length = 40
      character(max_field_name_length), allocatable :: field_names(:)

      call source%init("", eos, 0, 0, 0._dp, 1, 1, num_tracers)
      allocate(num_field_components(source%dof), field_dim(source%dof), &
           field_names(source%dof))
      call source%destroy()
      num_field_components = 1
      field_dim = 0
      field_names(1: num_source_scalar_variables) = &
           source_scalar_variable_names ! scalar fields
      i = num_source_scalar_variables + 1
      ! flow fields:
      do j = 1, eos%num_components
         field_names(i) = trim(eos%component_names(j)) // '_' // &
              trim(source_array_variable_names(1))
         i = i + 1
      end do
      if (.not. eos%isothermal) then
         field_names(i) = 'heat_' // &
              trim(source_array_variable_names(1))
         i = i + 1
      end if
      ! tracer injection rates:
      do j = 1, num_tracers
         field_names(i) = trim(tracer_names(j)) // '_' // &
              trim(source_array_variable_names(2))
         i = i + 1
      end do
      ! tracer flow rates:
      do j = 1, num_tracers
         field_names(i) = trim(tracer_names(j)) // '_' // &
              trim(source_array_variable_names(3))
         i = i + 1
      end do
      call dm_set_data_layout(dm_source, num_field_components, field_dim, &
           field_names)
      deallocate(num_field_components, field_dim, field_names)

    end subroutine setup_source_dm_data_layout

!........................................................................

    subroutine setup_source(source_spec_index, local_source_index, source_json, &
         ao, tracer_names, num_sources, separated_sources, thermo, &
         sources, source_dict, err)
      !! Sets up all cell sources for a source specification.

      PetscInt, intent(in) :: source_spec_index !! Index of source specification
      PetscInt, intent(in out) :: local_source_index !! Index of source
      type(fson_value), pointer, intent(in) :: source_json !! JSON input for specification
      AO, intent(in) :: ao !! Application ordering for natural to global cell indexing
      character(*), intent(in) :: tracer_names(:) !! Tracer names
      PetscInt, intent(in out) :: num_sources !! Current total number of sources (on all processes)
      type(list_type), intent(in out) :: separated_sources !! List of sources with separators
      class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics
      type(list_type), intent(in out) :: sources !! List of local sources
      type(dictionary_type), intent(in out) :: source_dict !! Dictionary of local named sources
      PetscErrorCode, intent(out) :: err
      ! Locals:
      type(source_type), pointer :: source
      character(len=64) :: srcstr
      character(len=12) :: istr
      PetscInt :: injection_component, production_component
      PetscReal :: initial_rate, initial_enthalpy, separator_pressure
      PetscInt :: num_cells, num_cells_all, i, source_offset
      PetscInt, allocatable :: natural_source_index(:), natural_cell_index(:)
      type(list_type) :: spec_sources
      IS :: cell_IS
      PetscInt, pointer, contiguous :: local_cell_index(:)
      ISLocalToGlobalMapping :: l2g
      PetscReal :: tracer_injection_rate(size(tracer_names))
      character(max_source_network_node_name_length) :: name

      err = 0
      call DMGetLocalToGlobalMapping(dm, l2g, ierr); CHKERRQ(ierr)
      write(istr, '(i0)') source_spec_index
      srcstr = 'source[' // trim(istr) // '].'

      call get_components(source_json, eos, &
           injection_component, production_component, logfile)
      call get_initial_rate(source_json, initial_rate)
      call get_initial_enthalpy(source_json, eos, &
           injection_component, initial_enthalpy)
      call get_separator_pressure(source_json, srcstr, separator_pressure, logfile)
      call get_tracer_injection_rate(source_json, tracer_names, &
           srcstr, tracer_injection_rate, logfile, err)

      if (err == 0) then

         if (fson_has_mpi(source_json, "name")) then
            call fson_get_mpi(source_json, "name", val = name)
         else
            name = ""
         end if

         call DMGetStratumSize(dm, source_label_name, source_spec_index, &
              num_cells, ierr); CHKERRQ(ierr)
         call MPI_allreduce(num_cells, num_cells_all, 1, MPI_INTEGER, MPI_SUM, &
              PETSC_COMM_WORLD, ierr)
         if (num_cells > 0) then
            call DMGetStratumIS(dm, source_label_name, source_spec_index, &
                 cell_IS, ierr); CHKERRQ(ierr)
            call ISGetIndicesF90(cell_IS, local_cell_index, ierr); CHKERRQ(ierr)
            allocate(natural_cell_index(num_cells))
            natural_cell_index = local_to_natural_cell_index(ao, l2g, local_cell_index)
         else
            allocate(natural_cell_index(1))
         end if
         natural_source_index = get_natural_source_indices(num_cells, num_cells_all, &
              natural_cell_index)

         call spec_sources%init(owner = PETSC_FALSE)
         if (num_cells > 0) then
            do i = 1, num_cells
               source_offset = global_section_offset(source_section, local_source_index, &
                    source_range_start)
               allocate(source)
               call source%init(name, eos, local_source_index, local_cell_index(i), &
                    initial_enthalpy, injection_component, production_component, &
                    num_tracers)
               call source%assign(source_data, source_offset)
               call source%setup(natural_source_index(i), natural_cell_index(i), &
                    initial_rate, tracer_injection_rate, separator_pressure, thermo)
               call sources%append(source)
               call spec_sources%append(source)
               if (separator_pressure > 0._dp) then
                  call separated_sources%append(source)
               end if
               local_source_index = local_source_index + 1
            end do
            call ISRestoreIndicesF90(cell_IS, local_cell_index, ierr); CHKERRQ(ierr)
            deallocate(natural_cell_index, natural_source_index)
         end if
         call setup_inline_source_controls(source_json, eos, thermo, &
              start_time, source_data, source_section, source_range_start, &
              fluid_data, fluid_section, fluid_range_start, srcstr, &
              tracer_names, spec_sources, source_controls, logfile, err)

         if ((name /= "") .and. (num_cells == 1) .and. (num_cells_all == 1)) then
            ! Uniquely named source- add to dictionary:
            call source_dict%add(name, spec_sources%head%data)
         end if

         call spec_sources%destroy()
         num_sources = num_sources + num_cells_all

      end if

    end subroutine setup_source

!........................................................................

    function get_natural_source_indices(num_cells, num_cells_all, &
         natural_cell_index) result(natural_source_index)
      !! Gets natural source indices for the current source.

      use utils_module, only: array_cumulative_sum, array_indices_in_int_array
      use mpi_utils_module, only: get_mpi_int_gather_array

      PetscInt, intent(in) :: num_cells !! Number of source cells on current process
      PetscInt, intent(in) :: num_cells_all !! Number of source cells on all processes
      PetscInt, intent(in) :: natural_cell_index(:) !! Natural cell indices on current process
      PetscInt, allocatable :: natural_source_index(:)
      ! Locals:
      PetscMPIInt :: rank, num_procs
      PetscInt, allocatable :: cell_counts(:), cell_displacements(:)
      PetscInt, allocatable :: natural_cell_index_all(:)
      PetscInt, allocatable :: ordered_natural_cell_index_all(:)
      PetscInt, allocatable :: natural_source_index_all(:)
      PetscErrorCode :: ierr

      call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
      call MPI_COMM_SIZE(PETSC_COMM_WORLD, num_procs, ierr)

      cell_counts = get_mpi_int_gather_array()
      cell_displacements = get_mpi_int_gather_array()
      call MPI_gather(num_cells, 1, MPI_INTEGER, cell_counts, 1, &
           MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      if (rank == 0) then
         allocate(natural_cell_index_all(num_cells_all), &
              natural_source_index_all(num_cells_all))
         cell_displacements = [[0], &
              array_cumulative_sum(cell_counts(1: num_procs - 1))]
      else
         allocate(natural_cell_index_all(1), natural_source_index_all(1))
      end if

      if (num_cells > 0) then
         allocate(natural_source_index(num_cells))
      else
         allocate(natural_source_index(1))
      end if

      call MPI_gatherv(natural_cell_index, num_cells, MPI_INTEGER, &
           natural_cell_index_all, cell_counts, cell_displacements, &
           MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      ordered_natural_cell_index_all = get_source_cell_indices(source_json, &
           natural_cell_index_all)

      if (rank == 0) then
         natural_source_index_all = array_indices_in_int_array( &
              ordered_natural_cell_index_all, natural_cell_index_all)
         natural_source_index_all = natural_source_index_all + num_sources - 1
      end if
      call MPI_scatterv(natural_source_index_all, cell_counts, cell_displacements, &
           MPI_INTEGER, natural_source_index, num_cells, MPI_INTEGER, &
           0, PETSC_COMM_WORLD, ierr)

      deallocate(cell_counts, cell_displacements, natural_source_index_all, &
           natural_cell_index_all, ordered_natural_cell_index_all)

    end function get_natural_source_indices

!........................................................................

    function get_source_cell_indices(source_json, cells) result(ordered_cells)

      !! Gets array of all natural cell indices for the source, in
      !! natural source order, on rank 0. For sources with the "cell"
      !! or "cells" properties this can simply be read from the
      !! input. For sources with cells defined by zones, cells from
      !! all ranks are ordered by natural cell index.

      type(fson_value), pointer, intent(in) :: source_json
      PetscInt, intent(in) :: cells(:)
      PetscInt, allocatable :: ordered_cells(:)
      ! Locals:
      type(fson_value), pointer :: cell_json, cells_json, zones_json
      PetscInt :: c
      PetscMPIInt :: rank
      PetscErrorCode :: ierr

      call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

      if (rank == 0) then

         call fson_get(source_json, "cell", cell_json)
         if (associated(cell_json)) then
            call fson_get(cell_json, ".", c)
            ordered_cells = [c]
         end if

         call fson_get(source_json, "cells", cells_json)
         if (associated(cells_json)) then
            select case (cells_json%value_type)
            case (TYPE_INTEGER)
               call fson_get(cells_json, ".", c)
               ordered_cells = [c]
            case (TYPE_ARRAY)
               call fson_get(cells_json, ".", ordered_cells)
            end select
         end if

         call fson_get(source_json, "zones", zones_json)
         if (associated(zones_json)) then
            ordered_cells = cells
            call PetscSortInt(size(cells), ordered_cells, ierr); CHKERRQ(ierr)
         end if

      else
         allocate(ordered_cells(1))
      end if

    end function get_source_cell_indices

!........................................................................

    subroutine source_indices_iterator(node, stopped)
      !! Gets indices from all local sources.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: s, source_offset

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         s = source%local_source_index
         source_offset = global_section_offset(source_section, s, &
              source_range_start)
         call source%assign(source_data, source_offset)
         indices(s + 1) = nint(source%source_index)
      end select

    end subroutine source_indices_iterator

!........................................................................

    subroutine setup_source_index(indices, source_index)
      !! Sets up natural-to-global source ordering IS. This gives the
      !! global index corresponding to a natural source index.

      use utils_module, only: array_cumulative_sum
      use mpi_utils_module, only: get_mpi_int_gather_array

      PetscInt, intent(in) :: indices(num_local_sources)
      IS, intent(in out) :: source_index
      ! Locals:
      PetscInt :: i
      PetscMPIInt :: rank, num_procs, num_all, is_count
      PetscInt, allocatable :: counts(:), displacements(:)
      PetscInt, allocatable :: indices_all(:), global_indices(:)
      PetscErrorCode :: ierr

      call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
      call MPI_COMM_SIZE(PETSC_COMM_WORLD, num_procs, ierr)

      counts = get_mpi_int_gather_array()
      displacements = get_mpi_int_gather_array()
      call MPI_gather(num_local_sources, 1, MPI_INTEGER, counts, 1, &
           MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      if (rank == 0) then
         displacements = [[0], &
              array_cumulative_sum(counts(1: num_procs - 1))]
         num_all = sum(counts)
         is_count = num_all
      else
         num_all = 1
         is_count = 0
      end if
      allocate(indices_all(0: num_all - 1), global_indices(0: num_all - 1))
      global_indices = -1
      call MPI_gatherv(indices, num_local_sources, MPI_INTEGER, &
           indices_all, counts, displacements, &
           MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr)
      if (rank == 0) then
         do i = 0, num_all - 1
            global_indices(indices_all(i)) = i
         end do
      end if
      deallocate(indices_all, counts, displacements)
      call ISCreateGeneral(PETSC_COMM_WORLD, is_count, &
           global_indices, PETSC_COPY_VALUES, source_index, ierr)
      CHKERRQ(ierr)
      call PetscObjectSetName(source_index, "source_index", ierr)
      deallocate(global_indices)

    end subroutine setup_source_index

!........................................................................

    subroutine setup_group_iterator(node, stopped)
      !! Sets up source group.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      character(max_source_network_node_name_length) :: name
      character(max_source_network_node_name_length), allocatable :: names(:)
      PetscInt :: i
      type(list_node_type), pointer :: dict_node
      type(source_group_type), pointer :: group

      stopped = PETSC_FALSE
      select type (group_json => node%data)
      type is (fson_value)

         call fson_get_mpi(group_json, "name", &
              "untitled", name)
         call fson_get_mpi(group_json, "source", &
              string_length = max_source_network_node_name_length, &
              val = names)

         allocate(group)
         call group%init(name)

         associate(num_names => size(names))
           do i = 1, num_names
              dict_node => source_dict%get(names(i))
              if (associated(dict_node)) then
                 select type (node => dict_node%data)
                 type is (source_type)
                    call group%nodes%append(node)
                 end select
              end if
           end do
         end associate

         call source_groups%append(group)

      end select

    end subroutine setup_group_iterator

!........................................................................

  end subroutine setup_sources

!------------------------------------------------------------------------

  PetscInt function get_component(source_json, eos, name) result(component)
    !! Gets a specified component type for the source.

    type(fson_value), pointer, intent(in) :: source_json
    class(eos_type), intent(in) :: eos
    character(len=*), intent(in) :: name
    ! Locals:
    character(max_component_name_length) :: component_str
    PetscInt :: component_type

    component_type = fson_type_mpi(source_json, name)
    if (component_type == TYPE_INTEGER) then
       call fson_get_mpi(source_json, name, val = component)
    else if (component_type == TYPE_STRING) then
       call fson_get_mpi(source_json, name, val = component_str)
       component = eos%component_index(component_str)
    end if

  end function get_component

!------------------------------------------------------------------------

  subroutine get_components(source_json, eos, &
       injection_component, production_component, logfile)
    !! Gets injection and production components for the source.

    type(fson_value), pointer, intent(in) :: source_json
    class(eos_type), intent(in) :: eos
    PetscInt, intent(out) :: injection_component
    PetscInt, intent(out) :: production_component
    type(logfile_type), intent(in out), optional :: logfile

    if (fson_has_mpi(source_json, "component")) then
       injection_component = get_component(source_json, eos, &
            "component")
    else
       injection_component = default_source_component
    end if

    if (fson_has_mpi(source_json, "production_component")) then
       production_component = get_component(source_json, eos, &
            "production_component")
    else
       associate(np => eos%num_primary_variables)
         if ((.not. eos%isothermal) .and. (injection_component == np)) then
            production_component = injection_component
         else
            production_component = default_source_production_component
         end if
       end associate
    end if

  end subroutine get_components

!------------------------------------------------------------------------

  subroutine get_initial_rate(source_json, initial_rate)
    !! Gets initial flow rate. This can only be determined for
    !! constant-rate sources. For other source types, a default
    !! initial rate is assigned, which may be modified by any source
    !! controls acting on the source.

    type(fson_value), pointer, intent(in) :: source_json
    PetscReal, intent(out) :: initial_rate
    ! Locals:
    PetscInt :: rate_type

    if (fson_has_mpi(source_json, "rate")) then

       rate_type = fson_type_mpi(source_json, "rate")

       select case (rate_type)
       case (TYPE_REAL, TYPE_INTEGER)
          call fson_get_mpi(source_json, "rate", val = initial_rate)
       case default
          initial_rate = default_source_rate
       end select

    else
       initial_rate = default_source_rate
    end if

  end subroutine get_initial_rate

!------------------------------------------------------------------------

  subroutine get_initial_enthalpy(source_json, eos, injection_component, &
       enthalpy)
    !! Gets initial injection enthalpy for the source, if needed.

    type(fson_value), pointer, intent(in) :: source_json
    class(eos_type), intent(in) :: eos
    PetscInt, intent(in) :: injection_component
    PetscReal, intent(out) :: enthalpy
    ! Locals:
    PetscInt :: enthalpy_type

    associate(np => eos%num_primary_variables)

      if (injection_component < np) then

         if (fson_has_mpi(source_json, "enthalpy")) then

            enthalpy_type = fson_type_mpi(source_json, "enthalpy")

            select case (enthalpy_type)
            case(TYPE_REAL, TYPE_INTEGER)
               call fson_get_mpi(source_json, "enthalpy", &
                    val = enthalpy)
            case default
               enthalpy = default_source_injection_enthalpy
            end select

         else
            enthalpy = default_source_injection_enthalpy
         end if

      else ! heat injection - no enthalpy needed
         enthalpy = 0._dp
      end if

    end associate

  end subroutine get_initial_enthalpy

!------------------------------------------------------------------------

  subroutine get_tracer_injection_rate(source_json, tracer_names, &
       srcstr, tracer_injection_rate, logfile, err)
    !! Gets constant tracer injection rates.

    use tracer_module, only: max_tracer_name_length
    use utils_module, only: str_array_index

    type(fson_value), pointer, intent(in) :: source_json
    character(*), intent(in) :: tracer_names(:) !! Tracer names
    character(*), intent(in) :: srcstr !! source identifier
    PetscReal, intent(out) :: tracer_injection_rate(:)
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal, allocatable :: injection_rate(:)
    PetscReal :: scalar_injection_rate
    PetscInt :: tracer_json_type, tracer_type, rank
    PetscInt :: num_tracers_specified, i, tracer_index
    type(fson_value), pointer :: tracers_json, tracer_json
    character(max_tracer_name_length) :: name
    PetscReal, parameter :: default_tracer_injection_rate = 0._dp

    err = 0

    if (size(tracer_names) > 0) then

       if (fson_has_mpi(source_json, "tracer")) then
          tracer_json_type = fson_type_mpi(source_json, "tracer")
          select case (tracer_json_type)
          case (TYPE_REAL, TYPE_INTEGER)
             ! apply same value to all tracers:
             call fson_get_mpi(source_json, "tracer", val = scalar_injection_rate)
             tracer_injection_rate = scalar_injection_rate
          case (TYPE_ARRAY)
            rank = fson_mpi_array_rank(source_json, "tracer")
            if (rank == 1) then
               ! array of values for different tracers:
               call fson_get_mpi(source_json, "tracer", val = injection_rate)
               tracer_injection_rate = default_tracer_injection_rate
               tracer_injection_rate(1: size(injection_rate)) = &
                    injection_rate
            end if
          case (TYPE_OBJECT)
             ! values specified by tracer name:
             tracer_injection_rate = default_tracer_injection_rate
             call fson_get_mpi(source_json, "tracer", tracers_json)
             num_tracers_specified = fson_value_count_mpi(tracers_json, ".")
             tracer_json => fson_value_children_mpi(tracers_json)
             do i = 1, num_tracers_specified
                tracer_type = fson_type_mpi(tracer_json, ".")
                select case (tracer_type)
                case (TYPE_REAL, TYPE_INTEGER)
                   name = fson_get_name_mpi(tracer_json)
                   call fson_get_mpi(tracer_json, ".", val = scalar_injection_rate)
                   tracer_index = str_array_index(name, tracer_names)
                   if (tracer_index > 0) then
                      tracer_injection_rate(tracer_index) = &
                           scalar_injection_rate
                   else
                      call logfile%write(LOG_LEVEL_ERR, "input", &
                           "unrecognised_tracer", &
                           str_key = trim(srcstr) // "name", &
                           str_value = name)
                      err = 1
                      exit
                   end if
                end select
                tracer_json => fson_value_next_mpi(tracer_json)
             end do
          case default
             tracer_injection_rate = default_tracer_injection_rate
          end select
       else
          tracer_injection_rate = default_tracer_injection_rate
       end if
    end if

  end subroutine get_tracer_injection_rate

!------------------------------------------------------------------------

  subroutine get_separator_pressure(source_json, srcstr, separator_pressure, &
       logfile)
    !! Gets separator pressure. Returns -1 if no separator is
    !! specified.

    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: source_json
    character(*), intent(in) :: srcstr !! source identifier
    PetscReal, intent(out) :: separator_pressure !! Separator pressure
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscInt :: separator_json_type
    PetscBool :: has_separator
    character(8) :: limiter_type_str

    separator_pressure = -1._dp

    if (fson_has_mpi(source_json, "separator")) then
       separator_json_type = fson_type_mpi(source_json, "separator")
       select case (separator_json_type)
       case (TYPE_LOGICAL)
          call fson_get_mpi(source_json, "separator", val = has_separator)
          if (has_separator) then
             separator_pressure = default_separator_pressure
             if (present(logfile) .and. logfile%active) then
                call logfile%write(LOG_LEVEL_INFO, 'input', 'default', real_keys = &
                     [trim(srcstr) // "separator.pressure"], &
                     real_values = [separator_pressure])
             end if
          end if
       case (TYPE_OBJECT)
          call fson_get_mpi(source_json, "separator.pressure", &
               default_separator_pressure, separator_pressure, logfile, srcstr)
       end select
    else if (fson_has_mpi(source_json, "limiter")) then
       call fson_get_mpi(source_json, "limiter.type", &
            default_source_control_limiter_type_str, limiter_type_str)
       limiter_type_str = str_to_lower(limiter_type_str)
       if (limiter_type_str /= "total") then
          call fson_get_mpi(source_json, "limiter.separator_pressure", &
               default_separator_pressure, separator_pressure, logfile, srcstr)
       end if
    end if

  end subroutine get_separator_pressure

!------------------------------------------------------------------------

  subroutine setup_inline_source_controls(source_json, eos, thermo, &
       start_time, source_data, source_section, source_range_start, &
       fluid_data, fluid_section, fluid_range_start, &
       srcstr, tracer_names, spec_sources, source_controls, &
       logfile, err)
    !! Sets up any 'inline' source controls for the source specification,
    !! i.e. controls defined implicitly in the specification.

    use interpolation_module, only: interpolation_type_from_str, &
         averaging_type_from_str, max_interpolation_str_length, &
         max_averaging_str_length, default_interpolation_str, &
         default_averaging_str

    type(fson_value), pointer, intent(in) :: source_json
    class(eos_type), intent(in) :: eos
    class(thermodynamics_type), intent(in) :: thermo
    PetscReal, intent(in) :: start_time
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start
    character(len = *), intent(in) :: srcstr
    character(len = *), intent(in) :: tracer_names(:)
    type(list_type), intent(in out) :: spec_sources
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: interpolation_type, averaging_type
    character(max_interpolation_str_length) :: interpolation_str
    character(max_averaging_str_length) :: averaging_str

    call fson_get_mpi(source_json, "interpolation", &
         default_interpolation_str, interpolation_str)
    interpolation_type = interpolation_type_from_str(interpolation_str)

    call fson_get_mpi(source_json, "averaging", &
         default_averaging_str, averaging_str)
    averaging_type = averaging_type_from_str(averaging_str)

    call setup_table_source_control(source_json, srcstr, interpolation_type, &
         averaging_type, tracer_names, spec_sources, source_controls, &
         logfile, err)

    if (err == 0) then

       call setup_deliverability_source_controls(source_json, srcstr, &
            start_time, source_data, source_section, source_range_start, &
            fluid_data, fluid_section, fluid_range_start, &
            interpolation_type, averaging_type, spec_sources, eos, &
            source_controls, logfile, err)

       if (err == 0) then

          call setup_recharge_source_controls(source_json, srcstr, &
               source_data, source_section, source_range_start, &
               fluid_data, fluid_section, fluid_range_start, &
               interpolation_type, averaging_type, spec_sources, eos, &
               source_controls, logfile, err)

          if (err == 0) then

             call setup_limiter_source_controls(source_json, srcstr, thermo, &
                  spec_sources, source_controls, logfile)

             call setup_direction_source_control(source_json, srcstr, thermo, &
                  spec_sources, source_controls, logfile)

             call setup_factor_source_control(source_json, srcstr, &
                  interpolation_type, averaging_type, spec_sources, &
                  source_controls, logfile, err)

          end if

       end if

    end if

  end subroutine setup_inline_source_controls

!------------------------------------------------------------------------

  subroutine setup_table_source_control(source_json, srcstr, &
       interpolation_type, averaging_type, tracer_names, &
       spec_sources, source_controls, logfile, err)
    !! Set up rate, enthalpy and tracer table source controls.

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    PetscInt, intent(in) :: interpolation_type, averaging_type
    character(*), intent(in) :: tracer_names(:)
    type(list_type), intent(in out) :: spec_sources
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err

    call setup_rate_table_control()
    if (err == 0) then
       call setup_enthalpy_table_control()
       if (err == 0) then
          call setup_tracer_table_controls()
       end if
    end if

  contains

!........................................................................

    subroutine setup_rate_table_control()

      PetscInt :: variable_type
      type(fson_value), pointer :: table
      type(source_control_rate_table_type), pointer :: control
      PetscReal, allocatable :: data_array(:,:)

      if (fson_has_mpi(source_json, "rate")) then
         variable_type = fson_type_mpi(source_json, "rate")
         if (variable_type == TYPE_ARRAY) then
            call fson_get_mpi(source_json, "rate", val = data_array)
         else if (variable_type == TYPE_OBJECT) then
            call fson_get_mpi(source_json, "rate", table)
            if (fson_has_mpi(table, "time")) then
               call fson_get_mpi(table, "time", val = data_array)
            end if
         end if
      end if

      if (allocated(data_array)) then
         if (spec_sources%count > 0) then
            allocate(control)
            call control%init(data_array, interpolation_type, &
                 averaging_type, spec_sources%copy())
            call source_controls%append(control)
         end if
         deallocate(data_array)
      end if
      
      if (err > 0) then
         call logfile%write(LOG_LEVEL_ERR, "input", "unsorted_array", &
              real_array_key = trim(srcstr) // "rate", &
              real_array_value = data_array(:, 1))
      end if

    end subroutine setup_rate_table_control

!........................................................................

    subroutine setup_enthalpy_table_control()

      PetscInt :: variable_type
      type(fson_value), pointer :: table
      type(source_control_enthalpy_table_type), pointer :: control
      PetscReal, allocatable :: data_array(:,:)

      if (fson_has_mpi(source_json, "enthalpy")) then
         variable_type = fson_type_mpi(source_json, "enthalpy")
         if (variable_type == TYPE_ARRAY) then
            call fson_get_mpi(source_json, "enthalpy", val = data_array)
         else if (variable_type == TYPE_OBJECT) then
            call fson_get_mpi(source_json, "enthalpy", table)
            if (fson_has_mpi(table, "time")) then
               call fson_get_mpi(table, "time", val = data_array)
            end if
         end if
      end if

      if (allocated(data_array)) then
         if (spec_sources%count > 0) then
            allocate(control)
            call control%init(data_array, interpolation_type, &
                 averaging_type, spec_sources%copy())
            call source_controls%append(control)
         end if
         deallocate(data_array)
      end if

      if (err > 0) then
         call logfile%write(LOG_LEVEL_ERR, "input", "unsorted_array", &
              real_array_key = trim(srcstr) // "enthalpy", &
              real_array_value = data_array(:, 1))
      end if

    end subroutine setup_enthalpy_table_control

!........................................................................

    subroutine setup_tracer_table_controls()

      use utils_module, only: str_array_index
      use tracer_module, only: max_tracer_name_length

      PetscInt :: variable_type, rank
      type(fson_value), pointer :: tracers_json, tracer_json
      PetscInt :: num_tracers_specified, i, tracer_type, tracer_index
      type(source_control_tracer_table_type), pointer :: control
      PetscReal, allocatable :: data_array(:,:)
      character(max_tracer_name_length) :: name

      if (fson_has_mpi(source_json, "tracer")) then

         variable_type = fson_type_mpi(source_json, "tracer")
         select case (variable_type)
         case (TYPE_ARRAY)

            rank = fson_mpi_array_rank(source_json, "tracer")
            if (rank == 2) then
               call fson_get_mpi(source_json, "tracer", val = data_array)
               if (allocated(data_array)) then
                  if (spec_sources%count > 0) then
                     do i = 1, size(tracer_names)
                        allocate(control)
                        call control%init(data_array, interpolation_type, &
                             averaging_type, spec_sources%copy())
                        control%tracer_index = i
                        call source_controls%append(control)
                     end do
                  end if
                  deallocate(data_array)
               end if
            end if

         case (TYPE_OBJECT)

            call fson_get_mpi(source_json, "tracer", tracers_json)
            num_tracers_specified = fson_value_count_mpi(tracers_json, ".")
            tracer_json => fson_value_children_mpi(tracers_json)
            do i = 1, num_tracers_specified
               tracer_type = fson_type_mpi(tracer_json, ".")
               select case (tracer_type)
               case (TYPE_ARRAY)
                  rank = fson_mpi_array_rank(tracer_json, ".")
                  if (rank == 2) then
                     name = fson_get_name_mpi(tracer_json)
                     call fson_get_mpi(tracer_json, ".", val = data_array)
                     if (allocated(data_array)) then
                        tracer_index = str_array_index(name, tracer_names)
                        if (tracer_index > 0) then
                           allocate(control)
                           call control%init(data_array, interpolation_type, &
                                averaging_type, spec_sources%copy())
                           control%tracer_index = tracer_index
                           call source_controls%append(control)
                        else
                           call logfile%write(LOG_LEVEL_ERR, "input", &
                                "unrecognised_tracer", &
                                str_key = trim(srcstr) // "name", &
                                str_value = name)
                           err = 1
                           exit
                        end if
                        deallocate(data_array)
                     end if
                  end if
               end select
               tracer_json => fson_value_next_mpi(tracer_json)
            end do
         end select
      end if

    end subroutine setup_tracer_table_controls

  end subroutine setup_table_source_control

!------------------------------------------------------------------------

  subroutine setup_factor_source_control(source_json, srcstr, &
       interpolation_type, averaging_type, spec_sources, &
       source_controls, logfile, err)
    !! Set up rate factor source controls.

    use interpolation_module, only: interpolation_type_from_str, &
         averaging_type_from_str, max_interpolation_str_length, &
         max_averaging_str_length

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    PetscInt, intent(in) :: interpolation_type, averaging_type
    type(list_type), intent(in out) :: spec_sources
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(source_control_rate_factor_type), pointer :: factor_control
    PetscInt :: variable_type
    PetscReal, allocatable :: factor_data_array(:,:)
    type(fson_value), pointer :: table
    PetscInt :: effective_interpolation_type, effective_averaging_type
    character(max_interpolation_str_length) :: interpolation_str
    character(max_averaging_str_length) :: averaging_str
    PetscReal :: const_factor

    ! Use interpolation/ averaging parameters from source itself by default:
    effective_interpolation_type = interpolation_type
    effective_averaging_type = averaging_type

    if (fson_has_mpi(source_json, "factor")) then
       variable_type = fson_type_mpi(source_json, "factor")
       select case (variable_type)
       case (TYPE_REAL, TYPE_INTEGER)
          call fson_get_mpi(source_json, "factor", val = const_factor)
          allocate(factor_data_array(1, 2))
          factor_data_array(1, :) = [0._dp, const_factor]
       case (TYPE_ARRAY)
          call fson_get_mpi(source_json, "factor", val = factor_data_array)
       case (TYPE_OBJECT)
          call fson_get_mpi(source_json, "factor", table)
          if (fson_has_mpi(table, "time")) then
             call fson_get_mpi(table, "time", val = factor_data_array)
          end if
          ! Check for overridden interpolation/ averaging parameters:
          if (fson_has_mpi(source_json, "factor.interpolation")) then
             call fson_get_mpi(source_json, "factor.interpolation", &
               val = interpolation_str)
             effective_interpolation_type = interpolation_type_from_str(interpolation_str)
          end if
          if (fson_has_mpi(source_json, "factor.averaging")) then
             call fson_get_mpi(source_json, "factor.averaging", &
               val = averaging_str)
             effective_averaging_type = averaging_type_from_str(averaging_str)
          end if
       end select
    end if

    if (allocated(factor_data_array) .and. (spec_sources%count > 0)) then
       allocate(factor_control)
       call factor_control%init(factor_data_array, effective_interpolation_type, &
            effective_averaging_type, spec_sources%copy())
       call source_controls%append(factor_control)
    end if

    if (allocated(factor_data_array)) deallocate(factor_data_array)

  end subroutine setup_factor_source_control

!------------------------------------------------------------------------

  subroutine get_reference_pressure(json, srcstr, source_type, &
       reference_pressure_array, calculate_reference_pressure, &
       pressure_table_coordinate, logfile)
    !! Get reference pressure for deliverability or recharge source control.

    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: json
    character(len = *), intent(in) :: srcstr, source_type
    PetscReal, allocatable, intent(out) :: reference_pressure_array(:,:)
    PetscBool, intent(out) :: calculate_reference_pressure
    PetscInt, intent(out) :: pressure_table_coordinate
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscReal :: reference_pressure
    PetscInt :: pressure_type
    PetscInt, parameter :: max_pressure_str_length = 8
    character(max_pressure_str_length) :: pressure_str
    type(fson_value), pointer :: pressure_json
    PetscReal, parameter :: default_time = 0._dp

    calculate_reference_pressure = PETSC_FALSE
    reference_pressure = default_deliverability_reference_pressure
    reference_pressure_array = reshape([default_time, &
         reference_pressure], [1,2])
    pressure_table_coordinate = default_source_pressure_table_coordinate

    if (fson_has_mpi(json, "pressure")) then
       pressure_type = fson_type_mpi(json, "pressure")
       select case (pressure_type)
       case (TYPE_REAL, TYPE_INTEGER)
          call fson_get_mpi(json, "pressure", &
               val = reference_pressure)
          reference_pressure_array(1,2) = reference_pressure
          pressure_table_coordinate = SRC_PRESSURE_TABLE_COORD_TIME
       case (TYPE_ARRAY)
          call fson_get_mpi(json, "pressure", &
               val = reference_pressure_array)
          pressure_table_coordinate = SRC_PRESSURE_TABLE_COORD_TIME
       case (TYPE_OBJECT)
          call fson_get_mpi(json, "pressure", pressure_json)
          if (fson_has_mpi(pressure_json, "time")) then
             call fson_get_mpi(pressure_json, "time", &
                  val = reference_pressure_array)
             pressure_table_coordinate = SRC_PRESSURE_TABLE_COORD_TIME
          else if (fson_has_mpi(pressure_json, "enthalpy")) then
             call fson_get_mpi(pressure_json, "enthalpy", &
                  val = reference_pressure_array)
             pressure_table_coordinate = SRC_PRESSURE_TABLE_COORD_ENTHALPY
          end if
       case (TYPE_STRING)
          call fson_get_mpi(json, "pressure", &
               val = pressure_str)
          if (trim(str_to_lower(pressure_str)) == 'initial') then
             calculate_reference_pressure = PETSC_TRUE
             pressure_table_coordinate = SRC_PRESSURE_TABLE_COORD_TIME
          else
             if (present(logfile) .and. logfile%active) then
                call logfile%write(LOG_LEVEL_WARN, 'input', 'unrecognised', &
                     str_key = trim(srcstr) // source_type // ".pressure", &
                     str_value = pressure_str)
             end if
          end if
       end select
    else
       if (present(logfile) .and. logfile%active) then
          call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
               real_keys = [trim(srcstr) // source_type // ".pressure"], &
               real_values = [reference_pressure])
       end if
    end if

  end subroutine get_reference_pressure

!------------------------------------------------------------------------
  
  subroutine get_deliverability_productivity(json, source_json, srcstr, &
       productivity_array, calculate_PI_from_rate, logfile)
    !! Gets productivity index for deliverability source control.

    type(fson_value), pointer, intent(in) :: json, source_json
    character(len = *), intent(in) :: srcstr
    PetscReal, allocatable, intent(out) :: productivity_array(:,:)
    PetscBool, intent(out) :: calculate_PI_from_rate
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscInt :: PI_type
    PetscReal :: productivity
    type(fson_value), pointer :: PI_json
    PetscReal, parameter :: default_time = 0._dp

    calculate_PI_from_rate = PETSC_FALSE
    productivity = default_deliverability_productivity
    productivity_array = reshape([default_time, productivity], [1,2])

    if (fson_has_mpi(json, "productivity")) then

       PI_type = fson_type_mpi(json, "productivity")

       select case(PI_type)
       case (TYPE_REAL, TYPE_INTEGER)
          call fson_get_mpi(json, "productivity", val = productivity)
          productivity_array = reshape([default_time, &
               productivity], [1,2])
       case (TYPE_ARRAY)
          call fson_get_mpi(json, "productivity", &
               val = productivity_array)
       case (TYPE_OBJECT)
          call fson_get_mpi(json, "productivity", PI_json)
          if (fson_has_mpi(PI_json, "time")) then
             call fson_get_mpi(PI_json, "time", &
                  val = productivity_array)
          end if
       case default
          if (present(logfile) .and. logfile%active) then
             call logfile%write(LOG_LEVEL_WARN, 'input', 'unrecognised', &
                  str_key = trim(srcstr) // "deliverability.productivity", &
                  str_value = "...")
          end if
       end select

    else if (fson_has_mpi(source_json, "rate")) then

       calculate_PI_from_rate = PETSC_TRUE

    else
       if (present(logfile) .and. logfile%active) then
          call logfile%write(LOG_LEVEL_INFO, 'input', 'default', real_keys = &
               [trim(srcstr) // "deliverability.productivity"], &
               real_values = [productivity])
       end if
    end if

  end subroutine get_deliverability_productivity

!------------------------------------------------------------------------

  subroutine setup_deliverability_source_controls(source_json, srcstr, &
       start_time, source_data, source_section, source_range_start, &
       fluid_data, fluid_section, fluid_range_start, &
       interpolation_type, averaging_type, spec_sources, eos, &
       source_controls, logfile, err)
    !! Set up deliverability source controls. Deliverability controls
    !! can control only one source, so if multiple cells are
    !! specified, multiple corresponding deliverability controls are
    !! created.

    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    PetscReal, intent(in) :: start_time
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start
    PetscInt, intent(in) :: interpolation_type, averaging_type
    type(list_type), intent(in out) :: spec_sources
    class(eos_type), intent(in) :: eos
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(fson_value), pointer :: deliv_json
    PetscReal, allocatable :: reference_pressure_array(:,:)
    type(source_control_deliverability_type), pointer :: deliv
    PetscBool :: calculate_reference_pressure
    PetscBool :: calculate_PI_from_rate
    PetscReal :: initial_rate, threshold
    PetscReal, allocatable :: productivity_array(:,:)
    PetscInt :: pressure_table_coordinate
    PetscReal, parameter :: default_rate = 0._dp
    PetscReal, parameter :: default_threshold = -1._dp

    if (fson_has_mpi(source_json, "deliverability")) then

       call fson_get_mpi(source_json, "deliverability", deliv_json)

       call fson_get_mpi(deliv_json, "threshold", default_threshold, threshold)

       if (threshold <= 0._dp) then
          call fson_get_mpi(source_json, "rate", default_rate, initial_rate)
       end if

       call get_reference_pressure(deliv_json, srcstr, "deliverability", &
            reference_pressure_array, calculate_reference_pressure, &
            pressure_table_coordinate, logfile)

       call get_deliverability_productivity(deliv_json, source_json, &
            srcstr, productivity_array, calculate_PI_from_rate, logfile)

       call spec_sources%traverse(setup_deliverability_iterator)

    end if

  contains

    subroutine setup_deliverability_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      type(list_type) :: single_source

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)

         call single_source%init(owner = PETSC_FALSE)
         call single_source%append(source)

         allocate(deliv)
          call deliv%init(productivity_array, interpolation_type, &
               averaging_type, reference_pressure_array, &
               pressure_table_coordinate, threshold, single_source)

          if (calculate_reference_pressure) then
             call deliv%set_reference_pressure_initial(source_data, &
                  source_section, source_range_start, fluid_data, &
                  fluid_section, fluid_range_start, eos)
          end if
          if (calculate_PI_from_rate) then
             call deliv%calculate_PI_from_rate(start_time, initial_rate, &
                  source_data, source_section, source_range_start, &
                  fluid_data, fluid_section, fluid_range_start, eos, &
                  deliv%productivity%val(1, 1))
          end if
          if (deliv%threshold > 0._dp) then
             deliv%threshold_productivity = &
                  deliv%productivity%interpolate(start_time, 1)
          end if
          call source_controls%append(deliv)

      end select

    end subroutine setup_deliverability_iterator

  end subroutine setup_deliverability_source_controls

!------------------------------------------------------------------------

  subroutine get_recharge_coefficient(json, source_json, srcstr, &
       recharge_array, logfile)
    !! Gets recharge coefficient for recharge source control.

    type(fson_value), pointer, intent(in) :: json, source_json
    character(len = *), intent(in) :: srcstr
    PetscReal, allocatable :: recharge_array(:,:)
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscInt :: coef_type
    PetscReal :: recharge_coefficient
    type(fson_value), pointer :: coef_json
    PetscReal, parameter :: default_time = 0._dp

    recharge_coefficient = default_recharge_coefficient
    recharge_array = reshape([default_time, &
         recharge_coefficient], [1,2])

    if (fson_has_mpi(json, "coefficient")) then

       coef_type = fson_type_mpi(json, "coefficient")

       select case(coef_type)
       case (TYPE_REAL, TYPE_INTEGER)
          call fson_get_mpi(json, "coefficient", &
               val = recharge_coefficient)
          recharge_array(1,2) = recharge_coefficient
       case (TYPE_ARRAY)
          call fson_get_mpi(json, "coefficient", &
               val = recharge_array)
       case (TYPE_OBJECT)
          call fson_get_mpi(json, "coefficient", coef_json)
          if (fson_has_mpi(coef_json, "time")) then
             call fson_get_mpi(coef_json, "time", &
                  val = recharge_array)
          end if
       case default
          if (present(logfile) .and. logfile%active) then
             call logfile%write(LOG_LEVEL_WARN, 'input', 'unrecognised', &
                  str_key = trim(srcstr) // "recharge.coefficient", &
                  str_value = "...")
          end if
       end select

    else
       if (present(logfile) .and. logfile%active) then
          call logfile%write(LOG_LEVEL_INFO, 'input', 'default', real_keys = &
               [trim(srcstr) // "recharge.coefficient"], &
               real_values = [recharge_coefficient])
       end if
    end if

  end subroutine get_recharge_coefficient

!------------------------------------------------------------------------

  subroutine setup_recharge_source_controls(source_json, srcstr, &
       source_data, source_section, source_range_start, &
       fluid_data, fluid_section, fluid_range_start, &
       interpolation_type, averaging_type, spec_sources, eos, &
       source_controls, logfile, err)
    !! Set up recharge source controls. Recharge controls
    !! can control only one source, so if multiple cells are
    !! specified, multiple corresponding recharge controls are
    !! created.

    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscInt, intent(in) :: source_range_start
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start
    PetscInt, intent(in) :: interpolation_type, averaging_type
    type(list_type), intent(in out) :: spec_sources
    class(eos_type), intent(in) :: eos
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(fson_value), pointer :: recharge_json
    PetscReal, allocatable :: reference_pressure_array(:,:)
    type(source_control_recharge_type), pointer :: recharge
    PetscBool :: calculate_reference_pressure
    PetscReal, allocatable :: recharge_array(:,:)
    PetscInt :: pressure_table_coordinate

    if (fson_has_mpi(source_json, "recharge")) then

       call fson_get_mpi(source_json, "recharge", recharge_json)

       call get_reference_pressure(recharge_json, srcstr, "recharge", &
            reference_pressure_array, calculate_reference_pressure, &
            pressure_table_coordinate, logfile)

       if (pressure_table_coordinate == SRC_PRESSURE_TABLE_COORD_TIME) then

          call get_recharge_coefficient(recharge_json, source_json, &
               srcstr, recharge_array, logfile)

          call spec_sources%traverse(setup_recharge_iterator)

       else
          if (present(logfile) .and. logfile%active) then
             call logfile%write(LOG_LEVEL_WARN, 'input', 'not_supported', &
                  str_key = trim(srcstr) // "recharge.pressure", &
                  str_value = "...")
          end if
       end if

    end if

  contains

    subroutine setup_recharge_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      type(list_type) :: single_source

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)

         call single_source%init(owner = PETSC_FALSE)
         call single_source%append(source)

         allocate(recharge)
         call recharge%init(recharge_array, interpolation_type, &
              averaging_type, reference_pressure_array, single_source)

         if (calculate_reference_pressure) then
            call recharge%set_reference_pressure_initial(source_data, &
                 source_section, source_range_start, fluid_data, &
                 fluid_section, fluid_range_start, eos)
         end if

         call source_controls%append(recharge)

      end select

    end subroutine setup_recharge_iterator

  end subroutine setup_recharge_source_controls

!------------------------------------------------------------------------

  subroutine setup_limiter_source_controls(source_json, srcstr, &
       thermo, spec_sources, source_controls, logfile)
    !! Set up limiter source control for each cell source.

    use utils_module, only: str_to_lower, str_array_index

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    class(thermodynamics_type), intent(in) :: thermo
    type(list_type), intent(in out) :: spec_sources
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    type(fson_value), pointer :: limiter_json
    character(8) :: limiter_type_str
    PetscReal :: limit
    class(source_control_limiter_type), pointer :: limiter

    if (fson_has_mpi(source_json, "limiter")) then

       call fson_get_mpi(source_json, "limiter", limiter_json)

       call fson_get_mpi(limiter_json, "type", &
            default_source_control_limiter_type_str, &
            limiter_type_str, logfile, srcstr)
       limiter_type_str = str_to_lower(limiter_type_str)

       call fson_get_mpi(limiter_json, "limit", &
            default_source_control_limiter_limit, limit, logfile, srcstr)

       call spec_sources%traverse(limiter_iterator)

    end if

  contains

    subroutine limiter_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      type(list_type) :: single_source

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         call single_source%init(owner = PETSC_FALSE)
         call single_source%append(source)
          select case (limiter_type_str)
          case ("total")
             allocate(source_control_total_limiter_type :: limiter)
          case ("water")
             allocate(source_control_water_limiter_type :: limiter)
          case ("steam")
             allocate(source_control_steam_limiter_type :: limiter)
          case default
             allocate(source_control_total_limiter_type :: limiter)
          end select
          call limiter%init(limit, single_source)
          call source_controls%append(limiter)
       end select

    end subroutine limiter_iterator

  end subroutine setup_limiter_source_controls

!------------------------------------------------------------------------

  subroutine setup_direction_source_control(source_json, srcstr, &
       thermo, spec_sources, source_controls, logfile)
    !! Set up direction source control. This can control multiple
    !! sources, so only one is created.

    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    class(thermodynamics_type), intent(in) :: thermo
    type(list_type), intent(in out) :: spec_sources
    type(list_type), intent(in out) :: source_controls
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscInt :: direction
    type(source_control_direction_type), pointer :: direction_control
    PetscInt, parameter :: max_direction_str_length = 16
    character(max_direction_str_length) :: direction_str

    if (fson_has_mpi(source_json, "direction")) then

       call fson_get_mpi(source_json, "direction", val = direction_str)
       select case (trim(str_to_lower(direction_str)))
       case ("production", "out")
          direction = SRC_DIRECTION_PRODUCTION
       case ("injection", "in")
          direction = SRC_DIRECTION_INJECTION
       case ("both")
          direction = SRC_DIRECTION_BOTH
       case default
          if (present(logfile) .and. logfile%active) then
             call logfile%write(LOG_LEVEL_WARN, 'input', 'unrecognised', &
                  str_key = trim(srcstr) // "direction", &
                  str_value = direction_str)
          end if
          direction = default_source_direction
       end select

       if ((direction /= SRC_DIRECTION_BOTH) .and. &
            spec_sources%count > 0) then
         allocate(direction_control)
         call direction_control%init(direction, spec_sources%copy())
         call source_controls%append(direction_control)
       end if

    end if

  end subroutine setup_direction_source_control

!------------------------------------------------------------------------

end module source_setup_module
