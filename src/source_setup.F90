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
       TYPE_STRING, TYPE_OBJECT, TYPE_LOGICAL, TYPE_NULL
  use fson_mpi_module
  use list_module
  use logfile_module
  use eos_module
  use thermodynamics_module, only: thermodynamics_type
  use source_network_module
  use source_module
  use source_network_node_module
  use source_network_control_module
  use source_network_group_module
  use source_network_reinjector_module
  use source_control_module
  use separator_module
  use tracer_module

  implicit none
  private

  public :: setup_source_network

contains

!------------------------------------------------------------------------

  subroutine setup_source_network(json, dm, ao, eos, tracers, thermo, start_time, &
       fluid_vector, fluid_range_start, source_network, logfile, err)
    !! Sets up source network, including sinks / sources, source
    !! controls and source groups.

    use dm_utils_module
    use dictionary_module
    use mpi_utils_module, only: invert_indices
    use dag_module
    use fson_utils_module, only: pfson_value_type

    type(fson_value), pointer, intent(in) :: json !! JSON file object
    DM, intent(in) :: dm !! Mesh DM
    AO, intent(in) :: ao !! Application ordering for natural to global cell indexing
    class(eos_type), intent(in) :: eos !! Equation of state
    type(tracer_type), intent(in) :: tracers(:) !! Tracers
    class(thermodynamics_type), intent(in out) :: thermo !! Thermodynamics formulation
    PetscReal, intent(in) :: start_time
    Vec, intent(in) :: fluid_vector !! Fluid vector
    PetscInt, intent(in) :: fluid_range_start !! Range start for global fluid vector
    type(source_network_type), intent(in out) :: source_network !! Source network
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    PetscInt, allocatable :: dependency_indices(:)
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscMPIInt :: rank
    DM :: dm_source, dm_group, dm_reinjector
    PetscInt :: num_local_sources, source_spec_index, local_source_index, idep
    PetscInt :: sorted_group_index, sorted_reinjector_index
    type(fson_value), pointer :: sources_json, source_json, groups_json, reinjectors_json
    PetscInt :: num_source_specs, num_tracers, num_local_root_groups, num_local_root_reinjectors
    PetscReal, pointer, contiguous :: fluid_data(:), source_data(:), &
         source_network_group_data(:), source_network_reinjector_data(:)
    PetscSection :: fluid_section, source_section, source_network_group_section, &
         source_network_reinjector_section
    type(pfson_value_type), allocatable :: group_specs_array(:), reinjector_specs_array(:)
    type(dag_type) :: group_dag, reinjector_dag
    type(dictionary_type) :: source_dict, source_dict_all
    type(dictionary_type) :: group_dict, group_index_dict
    type(dictionary_type) :: reinjector_dict, reinjector_index_dict
    PetscInt, allocatable :: indices(:), group_indices(:), reinjector_indices(:)
    PetscInt, allocatable :: group_order(:), reinjector_order(:)
    character(max_tracer_name_length), allocatable :: tracer_names(:)
    PetscErrorCode :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    source_network%num_sources = 0
    source_network%num_groups = 0
    source_network%num_reinjectors = 0
    num_tracers = size(tracers)
    tracer_names = tracers%name
    num_local_root_groups = 0
    num_local_root_reinjectors = 0
    err = 0

    call label_source_zones(json, num_local_sources, logfile, err)
    if (err == 0) then

       call source_network%init()

       call source_dict%init(owner = PETSC_FALSE)
       call source_dict_all%init(owner = PETSC_FALSE)
       call group_dict%init(owner = PETSC_FALSE)
       call group_index_dict%init(owner = PETSC_TRUE)
       call reinjector_dict%init(owner = PETSC_FALSE)
       call reinjector_index_dict%init(owner = PETSC_TRUE)

       call create_path_dm(num_local_sources, dm_source)
       call setup_source_dm_data_layout(dm_source)
       call DMCreateGlobalVector(dm_source, source_network%source, ierr); CHKERRQ(ierr)
       call PetscObjectSetName(source_network%source, "source", ierr); CHKERRQ(ierr)
       call global_vec_range_start(source_network%source, source_network%source_range_start)

       call global_vec_section(fluid_vector, fluid_section)
       call VecGetArrayReadF90(fluid_vector, fluid_data, ierr); CHKERRQ(ierr)
       call global_vec_section(source_network%source, source_section)
       call VecGetArrayF90(source_network%source, source_data, ierr); CHKERRQ(ierr)

       if (fson_has_mpi(json, "source")) then
          call fson_get_mpi(json, "source", sources_json)
          num_source_specs = fson_value_count_mpi(sources_json, ".")
          source_json => fson_value_children_mpi(sources_json)
          local_source_index = 0
          do source_spec_index = 0, num_source_specs - 1
             call setup_source(source_spec_index, local_source_index, &
                  source_json, ao, tracer_names, thermo, source_network, &
                  source_dict, source_dict_all, err)
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
          call source_network%sources%traverse(source_indices_iterator)
          source_network%source_index = invert_indices(indices, "source_index")
          deallocate(indices)

          if (fson_has_mpi(json, "network.group")) then
             call fson_get_mpi(json, "network.group", groups_json)
             source_network%num_groups = fson_value_count_mpi(groups_json, ".")
             call setup_item_index_dict(groups_json, &
                  source_network%num_groups, group_index_dict)
             call setup_group_dag(groups_json, source_network%num_groups, &
                  group_index_dict, group_dag, group_specs_array)
             call group_dag%sort(group_order, err)
             if (err == 0) then
                call init_source_network_groups(source_network, &
                     num_local_root_groups, logfile, err)
             end if
          end if
          
          if (err == 0) then
             call create_path_dm(num_local_root_groups, dm_group)
             call setup_group_dm_data_layout(dm_group)
             call DMCreateGlobalVector(dm_group, source_network%group, &
                  ierr); CHKERRQ(ierr)
             call PetscObjectSetName(source_network%group, "network_group", &
                  ierr); CHKERRQ(ierr)
             call global_vec_range_start(source_network%group, &
                  source_network%group_range_start)
             call global_vec_section(source_network%group, source_network_group_section)
             call VecGetArrayF90(source_network%group, source_network_group_data, &
                  ierr); CHKERRQ(ierr)
             sorted_group_index = 0
             call source_network%groups%traverse(group_init_data_iterator)
             allocate(group_indices(num_local_root_groups))
             call source_network%groups%traverse(source_network_group_indices_iterator)
             source_network%group_index = invert_indices(group_indices, "network_group_index")
             deallocate(group_indices)
             call VecRestoreArrayF90(source_network%group, source_network_group_data, ierr)
             CHKERRQ(ierr)

             if (fson_has_mpi(json, "network.reinject")) then
                call fson_get_mpi(json, "network.reinject", reinjectors_json)
                source_network%num_reinjectors = fson_value_count_mpi(reinjectors_json, ".")
                call setup_item_index_dict(reinjectors_json, &
                     source_network%num_reinjectors, reinjector_index_dict)
                call setup_reinjector_dag(reinjectors_json, source_network%num_reinjectors, &
                     reinjector_index_dict, reinjector_dag, reinjector_specs_array)
                call reinjector_dag%sort(reinjector_order, err)
                call init_source_network_reinjectors(source_network, &
                     num_local_root_reinjectors, logfile, err)
             end if
             if (err == 0) then
                call create_path_dm(num_local_root_reinjectors, dm_reinjector)
                call setup_reinjector_dm_data_layout(dm_reinjector)
                call DMCreateGlobalVector(dm_reinjector, source_network%reinjector, &
                     ierr); CHKERRQ(ierr)
                call PetscObjectSetName(source_network%reinjector, "network_reinject", &
                     ierr); CHKERRQ(ierr)
                call global_vec_range_start(source_network%reinjector, &
                     source_network%reinjector_range_start)
                call global_vec_section(source_network%reinjector, &
                     source_network_reinjector_section)
                call VecGetArrayF90(source_network%reinjector, source_network_reinjector_data, &
                     ierr); CHKERRQ(ierr)
                sorted_reinjector_index = 0
                call source_network%reinjectors%traverse(reinjector_init_data_iterator)
                allocate(reinjector_indices(num_local_root_reinjectors))
                call source_network%reinjectors%traverse(source_network_reinjector_indices_iterator)
                source_network%reinjector_index = invert_indices(reinjector_indices, &
                     "network_reinject_index")
                deallocate(reinjector_indices)
                call VecRestoreArrayF90(source_network%reinjector, &
                     source_network_reinjector_data, ierr); CHKERRQ(ierr)
             end if

          end if

       end if

       call VecRestoreArrayF90(source_network%source, source_data, ierr); CHKERRQ(ierr)
       call VecRestoreArrayReadF90(fluid_vector, fluid_data, ierr)
       CHKERRQ(ierr)
       call source_dict%destroy()
       call source_dict_all%destroy()
       call group_dict%destroy()
       call group_index_dict%destroy()
       call reinjector_dict%destroy()
       call reinjector_index_dict%destroy()

    end if

  contains

!........................................................................

    PetscBool function null_cell_source(source_json)
      !! Returns true if source specification does not specify any
      !! cells (or zones).

      type(fson_value), pointer, intent(in) :: source_json
      ! Locals:
      PetscInt :: cell_type
      PetscInt, allocatable :: cells(:)

      if (fson_has_mpi(source_json, "cell")) then
         cell_type = fson_type_mpi(source_json, "cell")
         null_cell_source = (cell_type == TYPE_NULL)
      else if (fson_has_mpi(source_json, "cells")) then
         cell_type = fson_type_mpi(source_json, "cells")
         if (cell_type == TYPE_NULL) then
            null_cell_source = PETSC_TRUE
         else if (cell_type == TYPE_INTEGER) then
            null_cell_source = PETSC_FALSE
         else if (cell_type == TYPE_ARRAY) then
            call fson_get_mpi(source_json, "cells", val = cells)
            null_cell_source = (.not. allocated(cells))
         end if
      else if (fson_has_mpi(source_json, "zones")) then
         null_cell_source = PETSC_FALSE
      else
         null_cell_source = PETSC_TRUE
      end if

    end function null_cell_source

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
      PetscBool :: has_label, null_cell
      IS :: cell_IS
      PetscInt, pointer, contiguous :: zone_cells(:), cells(:)
      PetscInt, allocatable :: labelled_cells(:)
      DMLabel :: ghost_label, source_label
      PetscMPIInt :: rank
      PetscErrorCode :: ierr

      num_local_sources = 0
      call DMGetLabel(dm, "ghost", ghost_label, ierr)
      call DMGetLabel(dm, source_label_name, source_label, ierr)
      call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)

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
            else ! check for source with no cell:
               null_cell = null_cell_source(source_json)
               if ((null_cell) .and. (rank == 0)) then
                  ! assign no-cell source to root rank:
                  num_local_sources = num_local_sources + 1
               end if
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
      PetscBool, parameter :: rate_specified = PETSC_FALSE, &
           enthalpy_specified = PETSC_FALSE
      PetscReal, parameter :: specified_rate = -1._dp, specified_enthalpy = 0._dp

      call source%init("", eos, 0, 0, 0._dp, 1, 1, rate_specified, &
           specified_rate, enthalpy_specified, specified_enthalpy, num_tracers)
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

    subroutine setup_group_dm_data_layout(dm_group)
      !! Sets up data layout on source group DM.

      DM, intent(in out) :: dm_group
      ! Locals:
      PetscInt, allocatable :: num_field_components(:), field_dim(:)
      character(max_source_network_variable_name_length), allocatable :: field_names(:)

      allocate(num_field_components(num_source_network_group_variables), &
           field_dim(num_source_network_group_variables), &
           field_names(num_source_network_group_variables))
      num_field_components = 1
      field_dim = 0
      field_names = source_network_group_variable_names

      call dm_set_data_layout(dm_group, num_field_components, field_dim, &
           field_names)
      deallocate(num_field_components, field_dim, field_names)

    end subroutine setup_group_dm_data_layout

!........................................................................

    subroutine setup_source(source_spec_index, local_source_index, source_json, &
         ao, tracer_names, thermo, source_network, source_dict, source_dict_all, err)
      !! Sets up all cell sources for a source specification.

      PetscInt, intent(in) :: source_spec_index !! Index of source specification
      PetscInt, intent(in out) :: local_source_index !! Index of source
      type(fson_value), pointer, intent(in) :: source_json !! JSON input for specification
      AO, intent(in) :: ao !! Application ordering for natural to global cell indexing
      character(*), intent(in) :: tracer_names(:) !! Tracer names
      class(thermodynamics_type), intent(in out) :: thermo !! Water thermodynamics
      type(source_network_type), intent(in out) :: source_network !! Source network
      type(dictionary_type), intent(in out) :: source_dict !! Dictionary of local named sources
      type(dictionary_type), intent(in out) :: source_dict_all !! Dictionary of all named sources
      PetscErrorCode, intent(out) :: err
      ! Locals:
      type(source_type), pointer :: source
      character(len=64) :: srcstr
      character(len=12) :: istr
      PetscInt :: injection_component, production_component
      PetscReal :: initial_rate, initial_enthalpy
      PetscReal, allocatable :: separator_pressure(:)
      PetscInt :: num_cells, num_cells_all, i, source_offset
      PetscInt, allocatable :: natural_source_index(:), natural_cell_index(:)
      type(list_type) :: spec_sources
      IS :: cell_IS
      PetscInt, pointer, contiguous :: local_cell_index(:)
      ISLocalToGlobalMapping :: l2g
      PetscReal :: tracer_injection_rate(size(tracer_names))
      character(max_source_network_node_name_length) :: name
      PetscBool :: rate_specified, enthalpy_specified, null_cell
      PetscReal :: specified_rate, specified_enthalpy

      err = 0
      call DMGetLocalToGlobalMapping(dm, l2g, ierr); CHKERRQ(ierr)
      write(istr, '(i0)') source_spec_index
      srcstr = 'source[' // trim(istr) // '].'

      null_cell = null_cell_source(source_json)
      call get_components(source_json, eos, &
           injection_component, production_component, logfile)
      call get_initial_rate(source_json, initial_rate, rate_specified, specified_rate)
      call get_initial_enthalpy(source_json, eos, &
           injection_component, initial_enthalpy, enthalpy_specified, specified_enthalpy)
      call get_separator_pressure(source_json, srcstr, separator_pressure, logfile)
      call get_tracer_injection_rate(source_json, tracer_names, &
           srcstr, tracer_injection_rate, logfile, err)

      if (err == 0) then

         if (fson_has_mpi(source_json, "name")) then
            call fson_get_mpi(source_json, "name", val = name)
         else
            name = ""
         end if

         if (null_cell) then
            if (rank == 0) then
               num_cells = 1
            else
               num_cells = 0
            end if
         else
            call DMGetStratumSize(dm, source_label_name, source_spec_index, &
                 num_cells, ierr); CHKERRQ(ierr)
         end if
         call MPI_allreduce(num_cells, num_cells_all, 1, MPI_INTEGER, MPI_SUM, &
              PETSC_COMM_WORLD, ierr)

         if (num_cells > 0) then
            allocate(natural_cell_index(num_cells))
            if (null_cell) then
               allocate(local_cell_index(num_cells))
               local_cell_index = -1
               natural_cell_index = -1
            else
               call DMGetStratumIS(dm, source_label_name, source_spec_index, &
                    cell_IS, ierr); CHKERRQ(ierr)
               call ISGetIndicesF90(cell_IS, local_cell_index, ierr); CHKERRQ(ierr)
               natural_cell_index = local_to_natural_cell_index(ao, l2g, local_cell_index)
            end if
         else
            allocate(natural_cell_index(1))
         end if
         natural_source_index = get_natural_source_indices(num_cells, num_cells_all, &
              natural_cell_index, null_cell)

         call spec_sources%init(owner = PETSC_FALSE)
         if (num_cells > 0) then
            do i = 1, num_cells
               source_offset = global_section_offset(source_section, local_source_index, &
                    source_network%source_range_start)
               allocate(source)
               call source%init(name, eos, local_source_index, local_cell_index(i), &
                    initial_enthalpy, injection_component, production_component, &
                    rate_specified, specified_rate, enthalpy_specified, specified_enthalpy, &
                    num_tracers)
               call source%assign(source_data, source_offset)
               call source%init_data(natural_source_index(i), natural_cell_index(i), &
                    initial_rate, tracer_injection_rate, separator_pressure, thermo)
               call source_network%sources%append(source)
               call spec_sources%append(source)
               if (source%separator%on) then
                  call source_network%separated_sources%append(source)
               end if
               local_source_index = local_source_index + 1
            end do
            if (.not. null_cell) then
               call ISRestoreIndicesF90(cell_IS, local_cell_index, ierr); CHKERRQ(ierr)
            end if
            deallocate(natural_cell_index, natural_source_index)
         end if
         call setup_inline_source_controls(source_json, eos, thermo, &
              start_time, source_data, source_section, &
              fluid_data, fluid_section, fluid_range_start, srcstr, &
              tracer_names, null_cell, spec_sources, source_network, logfile, err)

         if ((name /= "") .and. (num_cells_all == 1)) then
            ! Uniquely named source- add to dictionaries:
            call source_dict_all%add(name)
            if (num_cells == 1) then
               call source_dict%add(name, spec_sources%head%data)
            end if
         end if

         call spec_sources%destroy()
         source_network%num_sources = source_network%num_sources + num_cells_all

      end if

    end subroutine setup_source

!........................................................................

    function get_natural_source_indices(num_cells, num_cells_all, &
         natural_cell_index, null_cell) result(natural_source_index)
      !! Gets natural source indices for the current source.

      use utils_module, only: array_cumulative_sum, array_indices_in_int_array
      use mpi_utils_module, only: get_mpi_int_gather_array

      PetscInt, intent(in) :: num_cells !! Number of source cells on current process
      PetscInt, intent(in) :: num_cells_all !! Number of source cells on all processes
      PetscInt, intent(in) :: natural_cell_index(:) !! Natural cell indices on current process
      PetscBool, intent(in) :: null_cell !! If no cells specified in source
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
           natural_cell_index_all, null_cell)

      if (rank == 0) then
         natural_source_index_all = array_indices_in_int_array( &
              ordered_natural_cell_index_all, natural_cell_index_all)
         natural_source_index_all = natural_source_index_all + &
              source_network%num_sources - 1
      end if
      call MPI_scatterv(natural_source_index_all, cell_counts, cell_displacements, &
           MPI_INTEGER, natural_source_index, num_cells, MPI_INTEGER, &
           0, PETSC_COMM_WORLD, ierr)

      deallocate(cell_counts, cell_displacements, natural_source_index_all, &
           natural_cell_index_all, ordered_natural_cell_index_all)

    end function get_natural_source_indices

!........................................................................

    function get_source_cell_indices(source_json, cells, null_cell) &
         result(ordered_cells)

      !! Gets array of all natural cell indices for the source, in
      !! natural source order, on rank 0. For sources with the "cell"
      !! or "cells" properties this can simply be read from the
      !! input. For sources with cells defined by zones, cells from
      !! all ranks are ordered by natural cell index.

      type(fson_value), pointer, intent(in) :: source_json
      PetscInt, intent(in) :: cells(:)
      PetscBool, intent(in) :: null_cell
      PetscInt, allocatable :: ordered_cells(:)
      ! Locals:
      type(fson_value), pointer :: cell_json, cells_json, zones_json
      PetscInt :: c
      PetscMPIInt :: rank
      PetscErrorCode :: ierr

      call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

      if (rank == 0) then

         if (null_cell) then
            ordered_cells = [-1]
         else

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
              source_network%source_range_start)
         call source%assign(source_data, source_offset)
         indices(s + 1) = nint(source%source_index)
      end select

    end subroutine source_indices_iterator

!........................................................................

    subroutine setup_item_index_dict(items_json, num_items, item_index_dict)
      !! Forms dictionary mapping item names to their item
      !! indices. Items can be groups or reinjectors.

      type(fson_value), pointer, intent(in out) :: items_json
      PetscInt, intent(in) :: num_items
      type(dictionary_type), intent(in out) :: item_index_dict
      ! Locals:
      type(fson_value), pointer :: item_json
      PetscInt :: item_index
      character(max_source_network_node_name_length) :: name
      PetscInt, pointer :: i

      item_json => fson_value_children_mpi(items_json)
      do item_index = 0, num_items - 1
         call fson_get_mpi(item_json, "name", "", name)
         if (name /= "") then
            allocate(i)
            i = item_index
            call item_index_dict%add(name, i)
         end if
         item_json => fson_value_next_mpi(item_json)
      end do

    end subroutine setup_item_index_dict

!........................................................................

    subroutine setup_group_dag(groups_json, num_groups, &
         group_index_dict, group_dag, group_specs_array)
      !! Forms directed acyclic graph (DAG) representing group
      !! dependencies. (Also populates group_specs_array, for random
      !! access into group JSON specifications.)

      type(fson_value), pointer, intent(in out) :: groups_json
      PetscInt, intent(in) :: num_groups
      type(dictionary_type), intent(in out) :: group_index_dict
      type(dag_type), intent(in out) :: group_dag
      type(pfson_value_type), allocatable, intent(in out) :: group_specs_array(:)
      ! Locals:
      type(fson_value), pointer :: group_json
      character(max_source_network_node_name_length), allocatable :: node_names(:)
      type(list_node_type), pointer :: dict_node
      PetscInt :: group_index, i
      PetscMPIInt :: rank
      PetscErrorCode :: ierr

      call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
      call group_dag%init(num_groups)
      allocate(group_specs_array(num_groups))

      group_json => fson_value_children_mpi(groups_json)
      do group_index = 0, num_groups - 1
         if (fson_has_mpi(group_json, "in")) then
            call fson_get_mpi(group_json, "in", &
                 string_length = max_source_network_node_name_length, &
                 val = node_names)
            associate(num_nodes => size(node_names))
              allocate(dependency_indices(num_nodes))
              do i = 1, num_nodes
                 dict_node => group_index_dict%get(node_names(i))
                 if (associated(dict_node)) then
                    select type (idx => dict_node%data)
                    type is (PetscInt)
                       dependency_indices(i) = idx
                    end select
                 else
                    dependency_indices(i) = -1
                 end if
              end do
            end associate
            dependency_indices = pack(dependency_indices, &
                 dependency_indices > 0)
         end if
         call group_dag%set_edges(group_index, dependency_indices)
         if (rank == 0) then
            call group_specs_array(group_index + 1)%set(group_json)
         end if
         deallocate(node_names, dependency_indices)
         group_json => fson_value_next_mpi(group_json)
      end do

    end subroutine setup_group_dag

!........................................................................

    subroutine source_network_group_indices_iterator(node, stopped)
      !! Gets indices from all local source groups.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: g, source_network_group_offset

      stopped = PETSC_FALSE
      select type(group => node%data)
      class is (source_network_group_type)
         if (group%rank == 0) then
            g = group%local_group_index
            source_network_group_offset = global_section_offset( &
                 source_network_group_section, g, source_network%group_range_start)
            call group%assign(source_network_group_data, &
                 source_network_group_offset)
            group_indices(g + 1) = nint(group%group_index)
         end if
      end select

    end subroutine source_network_group_indices_iterator

!........................................................................

    PetscErrorCode function check_valid_scaling(group_json, scaling_type) &
         result(err)
      ! Returns true if the JSON spec for a group has valid
      ! scaling. Progressive scaling is invalid if the group has its
      ! own separator and a separated water or steam limiter.

      type(fson_value), pointer, intent(in out) :: group_json
      character(*), intent(in) :: scaling_type

      err = 0
      if ((scaling_type == 'progressive') .and. &
           (fson_has_mpi(group_json, 'separator')) .and. &
           (fson_has_mpi(group_json, 'limiter.water') .or. &
           fson_has_mpi(group_json, 'limiter.steam'))) then
            err = 1
      end if

    end function check_valid_scaling

!........................................................................

    subroutine init_source_network_groups(source_network, &
         num_local_root_groups, logfile, err)
      !! Initialise source network groups and controls and return
      !! number of local root groups. An error is returned if any
      !! unrecognised group input nodes are specified.

      use utils_module, only: str_to_lower

      type(source_network_type), intent(in out) :: source_network
      PetscInt, intent(in out) :: num_local_root_groups
      type(logfile_type), intent(in out), optional :: logfile
      PetscErrorCode, intent(out) :: err
      ! Locals:
      PetscInt :: ig, i, g, group_index, num_nodes, comm_size
      PetscInt, allocatable :: source_cell_indices(:)
      type(fson_value), pointer :: group_json
      character(max_source_network_node_name_length) :: name
      character(max_source_network_node_name_length), allocatable :: node_names(:)
      class(source_network_group_type), pointer :: group
      type(list_node_type), pointer :: source_dict_node, group_dict_node
      type(dictionary_type) :: group_source_dict
      character(len=64) :: grpstr
      character(len=12) :: igstr
      character(len=16) :: scaling_type
      character(len=8), parameter :: default_group_scaling_type = "uniform"

      err = 0
      call group_source_dict%init(PETSC_FALSE)

      do ig = 0, source_network%num_groups - 1

         group_index = group_order(ig) + 1
         group_json => group_specs_array(group_index)%ptr
         write(igstr, '(i0)') group_index - 1
         grpstr = 'network.group[' // trim(igstr) // '].'

         call fson_get_mpi(group_json, "name", "", name)
         call fson_get_mpi(group_json, "in", &
              string_length = max_source_network_node_name_length, &
              val = node_names)
         call fson_get_mpi(group_json, "scaling", default_group_scaling_type, &
              scaling_type, logfile, trim(grpstr) // "scaling")
         scaling_type = str_to_lower(scaling_type)

         err = check_valid_scaling(group_json, scaling_type)
         if (err == 0) then

            select case (scaling_type)
            case ("progressive")
               allocate(progressive_scaling_source_network_group_type :: group)
            case default
               allocate(uniform_scaling_source_network_group_type :: group)
            end select

            call group%init(name)
            call setup_inline_source_group_controls(group_json, group, &
                 source_network, logfile)
            num_nodes = size(node_names)
            allocate(source_cell_indices(num_nodes))
            source_cell_indices = -1

            do i = 1, num_nodes
               associate(node_name => node_names(i))
                 if (source_dict_all%has(node_name)) then
                    if (group_source_dict%has(node_name)) then
                       if (present(logfile)) then
                          call logfile%write(LOG_LEVEL_ERR, "input", &
                               "source " // trim(node_name) // &
                               " outputs to more than one group.")
                       end if
                       err = 1
                       deallocate(group)
                       exit
                    else
                       source_dict_node => source_dict%get(node_name)
                       if (associated(source_dict_node)) then
                          select type (source => source_dict_node%data)
                          type is (source_type)
                             source%link_index = i
                             call group%in%append(source)
                             source_cell_indices(i) = source%natural_cell_index
                          end select
                       end if
                       call group_source_dict%add(node_name)
                    end if
                 else
                    group_dict_node => group_dict%get(node_name)
                    if (associated(group_dict_node)) then
                       select type (in_group => group_dict_node%data)
                       class is (source_network_group_type)
                          if (associated(in_group%out)) then
                             if (present(logfile)) then
                                call logfile%write(LOG_LEVEL_ERR, "input", &
                                     "group " // trim(in_group%name) // &
                                     " outputs to more than one group.")
                             end if
                             err = 1
                             deallocate(group)
                             exit
                          else
                             in_group%out => group
                             if (in_group%rank == 0) then
                                in_group%link_index = i
                                source_cell_indices = [source_cell_indices, &
                                     in_group%source_cell_indices]
                             end if
                             call group%in%append(in_group)
                          end if
                       end select
                    else
                       if (present(logfile)) then
                          call logfile%write(LOG_LEVEL_ERR, "input", &
                               "unrecognised group input: " // trim(node_name))
                       end if
                       err = 1
                       deallocate(group)
                       exit
                    end if
                 end if
               end associate
            end do

            if (err == 0) then

               call group%init_comm()

               if (group%rank == 0) then
                  g = num_local_root_groups
                  num_local_root_groups = num_local_root_groups + 1
               else
                  g = -1
               end if
               group%local_group_index = g

               source_cell_indices = pack(source_cell_indices, &
                    source_cell_indices >= 0)
               group%source_cell_indices = group%gatherv(source_cell_indices)

               call source_network%groups%append(group)
               if (name /= "") then
                  call group_dict%add(name, group)
               end if
            else
               exit
            end if

         else
            if (present(logfile)) then
               call logfile%write(LOG_LEVEL_ERR, "input", &
                    "invalid_scaling", str_key = "group", &
                    str_value = trim(name))
            end if
            exit
         end if

         deallocate(node_names, source_cell_indices)

      end do

      call group_source_dict%destroy()

    end subroutine init_source_network_groups

!........................................................................

    subroutine group_init_data_iterator(node, stopped)
      !! Initialises data in a source group.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: g, source_network_group_offset, group_index
      character(len=64) :: grpstr
      character(len=12) :: istr
      type(fson_value), pointer :: group_json
      PetscReal, allocatable :: separator_pressure(:)

      stopped = PETSC_FALSE
      select type (group => node%data)
      class is (source_network_group_type)
         group_index = group_order(sorted_group_index)
         write(istr, '(i0)') group_index
         grpstr = 'network.group[' // trim(istr) // '].'
         group_json => group_specs_array(group_index + 1)%ptr
         call get_separator_pressure(group_json, grpstr, &
              separator_pressure, logfile)
         if (group%rank == 0) then
            g = group%local_group_index
            source_network_group_offset = global_section_offset( &
                 source_network_group_section, g, source_network%group_range_start)
            call group%assign(source_network_group_data, source_network_group_offset)
            call group%init_data(group_index, separator_pressure, thermo)
         end if
         sorted_group_index = sorted_group_index + 1
      end select

    end subroutine group_init_data_iterator

!........................................................................

    subroutine reinjector_init_data_iterator(node, stopped)
      !! Initialises data in a source reinjector.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: r, source_network_reinjector_offset, reinjector_index
      character(len=64) :: rstr
      character(len=12) :: istr
      type(fson_value), pointer :: reinjector_json

      stopped = PETSC_FALSE
      select type (reinjector => node%data)
      class is (source_network_reinjector_type)
         reinjector_index = reinjector_order(sorted_reinjector_index)
         write(istr, '(i0)') reinjector_index
         rstr = 'network.reinjector[' // trim(istr) // '].'
         reinjector_json => reinjector_specs_array(reinjector_index + 1)%ptr
         if (reinjector%rank == 0) then
            r = reinjector%local_reinjector_index
            source_network_reinjector_offset = global_section_offset( &
                 source_network_reinjector_section, r, source_network%reinjector_range_start)
            call reinjector%assign(source_network_reinjector_data, &
                 source_network_reinjector_offset)
            call reinjector%init_data(reinjector_index)
         end if
         sorted_reinjector_index = sorted_reinjector_index + 1
      end select

    end subroutine reinjector_init_data_iterator

!........................................................................

    subroutine source_network_reinjector_indices_iterator(node, stopped)
      !! Gets indices from all local source reinjectors.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: r, source_network_reinjector_offset

      stopped = PETSC_FALSE
      select type(reinjector => node%data)
      class is (source_network_reinjector_type)
         if (reinjector%rank == 0) then
            r = reinjector%local_reinjector_index
            source_network_reinjector_offset = global_section_offset( &
                 source_network_reinjector_section, r, &
                 source_network%reinjector_range_start)
            call reinjector%assign(source_network_reinjector_data, &
                 source_network_reinjector_offset)
            reinjector_indices(r + 1) = nint(reinjector%reinjector_index)
         end if
      end select

    end subroutine source_network_reinjector_indices_iterator

!........................................................................

    subroutine setup_reinjector_dm_data_layout(dm_reinjector)
      !! Sets up data layout on source reinjector DM.

      DM, intent(in out) :: dm_reinjector
      ! Locals:
      PetscInt, allocatable :: num_field_components(:), field_dim(:)
      character(max_source_network_variable_name_length), allocatable :: field_names(:)

      allocate(num_field_components(num_source_network_reinjector_variables), &
           field_dim(num_source_network_reinjector_variables), &
           field_names(num_source_network_reinjector_variables))
      num_field_components = 1
      field_dim = 0
      field_names = source_network_reinjector_variable_names

      call dm_set_data_layout(dm_reinjector, num_field_components, field_dim, &
           field_names)
      deallocate(num_field_components, field_dim, field_names)

    end subroutine setup_reinjector_dm_data_layout

!........................................................................

    subroutine setup_reinjector_dag(reinjectors_json, num_reinjectors, &
         reinjector_index_dict, reinjector_dag, reinjector_specs_array)
      !! Forms directed acyclic graph (DAG) representing reinjectors
      !! dependencies. (Also populates reinjector_specs_array, for random
      !! access into reinjector JSON specifications.)

      use fson_value_m, only : TYPE_STRING

      type(fson_value), pointer, intent(in out) :: reinjectors_json
      PetscInt, intent(in) :: num_reinjectors
      type(dictionary_type), intent(in out) :: reinjector_index_dict
      type(dag_type), intent(in out) :: reinjector_dag
      type(pfson_value_type), allocatable, intent(in out) :: reinjector_specs_array(:)
      ! Locals:
      type(fson_value), pointer :: reinjector_json, outputs_json, output_json
      character(5) :: key
      type(list_type) :: dep_list
      type(list_node_type), pointer :: dict_node
      PetscInt :: reinjector_index, i, k, num_outputs, overflow_type, out_type
      character(max_source_network_node_name_length) :: node_name
      PetscInt, parameter :: num_keys = 2
      character(5), parameter :: keys(num_keys) = ["water", "steam"]
      PetscMPIInt :: rank
      PetscErrorCode :: ierr

      call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

      call reinjector_dag%init(num_reinjectors)
      allocate(reinjector_specs_array(num_reinjectors))

      reinjector_json => fson_value_children_mpi(reinjectors_json)
      do reinjector_index = 0, num_reinjectors - 1

         call dep_list%init(owner = PETSC_FALSE)
         do k = 1, num_keys
            key = keys(k)
            if (fson_has_mpi(reinjector_json, trim(key))) then
               call fson_get_mpi(reinjector_json, trim(key), &
                    outputs_json)
               num_outputs = fson_value_count_mpi(outputs_json, ".")
               output_json => fson_value_children_mpi(outputs_json)
               do i = 1, num_outputs
                  if (fson_has_mpi(output_json, "out")) then
                     out_type = fson_type_mpi(output_json, "out")
                     if (out_type == TYPE_STRING) then
                        call fson_get_mpi(output_json, "out", &
                             val = node_name)
                        dict_node => reinjector_index_dict%get(node_name)
                        if (associated(dict_node)) then
                           select type (idx => dict_node%data)
                           type is (PetscInt)
                              call dep_list%append(idx)
                           end select
                        end if
                     end if
                  end if
                  output_json => fson_value_next_mpi(output_json)
               end do
            end if
         end do

         node_name = ""
         if (fson_has_mpi(reinjector_json, "overflow")) then
            overflow_type = fson_type_mpi(reinjector_json, "overflow")
            select case (overflow_type)
            case (TYPE_STRING)
               call fson_get_mpi(reinjector_json, "overflow", val = node_name)
            case (TYPE_OBJECT)
               if (fson_has_mpi(reinjector_json, "overflow.out")) then
                  call fson_get_mpi(reinjector_json, "overflow.out", &
                       val = node_name)
               end if
            end select
         end if
         if (node_name /= "") then
            dict_node => reinjector_index_dict%get(node_name)
            if (associated(dict_node)) then
               select type (idx => dict_node%data)
               type is (PetscInt)
                  call dep_list%append(idx)
               end select
            end if
         end if

         allocate(dependency_indices(dep_list%count))
         idep = 1
         call dep_list%traverse(dag_dependency_iterator)
         call dep_list%destroy()
         call reinjector_dag%set_edges(reinjector_index, dependency_indices)
         deallocate(dependency_indices)
         if (rank == 0) then
            call reinjector_specs_array(reinjector_index + 1)%set(reinjector_json)
         end if

         reinjector_json => fson_value_next_mpi(reinjector_json)
      end do

    end subroutine setup_reinjector_dag

!........................................................................

    subroutine dag_dependency_iterator(node, stopped)
      !! Converts dependency list into an integer array.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type(idx => node%data)
      type is (PetscInt)
         dependency_indices(idep) = idx
         idep = idep + 1
      end select

    end subroutine dag_dependency_iterator

!........................................................................

    subroutine init_reinjector_input_source(name, reinjector_input_dict, &
         source_dict, reinjector, err)
      !! Initialises reinjector input from a source.

      use mpi_utils_module, only: mpi_broadcast_error_flag

      character(*), intent(in) :: name
      type(dictionary_type), intent(in out) :: reinjector_input_dict, source_dict
      type(source_network_reinjector_type), pointer, intent(in out) :: reinjector
      PetscErrorCode, intent(in out) :: err
      ! Locals:
      type(list_node_type), pointer :: source_dict_node

      if (reinjector_input_dict%has(name)) then
         if (present(logfile)) then
            call logfile%write(LOG_LEVEL_ERR, "input", &
                 "duplicate reinjector input: " // trim(name))
         end if
         err = 1
         deallocate(reinjector)
      else
         source_dict_node => source_dict%get(name)
         if (associated(source_dict_node)) then
            select type (source => source_dict_node%data)
            type is (source_type)
               reinjector%in => source
            end select
         end if
         call reinjector_input_dict%add(name)
      end if

      call mpi_broadcast_error_flag(err)

    end subroutine init_reinjector_input_source

!........................................................................

    subroutine init_reinjector_input_group(group, reinjector_input_dict, &
         reinjector, err)
      !! Initialises reinjector input from a group.

      use mpi_utils_module, only: mpi_broadcast_error_flag

      type(source_network_group_type), target, intent(in out) :: group
      type(dictionary_type), intent(in out) :: reinjector_input_dict
      type(source_network_reinjector_type), pointer, intent(in out) :: reinjector
      PetscErrorCode, intent(in out) :: err

      if (associated(group%out)) then
         if (present(logfile)) then
            call logfile%write(LOG_LEVEL_ERR, "input", &
                 "duplicate reinjector input: " // trim(group%name))
         end if
         err = 1
         deallocate(reinjector)
      else
         group%out => reinjector
         if (group%rank == 0) reinjector%in => group
         call reinjector_input_dict%add(group%name)
      end if

      call mpi_broadcast_error_flag(err)

    end subroutine init_reinjector_input_group

!........................................................................

    subroutine init_reinjector_outputs(reinjector_json, reinjector_str, &
         flow_type_str, source_dict, source_dict_all, reinjector_output_dict, &
         source_network, output_index, reinjector, err)
      !! Initialises reinjector outputs of the given flow type.

      use mpi_utils_module, only: mpi_broadcast_error_flag

      type(fson_value), pointer, intent(in out) :: reinjector_json
      character(*), intent(in) :: reinjector_str, flow_type_str
      type(dictionary_type), intent(in out) :: source_dict, source_dict_all, &
           reinjector_output_dict
      type(source_network_type), intent(in out) :: source_network
      type(source_network_reinjector_type), intent(in out) :: reinjector
      PetscInt, intent(in out) :: output_index
      PetscInt, allocatable :: source_cell_indices(:)
      PetscErrorCode, intent(in out) :: err
      ! Locals:
      PetscInt :: flow_type, num_outputs, i, out_type
      type(fson_value), pointer :: outputs_json, output_json
      class(specified_reinjector_output_type), pointer :: output
      character(max_source_network_node_name_length) :: out_name
      type(list_node_type), pointer :: source_dict_node, reinjector_dict_node
      PetscMPIInt :: rank

      call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

      if (fson_has_mpi(reinjector_json, flow_type_str)) then

         flow_type = separated_flow_type_from_str(flow_type_str)

         call fson_get_mpi(reinjector_json, trim(flow_type_str), outputs_json)
         num_outputs = fson_value_count_mpi(outputs_json, ".")
         allocate(source_cell_indices(num_outputs))
         source_cell_indices = -1
         output_json => fson_value_children_mpi(outputs_json)

         do i = 1, num_outputs

            if (fson_has_mpi(output_json, "proportion")) then
               allocate(proportion_reinjector_output_type :: output)
               call get_initial_reinjector_output_proportion(output_json, output)
            else
               allocate(rate_reinjector_output_type :: output)
               call get_initial_reinjector_output_rate(output_json, output)
            end if
            call output%init(reinjector, flow_type)
            call get_initial_reinjector_output_enthalpy(output_json, output)

            if (fson_has_mpi(output_json, "out")) then
               out_type = fson_type_mpi(output_json, "out")
               if (out_type == TYPE_NULL) then
                  ! Reinjector output with null source/reinjector:
                  ! assign to global root rank and store index in output
                  if (rank == 0) then
                     output%link_index = output_index
                     call reinjector%out%append(output)
                  end if
               else
                  call fson_get_mpi(output_json, "out", val = out_name)
                  if (reinjector_output_dict%has(out_name)) then
                     if (present(logfile)) then
                        call logfile%write(LOG_LEVEL_ERR, "input", &
                             "duplicate reinjector output: " // trim(out_name))
                     end if
                     err = 1
                     deallocate(output)
                     exit
                  else
                     if (source_dict_all%has(out_name)) then
                        source_dict_node => source_dict%get(out_name)
                        if (associated(source_dict_node)) then
                           ! source is on this process:
                           select type (source => source_dict_node%data)
                           type is (source_type)
                              output%out => source
                              source%link_index = output_index
                              call reinjector%out%append(output)
                              source_cell_indices(i) = source%natural_cell_index
                           end select
                        else ! source is not on this process:
                           deallocate(output)
                        end if
                        call reinjector_output_dict%add(out_name)
                     else
                        reinjector_dict_node => reinjector_dict%get(out_name)
                        if (associated(reinjector_dict_node)) then
                           select type (out_reinjector => reinjector_dict_node%data)
                           class is (source_network_reinjector_type)
                              if (out_reinjector%rank == 0) then
                                 output%out => out_reinjector
                                 out_reinjector%in => output
                                 out_reinjector%link_index = output_index
                                 source_cell_indices = [source_cell_indices, &
                                     out_reinjector%source_cell_indices]
                                 call reinjector%out%append(output)
                              else
                                 deallocate(output)
                              end if
                              call reinjector_output_dict%add(out_name)
                           end select
                        else
                           if (present(logfile)) then
                              call logfile%write(LOG_LEVEL_ERR, "input", &
                                   "unrecognised reinjector output: " // trim(out_name))
                           end if
                           err = 1
                           deallocate(output)
                           exit
                        end if
                     end if
                  end if
               end if
            else
               ! Reinjector output with no source/reinjector assigned:
               ! assign to global root rank and store index in output
               if (rank == 0) then
                  output%link_index = output_index
                  call reinjector%out%append(output)
               end if
            end if

            call setup_source_reinjector_output_controls(output_json, output, &
                 source_network)

            output_json => fson_value_next_mpi(output_json)
            output_index = output_index + 1
         end do

         source_cell_indices = pack(source_cell_indices, &
              source_cell_indices >= 0)
         reinjector%source_cell_indices = [reinjector%source_cell_indices, &
              source_cell_indices]
         deallocate(source_cell_indices)

      end if

      call mpi_broadcast_error_flag(err)

    end subroutine init_reinjector_outputs

!------------------------------------------------------------------------

    subroutine get_initial_reinjector_output_proportion(output_json, output)
      !! Gets initial flow proportion for proportion reinjector output.

      type(fson_value), pointer, intent(in out) :: output_json
      class(specified_reinjector_output_type), intent(in out) :: output
      ! Locals:
      PetscInt :: proportion_type

      select type (prop_output => output)
      class is (proportion_reinjector_output_type)
         proportion_type = fson_type_mpi(output_json, "proportion")
         select case (proportion_type)
         case (TYPE_REAL, TYPE_INTEGER)
            call fson_get_mpi(output_json, "proportion", &
                 val = prop_output%proportion)
         case default
            prop_output%proportion = default_reinjector_output_proportion
         end select
      end select

    end subroutine get_initial_reinjector_output_proportion

!------------------------------------------------------------------------

    subroutine get_initial_reinjector_output_rate(output_json, output)
      !! Gets initial flow rate for rate reinjector output.

      type(fson_value), pointer, intent(in out) :: output_json
      class(specified_reinjector_output_type), intent(in out) :: output
      ! Locals:
      PetscInt :: rate_type

      select type (rate_output => output)
      class is (rate_reinjector_output_type)
         if (fson_has_mpi(output_json, "rate")) then
            rate_type = fson_type_mpi(output_json, "rate")
            select case (rate_type)
            case (TYPE_REAL, TYPE_INTEGER)
               call fson_get_mpi(output_json, "rate", &
                    val = rate_output%specified_rate)
            case default
               rate_output%specified_rate = default_reinjector_output_rate
            end select
         else
            rate_output%specified_rate = default_reinjector_output_rate
         end if
      end select

    end subroutine get_initial_reinjector_output_rate

!------------------------------------------------------------------------

    subroutine get_initial_reinjector_output_enthalpy(output_json, output)
      !! Gets initial enthalpy for rate reinjector output.

      type(fson_value), pointer, intent(in out) :: output_json
      class(specified_reinjector_output_type), intent(in out) :: output
      ! Locals:
      PetscInt :: enthalpy_type

      if (fson_has_mpi(output_json, "enthalpy")) then
         enthalpy_type = fson_type_mpi(output_json, "enthalpy")
         select case (enthalpy_type)
         case (TYPE_REAL, TYPE_INTEGER)
            call fson_get_mpi(output_json, "enthalpy", &
                 val = output%specified_enthalpy)
         case default
            output%specified_enthalpy = default_reinjector_output_enthalpy
         end select
      else
         output%specified_enthalpy = default_reinjector_output_enthalpy
      end if

    end subroutine get_initial_reinjector_output_enthalpy

!------------------------------------------------------------------------

    subroutine setup_rate_reinjector_output_control(output_json, &
         interpolation_type, averaging_type, output, source_network)
      !! Sets up control for rate reinjector output.

      type(fson_value), pointer, intent(in out) :: output_json
      PetscInt, intent(in) :: interpolation_type, averaging_type
      class(specified_reinjector_output_type), pointer, intent(in out) :: output
      class(source_network_type), intent(in out) :: source_network
      ! Locals:
      PetscInt :: rate_type
      PetscReal, allocatable :: data_array(:,:)
      type(reinjector_rate_table_source_network_control_type), pointer :: control
      type(list_type) :: outputs

      if (fson_has_mpi(output_json, "rate")) then
         rate_type = fson_type_mpi(output_json, "rate")
         if (rate_type == TYPE_ARRAY) then
            call fson_get_mpi(output_json, "rate", val = data_array)
            if (associated(output)) then
               allocate(control)
               call outputs%init(owner = PETSC_FALSE)
               call outputs%append(output)
               call control%init(outputs, data_array, interpolation_type, &
                    averaging_type)
               call source_network%network_controls%append(control)
            end if
         end if
      end if

    end subroutine setup_rate_reinjector_output_control

!------------------------------------------------------------------------

    subroutine setup_proportion_reinjector_output_control(output_json, &
         interpolation_type, averaging_type, output, source_network)
      !! Sets up control for proportion reinjector output.

      type(fson_value), pointer, intent(in out) :: output_json
      PetscInt, intent(in) :: interpolation_type, averaging_type
      class(specified_reinjector_output_type), pointer, intent(in out) :: output
      class(source_network_type), intent(in out) :: source_network
      ! Locals:
      PetscInt :: proportion_type
      PetscReal, allocatable :: data_array(:,:)
      type(reinjector_proportion_table_source_network_control_type), pointer :: control
      type(list_type) :: outputs

      if (fson_has_mpi(output_json, "proportion")) then
         proportion_type = fson_type_mpi(output_json, "proportion")
         if (proportion_type == TYPE_ARRAY) then
            call fson_get_mpi(output_json, "proportion", val = data_array)
            if (associated(output)) then
               allocate(control)
               call outputs%init(owner = PETSC_FALSE)
               call outputs%append(output)
               call control%init(outputs, data_array, interpolation_type, &
                    averaging_type)
               call source_network%network_controls%append(control)
            end if
         end if
      end if

    end subroutine setup_proportion_reinjector_output_control

!------------------------------------------------------------------------

    subroutine setup_enthalpy_reinjector_output_control(output_json, &
         interpolation_type, averaging_type, output, source_network)
      !! Sets up enthalpy control for reinjector output.

      type(fson_value), pointer, intent(in out) :: output_json
      PetscInt, intent(in) :: interpolation_type, averaging_type
      class(specified_reinjector_output_type), pointer, intent(in out) :: output
      class(source_network_type), intent(in out) :: source_network
      ! Locals:
      PetscInt :: enthalpy_type
      PetscReal, allocatable :: data_array(:,:)
      type(reinjector_enthalpy_table_source_network_control_type), pointer :: control
      type(list_type) :: outputs

      if (fson_has_mpi(output_json, "enthalpy")) then
         enthalpy_type = fson_type_mpi(output_json, "enthalpy")
         if (enthalpy_type == TYPE_ARRAY) then
            call fson_get_mpi(output_json, "enthalpy", val = data_array)
            if (associated(output)) then
               allocate(control)
               call outputs%init(owner = PETSC_FALSE)
               call outputs%append(output)
               call control%init(outputs, data_array, interpolation_type, &
                    averaging_type)
               call source_network%network_controls%append(control)
            end if
         end if
      end if

    end subroutine setup_enthalpy_reinjector_output_control

!------------------------------------------------------------------------

    subroutine setup_source_reinjector_output_controls(output_json, output, &
         source_network)
      !! Sets up controls on reinjector outputs.

      use interpolation_module, only: interpolation_type_from_str, &
           averaging_type_from_str, max_interpolation_str_length, &
           max_averaging_str_length, default_interpolation_str, &
           default_averaging_str

      type(fson_value), pointer, intent(in out) :: output_json
      class(specified_reinjector_output_type), pointer, intent(in out) :: output
      class(source_network_type), intent(in out) :: source_network
      ! Locals:
      PetscInt :: interpolation_type, averaging_type
      character(max_interpolation_str_length) :: interpolation_str
      character(max_averaging_str_length) :: averaging_str

      call fson_get_mpi(output_json, "interpolation", &
           default_interpolation_str, interpolation_str)
      interpolation_type = interpolation_type_from_str(interpolation_str)

      call fson_get_mpi(output_json, "averaging", &
           default_averaging_str, averaging_str)
      averaging_type = averaging_type_from_str(averaging_str)

      call setup_rate_reinjector_output_control(output_json, &
           interpolation_type, averaging_type, output, source_network)

      call setup_proportion_reinjector_output_control(output_json, &
           interpolation_type, averaging_type, output, source_network)

      call setup_enthalpy_reinjector_output_control(output_json, &
           interpolation_type, averaging_type, output, source_network)

    end subroutine setup_source_reinjector_output_controls

!------------------------------------------------------------------------

    subroutine init_reinjector_overflow(reinjector_json, reinjector_str, &
         source_dict, source_dict_all, reinjector_output_dict, &
         reinjector, err)
      !! Initialises reinjector overflow.

      use mpi_utils_module, only: mpi_broadcast_error_flag

      type(fson_value), pointer, intent(in out) :: reinjector_json
      character(*), intent(in) :: reinjector_str
      type(dictionary_type), intent(in out) :: source_dict, source_dict_all, &
           reinjector_output_dict
      type(source_network_reinjector_type), intent(in out) :: reinjector
      PetscErrorCode, intent(out) :: err
      ! Locals:
      PetscInt :: overflow_type
      PetscInt, allocatable :: source_cell_indices(:)
      character(max_source_network_node_name_length) :: node_name
      type(list_node_type), pointer :: source_dict_node, reinjector_dict_node

      err = 0
      allocate(source_cell_indices(0))
      if (fson_has_mpi(reinjector_json, "overflow")) then
         overflow_type = fson_type_mpi(reinjector_json, "overflow")
         node_name = ""
         select case (overflow_type)
         case (TYPE_STRING)
            call fson_get_mpi(reinjector_json, "overflow", val = node_name)
         case (TYPE_OBJECT)
            if (fson_has_mpi(reinjector_json, "overflow.out")) then
               call fson_get_mpi(reinjector_json, "overflow.out", &
                    val = node_name)
            end if
         end select
         if (node_name /= "") then
            if (reinjector_output_dict%has(node_name)) then
               if (present(logfile)) then
                  call logfile%write(LOG_LEVEL_ERR, "input", &
                       "duplicate reinjector output: " // trim(node_name))
               end if
               err = 1
            else
               if (source_dict_all%has(node_name)) then
                  source_dict_node => source_dict%get(node_name)
                  if (associated(source_dict_node)) then
                     ! source is on this process:
                     select type (source => source_dict_node%data)
                     type is (source_type)
                        reinjector%overflow%out => source
                        source%link_index = 0 ! not used
                        source_cell_indices = [source%natural_cell_index]
                     end select
                  end if
                  call reinjector_output_dict%add(node_name)
               else
                  reinjector_dict_node => reinjector_dict%get(node_name)
                  if (associated(reinjector_dict_node)) then
                     select type (out_reinjector => reinjector_dict_node%data)
                     class is (source_network_reinjector_type)
                        if (out_reinjector%rank == 0) then
                           reinjector%overflow%out => out_reinjector
                           out_reinjector%in => reinjector%overflow
                           out_reinjector%link_index = 0 ! not used
                           source_cell_indices = out_reinjector%source_cell_indices
                        end if
                        call reinjector_output_dict%add(node_name)
                     end select
                  else
                     if (present(logfile)) then
                        call logfile%write(LOG_LEVEL_ERR, "input", &
                             "unrecognised reinjector output: " // trim(node_name))
                     end if
                     err = 1
                  end if
               end if
            end if
         end if
      end if

      call mpi_broadcast_error_flag(err)

      if (err == 0) then
         call reinjector%overflow%get_world_rank()
         reinjector%source_cell_indices = [reinjector%source_cell_indices, &
              source_cell_indices]
      end if

    end subroutine init_reinjector_overflow

!........................................................................

    subroutine init_source_network_reinjectors(source_network, &
         num_local_root_reinjectors, logfile, err)
      !! Initialise source network reinjectors and return number of
      !! local root reinjectors. An error is returned if any
      !! unrecognised reinjector output nodes are specified, a
      !! reinjector is not assigned an input node name, or a
      !! reinjector input is used by more than one reinjector.

      use utils_module, only: str_to_lower

      type(source_network_type), intent(in out) :: source_network
      PetscInt, intent(in out) :: num_local_root_reinjectors
      type(logfile_type), intent(in out), optional :: logfile
      PetscErrorCode, intent(out) :: err
      ! Locals:
      PetscInt :: ir, reinjector_index, r, output_index
      type(fson_value), pointer :: reinjector_json
      character(max_source_network_node_name_length) :: name, in_name
      character(len=64) :: rstr
      character(len=12) :: irstr
      type(dictionary_type) :: reinjector_input_dict, reinjector_output_dict
      type(source_network_reinjector_type), pointer :: reinjector
      type(list_node_type), pointer :: group_dict_node

      err = 0
      call reinjector_input_dict%init(PETSC_FALSE)
      call reinjector_output_dict%init(PETSC_FALSE)

      do ir = 0, source_network%num_reinjectors - 1

         reinjector_index = reinjector_order(ir) + 1
         reinjector_json => reinjector_specs_array(reinjector_index)%ptr
         write(irstr, '(i0)') reinjector_index - 1
         rstr = 'network.reinject[' // trim(irstr) // '].'

         call fson_get_mpi(reinjector_json, "name", "", name)

         allocate(reinjector)
         call reinjector%init(name)

         if (fson_has_mpi(reinjector_json, "in")) then

            call fson_get_mpi(reinjector_json, "in", val = in_name)
            if (source_dict_all%has(in_name)) then
               call init_reinjector_input_source(in_name, reinjector_input_dict, &
                    source_dict, reinjector, err)
               if (err > 0) exit
            else
               group_dict_node => group_dict%get(in_name)
               if (associated(group_dict_node)) then
                  select type (group => group_dict_node%data)
                  class is (source_network_group_type)
                     call init_reinjector_input_group(group, reinjector_input_dict, &
                          reinjector, err)
                     if (err > 0) exit
                  end select
               else
                  if (present(logfile)) then
                     call logfile%write(LOG_LEVEL_ERR, "input", &
                          "unrecognised reinjector input: " // trim(in_name))
                  end if
                  err = 1
                  deallocate(reinjector)
                  exit
               end if
            end if

         end if

         output_index = 1
         call init_reinjector_outputs(reinjector_json, rstr, "water", &
              source_dict, source_dict_all, reinjector_output_dict, &
              source_network, output_index, reinjector, err)
         if (err == 0) then
            call init_reinjector_outputs(reinjector_json, rstr, "steam", &
                 source_dict, source_dict_all, reinjector_output_dict, &
                 source_network, output_index, reinjector, err)
            if (err == 0) then
               call init_reinjector_overflow(reinjector_json, rstr, &
                 source_dict, source_dict_all, reinjector_output_dict, &
                 reinjector, err)
               if (err == 0) then
                  call reinjector%init_comm()
                  call reinjector%overflow%allocate_variables()
                  if (reinjector%rank == 0) then
                     r = num_local_root_reinjectors
                     num_local_root_reinjectors = num_local_root_reinjectors + 1
                  else
                     r = -1
                  end if
                  reinjector%local_reinjector_index = r
                  reinjector%source_cell_indices = reinjector%gatherv( &
                       reinjector%source_cell_indices)
                  call source_network%reinjectors%append(reinjector)
                  if (name /= "") then
                     call reinjector_dict%add(name, reinjector)
                  end if
               else
                  deallocate(reinjector)
                  exit
               end if
            else
               deallocate(reinjector)
               exit
            end if
         else
            deallocate(reinjector)
            exit
         end if

      end do

      call source_network%reinjectors%traverse(init_reinjector_in_comms_iterator)
      call reinjector_input_dict%destroy()
      call reinjector_output_dict%destroy()

    end subroutine init_source_network_reinjectors

!........................................................................

    subroutine init_reinjector_in_comms_iterator(node, stopped)
      !! Initialises MPI communicators for all reinjectors. Note this
      !! cannot be done until after all reinjectors have been created
      !! (the inputs for some reinjectors are assigned only during the
      !! creation of the outputs of other reinjectors).

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type(reinjector => node%data)
      class is (source_network_reinjector_type)
         call reinjector%init_in_comm()
      end select

    end subroutine init_reinjector_in_comms_iterator

  end subroutine setup_source_network

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

  subroutine get_initial_rate(source_json, initial_rate, rate_specified, &
       specified_rate)
    !! Gets initial flow rate. This can only be determined for
    !! constant-rate sources. For other source types, a default
    !! initial rate is assigned, which may be modified by any source
    !! controls acting on the source.

    type(fson_value), pointer, intent(in) :: source_json
    PetscReal, intent(out) :: initial_rate
    PetscBool, intent(out) :: rate_specified
    PetscReal, intent(out) :: specified_rate
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

       rate_specified = PETSC_TRUE

    else
       initial_rate = default_source_rate
       rate_specified = PETSC_FALSE
    end if

    specified_rate = initial_rate

  end subroutine get_initial_rate

!------------------------------------------------------------------------

  subroutine get_initial_enthalpy(source_json, eos, injection_component, &
       enthalpy, enthalpy_specified, specified_enthalpy)
    !! Gets initial injection enthalpy for the source, if needed.

    type(fson_value), pointer, intent(in) :: source_json
    class(eos_type), intent(in) :: eos
    PetscInt, intent(in) :: injection_component
    PetscReal, intent(out) :: enthalpy
    PetscBool, intent(out) :: enthalpy_specified
    PetscReal, intent(out) :: specified_enthalpy
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

            enthalpy_specified = PETSC_TRUE

         else
            enthalpy = default_source_injection_enthalpy
            enthalpy_specified = PETSC_FALSE
         end if

      else ! heat injection - no enthalpy needed
         enthalpy = 0._dp
         enthalpy_specified = PETSC_FALSE
      end if

    end associate

    specified_enthalpy = enthalpy

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
    !! Gets separator pressure array. Returns [-1] if no separator is
    !! specified.

    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: source_json
    character(*), intent(in) :: srcstr !! source identifier
    PetscReal, allocatable, intent(out) :: separator_pressure(:) !! Separator pressures
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscInt :: separator_json_type, separator_pressure_type
    PetscBool :: has_separator
    character(8) :: limiter_type_str
    PetscReal :: scalar_separator_pressure

    if (fson_has_mpi(source_json, "separator")) then
       separator_json_type = fson_type_mpi(source_json, "separator")
       select case (separator_json_type)
       case (TYPE_LOGICAL)
          call fson_get_mpi(source_json, "separator", val = has_separator)
          if (has_separator) then
             separator_pressure = [default_separator_pressure]
             if (present(logfile) .and. logfile%active) then
                call logfile%write(LOG_LEVEL_INFO, 'input', 'default', real_keys = &
                     [trim(srcstr) // "separator.pressure"], &
                     real_values = separator_pressure)
             end if
          end if
       case (TYPE_OBJECT)
          if (fson_has_mpi(source_json, "separator.pressure")) then
             separator_pressure_type = fson_type_mpi(source_json, "separator.pressure")
             select case (separator_pressure_type)
             case (TYPE_REAL, TYPE_INTEGER)
                call fson_get_mpi(source_json, "separator.pressure", &
                     val = scalar_separator_pressure)
                separator_pressure = [scalar_separator_pressure]
             case (TYPE_ARRAY)
                call fson_get_mpi(source_json, "separator.pressure", &
                     val = separator_pressure)
             end select
          else
             separator_pressure = [default_separator_pressure]
             if (present(logfile) .and. logfile%active) then
                call logfile%write(LOG_LEVEL_INFO, 'input', 'default', real_keys = &
                     [trim(srcstr) // "separator.pressure"], &
                     real_values = separator_pressure)
             end if
          end if
       end select
    else if (fson_has_mpi(source_json, "limiter")) then
       call fson_get_mpi(source_json, "limiter.type", &
            default_source_control_limiter_type_str, limiter_type_str)
       limiter_type_str = str_to_lower(limiter_type_str)
       if (limiter_type_str /= "total") then
          separator_pressure_type = fson_type_mpi(source_json, &
               "limiter.separator_pressure")
          select case (separator_pressure_type)
          case (TYPE_REAL, TYPE_INTEGER)
             call fson_get_mpi(source_json, "limiter.separator_pressure", &
                  val = scalar_separator_pressure)
             separator_pressure = [scalar_separator_pressure]
          case (TYPE_ARRAY)
             call fson_get_mpi(source_json, "limiter.separator_pressure", &
                  val = separator_pressure)
          end select
       end if
    end if

    if (.not. allocated(separator_pressure)) then
       ! no separator:
       separator_pressure = [-1._dp]
    end if

  end subroutine get_separator_pressure

!------------------------------------------------------------------------

  subroutine setup_inline_source_controls(source_json, eos, thermo, &
       start_time, source_data, source_section, &
       fluid_data, fluid_section, fluid_range_start, &
       srcstr, tracer_names, null_cell, spec_sources, source_network, &
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
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start
    character(len = *), intent(in) :: srcstr
    character(len = *), intent(in) :: tracer_names(:)
    PetscBool, intent(in) :: null_cell
    type(list_type), intent(in out) :: spec_sources
    type(source_network_type), intent(in out) :: source_network
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
         averaging_type, tracer_names, spec_sources, source_network, &
         logfile, err)

    if (err == 0) then

       call setup_deliverability_source_controls(source_json, srcstr, &
            start_time, source_data, source_section, &
            fluid_data, fluid_section, fluid_range_start, &
            interpolation_type, averaging_type, spec_sources, eos, &
            null_cell, source_network, logfile, err)

       if (err == 0) then

          call setup_recharge_source_controls(source_json, srcstr, &
               source_data, source_section, &
               fluid_data, fluid_section, fluid_range_start, &
               interpolation_type, averaging_type, spec_sources, eos, &
               null_cell, source_network, logfile, err)

          if (err == 0) then

             call setup_limiter_source_controls(source_json, srcstr, thermo, &
                  interpolation_type, averaging_type, spec_sources, &
                  source_network, logfile)

             call setup_direction_source_control(source_json, srcstr, thermo, &
                  null_cell, spec_sources, source_network, logfile, err)

             if (err == 0) then
                call setup_factor_source_control(source_json, srcstr, &
                     interpolation_type, averaging_type, spec_sources, &
                     source_network, logfile, err)
             end if
          end if

       end if

    end if

  end subroutine setup_inline_source_controls

!------------------------------------------------------------------------

    subroutine setup_table_source_control(source_json, srcstr, &
       interpolation_type, averaging_type, tracer_names, &
       spec_sources, source_network, logfile, err)
    !! Set up rate, enthalpy and tracer table source controls.

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    PetscInt, intent(in) :: interpolation_type, averaging_type
    character(*), intent(in) :: tracer_names(:)
    type(list_type), intent(in out) :: spec_sources
    type(source_network_type), intent(in out) :: source_network
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err

    err = 0
    call setup_rate_table_control()
    call setup_enthalpy_table_control()
    call setup_tracer_table_controls()

  contains

!........................................................................

    subroutine setup_rate_table_control()

      PetscInt :: variable_type
      type(fson_value), pointer :: table
      type(rate_table_source_control_type), pointer :: control
      PetscReal :: fixed_rate
      PetscReal, allocatable :: data_array(:,:)

      if (fson_has_mpi(source_json, "rate")) then
         variable_type = fson_type_mpi(source_json, "rate")
         select case (variable_type)
         case (TYPE_REAL, TYPE_INTEGER)
            call fson_get_mpi(source_json, "rate", val = fixed_rate)
            data_array = reshape([0._dp, fixed_rate], [1, 2])
         case (TYPE_ARRAY)
            call fson_get_mpi(source_json, "rate", val = data_array)
         case (TYPE_OBJECT)
            call fson_get_mpi(source_json, "rate", table)
            if (fson_has_mpi(table, "time")) then
               call fson_get_mpi(table, "time", val = data_array)
            end if
         end select
      end if

      if (allocated(data_array)) then
         if (spec_sources%count > 0) then
            allocate(control)
            call control%init(spec_sources%copy(), data_array, &
                 interpolation_type, averaging_type)
            call source_network%source_controls%append(control)
         end if
         deallocate(data_array)
      end if

    end subroutine setup_rate_table_control

!........................................................................

    subroutine setup_enthalpy_table_control()

      PetscInt :: variable_type
      type(fson_value), pointer :: table
      type(enthalpy_table_source_control_type), pointer :: control
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
            call control%init(spec_sources%copy(), data_array, &
                 interpolation_type, averaging_type)
            call source_network%source_controls%append(control)
         end if
         deallocate(data_array)
      end if

    end subroutine setup_enthalpy_table_control

!........................................................................

    subroutine setup_tracer_table_controls()

      use utils_module, only: str_array_index
      use tracer_module, only: max_tracer_name_length

      PetscInt :: variable_type, rank
      type(fson_value), pointer :: tracers_json, tracer_json
      PetscInt :: num_tracers_specified, i, tracer_type, tracer_index
      type(tracer_table_source_control_type), pointer :: control
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
                        call control%init(spec_sources%copy(), data_array, &
                             interpolation_type, averaging_type)
                        control%tracer_index = i
                        call source_network%source_controls%append(control)
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
                           call control%init(spec_sources%copy(), data_array, &
                                interpolation_type, averaging_type)
                           control%tracer_index = tracer_index
                           call source_network%source_controls%append(control)
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
       source_network, logfile, err)
    !! Set up rate factor source controls.

    use interpolation_module, only: interpolation_type_from_str, &
         averaging_type_from_str, max_interpolation_str_length, &
         max_averaging_str_length

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    PetscInt, intent(in) :: interpolation_type, averaging_type
    type(list_type), intent(in out) :: spec_sources
    type(source_network_type), intent(in out) :: source_network
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(rate_factor_table_source_control_type), pointer :: factor_control
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
       call factor_control%init(spec_sources%copy(), factor_data_array, &
            effective_interpolation_type, effective_averaging_type)
       call source_network%source_controls%append(factor_control)
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
       start_time, source_data, source_section, &
       fluid_data, fluid_section, fluid_range_start, &
       interpolation_type, averaging_type, spec_sources, eos, null_cell, &
       source_network, logfile, err)
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
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start
    PetscInt, intent(in) :: interpolation_type, averaging_type
    type(list_type), intent(in out) :: spec_sources
    class(eos_type), intent(in) :: eos
    PetscBool, intent(in) :: null_cell
    type(source_network_type), intent(in out) :: source_network
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(fson_value), pointer :: deliv_json
    PetscReal, allocatable :: reference_pressure_array(:,:)
    type(deliverability_source_control_type), pointer :: deliv
    PetscBool :: calculate_reference_pressure
    PetscBool :: calculate_PI_from_rate
    PetscReal :: initial_rate, threshold
    PetscReal, allocatable :: productivity_array(:,:)
    PetscInt :: pressure_table_coordinate
    PetscReal, parameter :: default_rate = 0._dp
    PetscReal, parameter :: default_threshold = -1._dp

    err = 0

    if (fson_has_mpi(source_json, "deliverability")) then

       if (null_cell) then
          if (present(logfile) .and. logfile%active) then
             call logfile%write(LOG_LEVEL_ERR, 'input', 'invalid', &
                  str_key = trim(srcstr) // "deliverability", &
                  str_value = "... (null cell)")
          end if
          err = 1
       else

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
         call deliv%init(single_source, productivity_array, &
              interpolation_type, averaging_type, reference_pressure_array, &
               pressure_table_coordinate, threshold)

          if (calculate_reference_pressure) then
             call deliv%set_reference_pressure_initial(fluid_data, &
                  fluid_section, fluid_range_start)
          end if
          if (calculate_PI_from_rate) then
             call deliv%calculate_PI_from_rate(start_time, initial_rate, &
                  fluid_data, fluid_section, fluid_range_start, &
                  deliv%productivity%val(1, 1))
          end if
          if (deliv%threshold > 0._dp) then
             deliv%threshold_productivity = &
                  deliv%productivity%interpolate(start_time, 1)
          end if
          source%rate_specified = PETSC_TRUE
          call source_network%source_controls%append(deliv)

      end select

    end subroutine setup_deliverability_iterator

  end subroutine setup_deliverability_source_controls

!------------------------------------------------------------------------

  subroutine get_recharge_coefficient(json, srcstr, key, source_json, &
       recharge_array, logfile)
    !! Gets recharge/injectivity coefficient for recharge source control.

    type(fson_value), pointer, intent(in) :: json, source_json
    character(len = *), intent(in) :: srcstr, key
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
                  str_key = trim(srcstr) // trim(key) // ".coefficient", &
                  str_value = "...")
          end if
       end select

    else
       if (present(logfile) .and. logfile%active) then
          call logfile%write(LOG_LEVEL_INFO, 'input', 'default', real_keys = &
               [trim(srcstr) // trim(key) // ".coefficient"], &
               real_values = [recharge_coefficient])
       end if
    end if

  end subroutine get_recharge_coefficient

!------------------------------------------------------------------------

  subroutine setup_recharge_source_controls(source_json, srcstr, &
       source_data, source_section, &
       fluid_data, fluid_section, fluid_range_start, &
       interpolation_type, averaging_type, spec_sources, eos, &
       null_cell, source_network, logfile, err)
    !! Set up recharge/injectivity source controls. These controls
    !! can control only one source, so if multiple cells are
    !! specified, multiple corresponding recharge controls are
    !! created.

    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    PetscReal, pointer, contiguous, intent(in) :: source_data(:)
    PetscSection, intent(in) :: source_section
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    PetscInt, intent(in) :: fluid_range_start
    PetscInt, intent(in) :: interpolation_type, averaging_type
    type(list_type), intent(in out) :: spec_sources
    class(eos_type), intent(in) :: eos
    PetscBool, intent(in) :: null_cell
    type(source_network_type), intent(in out) :: source_network
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(fson_value), pointer :: recharge_json
    PetscReal, allocatable :: reference_pressure_array(:,:)
    type(recharge_source_control_type), pointer :: recharge
    PetscBool :: calculate_reference_pressure
    PetscReal, allocatable :: recharge_array(:,:)
    PetscInt :: pressure_table_coordinate, k
    PetscInt, parameter :: num_keys = 2
    character(12), parameter :: keys(num_keys) = ["recharge   ", "injectivity"]
    character(12) :: key

    do k = 1, num_keys
       key = keys(k)
       if (fson_has_mpi(source_json, key)) then

          if (null_cell) then
             if (present(logfile) .and. logfile%active) then
                call logfile%write(LOG_LEVEL_ERR, 'input', 'invalid', &
                     str_key = trim(srcstr) // key, &
                     str_value = "... (null cell)")
             end if
             err = 1
             exit
          else

             call fson_get_mpi(source_json, key, recharge_json)

             call get_reference_pressure(recharge_json, srcstr, key, &
                  reference_pressure_array, calculate_reference_pressure, &
                  pressure_table_coordinate, logfile)

             if (pressure_table_coordinate == SRC_PRESSURE_TABLE_COORD_TIME) then

                call get_recharge_coefficient(recharge_json, srcstr, key, &
                     source_json, recharge_array, logfile)

                call spec_sources%traverse(setup_recharge_iterator)

             else
                if (present(logfile) .and. logfile%active) then
                   call logfile%write(LOG_LEVEL_WARN, 'input', 'not_supported', &
                        str_key = trim(srcstr) // trim(key) // ".pressure", &
                        str_value = "...")
                end if
             end if

             exit
          end if

       end if
    end do

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
         call recharge%init(single_source, recharge_array, interpolation_type, &
              averaging_type, reference_pressure_array)

         if (calculate_reference_pressure) then
            call recharge%set_reference_pressure_initial(fluid_data, &
                 fluid_section, fluid_range_start)
         end if
         source%rate_specified = PETSC_TRUE

         call source_network%source_controls%append(recharge)

      end select

    end subroutine setup_recharge_iterator

  end subroutine setup_recharge_source_controls

!------------------------------------------------------------------------

  subroutine setup_limiter_source_controls(source_json, srcstr, &
       thermo, interpolation_type, averaging_type, spec_sources, &
       source_network, logfile)
    !! Set up limiter source control for each cell source.

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    class(thermodynamics_type), intent(in) :: thermo
    PetscInt, intent(in) :: interpolation_type, averaging_type
    type(list_type), intent(in out) :: spec_sources
    type(source_network_type), intent(in out) :: source_network
    type(logfile_type), intent(in out), optional :: logfile

    if (fson_has_mpi(source_json, "limiter")) then
       call add_limiter(source_json, srcstr, spec_sources, &
            interpolation_type, averaging_type, PETSC_TRUE, &
            source_network%source_controls, logfile)
    end if

  end subroutine setup_limiter_source_controls

!------------------------------------------------------------------------

  subroutine add_limiter(node_json, name, network_nodes, &
       default_interpolation_type, default_averaging_type, &
       require_nodes, controls, logfile)
    !! Creates a limiter control with the specified network nodes and
    !! parameters from JSON input and adds it to the specified list of
    !! controls.

    use interpolation_module, only: interpolation_type_from_str, &
         averaging_type_from_str, max_interpolation_str_length, &
         max_averaging_str_length

    type(fson_value), pointer, intent(in) :: node_json
    character(*), intent(in) :: name
    type(list_type), intent(in out) :: network_nodes
    PetscInt, intent(in) :: default_interpolation_type, default_averaging_type
    PetscBool, intent(in) :: require_nodes
    type(list_type), intent(in out) :: controls
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    type(fson_value), pointer :: limiter_json
    PetscInt :: flow_type, variable_type
    PetscReal, allocatable :: limit_data_array(:, :)
    PetscReal :: const_limit
    character(max_limiter_type_length) :: limiter_type_str
    character(max_averaging_str_length) :: averaging_str
    character(max_interpolation_str_length) :: interpolation_str
    PetscInt :: interpolation_type, averaging_type
    PetscBool :: create
    type(limiter_table_source_network_control_type), pointer :: limiter

    call fson_get_mpi(node_json, "limiter", limiter_json)
    call get_table_parameters()
    if (fson_has_mpi(limiter_json, "limit") .or. &
         (fson_has_mpi(limiter_json, "type"))) then
       call add_limiter_single()
    else
       call add_limiter_multi()
    end if

  contains

    subroutine get_table_parameters()

      interpolation_type = default_interpolation_type
      averaging_type = default_averaging_type
      if (fson_has_mpi(limiter_json, "interpolation")) then
         call fson_get_mpi(limiter_json, "interpolation", val = interpolation_str)
         interpolation_type = interpolation_type_from_str(interpolation_str)
      end if
      if (fson_has_mpi(limiter_json, "averaging")) then
         call fson_get_mpi(limiter_json, "averaging", val = averaging_str)
         averaging_type = averaging_type_from_str(averaging_str)
      end if

    end subroutine get_table_parameters

!........................................................................

    function get_limit_data(name) result(data)

      character(*), intent(in) :: name !! Limiter variable key
      PetscReal, allocatable :: data(:, :)

      if (fson_has_mpi(limiter_json, name)) then
         variable_type = fson_type_mpi(limiter_json, name)
         select case (variable_type)
         case (TYPE_REAL, TYPE_INTEGER)
            call fson_get_mpi(limiter_json, name, &
                 default_source_control_limiter_limit, val = const_limit)
            allocate(data(1, 2))
            data(1, :) = [0._dp, const_limit]
         case (TYPE_ARRAY)
            call fson_get_mpi(limiter_json, name, val = data)
         end select
      else
         allocate(data(1, 2))
         data(1, :) = [0._dp, const_limit]
         call logfile%write(LOG_LEVEL_INFO, 'input', 'default', real_keys = &
              [trim(name) // "limiter." // trim(name)], &
              real_values = [const_limit])
      end if

    end function get_limit_data

!........................................................................

    subroutine add_limiter_single()
      !! Adds a limiter for a single flow type, using original limiter
      !! JSON syntax in which a flow type and limit are specified.

      call fson_get_mpi(limiter_json, "type", &
           default_source_control_limiter_type_str, &
           limiter_type_str, logfile, name)
      flow_type = separated_flow_type_from_str(limiter_type_str)
      limit_data_array = get_limit_data("limit")
      create = ((require_nodes .and. (network_nodes%count > 0)) .or. &
           (.not. require_nodes))

      if (create) then
         allocate(limiter)
         call limiter%init(network_nodes%copy(), 1)
         limiter%flow_type = [flow_type]
         call limiter%init_table(1, limit_data_array, interpolation_type, &
              averaging_type)
         call controls%append(limiter)
      end if

    end subroutine add_limiter_single

!........................................................................

    subroutine add_limiter_multi()
      !! Adds a limiter for multiple flow types, with JSON syntax of
      !! (key, value) pairs with the keys as flow type strings and
      !! values as limits.

      ! Locals:
      PetscInt, parameter :: num_limiter_types = 3
      character(max_limiter_type_length), parameter :: &
           possible_limiter_types(num_limiter_types) = ["total", "water", "steam"]
      character(max_limiter_type_length) :: limiter_type_str
      PetscInt :: i, num_types
      character(max_limiter_type_length), allocatable :: limiter_types(:)
      PetscInt, allocatable :: flow_types(:)

      num_types = 0
      allocate(flow_types(num_limiter_types), limiter_types(num_limiter_types))
      flow_types = -1
      do i = 1, num_limiter_types
         limiter_type_str = possible_limiter_types(i)
         if (fson_has_mpi(limiter_json, limiter_type_str)) then
            num_types = num_types + 1
            limiter_types(num_types) = limiter_type_str
            flow_types(num_types) = separated_flow_type_from_str( &
                 limiter_type_str)
         end if
      end do
      flow_types = flow_types(1: num_types)
      limiter_types = limiter_types(1: num_types)

      create = (((require_nodes .and. (network_nodes%count > 0)) .or. &
           (.not. require_nodes)) .and. (num_types > 0))

      if (create) then
         allocate(limiter)
         call limiter%init(network_nodes%copy(), num_types)
         limiter%flow_type = flow_types
      end if

      do i = 1, num_types
         limit_data_array = get_limit_data(limiter_types(i))
         if (create) then
            call limiter%init_table(i, limit_data_array, &
                 interpolation_type, averaging_type)
         end if
      end do

      if (create) then
         call controls%append(limiter)
      end if

    end subroutine add_limiter_multi

  end subroutine add_limiter

!------------------------------------------------------------------------

  subroutine setup_direction_source_control(source_json, srcstr, &
       thermo, null_cell, spec_sources, source_network, logfile, err)
    !! Set up direction source control. This can control multiple
    !! sources, so only one is created.

    use utils_module, only: str_to_lower

    type(fson_value), pointer, intent(in) :: source_json
    character(len=*) :: srcstr
    class(thermodynamics_type), intent(in) :: thermo
    PetscBool, intent(in) :: null_cell
    type(list_type), intent(in out) :: spec_sources
    type(source_network_type), intent(in out) :: source_network
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: direction
    type(direction_source_control_type), pointer :: direction_control
    PetscInt, parameter :: max_direction_str_length = 16
    character(max_direction_str_length) :: direction_str

    err = 0
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

       if (null_cell .and. (direction /= SRC_DIRECTION_INJECTION)) then
          if (present(logfile) .and. logfile%active) then
             call logfile%write(LOG_LEVEL_ERR, 'input', 'invalid', &
                  str_key = trim(srcstr) // "direction", &
                  str_value = trim(direction_str) // " (null cell)")
          end if
          err = 1
       else
          if ((direction /= SRC_DIRECTION_BOTH) .and. &
               spec_sources%count > 0) then
             allocate(direction_control)
             call direction_control%init(spec_sources%copy(), direction)
             call source_network%source_controls%append(direction_control)
          end if
       end if
    end if

  end subroutine setup_direction_source_control

!------------------------------------------------------------------------

  subroutine setup_inline_source_group_controls(group_json, group, &
       source_network, logfile)
    !! Sets up any 'inline' source group controls for the source group
    !! specification, i.e. controls defined implicitly in the
    !! specification, and adds them to the source_network network_controls
    !! list.

    use interpolation_module, only: interpolation_type_from_str, &
         averaging_type_from_str, max_interpolation_str_length, &
         max_averaging_str_length, default_interpolation_str, &
         default_averaging_str

    type(fson_value), pointer, intent(in) :: group_json
    class(source_network_group_type), intent(in) :: group
    type(source_network_type), intent(in out) :: source_network
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    PetscInt :: interpolation_type, averaging_type
    character(max_interpolation_str_length) :: interpolation_str
    character(max_averaging_str_length) :: averaging_str

    call fson_get_mpi(group_json, "interpolation", &
         default_interpolation_str, interpolation_str)
    interpolation_type = interpolation_type_from_str(interpolation_str)

    call fson_get_mpi(group_json, "averaging", &
         default_averaging_str, averaging_str)
    averaging_type = averaging_type_from_str(averaging_str)

    call setup_group_limiter_source_controls(group_json, group, &
         interpolation_type, averaging_type, source_network, &
         logfile)

  end subroutine setup_inline_source_group_controls

!------------------------------------------------------------------------

  subroutine setup_group_limiter_source_controls(group_json, group, &
       interpolation_type, averaging_type, source_network, &
       logfile)
    !! Set up limiter control on source group.

    type(fson_value), pointer, intent(in) :: group_json
    class(source_network_group_type), intent(in) :: group
    PetscInt, intent(in) :: interpolation_type, averaging_type
    type(source_network_type), intent(in out) :: source_network
    type(logfile_type), intent(in out), optional :: logfile
    ! Locals:
    type(list_type) :: group_list

    if (fson_has_mpi(group_json, "limiter")) then

       call group_list%init(owner = PETSC_FALSE)
       call group_list%append(group)

       call add_limiter(group_json, group%name, group_list, &
            interpolation_type, averaging_type, PETSC_FALSE, &
            source_network%network_controls, logfile)

       call group_list%destroy()

    end if

  end subroutine setup_group_limiter_source_controls

!------------------------------------------------------------------------

end module source_setup_module
