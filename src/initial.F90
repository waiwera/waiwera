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

module initial_module
  !! Module for setting up initial conditions.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use fson
  use fson_mpi_module

  implicit none
  private

  public :: setup_initial, scale_initial_primary

contains

!------------------------------------------------------------------------

  subroutine setup_initial_primary_constant(mesh, primary, &
       eos, y, y_range_start, logfile, err)

    !! Initializes solution vector y with constant values over the mesh.

    use dm_utils_module, only: global_vec_section, global_section_offset, &
         dm_get_end_interior_cell
    use mesh_module, only: mesh_type
    use eos_module, only: eos_type
    use logfile_module

    type(mesh_type), intent(in) :: mesh
    PetscReal, intent(in) :: primary(:)
    class(eos_type), intent(in) :: eos
    Vec, intent(in out) :: y
    PetscInt, intent(in) :: y_range_start
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: np, c, ghost, primary_size
    PetscInt :: start_cell, end_cell, end_interior_cell, y_offset
    PetscErrorCode :: ierr
    PetscReal, pointer, contiguous :: cell_primary(:), y_array(:)
    PetscSection :: y_section
    DMLabel :: ghost_label

    err = 0
    np = eos%num_primary_variables
    primary_size = size(primary)
    if (primary_size /= np) then
       err = 1
       if (present(logfile)) then
          call logfile%write(LOG_LEVEL_ERR, 'initialize', 'initial.primary', &
               int_keys = ['size'], int_values = [primary_size], rank = 0)
       end if
    end if

    if (err == 0) then

       call global_vec_section(y, y_section)
       call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)

       call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr)
       CHKERRQ(ierr)
       call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
       CHKERRQ(ierr)
       end_interior_cell = dm_get_end_interior_cell(mesh%dm, end_cell)

       do c = start_cell, end_interior_cell - 1
          call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
          if (ghost < 0) then
             y_offset = global_section_offset(y_section, c, y_range_start)
             cell_primary => y_array(y_offset : y_offset + np - 1)
             cell_primary = primary
          end if
       end do

       call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)

    end if

  end subroutine setup_initial_primary_constant

!------------------------------------------------------------------------

  subroutine setup_initial_primary_array(json, mesh, eos, y, range_start, &
       minc_specified, logfile, err)

    !! Initializes solution vector y from a rank-2 array of values for
    !! all mesh cells in the JSON input.

    use fson_mpi_module
    use mesh_module, only: mesh_type
    use eos_module
    use dm_utils_module, only: global_vec_section, global_section_offset, &
         natural_to_local_cell_index, dm_distribute_local_vec, dm_cell_counts
    use mpi_utils_module, only: mpi_broadcast_error_flag
    use logfile_module

    type(fson_value), pointer, intent(in) :: json
    type(mesh_type), intent(in) :: mesh
    class(eos_type), intent(in) :: eos
    Vec, intent(in out) :: y
    PetscInt, intent(in) :: range_start
    PetscBool, intent(in) :: minc_specified
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscMPIInt :: nprocs, rank
    PetscInt :: num_data
    PetscReal, allocatable :: primary_array(:,:)
    PetscReal, pointer, contiguous :: primary_data(:)
    Vec :: primary
    PetscSection :: section
    IS :: interior
    PetscErrorCode :: ierr
    PetscReal, pointer, contiguous :: y_array(:)
    DMLabel :: ghost_label
    ISLocalToGlobalMapping :: l2g
    PetscInt :: i, c, offset, ghost, cells_total, cells_min, cells_max, np
    PetscInt, allocatable :: primary_shape(:)
    PetscReal, pointer, contiguous :: cell_primary(:)

    err = 0
    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call MPI_comm_size(PETSC_COMM_WORLD, nprocs, ierr)

    np = eos%num_primary_variables

    if ((.not. mesh%has_minc) .or. (.not. minc_specified)) then

       call dm_cell_counts(mesh%serial_dm, cells_total, cells_min, cells_max)
       if (rank == 0) then
          call fson_get(json, "initial.primary", primary_array)
          primary_shape = shape(primary_array)
          if (any(primary_shape /= [cells_total, np])) then
             err = 1
             if (present(logfile)) then
                call logfile%write(LOG_LEVEL_ERR, 'initialize', 'initial.primary', &
                     int_array_key = 'shape', int_array_value = primary_shape, rank = 0)
             end if
          end if
       else
          allocate(primary_array(0, 0))
       end if
       call mpi_broadcast_error_flag(err)
       if (err == 0) then
          call DMCreateLocalVector(mesh%serial_dm, primary, ierr); CHKERRQ(ierr)
          call VecGetArrayF90(primary, primary_data, ierr); CHKERRQ(ierr)
          primary_data = pack(transpose(primary_array), PETSC_TRUE)
          call VecRestoreArrayF90(primary, primary_data, ierr); CHKERRQ(ierr)
          deallocate(primary_array)

          if (nprocs > 1) then
             call dm_distribute_local_vec(mesh%dm, mesh%dist_sf, primary)
             call DMLocalToGlobal(mesh%dm, primary, INSERT_VALUES, y, &
                  ierr); CHKERRQ(ierr)
          else
             call VecGetLocalSize(primary, num_data, ierr); CHKERRQ(ierr)
             call ISCreateStride(PETSC_COMM_WORLD, num_data, &
                  0, 1, interior, ierr); CHKERRQ(ierr)
             call VecISCopy(y, interior, SCATTER_FORWARD, primary, &
                  ierr); CHKERRQ(ierr)
             call ISDestroy(interior, ierr); CHKERRQ(ierr)
          end if
          call VecDestroy(primary, ierr); CHKERRQ(ierr)
       end if
    else

       ! Initialising entire MINC solution from JSON array - note this is not scalable:
       call fson_get_mpi(json, "initial.primary", val = primary_array)
       num_data = size(primary_array, 1)
       call global_vec_section(y, section)
       call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)
       call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
       call DMGetLocalToGlobalMapping(mesh%dm, l2g, ierr); CHKERRQ(ierr)
       do i = 1, num_data
          associate(natural_cell_index => i - 1)
            c = natural_to_local_cell_index(mesh%cell_natural_global, l2g, &
                 natural_cell_index)
            if (c >= 0) then
               call DMLabelGetValue(ghost_label, c, ghost, ierr)
               if (ghost < 0) then
                  offset = global_section_offset(section, c, range_start)
                  cell_primary => y_array(offset : offset + &
                       eos%num_primary_variables - 1)
                  cell_primary = primary_array(i, :)
               end if
            end if
          end associate
       end do
       call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)

    end if

  end subroutine setup_initial_primary_array

!------------------------------------------------------------------------

  subroutine setup_initial_region_constant(mesh, region, eos, &
       fluid_vector, fluid_range_start)

    !! Initializes fluid regions with constant values over the mesh.

    use dm_utils_module, only: global_vec_section, global_section_offset, &
         dm_get_end_interior_cell
    use mesh_module, only: mesh_type
    use fluid_module, only: fluid_type
    use eos_module, only: eos_type

    type(mesh_type), intent(in) :: mesh
    PetscInt, intent(in) :: region
    class(eos_type), intent(in) :: eos
    Vec, intent(in out) :: fluid_vector
    PetscInt, intent(in) :: fluid_range_start
    ! Locals:
    PetscInt :: c, ghost, start_cell, end_cell, end_interior_cell
    PetscInt :: fluid_offset
    PetscErrorCode :: ierr
    type(fluid_type) :: fluid
    PetscReal, pointer, contiguous :: fluid_array(:)
    PetscSection :: fluid_section
    DMLabel :: ghost_label

    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call fluid%init(eos%num_components, eos%num_phases)

    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)
    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    end_interior_cell = dm_get_end_interior_cell(mesh%dm, end_cell)

    do c = start_cell, end_interior_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          fluid_offset = global_section_offset(fluid_section, c, fluid_range_start)
          call fluid%assign(fluid_array, fluid_offset)
          fluid%region = dble(region)
       end if
    end do

    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call fluid%destroy()

  end subroutine setup_initial_region_constant

!------------------------------------------------------------------------

  subroutine setup_initial_region_array(json, mesh, eos, fluid_vector, &
       fluid_range_start, minc_specified, logfile, err)

    !! Initializes fluid regions from a rank-1 array of values for
    !! all mesh cells in the JSON input.

    use fson_mpi_module
    use mesh_module, only: mesh_type
    use eos_module
    use dm_utils_module
    use fluid_module, only: fluid_type
    use minc_module, only: minc_level_label_name
    use logfile_module
    use mpi_utils_module, only: mpi_broadcast_error_flag

    type(fson_value), pointer, intent(in) :: json
    type(mesh_type), intent(in) :: mesh
    class(eos_type), intent(in) :: eos
    Vec, intent(in out) :: fluid_vector
    PetscInt, intent(in) :: fluid_range_start
    PetscBool, intent(in) :: minc_specified
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscMPIInt :: nprocs, rank
    PetscInt :: start_cell, end_cell, num_cells, end_interior_cell
    PetscInt, allocatable :: region_array(:)
    IS :: serial_region
    PetscSection :: serial_section, section, fluid_section
    DM :: dm_is
    IS :: region
    PetscInt :: i, c, ghost, offset, minc_level
    PetscInt :: cells_total, cells_min, cells_max, region_size
    DMLabel :: ghost_label, minc_label
    type(fluid_type) :: fluid
    PetscInt, pointer, contiguous :: region_indices(:)
    PetscReal, pointer, contiguous :: fluid_array(:)
    ISLocalToGlobalMapping :: l2g
    PetscErrorCode :: ierr

    err = 0
    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call MPI_comm_size(PETSC_COMM_WORLD, nprocs, ierr)
    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    if (mesh%has_minc) then
       call DMGetLabel(mesh%dm, minc_level_label_name, minc_label, &
            ierr); CHKERRQ(ierr)
    end if

    call fluid%init(eos%num_components, eos%num_phases)

    if ((.not. mesh%has_minc) .or. (.not. minc_specified)) then

       call dm_cell_counts(mesh%serial_dm, cells_total, cells_min, cells_max)
       if (rank == 0) then
          call fson_get(json, "initial.region", region_array)
          region_size = size(region_array)
          if (region_size /= cells_total) then
             err = 1
             if (present(logfile)) then
                call logfile%write(LOG_LEVEL_ERR, 'initialize', 'initial.region', &
                     int_keys = ['size'], int_values = [region_size], rank = 0)
             end if
          end if
       else
          allocate(region_array(0))
       end if

       call mpi_broadcast_error_flag(err)
       if (err == 0) then

          num_cells = size(region_array)
          call ISCreateGeneral(PETSC_COMM_WORLD, num_cells, region_array, &
               PETSC_COPY_VALUES, serial_region, ierr); CHKERRQ(ierr)

          if (nprocs > 1) then
             call DMClone(mesh%serial_dm, dm_is, ierr); CHKERRQ(ierr)
             call dm_set_default_data_layout(dm_is, 1)
             call DMGetSection(dm_is, serial_section, ierr); CHKERRQ(ierr)
             call PetscSectionCreate(PETSC_COMM_WORLD, section, ierr); CHKERRQ(ierr)
             call DMPlexDistributeFieldIS(mesh%dm, mesh%dist_sf, serial_section, &
                  serial_region, section, region, ierr); CHKERRQ(ierr)
             call DMDestroy(dm_is, ierr); CHKERRQ(ierr)
             call PetscSectionDestroy(section, ierr); CHKERRQ(ierr)
          else
             call ISDuplicate(serial_region, region, ierr); CHKERRQ(ierr)
          end if

          call global_vec_section(fluid_vector, fluid_section)
          call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
          call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
          CHKERRQ(ierr)
          end_interior_cell = dm_get_end_interior_cell(mesh%dm, end_cell)
          call ISGetIndicesF90(region, region_indices, ierr); CHKERRQ(ierr)
          i = 1
          do c = start_cell, end_interior_cell - 1
             call DMLabelGetValue(ghost_label, c, ghost, ierr)
             if (ghost < 0) then
                if (mesh%has_minc) then
                   call DMLabelGetValue(minc_label, c, minc_level, ierr)
                   CHKERRQ(ierr)
                else
                   minc_level = -1
                end if
                if (minc_level <= 0) then
                   offset =  global_section_offset(fluid_section, c, fluid_range_start)
                   call fluid%assign(fluid_array, offset)
                   fluid%region = dble(region_indices(i))
                end if
             end if
             i = i + 1
          end do

       end if

       deallocate(region_array)
       call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
       call ISRestoreIndicesF90(region, region_indices, ierr); CHKERRQ(ierr)
       call ISDestroy(region, ierr); CHKERRQ(ierr)
       call ISDestroy(serial_region, ierr); CHKERRQ(ierr)

    else

       ! Initialising entire MINC solution from JSON array - note this is not scalable:
       call fson_get_mpi(json, "initial.region", val = region_array)
       num_cells = size(region_array, 1)
       call global_vec_section(fluid_vector, fluid_section)
       call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
       call DMGetLocalToGlobalMapping(mesh%dm, l2g, ierr); CHKERRQ(ierr)

       do i = 1, num_cells
          associate(natural_cell_index => i - 1)
          c = natural_to_local_cell_index(mesh%cell_natural_global, l2g, &
               natural_cell_index)
          if (c >= 0) then
             call DMLabelGetValue(ghost_label, c, ghost, ierr)
             if (ghost < 0) then
                offset = global_section_offset(fluid_section, c, fluid_range_start)
                call fluid%assign(fluid_array, offset)
                fluid%region = dble(region_array(i))
             end if
          end if
        end associate
     end do

     call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
  end if

  call fluid%destroy()

  end subroutine setup_initial_region_array

!------------------------------------------------------------------------

  subroutine setup_initial_file(filename, mesh, eos, t, y, fluid_vector, &
       tracer_vector, y_range_start, fluid_range_start, tracer_range_start, &
       index, use_original_dm, tracers)
    !! Initializes fluid vector, solution vector y and tracer vector
    !! (if tracers are being simulated) from HDF5 file.

    use mesh_module
    use dm_utils_module, only: global_vec_section, global_section_offset, &
         dm_get_cell_index, vec_reorder, section_get_field_names, &
         vec_copy_subvector
    use eos_module, only: eos_type, max_component_name_length, &
         max_phase_name_length
    use fluid_module, only: fluid_type, create_fluid_vector
    use utils_module, only: str_array_index, str_to_lower
    use hdf5io_module, only: max_field_name_length, vec_load_fields_hdf5
    use tracer_module, only: tracer_type

    character(len = *), intent(in) :: filename
    type(mesh_type), intent(in) :: mesh
    class(eos_type), intent(in) :: eos
    PetscReal, intent(in) :: t
    Vec, intent(in out) :: y, fluid_vector, tracer_vector
    PetscInt, intent(in) :: y_range_start, fluid_range_start, tracer_range_start
    PetscInt, intent(in) :: index !! time index to fetch initial conditions from
    PetscBool, intent(in) :: use_original_dm !! Whether file results correspond to original_dm
    type(tracer_type), intent(in) :: tracers(:)
    ! Locals:
    PetscViewer :: viewer
    PetscInt, allocatable :: field_indices(:)
    PetscInt :: tracer_field_indices(size(tracers))
    PetscInt :: i, num_tracers
    IS :: original_cell_index, output_cell_index
    PetscErrorCode :: ierr

    call PetscViewerHDF5Open(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &
         viewer, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5PushGroup(viewer, "/", ierr); CHKERRQ(ierr)

    call get_required_field_indices()
    num_tracers = size(tracers)
    if (num_tracers > 0) then
       tracer_field_indices = [(i - 1, i = 1, num_tracers)]
    end if

    if (use_original_dm) then
       call dm_get_cell_index(mesh%original_dm, mesh%original_cell_natural_global, &
            original_cell_index)
       call ISDuplicate(original_cell_index, output_cell_index, ierr); CHKERRQ(ierr)
    else
       call ISDuplicate(mesh%cell_index, output_cell_index, ierr); CHKERRQ(ierr)
    end if
    call PetscObjectSetName(output_cell_index, "cell_index", ierr); CHKERRQ(ierr)
    call ISLoad(output_cell_index, viewer, ierr); CHKERRQ(ierr)

    if (use_original_dm) then
       call load_fluid_original_dm()
       if (num_tracers > 0) call load_tracers_original_dm()
       call ISDestroy(original_cell_index, ierr); CHKERRQ(ierr)
    else
       call load_fluid()
       if (num_tracers > 0) call load_tracers()
    end if

    call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
    call ISDestroy(output_cell_index, ierr); CHKERRQ(ierr)

    deallocate(field_indices)

    call assign_solution_vector()

  contains

!........................................................................

    subroutine get_required_field_indices()
      !! Gets field indices to be read in, corresponding to the
      !! required output fluid fields of the EOS.

      ! Locals:
      DM :: fluid_dm
      PetscSection :: section
      character(max_field_name_length), allocatable :: fields(:)
      PetscInt :: i
      PetscErrorCode :: ierr

      associate(num_required => size(eos%required_output_fluid_fields))
        allocate(field_indices(num_required))
        call VecGetDM(fluid_vector, fluid_dm, ierr); CHKERRQ(ierr)
        call DMGetSection(fluid_dm, section, ierr); CHKERRQ(ierr)
        call section_get_field_names(section, PETSC_TRUE, fields)
        do i = 1, num_required
           field_indices(i) = str_array_index( &
                str_to_lower(eos%required_output_fluid_fields(i)), &
                fields) - 1
        end do
      end associate
      deallocate(fields)

    end subroutine get_required_field_indices

!........................................................................

    subroutine load_fluid_original_dm()
      !! Loads fluid vector corresponding to original DM from HDF5 file.

      ! Locals:
      DM :: fluid_dm
      Vec :: original_fluid_vector
      PetscInt :: original_fluid_range_start
      PetscErrorCode :: ierr

      call create_fluid_vector(mesh%original_dm, max_component_name_length, &
           eos%component_names, max_phase_name_length, &
           eos%phase_names, original_fluid_vector, original_fluid_range_start)

      call VecGetDM(original_fluid_vector, fluid_dm, ierr); CHKERRQ(ierr)
      call DMSetOutputSequenceNumber(fluid_dm, index, t, ierr); CHKERRQ(ierr)
      call vec_load_fields_hdf5(original_fluid_vector, field_indices, &
           "/cell_fields", viewer, original_cell_index, output_cell_index)
      call vec_copy_subvector(original_fluid_vector, fluid_vector)
      call VecDestroy(original_fluid_vector, ierr); CHKERRQ(ierr)

    end subroutine load_fluid_original_dm

!........................................................................

    subroutine load_fluid()
      !! Loads fluid vector from HDF5 file.

      ! Locals:
      DM :: fluid_dm
      PetscErrorCode :: ierr

      call VecGetDM(fluid_vector, fluid_dm, ierr); CHKERRQ(ierr)
      call DMSetOutputSequenceNumber(fluid_dm, index, t, ierr); CHKERRQ(ierr)
      call vec_load_fields_hdf5(fluid_vector, field_indices, &
           "/cell_fields", viewer, mesh%cell_index, output_cell_index)

    end subroutine load_fluid

!........................................................................

    subroutine load_tracers_original_dm()
      !! Loads tracer vector corresponding to original DM from HDF5 file.

      use tracer_module, only: create_tracer_vector

      ! Locals:
      DM :: tracer_dm
      Vec :: original_tracer_vector
      PetscInt :: original_tracer_range_start
      PetscErrorCode :: ierr

      call create_tracer_vector(mesh%original_dm, tracers, &
           original_tracer_vector, original_tracer_range_start)

      call VecGetDM(original_tracer_vector, tracer_dm, ierr); CHKERRQ(ierr)
      call DMSetOutputSequenceNumber(tracer_dm, index, t, ierr); CHKERRQ(ierr)
      call vec_load_fields_hdf5(original_tracer_vector, tracer_field_indices, &
           "/cell_fields", viewer, original_cell_index, output_cell_index)
      call vec_copy_subvector(original_tracer_vector, tracer_vector)
      call VecDestroy(original_tracer_vector, ierr); CHKERRQ(ierr)

    end subroutine load_tracers_original_dm

!........................................................................

    subroutine load_tracers()
      !! Loads tracer vector from HDF5 file.

      ! Locals:
      DM :: tracer_dm
      PetscErrorCode :: ierr

      call VecGetDM(tracer_vector, tracer_dm, ierr); CHKERRQ(ierr)
      call DMSetOutputSequenceNumber(tracer_dm, index, t, ierr); CHKERRQ(ierr)
      call vec_load_fields_hdf5(tracer_vector, tracer_field_indices, &
           "/cell_fields", viewer, mesh%cell_index, output_cell_index)

    end subroutine load_tracers

!........................................................................

    subroutine assign_solution_vector()
      !! Determine solution vector from fluid vector.

      ! Locals:
      PetscSection :: y_section, fluid_section
      PetscInt :: start_cell, end_cell, c
      PetscInt :: ghost, np
      PetscInt :: y_offset, fluid_offset
      PetscReal, pointer, contiguous :: y_array(:), fluid_array(:)
      PetscReal, pointer, contiguous :: cell_primary(:)
      DMLabel :: ghost_label
      type(fluid_type) :: fluid
      PetscErrorCode :: ierr
      
      call fluid%init(eos%num_components, eos%num_phases)
      np = eos%num_primary_variables
      call global_vec_section(y, y_section)
      call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)

      call global_vec_section(fluid_vector, fluid_section)
      call VecGetArrayReadF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

      call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr)
      CHKERRQ(ierr)
      call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
      CHKERRQ(ierr)

      do c = start_cell, end_cell - 1
         call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
         if (ghost < 0) then
            fluid_offset =  global_section_offset(fluid_section, c, &
                 fluid_range_start)
            if (fluid_array(fluid_offset) > 0._dp) then
               y_offset = global_section_offset(y_section, c, y_range_start)
               cell_primary => y_array(y_offset : y_offset + np - 1)
               call fluid%assign(fluid_array, fluid_offset)
               call eos%primary_variables(fluid, cell_primary)
            end if
         end if
      end do

      call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)
      call VecRestoreArrayReadF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
      call fluid%destroy()

    end subroutine assign_solution_vector

  end subroutine setup_initial_file

!------------------------------------------------------------------------

  subroutine scale_initial_primary(mesh, eos, y, fluid_vector, &
      y_range_start, fluid_range_start)
    !! Scales primary variables y in each cell according to the eos
    !! primary variable scaling.

    use mesh_module, only: mesh_type
    use eos_module, only: eos_type
    use fluid_module, only: fluid_type
    use dm_utils_module, only: global_vec_section, global_section_offset, &
         dm_get_end_interior_cell

    type(mesh_type), intent(in) :: mesh
    class(eos_type), intent(in) :: eos
    Vec, intent(in out) :: y, fluid_vector
    PetscInt, intent(in) :: y_range_start, fluid_range_start
    ! Locals:
    PetscSection :: y_section, fluid_section
    PetscReal, pointer, contiguous :: y_array(:), fluid_array(:)
    PetscReal, pointer, contiguous :: cell_primary(:)
    type(fluid_type) :: fluid
    DMLabel :: ghost_label
    PetscInt :: np, c, ghost, y_offset, fluid_offset
    PetscInt :: start_cell, end_cell, end_interior_cell
    PetscErrorCode :: ierr

    call global_vec_section(y, y_section)
    call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)

    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayReadF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call fluid%init(eos%num_components, eos%num_phases)
    np = eos%num_primary_variables

    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)
    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    end_interior_cell = dm_get_end_interior_cell(mesh%dm, end_cell)

    do c = start_cell, end_interior_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          y_offset = global_section_offset(y_section, c, y_range_start)
          cell_primary => y_array(y_offset : y_offset + np - 1)
          fluid_offset = global_section_offset(fluid_section, c, &
               fluid_range_start)
          call fluid%assign(fluid_array, fluid_offset)
          cell_primary = eos%scale(cell_primary, nint(fluid%region))
       end if
    end do

    call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call fluid%destroy()

  end subroutine scale_initial_primary

!------------------------------------------------------------------------

  subroutine setup_initial(json, mesh, eos, t, y, fluid_vector, &
       tracer_vector, y_range_start, fluid_range_start, tracer_range_start, &
       tracers, logfile, err)
    !! Initializes time t and a Vec y with initial conditions read
    !! from JSON input 'initial'.  A full set of initial conditions
    !! may be read in from an HDF5 file, or a constant set of primary
    !! variables can be read from the JSON input and applied to all
    !! cells. The fluid thermodynamic region for each cell is also
    !! initialised (in the fluid vector), and the tracer mass
    !! fractions (if tracers are being simulated).

    use mesh_module
    use eos_module
    use dm_utils_module, only: global_vec_section, global_section_offset
    use tracer_module, only: tracer_type
    use logfile_module

    type(fson_value), pointer, intent(in) :: json
    type(mesh_type), intent(in) :: mesh
    class(eos_type), intent(in) :: eos
    PetscReal, intent(out) :: t
    Vec, intent(in out) :: y, fluid_vector, tracer_vector
    PetscInt, intent(in) :: y_range_start, fluid_range_start, tracer_range_start
    type(tracer_type), intent(in) :: tracers(:)
    type(logfile_type), intent(in out), optional :: logfile
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscReal, parameter :: default_start_time = 0.0_dp
    PetscReal, allocatable :: primary(:)
    PetscInt :: region, num_tracers
    PetscReal :: primary_scalar, tracer_scalar
    PetscInt :: primary_rank, region_rank, tracer_rank
    PetscReal, allocatable :: tracer_values(:)
    PetscInt, parameter :: max_filename_length = 240
    character(len = max_filename_length) :: filename
    PetscInt, parameter :: default_index = 0
    PetscInt :: index, tracer_array_size
    PetscBool :: minc_specified, use_original_dm, has_file
    PetscBool, parameter :: default_minc_specified = PETSC_FALSE
    PetscReal, parameter :: default_tracer = 0._dp

    err = 0
    num_tracers = size(tracers)
    call fson_get_mpi(json, "time.start", default_start_time, t, logfile)

    if (fson_has_mpi(json, "initial")) then

       if (mesh%has_minc) then
          call fson_get_mpi(json, "initial.minc", default_minc_specified, &
               minc_specified, logfile)
          use_original_dm = .not. minc_specified
       else
          minc_specified = PETSC_FALSE
          use_original_dm = PETSC_FALSE
       end if

       has_file = fson_has_mpi(json, "initial.filename")

       primary_rank = fson_mpi_array_rank(json, "initial.primary")
       select case (primary_rank)
       case (-1)
          call setup_initial_primary_constant(mesh, eos%default_primary, &
               eos, y, y_range_start, logfile, err)
          if (.not. has_file) then
             if (present(logfile)) then
                call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                     real_array_key = 'initial.primary', &
                     real_array_value = eos%default_primary)
             end if
          end if
       case (0, 1)
          if (primary_rank == 0) then
             call fson_get_mpi(json, "initial.primary", val = primary_scalar)
             primary = [primary_scalar]
          else
             call fson_get_mpi(json, "initial.primary", val = primary)
          end if
          call setup_initial_primary_constant(mesh, primary, eos, y, &
               y_range_start, logfile, err)
          deallocate(primary)
       case (2)
          call setup_initial_primary_array(json, mesh, eos, y, y_range_start, &
               minc_specified, logfile, err)
       case default
          if (present(logfile)) then
             call logfile%write(LOG_LEVEL_WARN, 'input', &
                  '"unrecognised initial.primary"')
          end if
       end select

       if (err == 0) then

          region_rank = fson_mpi_array_rank(json, "initial.region")
          select case (region_rank)
          case (-1)
             call setup_initial_region_constant(mesh, eos%default_region, &
                  eos, fluid_vector, fluid_range_start)
             if (.not. has_file) then
                if (present(logfile)) then
                   call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                        int_keys = ['initial.region'], &
                        int_values = [eos%default_region])
                end if
             end if
          case (0)
             call fson_get_mpi(json, "initial.region", val = region)
             call setup_initial_region_constant(mesh, region, eos, fluid_vector, &
                  fluid_range_start)
          case (1)
             call setup_initial_region_array(json, mesh, eos, fluid_vector, &
                  fluid_range_start, minc_specified, logfile, err)
          case default
             if (present(logfile)) then
                call logfile%write(LOG_LEVEL_WARN, 'input', &
                     '"unrecognised initial.region"')
             end if
          end select

          if (err == 0) then

             if (num_tracers > 0) then
                tracer_rank = fson_mpi_array_rank(json, "initial.tracer")
                select case (tracer_rank)
                case (-1)
                   call setup_initial_tracer_scalar(mesh, default_tracer, &
                        tracer_vector, num_tracers, tracer_range_start)
                   if (.not. has_file) then
                      if (present(logfile)) then
                         call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                              real_keys = ['initial.tracer'], &
                              real_values = [default_tracer])
                      end if
                   end if
                case (0)
                   call fson_get_mpi(json, "initial.tracer", val = tracer_scalar)
                   call setup_initial_tracer_scalar(mesh, tracer_scalar, &
                        tracer_vector, num_tracers, tracer_range_start)
                case (1)
                   tracer_array_size = fson_value_count_mpi(json, "initial.tracer")
                   call fson_get_mpi(json, "initial.tracer", val = tracer_values)
                   call setup_initial_tracer_constant(mesh, tracer_values, &
                        tracer_vector, num_tracers, tracer_range_start)
                   ! could add rank 2 option here for initialising tracer
                   ! from array of values over entire mesh
                case default
                   if (present(logfile)) then
                      call logfile%write(LOG_LEVEL_WARN, 'input', &
                           '"unrecognised initial.tracer"')
                   end if
                end select
             end if

             if (has_file) then
                call fson_get_mpi(json, "initial.filename", val = filename)
                call fson_get_mpi(json, "initial.index", default_index, index, logfile)
                call setup_initial_file(filename, mesh, eos, t, y, fluid_vector, &
                     tracer_vector, y_range_start, fluid_range_start, tracer_range_start, &
                     index, use_original_dm, tracers)
             end if

          end if
       end if

    else

       call setup_initial_primary_constant(mesh, eos%default_primary, &
            eos, y, y_range_start, logfile, err)
       call setup_initial_region_constant(mesh, eos%default_region, &
            eos, fluid_vector, fluid_range_start)
       if (present(logfile)) then
          call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
               real_array_key = 'initial.primary', &
               real_array_value = eos%default_primary)
          call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
               int_keys = ['initial.region'], &
               int_values = [eos%default_region])
       end if
       if (num_tracers > 0) then
          call setup_initial_tracer_scalar(mesh, default_tracer, &
               tracer_vector, num_tracers, tracer_range_start)
       end if
       minc_specified = PETSC_TRUE

    end if

    if (mesh%has_minc .and. (.not. minc_specified)) then
       call setup_minc_initial(mesh, eos, y, fluid_vector, y_range_start, &
            fluid_range_start)
       if (num_tracers > 0) then
          call setup_minc_initial_tracer(mesh, tracer_vector, &
               tracer_range_start, num_tracers)
       end if
    end if

  end subroutine setup_initial

!------------------------------------------------------------------------

  subroutine setup_minc_initial(mesh, eos, y, fluid_vector, &
       y_range_start, fluid_range_start)

    !! Sets up initial conditions in MINC matrix cells. These are
    !! simply copied from the corresponding fracture cells, which are
    !! assumed already initialised.

    use mesh_module
    use eos_module
    use dm_utils_module, only: global_vec_section, global_section_offset
    use fluid_module, only: fluid_type
    use minc_module, only: minc_zone_label_name

    class(mesh_type), intent(in) :: mesh !! Mesh object
    class(eos_type), intent(in) :: eos !! Equation of state module
    Vec, intent(in out) :: y, fluid_vector !! Solution and fluid vectors
    PetscInt, intent(in) :: y_range_start, fluid_range_start !! Range starts for vectors
    ! Locals:
    PetscSection :: y_section, fluid_section
    PetscReal, pointer, contiguous :: y_array(:), fluid_array(:)
    PetscInt :: iminc, c, h, i, m, minc_p, num_minc_zone_cells, max_num_levels
    PetscInt, allocatable :: ic(:)
    PetscInt :: y_offset, fluid_offset, y_minc_offset, fluid_minc_offset, np
    type(fluid_type) :: fluid
    PetscInt :: frac_region
    IS :: minc_IS
    PetscInt, pointer :: minc_cells(:)
    PetscErrorCode :: ierr

    call global_vec_section(y, y_section)
    call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call fluid%init(eos%num_components, eos%num_phases)
    np = eos%num_primary_variables

    max_num_levels = maxval(mesh%minc%num_levels)
    allocate(ic(max_num_levels))
    ic = 0
    h = 0

    do iminc = 1, size(mesh%minc)
       associate(minc => mesh%minc(iminc))
         call DMGetStratumSize(mesh%dm, minc_zone_label_name, iminc, &
              num_minc_zone_cells, ierr); CHKERRQ(ierr)
         if (num_minc_zone_cells > 0) then
            call DMGetStratumIS(mesh%dm, minc_zone_label_name, &
                 iminc, minc_IS, ierr); CHKERRQ(ierr)
            call ISGetIndicesF90(minc_IS, minc_cells, ierr); CHKERRQ(ierr)

            do i = 1, num_minc_zone_cells

               c = minc_cells(i)
               if (mesh%ghost_cell(c) < 0) then

                  minc_p = mesh%strata(h)%minc_point(c, 0)
                  y_offset = global_section_offset(y_section, minc_p, y_range_start)
                  associate(frac_primary => y_array(y_offset : y_offset + np - 1))
                    fluid_offset = global_section_offset(fluid_section, minc_p, fluid_range_start)
                    call fluid%assign(fluid_array, fluid_offset)
                    frac_region = nint(fluid%region)

                    do m = 1, minc%num_levels
                      minc_p = mesh%strata(h)%minc_point(ic(m), m)
                       y_minc_offset = global_section_offset(y_section, minc_p, y_range_start)
                       associate (primary => y_array(y_minc_offset : &
                            y_minc_offset + np - 1))
                         primary = frac_primary
                       end associate
                       fluid_minc_offset = global_section_offset(fluid_section, minc_p, &
                            fluid_range_start)
                       call fluid%assign(fluid_array, fluid_minc_offset)
                       fluid%region = dble(frac_region)
                       ic(m) = ic(m) + 1
                    end do

                  end associate
               end if
            end do

            call ISRestoreIndicesF90(minc_IS, minc_cells, ierr); CHKERRQ(ierr)
         end if
       end associate
    end do

    call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call fluid%destroy()
    deallocate(ic)

  end subroutine setup_minc_initial

!------------------------------------------------------------------------

  subroutine setup_initial_tracer_constant(mesh, tracer_values, &
       tracer_vector, num_tracers, tracer_range_start)
    !! Initialise tracer vector with constant values (one per tracer)
    !! over the mesh.

    use dm_utils_module, only: global_vec_section, global_section_offset, &
         dm_get_end_interior_cell
    use mesh_module, only: mesh_type

    type(mesh_type), intent(in) :: mesh
    PetscReal, intent(in) :: tracer_values(:)
    Vec, intent(in out) :: tracer_vector
    PetscInt, intent(in) :: num_tracers
    PetscInt, intent(in) :: tracer_range_start
    ! Locals:
    PetscSection :: section
    PetscReal, pointer, contiguous :: tracer_array(:), cell_tracer(:)
    DMLabel :: ghost_label
    PetscInt :: c, start_cell, end_cell, end_interior_cell, ghost, offset
    PetscErrorCode :: ierr

    call global_vec_section(tracer_vector, section)
    call VecGetArrayF90(tracer_vector, tracer_array, ierr); CHKERRQ(ierr)

    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)
    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    end_interior_cell = dm_get_end_interior_cell(mesh%dm, end_cell)

    do c = start_cell, end_interior_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          offset = global_section_offset(section, c, tracer_range_start)
          cell_tracer => tracer_array(offset : offset + num_tracers - 1)
          cell_tracer = tracer_values
       end if
    end do

    call VecRestoreArrayF90(tracer_vector, tracer_array, ierr); CHKERRQ(ierr)

  end subroutine setup_initial_tracer_constant

!------------------------------------------------------------------------

  subroutine setup_initial_tracer_scalar(mesh, tracer_scalar, &
       tracer_vector, num_tracers, tracer_range_start)
    !! Initialise tracer vector with the same scalar value for each tracer.

    use mesh_module, only: mesh_type

    type(mesh_type), intent(in) :: mesh
    PetscReal, intent(in) :: tracer_scalar
    Vec, intent(in out) :: tracer_vector
    PetscInt, intent(in) :: num_tracers
    PetscInt, intent(in) :: tracer_range_start
    ! Locals:
    PetscReal :: tracer_values(num_tracers)

    tracer_values = tracer_scalar
    call setup_initial_tracer_constant(mesh, tracer_values, &
         tracer_vector, num_tracers, tracer_range_start)

  end subroutine setup_initial_tracer_scalar

!------------------------------------------------------------------------

  subroutine setup_minc_initial_tracer(mesh, tracer_vector, &
       tracer_range_start, num_tracers)

    !! Sets up tracer initial conditions in MINC matrix cells. These
    !! are simply copied from the corresponding fracture cells, which
    !! are assumed already initialised.

    use mesh_module
    use dm_utils_module, only: global_vec_section, global_section_offset
    use minc_module, only: minc_zone_label_name

    class(mesh_type), intent(in) :: mesh !! Mesh object
    Vec, intent(in out) :: tracer_vector !! Tracer vector
    PetscInt, intent(in) :: tracer_range_start !! Range start for tracer vector
    PetscInt, intent(in) :: num_tracers !! Number of tracers (> 0)
    ! Locals:
    PetscSection :: tracer_section
    PetscReal, pointer, contiguous :: tracer_array(:)
    PetscInt :: iminc, c, h, i, m, minc_p, num_minc_zone_cells, max_num_levels
    PetscInt, allocatable :: ic(:)
    PetscInt :: tracer_offset, tracer_minc_offset
    IS :: minc_IS
    PetscInt, pointer :: minc_cells(:)
    PetscErrorCode :: ierr

    call global_vec_section(tracer_vector, tracer_section)
    call VecGetArrayF90(tracer_vector, tracer_array, ierr); CHKERRQ(ierr)

    max_num_levels = maxval(mesh%minc%num_levels)
    allocate(ic(max_num_levels))
    ic = 0
    h = 0

    do iminc = 1, size(mesh%minc)
       associate(minc => mesh%minc(iminc))
         call DMGetStratumSize(mesh%dm, minc_zone_label_name, iminc, &
              num_minc_zone_cells, ierr); CHKERRQ(ierr)
         if (num_minc_zone_cells > 0) then
            call DMGetStratumIS(mesh%dm, minc_zone_label_name, &
                 iminc, minc_IS, ierr); CHKERRQ(ierr)
            call ISGetIndicesF90(minc_IS, minc_cells, ierr); CHKERRQ(ierr)

            do i = 1, num_minc_zone_cells

               c = minc_cells(i)
               if (mesh%ghost_cell(c) < 0) then

                  minc_p = mesh%strata(h)%minc_point(c, 0)
                  tracer_offset = global_section_offset(tracer_section, &
                       minc_p, tracer_range_start)
                  associate(frac_tracer => tracer_array(tracer_offset : &
                       tracer_offset + num_tracers - 1))

                    do m = 1, minc%num_levels
                      minc_p = mesh%strata(h)%minc_point(ic(m), m)
                      tracer_minc_offset = global_section_offset(tracer_section, &
                           minc_p, tracer_range_start)
                       associate (tracer => tracer_array(tracer_minc_offset : &
                            tracer_minc_offset + num_tracers - 1))
                         tracer = frac_tracer
                       end associate
                       ic(m) = ic(m) + 1
                    end do

                  end associate
               end if
            end do

            call ISRestoreIndicesF90(minc_IS, minc_cells, ierr); CHKERRQ(ierr)
         end if
       end associate
    end do

    call VecRestoreArrayF90(tracer_vector, tracer_array, ierr); CHKERRQ(ierr)
    deallocate(ic)

  end subroutine setup_minc_initial_tracer

!------------------------------------------------------------------------

end module initial_module
