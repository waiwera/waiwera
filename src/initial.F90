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

  public :: setup_initial, scale_initial_primary, setup_minc_initial

contains

!------------------------------------------------------------------------

  subroutine setup_initial_primary_constant(mesh, primary, &
       eos, y, y_range_start)

    !! Initializes solution vector y with constant values over the mesh.

    use dm_utils_module, only: global_vec_section, global_section_offset
    use mesh_module, only: mesh_type
    use eos_module, only: eos_type

    type(mesh_type), intent(in) :: mesh
    PetscReal, intent(in) :: primary(:)
    class(eos_type), intent(in) :: eos
    Vec, intent(in out) :: y
    PetscInt, intent(in) :: y_range_start
    ! Locals:
    PetscInt :: np, c, ghost
    PetscInt :: y_offset
    PetscErrorCode :: ierr
    PetscReal, pointer, contiguous :: cell_primary(:), y_array(:)
    PetscSection :: y_section
    DMLabel :: ghost_label

    np = eos%num_primary_variables

    call global_vec_section(y, y_section)
    call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)

    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    do c = mesh%start_cell, mesh%end_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          call global_section_offset(y_section, c, &
               y_range_start, y_offset, ierr); CHKERRQ(ierr)
          cell_primary => y_array(y_offset : y_offset + np - 1)
          cell_primary = primary
       end if
    end do

    call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)

  end subroutine setup_initial_primary_constant

!------------------------------------------------------------------------

  subroutine setup_initial_primary_array(mesh, primary, eos, y, &
       range_start)

    !! Initializes solution vector y from a rank-2 array of values for
    !! all mesh cells.

    use dm_utils_module, only: global_vec_section, global_section_offset
    use mesh_module, only: mesh_type
    use cell_order_module, only: cell_order_label_name
    use eos_module, only: eos_type

    type(mesh_type), intent(in) :: mesh
    PetscReal, intent(in) :: primary(:,:)
    class(eos_type), intent(in) :: eos
    Vec, intent(in out) :: y
    PetscInt, intent(in) :: range_start
    ! Locals:
    PetscInt :: num_cells, np, global_cell_index, i
    PetscInt :: offset, ghost, num_matching, c
    PetscErrorCode :: ierr
    PetscReal, pointer, contiguous :: cell_primary(:), y_array(:)
    PetscSection :: section
    IS :: cell_IS
    PetscInt, pointer :: cells(:)
    DMLabel :: ghost_label

    num_cells = size(primary, 1)
    np = eos%num_primary_variables

    call global_vec_section(y, section)
    call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)

    do i = 1, num_cells
       global_cell_index = i - 1
       call DMGetStratumSize(mesh%dm, cell_order_label_name, &
            global_cell_index, num_matching, ierr); CHKERRQ(ierr)
       if (num_matching > 0) then
          call DMGetStratumIS(mesh%dm, cell_order_label_name, &
               global_cell_index, cell_IS, ierr); CHKERRQ(ierr)
          call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
          do c = 1, num_matching
             call DMLabelGetValue(ghost_label, cells(c), ghost, ierr)
             if (ghost < 0) then
                call global_section_offset(section, cells(c), range_start, &
                     offset, ierr)
                cell_primary => y_array(offset : offset + np - 1)
                cell_primary = primary(i, :)
             end if
          end do
          call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
          call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)
       end if
    end do

    call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)

  end subroutine setup_initial_primary_array

!------------------------------------------------------------------------

  subroutine setup_initial_region_constant(mesh, region, eos, &
       fluid_vector, fluid_range_start)

    !! Initializes fluid regions with constant values over the mesh.

    use dm_utils_module, only: global_vec_section, global_section_offset
    use mesh_module, only: mesh_type
    use fluid_module, only: fluid_type
    use eos_module, only: eos_type

    type(mesh_type), intent(in) :: mesh
    PetscInt, intent(in) :: region
    class(eos_type), intent(in) :: eos
    Vec, intent(in out) :: fluid_vector
    PetscInt, intent(in) :: fluid_range_start
    ! Locals:
    PetscInt :: c, ghost
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

    do c = mesh%start_cell, mesh%end_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          call global_section_offset(fluid_section, c, &
               fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)
          call fluid%assign(fluid_array, fluid_offset)
          fluid%region = dble(region)
       end if
    end do

    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call fluid%destroy()

  end subroutine setup_initial_region_constant

!------------------------------------------------------------------------

  subroutine setup_initial_region_array(mesh, region, eos, fluid_vector, &
       fluid_range_start)

    !! Initializes fluid regions from a rank-1 array of values for
    !! all mesh cells.

    use dm_utils_module, only: global_vec_section, global_section_offset
    use mesh_module, only: mesh_type
    use cell_order_module, only: cell_order_label_name
    use eos_module, only: eos_type
    use fluid_module, only: fluid_type

    type(mesh_type), intent(in) :: mesh
    PetscInt, intent(in) :: region(:)
    class(eos_type), intent(in) :: eos
    Vec, intent(in out) :: fluid_vector
    PetscInt, intent(in) :: fluid_range_start
    ! Locals:
    PetscInt :: num_cells, np, global_cell_index, i, ghost
    PetscInt :: offset, num_matching, c
    PetscErrorCode :: ierr
    type(fluid_type) :: fluid
    PetscReal, pointer, contiguous :: fluid_array(:)
    PetscSection :: section
    IS :: cell_IS
    PetscInt, pointer :: cells(:)
    DMLabel :: ghost_label

    num_cells = size(region, 1)
    np = eos%num_primary_variables

    call global_vec_section(fluid_vector, section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call fluid%init(eos%num_components, eos%num_phases)
    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)

    do i = 1, num_cells
       global_cell_index = i - 1
       call DMGetStratumSize(mesh%dm, cell_order_label_name, &
            global_cell_index, num_matching, ierr); CHKERRQ(ierr)
       if (num_matching > 0) then
          call DMGetStratumIS(mesh%dm, cell_order_label_name, &
               global_cell_index, cell_IS, ierr); CHKERRQ(ierr)
          call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
          do c = 1, num_matching
             call DMLabelGetValue(ghost_label, cells(c), ghost, ierr)
             if (ghost < 0) then
                call global_section_offset(section, cells(c), &
                     fluid_range_start, offset, ierr); CHKERRQ(ierr)
                call fluid%assign(fluid_array, offset)
                fluid%region = dble(region(i))
             end if
          end do
          call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
          call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)
       end if
    end do

    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

  end subroutine setup_initial_region_array

!------------------------------------------------------------------------

  subroutine setup_initial_file(filename, mesh, eos, t, y, fluid_vector, &
       y_range_start, fluid_range_start, index)
    !! Initializes fluid vector and solution vector y from HDF5 file.

    use mesh_module
    use dm_utils_module, only: global_vec_section, global_section_offset
    use eos_module, only: eos_type
    use fluid_module, only: fluid_type

    character(len = *), intent(in) :: filename
    type(mesh_type), intent(in) :: mesh
    class(eos_type), intent(in) :: eos
    PetscReal, intent(in) :: t
    Vec, intent(in out) :: y, fluid_vector
    PetscInt, intent(in) :: y_range_start, fluid_range_start
    PetscInt, intent(in) :: index !! time index to fetch initial conditions from
    ! Locals:
    PetscViewer :: viewer
    DM :: fluid_dm
    IS :: output_cell_index
    PetscSection :: y_section, fluid_section
    PetscReal, pointer, contiguous :: y_array(:), fluid_array(:)
    PetscReal, pointer, contiguous :: cell_primary(:)
    type(fluid_type) :: fluid
    DMLabel :: ghost_label
    PetscInt :: np, c, ghost, y_offset, fluid_offset
    PetscErrorCode :: ierr

    call PetscViewerHDF5Open(PETSC_COMM_WORLD, filename, FILE_MODE_READ, &
         viewer, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5PushGroup(viewer, "/", ierr); CHKERRQ(ierr)

    call ISDuplicate(mesh%cell_index, output_cell_index, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(output_cell_index, "cell_index", ierr)
    CHKERRQ(ierr)
    call ISLoad(output_cell_index, viewer, ierr); CHKERRQ(ierr)

    call VecGetDM(fluid_vector, fluid_dm, ierr); CHKERRQ(ierr)
    call DMSetOutputSequenceNumber(fluid_dm, index, t, ierr); CHKERRQ(ierr)
    call VecLoad(fluid_vector, viewer, ierr); CHKERRQ(ierr)

    call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

    call mesh%order_vector(fluid_vector, output_cell_index)
    call ISDestroy(output_cell_index, ierr); CHKERRQ(ierr)

    ! Determine solution vector from fluid vector:

    call global_vec_section(y, y_section)
    call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)

    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayReadF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call fluid%init(eos%num_components, eos%num_phases)
    np = eos%num_primary_variables

    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    do c = mesh%start_cell, mesh%end_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          call global_section_offset(y_section, c, &
               y_range_start, y_offset, ierr); CHKERRQ(ierr)
          cell_primary => y_array(y_offset : y_offset + np - 1)
          call global_section_offset(fluid_section, c, &
               fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)
          call fluid%assign(fluid_array, fluid_offset)
          call eos%primary_variables(fluid, cell_primary)
       end if
    end do

    call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call fluid%destroy()

  end subroutine setup_initial_file

!------------------------------------------------------------------------

  subroutine scale_initial_primary(mesh, eos, y, fluid_vector, &
      y_range_start, fluid_range_start)
    !! Scales primary variables y in each cell according to the eos
    !! primary variable scaling.

    use mesh_module, only: mesh_type
    use eos_module, only: eos_type
    use fluid_module, only: fluid_type
    use dm_utils_module, only: global_vec_section, global_section_offset

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
    PetscErrorCode :: ierr

    call global_vec_section(y, y_section)
    call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)

    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayReadF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call fluid%init(eos%num_components, eos%num_phases)
    np = eos%num_primary_variables

    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    do c = mesh%start_cell, mesh%end_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          call global_section_offset(y_section, c, &
               y_range_start, y_offset, ierr); CHKERRQ(ierr)
          cell_primary => y_array(y_offset : y_offset + np - 1)
          call global_section_offset(fluid_section, c, &
               fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)
          call fluid%assign(fluid_array, fluid_offset)
          cell_primary = cell_primary / &
                     eos%primary_scale(:, nint(fluid%region))
       end if
    end do

    call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call fluid%destroy()

  end subroutine scale_initial_primary

!------------------------------------------------------------------------

  subroutine setup_initial(json, mesh, eos, t, y, fluid_vector, &
       y_range_start, fluid_range_start, logfile)
    !! Initializes time t and a Vec y with initial conditions read
    !! from JSON input 'initial'.  A full set of initial conditions
    !! may be read in from an HDF5 file, or a constant set of primary
    !! variables can be read from the JSON input and applied to all
    !! cells. The fluid thermodynamic region for each cell is also
    !! initialised (in the fluid vector).

    use mesh_module
    use eos_module
    use dm_utils_module, only: global_vec_section, global_section_offset
    use logfile_module
    use fson_value_m, only: TYPE_ARRAY

    type(fson_value), pointer, intent(in) :: json
    type(mesh_type), intent(in) :: mesh
    class(eos_type), intent(in) :: eos
    PetscReal, intent(out) :: t
    Vec, intent(in out) :: y, fluid_vector
    PetscInt, intent(in) :: y_range_start, fluid_range_start
    type(logfile_type), intent(in out) :: logfile
    ! Locals:
    PetscReal, parameter :: default_start_time = 0.0_dp
    PetscReal, allocatable :: primary(:), primary_array(:,:)
    PetscInt :: region
    PetscInt, allocatable :: region_array(:)
    PetscReal :: primary_scalar
    PetscInt :: primary_rank, region_rank
    PetscInt, parameter :: max_filename_length = 240
    character(len = max_filename_length) :: filename
    PetscInt, parameter :: default_index = 0
    PetscInt :: index

    call fson_get_mpi(json, "time.start", default_start_time, t, logfile)

    if (fson_has_mpi(json, "initial")) then

       if (fson_has_mpi(json, "initial.filename")) then

          call fson_get_mpi(json, "initial.filename", val = filename)
          call fson_get_mpi(json, "initial.index", default_index, index, logfile)
          call setup_initial_file(filename, mesh, eos, t, y, fluid_vector, &
               y_range_start, fluid_range_start, index)

       else if (fson_has_mpi(json, "initial.primary")) then

          primary_rank = fson_mpi_array_rank(json, "initial.primary")
          select case (primary_rank)
          case(0)
             call fson_get_mpi(json, "initial.primary", val = primary_scalar)
             primary = [primary_scalar]
          case(1)
             call fson_get_mpi(json, "initial.primary", val = primary)
          case(2)
             call fson_get_mpi(json, "initial.primary", val = primary_array)
          end select
          
          region_rank = fson_mpi_array_rank(json, "initial.region")
          select case (region_rank)
          case(-1)
             region = eos%default_region
             call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
                  int_keys = ['initial.region'], &
                  int_values = [eos%default_region])
          case(0)
             call fson_get_mpi(json, "initial.region", val = region)
          case(1)
             call fson_get_mpi(json, "initial.region", val = region_array)
          end select

          if (primary_rank < 2) then
             call setup_initial_primary_constant(mesh, primary, eos, y, &
                  y_range_start)
             deallocate(primary)
          else
             call setup_initial_primary_array(mesh, primary_array, eos, y, &
                  y_range_start)
             deallocate(primary_array)
          end if

          if (region_rank < 1) then
             call setup_initial_region_constant(mesh, region, eos, fluid_vector, &
                  fluid_range_start)
          else
             call setup_initial_region_array(mesh, region_array, eos, &
                  fluid_vector, fluid_range_start)
             deallocate(region_array)
          end if

       else
          call logfile%write(LOG_LEVEL_WARN, 'input', &
               '"unrecognised initial.primary"')
       end if

    else

       call setup_initial_primary_constant(mesh, eos%default_primary, &
            eos, y, y_range_start)
       call setup_initial_region_constant(mesh, eos%default_region, &
            eos, fluid_vector, fluid_range_start)
       call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
            real_array_key = 'initial.primary', &
            real_array_value = eos%default_primary)
       call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
            int_keys = ['initial.region'], &
            int_values = [eos%default_region])

    end if

  end subroutine setup_initial

!------------------------------------------------------------------------

  subroutine setup_minc_initial(mesh, y, fluid_vector, &
       y_range_start, fluid_range_start, eos, logfile)

    !! Sets up initial conditions in MINC matrix cells. These are
    !! simply copied from the corresponding fracture cells, which are
    !! assumed already initialised.

    use mesh_module
    use eos_module
    use dm_utils_module, only: global_vec_section, global_section_offset
    use fluid_module, only: fluid_type
    use minc_module, only: minc_zone_label_name
    use logfile_module

    class(mesh_type), intent(in out) :: mesh
    Vec, intent(in out) :: y, fluid_vector !! Solution and fluid vectors
    PetscInt, intent(in) :: y_range_start, fluid_range_start !! Range starts for vectors
    class(eos_type), intent(in) :: eos !!
    type(logfile_type), intent(in out), optional :: logfile !! Log file
    ! Locals:
    PetscSection :: y_section, fluid_section
    PetscReal, pointer, contiguous :: y_array(:), fluid_array(:)
    PetscInt :: iminc, c, h, i, m, p, num_minc_zone_cells, max_num_levels
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

    max_num_levels = ubound(mesh%minc_shift, 2)
    allocate(ic(max_num_levels))
    ic = 0
    h = 0

    do iminc = 1, size(mesh%minc)
       associate(minc => mesh%minc(iminc))
         call DMGetStratumSize(mesh%original_dm, minc_zone_label_name, iminc, &
              num_minc_zone_cells, ierr); CHKERRQ(ierr)
         if (num_minc_zone_cells > 0) then
            call DMGetStratumIS(mesh%original_dm, minc_zone_label_name, &
                 iminc, minc_IS, ierr); CHKERRQ(ierr)
            call ISGetIndicesF90(minc_IS, minc_cells, ierr); CHKERRQ(ierr)

            do i = 1, num_minc_zone_cells

               c = minc_cells(i)
               if (mesh%ghost_cell(c) < 0) then

                  call global_section_offset(y_section, c, y_range_start, &
                       y_offset, ierr); CHKERRQ(ierr)
                  associate(frac_primary => y_array(y_offset : y_offset + np - 1))
                    call global_section_offset(fluid_section, c, fluid_range_start, &
                         fluid_offset, ierr); CHKERRQ(ierr)
                    call fluid%assign(fluid_array, fluid_offset)
                    frac_region = nint(fluid%region)

                    do m = 1, minc%num_levels
                       p = ic(m) + mesh%minc_shift(h, m)
                       call global_section_offset(y_section, p, y_range_start, &
                            y_minc_offset, ierr); CHKERRQ(ierr)
                       associate (primary => y_array(y_minc_offset : &
                            y_minc_offset + np - 1))
                         primary = frac_primary
                       end associate
                       call global_section_offset(fluid_section, p, fluid_range_start, &
                            fluid_minc_offset, ierr); CHKERRQ(ierr)
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

  end subroutine setup_minc_initial

!------------------------------------------------------------------------

end module
