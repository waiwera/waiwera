module initial_module
  !! Module for setting up initial conditions.

  use kinds_module
  use fson
  use fson_mpi_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  public :: setup_initial

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
    PetscReal, pointer :: cell_primary(:), y_array(:)
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
    use mesh_module, only: mesh_type, cell_order_label_name
    use eos_module, only: eos_type

    type(mesh_type), intent(in) :: mesh
    PetscReal, intent(in) :: primary(:,:)
    class(eos_type), intent(in) :: eos
    Vec, intent(in out) :: y
    PetscInt, intent(in) :: range_start
    ! Locals:
    PetscInt :: num_cells, np, global_cell_index, i, c
    PetscInt :: offset
    PetscErrorCode :: ierr
    PetscReal, pointer :: cell_primary(:), y_array(:)
    PetscSection :: section
    IS :: cell_IS
    PetscInt, pointer :: cells(:)

    num_cells = size(primary, 1)
    np = eos%num_primary_variables

    call global_vec_section(y, section)
    call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)

    do i = 1, num_cells
       global_cell_index = i - 1
       call DMGetStratumIS(mesh%dm, cell_order_label_name, &
            global_cell_index, cell_IS, ierr); CHKERRQ(ierr)
       if (cell_IS /= 0) then
          call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
          if (size(cells) > 0) then
             c = cells(1)
             call global_section_offset(section, c, range_start, &
                  offset, ierr)
             cell_primary => y_array(offset : offset + np - 1)
             cell_primary = primary(i, :)
          end if
          call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
       end if
       call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)
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
    PetscReal, pointer :: fluid_array(:)
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
    use mesh_module, only: mesh_type, cell_order_label_name
    use eos_module, only: eos_type
    use fluid_module, only: fluid_type

    type(mesh_type), intent(in) :: mesh
    PetscInt, intent(in) :: region(:)
    class(eos_type), intent(in) :: eos
    Vec, intent(in out) :: fluid_vector
    PetscInt, intent(in) :: fluid_range_start
    ! Locals:
    PetscInt :: num_cells, np, global_cell_index, i, c
    PetscInt :: offset
    PetscErrorCode :: ierr
    type(fluid_type) :: fluid
    PetscReal, pointer :: fluid_array(:)
    PetscSection :: section
    IS :: cell_IS
    PetscInt, pointer :: cells(:)

    num_cells = size(region, 1)
    np = eos%num_primary_variables

    call global_vec_section(fluid_vector, section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call fluid%init(eos%num_components, eos%num_phases)

    do i = 1, num_cells
       global_cell_index = i - 1
       call DMGetStratumIS(mesh%dm, cell_order_label_name, &
            global_cell_index, cell_IS, ierr); CHKERRQ(ierr)
       if (cell_IS /= 0) then
          call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
          if (size(cells) > 0) then
             c = cells(1)
             call global_section_offset(section, c, &
                  fluid_range_start, offset, ierr); CHKERRQ(ierr)
             call fluid%assign(fluid_array, offset)
             fluid%region = dble(region(i))
          end if
          call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
       end if
       call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)
    end do

    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

  end subroutine setup_initial_region_array

!------------------------------------------------------------------------

  subroutine setup_initial_file(filename, mesh, eos, t, y, fluid_vector, &
       y_range_start, fluid_range_start, fluid_initialized)
    !! Initializes fluid vector and solution vector y from HDF5 file.

    use mpi_module
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
    PetscBool, intent(out) :: fluid_initialized
    ! Locals:
    PetscViewer :: viewer
    DM :: fluid_dm
    IS :: output_cell_index
    PetscSection :: y_section, fluid_section
    PetscReal, pointer :: y_array(:), fluid_array(:)
    PetscReal, pointer :: cell_primary(:)
    type(fluid_type) :: fluid
    DMLabel :: ghost_label
    PetscInt :: np, c, ghost, y_offset, fluid_offset
    PetscErrorCode :: ierr

    call PetscViewerHDF5Open(mpi%comm, filename, FILE_MODE_READ, &
         viewer, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5PushGroup(viewer, "/", ierr); CHKERRQ(ierr)

    call ISDuplicate(mesh%cell_index, output_cell_index, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(output_cell_index, "cell_index", ierr)
    CHKERRQ(ierr)
    call ISLoad(output_cell_index, viewer, ierr); CHKERRQ(ierr)

    ! TODO :: navigate to last timestep- currently this will work only
    ! if the file has only results for one timestep in it.
    ! call PetscViewerHDF5SetTimestep(viewer, ?, ierr); CHKERRQ(ierr)

    call VecGetDM(fluid_vector, fluid_dm, ierr); CHKERRQ(ierr)
    call DMSetOutputSequenceNumber(fluid_dm, 0, t,ierr); CHKERRQ(ierr)
    call VecLoad(fluid_vector, viewer, ierr); CHKERRQ(ierr)

    call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

    call mesh%order_vector(fluid_vector, output_cell_index)
    call ISDestroy(output_cell_index, ierr); CHKERRQ(ierr)
    fluid_initialized = PETSC_TRUE

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

  subroutine setup_initial(json, mesh, eos, t, y, rock_vector, &
       fluid_vector, y_range_start, rock_range_start, fluid_range_start, &
       logfile, fluid_initialized)
    !! Initializes time t and a Vec y with initial conditions read
    !! from JSON input 'initial'.  A full set of initial conditions
    !! may be read in from an HDF5 file, or a constant set of primary
    !! variables can be read from the JSON input and applied to all
    !! cells.

    use mesh_module
    use eos_module
    use rock_module
    use dm_utils_module, only: global_vec_section, global_section_offset
    use logfile_module
    use fson_value_m, only: TYPE_ARRAY

    type(fson_value), pointer, intent(in) :: json
    type(mesh_type), intent(in) :: mesh
    class(eos_type), intent(in) :: eos
    PetscReal, intent(out) :: t
    Vec, intent(in out) :: y, rock_vector, fluid_vector
    PetscInt, intent(in) :: y_range_start, rock_range_start, fluid_range_start
    type(logfile_type), intent(in out) :: logfile
    PetscBool, intent(out) :: fluid_initialized
    ! Locals:
    PetscReal, parameter :: default_start_time = 0.0_dp
    PetscReal, allocatable :: primary(:), primary_array(:,:)
    PetscInt :: region
    PetscInt, allocatable :: region_array(:)
    PetscReal :: primary_scalar
    PetscInt :: primary_rank, region_rank
    PetscInt, parameter :: max_filename_length = 240
    character(len = max_filename_length) :: filename

    call fson_get_mpi(json, "time.start", default_start_time, t, logfile)

    if (fson_has_mpi(json, "initial")) then

       if (fson_has_mpi(json, "initial.filename")) then

          call fson_get_mpi(json, "initial.filename", val = filename)
          call setup_initial_file(filename, mesh, eos, t, y, fluid_vector, &
               y_range_start, fluid_range_start, fluid_initialized)

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

end module
