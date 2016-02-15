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

  subroutine setup_initial_constant(mesh, primary, region, &
       eos, y, fluid_vector, y_range_start, fluid_range_start)

    !! Initializes solution vector y and fluid regions in fluid vector
    !! with constant values over the mesh.

    use dm_utils_module, only: global_vec_section, global_section_offset
    use mesh_module, only: mesh_type
    use fluid_module, only: fluid_type
    use eos_module, only: eos_type

    type(mesh_type), intent(in) :: mesh
    PetscReal, intent(in) :: primary(:)
    PetscInt, intent(in) :: region
    class(eos_type), intent(in) :: eos
    Vec, intent(in out) :: y, fluid_vector
    PetscInt, intent(in) :: y_range_start, fluid_range_start
    ! Locals:
    PetscInt :: np, c, ghost
    PetscInt :: y_offset, fluid_offset
    PetscErrorCode :: ierr
    type(fluid_type) :: fluid
    PetscReal, pointer :: cell_primary(:), y_array(:), fluid_array(:)
    PetscSection :: y_section, fluid_section
    DMLabel :: ghost_label

    np = eos%num_primary_variables

    call global_vec_section(y, y_section)
    call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)

    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call fluid%init(eos%num_components, eos%num_phases)

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
          cell_primary = primary
          fluid%region = dble(region)
       end if
    end do

    call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call fluid%destroy()

  end subroutine setup_initial_constant

!------------------------------------------------------------------------

  subroutine setup_initial_file(filename, mesh, eos, t, y, fluid_vector, &
       y_range_start, fluid_range_start)
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
    ! Locals:
    PetscViewer :: viewer
    DM :: geom_dm, fluid_dm
    Vec :: output_cell_geom
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

    call VecGetDM(mesh%cell_geom, geom_dm, ierr); CHKERRQ(ierr)
    call DMGetGlobalVector(geom_dm, output_cell_geom, ierr)
    CHKERRQ(ierr)
    call PetscObjectSetName(output_cell_geom, "cell_geometry", ierr)
    CHKERRQ(ierr)
    call PetscViewerHDF5PushGroup(viewer, "geometry", ierr); CHKERRQ(ierr)
    call VecLoad(output_cell_geom, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)

    ! TODO :: navigate to last timestep- currently this will work only
    ! if the file has only results for one timestep in it.
    ! call PetscViewerHDF5SetTimestep(viewer, ?, ierr); CHKERRQ(ierr)

    call VecGetDM(fluid_vector, fluid_dm, ierr); CHKERRQ(ierr)
    call DMSetOutputSequenceNumber(fluid_dm, 0, t,ierr); CHKERRQ(ierr)
    call VecLoad(fluid_vector, viewer, ierr); CHKERRQ(ierr)

    call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

    call mesh%order_vector(output_cell_geom, fluid_vector)
    call DMRestoreGlobalVector(mesh%dm, output_cell_geom, ierr)
    CHKERRQ(ierr)

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
       logfile)
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

    type(fson_value), pointer, intent(in) :: json
    type(mesh_type), intent(in) :: mesh
    class(eos_type), intent(in) :: eos
    PetscReal, intent(out) :: t
    Vec, intent(in out) :: y, rock_vector, fluid_vector
    PetscInt, intent(in) :: y_range_start, rock_range_start, fluid_range_start
    type(logfile_type), intent(in out) :: logfile
    ! Locals:
    PetscReal, parameter :: default_start_time = 0.0_dp
    PetscReal, allocatable :: primary(:)
    PetscInt :: region
    PetscInt, parameter :: max_filename_length = 240
    character(len = max_filename_length) :: filename

    call fson_get_mpi(json, "time.start", default_start_time, t, logfile)

    if (fson_has_mpi(json, "initial")) then

       if (fson_has_mpi(json, "initial.filename")) then

          call fson_get_mpi(json, "initial.filename", val = filename)
          call setup_initial_file(filename, mesh, eos, t, y, fluid_vector, &
               y_range_start, fluid_range_start)

       else

          call fson_get_mpi(json, "initial.primary", eos%default_primary, &
               primary, logfile)
          call fson_get_mpi(json, "initial.region", eos%default_region, &
               region, logfile)

          call setup_initial_constant(mesh, primary, region, eos, y, &
               fluid_vector, y_range_start, fluid_range_start)
          deallocate(primary)

       end if
    else

       call setup_initial_constant(mesh, eos%default_primary, &
            eos%default_region, eos, y, fluid_vector, y_range_start, &
            fluid_range_start)
       call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
            real_array_key = 'initial.primary', &
            real_array_value = eos%default_primary, &
            output_rank_only = .true.)
       call logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
            int_keys = ['initial.region'], &
            int_values = [eos%default_region], &
            output_rank_only = .true.)

    end if

  end subroutine setup_initial

!------------------------------------------------------------------------

end module
