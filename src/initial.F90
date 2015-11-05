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

    call DMPlexGetLabel(mesh%dm, "ghost", ghost_label, ierr)
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

  subroutine setup_initial_file(filename, mesh, eos, y, fluid_vector, &
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
    Vec, intent(in out) :: y, fluid_vector
    PetscInt, intent(in) :: y_range_start, fluid_range_start
    ! Locals:
    PetscViewer :: viewer
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

    ! TODO :: navigate to last timestep:
    ! call PetscViewerHDF5SetTimestep(viewer, -1, ierr); CHKERRQ(ierr)

    call VecLoad(fluid_vector, viewer, ierr); CHKERRQ(ierr)

    call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

    ! Determine solution vector from fluid vector:

    call global_vec_section(y, y_section)
    call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)

    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayReadF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call fluid%init(eos%num_components, eos%num_phases)
    np = eos%num_primary_variables

    call DMPlexGetLabel(mesh%dm, "ghost", ghost_label, ierr)
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
       fluid_vector, y_range_start, rock_range_start, fluid_range_start)
    !! Initializes time t and a Vec y with initial conditions read
    !! from JSON input 'initial'.  Conditions may be specified as a
    !! constant value or as an array. The array may contain a complete
    !! set of initial conditions for all cells, or if a shorter array is
    !! given, this is repeated over initial conditions vector.

    use mesh_module
    use eos_module
    use rock_module
    use dm_utils_module, only: global_vec_section, global_section_offset

    type(fson_value), pointer, intent(in) :: json
    type(mesh_type), intent(in) :: mesh
    class(eos_type), intent(in) :: eos
    PetscReal, intent(out) :: t
    Vec, intent(in out) :: y, rock_vector, fluid_vector
    PetscInt, intent(in) :: y_range_start, rock_range_start, fluid_range_start
    ! Locals:
    PetscReal, parameter :: default_start_time = 0.0_dp
    PetscReal, allocatable :: primary(:)
    PetscInt :: region
    PetscInt, parameter :: max_filename_length = 240
    character(len = max_filename_length) :: filename

    call fson_get_mpi(json, "time.start", default_start_time, t)

    if (fson_has_mpi(json, "initial")) then

       if (fson_has_mpi(json, "initial.filename")) then

          call fson_get_mpi(json, "initial.filename", val = filename)
          call setup_initial_file(filename, mesh, eos, y, fluid_vector, &
               y_range_start, fluid_range_start)

       else

          call fson_get_mpi(json, "initial.primary", eos%default_primary, primary)
          call fson_get_mpi(json, "initial.region", eos%default_region, region)

          call setup_initial_constant(mesh, primary, region, eos, y, &
               fluid_vector, y_range_start, fluid_range_start)
          deallocate(primary)

       end if
    else

       call setup_initial_constant(mesh, eos%default_primary, &
            eos%default_region, eos, y, fluid_vector, y_range_start, &
            fluid_range_start)

    end if

  end subroutine setup_initial

!------------------------------------------------------------------------

end module
