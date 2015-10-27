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
       num_components, num_phases, y, fluid_vector, y_range_start, &
       fluid_range_start)

    !! Initializes solution vector y and fluid regions in fluid vector
    !! with constant values over the mesh.

    use dm_utils_module, only: global_vec_section, global_section_offset
    use mesh_module, only: mesh_type
    use fluid_module, only: fluid_type

    type(mesh_type), intent(in) :: mesh
    PetscReal, intent(in) :: primary(:)
    PetscInt, intent(in) :: region, num_components, num_phases
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

    np = size(primary)

    call global_vec_section(y, y_section)
    call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)

    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call fluid%init(num_components, num_phases)

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

  subroutine setup_initial_file(filename)
    !! Initializes solution vector y and fluid regions in fluid vector
    !! from HDF5 file.

    character(len = *), intent(in) :: filename

    ! ** TODO **

  end subroutine setup_initial_file

!------------------------------------------------------------------------

  subroutine setup_initial(json, mesh, eos, t, y, rock_vector, &
       fluid_vector, y_range_start, rock_range_start, fluid_range_start)
    !! Initializes time t and a Vec y with initial conditions read
    !! from JSON input 'initial'.  Conditions may be specified as a
    !! constant value or as an array. The array may contain a complete
    !! set of initial conditions for all cells, or if a shorter array is
    !! given, this is repeated over initial conditions vector.

    use fson_value_m, only : TYPE_OBJECT, TYPE_STRING
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

       select case (fson_type_mpi(json, "initial"))

       case (TYPE_OBJECT)

          call fson_get_mpi(json, "initial.primary", eos%default_primary, primary)
          call fson_get_mpi(json, "initial.region", eos%default_region, region)

          call setup_initial_constant(mesh, primary, region, eos%num_components, &
               eos%num_phases, y, fluid_vector, y_range_start, fluid_range_start)
          deallocate(primary)

       case (TYPE_STRING)

          call fson_get_mpi(json, "initial", val = filename)
          call setup_initial_file(filename)

       end select

    else
       call setup_initial_constant(mesh, eos%default_primary, eos%default_region, &
            eos%num_components, eos%num_phases, y, fluid_vector, y_range_start, &
            fluid_range_start)
    end if

  end subroutine setup_initial

!------------------------------------------------------------------------

end module
