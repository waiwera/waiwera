module initial_test

  ! Tests for initial module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use fruit
  use fson
  use initial_module
  use fson_mpi_module

  implicit none
  private

  public :: test_initial_hdf5

contains

!------------------------------------------------------------------------

  subroutine test_initial_hdf5
    ! HDF5 initial conditions

    use IAPWS_module
    use eos_we_module
    use mesh_module
    use eos_module, only: max_component_name_length, &
         max_phase_name_length
    use fluid_module, only: fluid_type, setup_fluid_vector
    use dm_utils_module, only: global_vec_section, global_section_offset, &
         local_vec_section, section_offset, global_vec_range_start
    use rock_module, only: rock_type
    use relative_permeability_module, only: relative_permeability_corey_type
    use capillary_pressure_module, only: capillary_pressure_zero_type
    use cell_module, only: cell_type
    use flow_simulation_test, only: vec_write
    use logfile_module

    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    type(mesh_type) :: mesh
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscMPIInt :: rank
    PetscErrorCode :: ierr, err
    PetscViewer :: viewer
    Vec :: fluid_vector, y
    PetscInt :: y_range_start, fluid_range_start, start_cell, end_cell
    PetscInt :: c, ghost, fluid_offset, cell_geom_offset
    PetscSection :: fluid_section, cell_geom_section
    PetscReal, pointer, contiguous :: fluid_array(:)
    DMLabel :: ghost_label
    type(fluid_type) :: fluid
    PetscReal, allocatable :: primary(:), expected_primary(:)
    PetscReal, pointer, contiguous :: rock_data(:)
    type(rock_type) :: rock
    type(relative_permeability_corey_type) :: relative_permeability
    type(capillary_pressure_zero_type) :: capillary_pressure
    PetscReal, pointer, contiguous :: cell_geom_array(:)
    type(cell_type) :: cell
    PetscReal :: t
    type(logfile_type) :: logfile
    PetscReal, parameter :: permeability(3) = [1.e-13_dp, 1.e-13_dp, 1.e-14_dp]
    PetscReal, parameter :: porosity = 0.1_dp
    PetscReal, parameter :: wet_conductivity = 2.5_dp
    PetscReal, parameter :: dry_conductivity = 1.5_dp
    PetscReal, parameter :: rock_density = 2200._dp
    PetscReal, parameter :: specific_heat = 1000._dp
    PetscInt, parameter :: rock_offset = 1
    PetscReal, parameter :: tol = 1.e-3

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)
    viewer = PETSC_NULL_VIEWER
    call logfile%init('', echo = PETSC_FALSE)

    json => fson_parse_mpi(str = &
         '{"mesh": {"filename": "data/mesh/col100.exo"},' // &
         ' "eos": {"name": "we"}, ' // &
         ' "initial": {"filename": "data/initial/fluid.h5"}}')
    call thermo%init()
    call eos%init(json, thermo)
    call mesh%init(json)
    call mesh%configure(eos, gravity, json, viewer = viewer, err = err)

    call DMCreateGlobalVector(mesh%dm, y, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(y, "primary", ierr); CHKERRQ(ierr)
    call global_vec_range_start(y, y_range_start)
    
    call setup_fluid_vector(mesh%dm, max_component_name_length, &
         eos%component_names, max_phase_name_length, &
         eos%phase_names, fluid_vector, fluid_range_start)

    t = 0._dp
    call setup_initial(json, mesh, eos, t, y, fluid_vector, &
         y_range_start, fluid_range_start, logfile)
    
    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call rock%init()
    allocate(rock_data(rock%dof))
    rock_data = [permeability, wet_conductivity, &
         dry_conductivity, porosity, rock_density, specific_heat]
    call rock%assign(rock_data, rock_offset)
    call relative_permeability%init(json)
    call rock%assign_relative_permeability(relative_permeability)
    call capillary_pressure%init(json)
    call rock%assign_capillary_pressure(capillary_pressure)
    call fson_destroy_mpi(json)
    
    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)
    call fluid%init(eos%num_components, eos%num_phases)
    allocate(primary(eos%num_primary_variables), &
         expected_primary(eos%num_primary_variables))

    call local_vec_section(mesh%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)
    call cell%init(eos%num_components, eos%num_phases)

    do c = start_cell, end_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          call global_section_offset(fluid_section, c, &
               fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)
          call fluid%assign(fluid_array, fluid_offset)
          call section_offset(cell_geom_section, c, cell_geom_offset, ierr)
          CHKERRQ(ierr)
          call cell%assign_geometry(cell_geom_array, cell_geom_offset)
          associate(z => cell%centroid(3))
            call primary_variables(z, expected_primary)
          end associate
          call eos%primary_variables(fluid, primary)
          call assert_equals(expected_primary, primary, &
               eos%num_primary_variables, tol, 'primary')
       end if
    end do
    
    call VecRestoreArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)
    call fluid%destroy()
    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call VecDestroy(fluid_vector, ierr); CHKERRQ(ierr)
    call VecDestroy(y, ierr); CHKERRQ(ierr)
    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()
    call rock%destroy()
    call relative_permeability%destroy()
    call capillary_pressure%destroy()
    call cell%destroy()
    call logfile%destroy()
    deallocate(primary, expected_primary, rock_data)

  contains

    subroutine primary_variables(z, primary)
      !! Return primary variables as function of elevation z.

      PetscReal, intent(in) :: z
      PetscReal, intent(out) :: primary(:)
      ! Locals:
      PetscReal, parameter :: P0 = 1.e5_dp, rho = 997._dp
      PetscReal, parameter :: T0 = 20._dp, dtdz = 25._dp / 1.e3_dp

      associate(pressure => primary(1), temperature => primary(2))
        pressure = P0 + rho * gravity(3) * z
        temperature = T0 - dtdz * z
      end associate

    end subroutine primary_variables

  end subroutine test_initial_hdf5

!------------------------------------------------------------------------

end module initial_test
