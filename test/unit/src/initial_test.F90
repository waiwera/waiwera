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
    use relative_permeability_module, only: relative_permeability_corey_type
    use capillary_pressure_module, only: capillary_pressure_zero_type
    use cell_module, only: cell_type
    use flow_simulation_test, only: vec_write
    use logfile_module

    character(:), allocatable :: json_str
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)

    json_str = &
         '{"mesh": {"filename": "data/mesh/col100.exo"},' // &
         ' "eos": {"name": "we"}, ' // &
         ' "initial": {"filename": "data/initial/fluid.h5"}}'
    call initial_hdf5_test('single porosity', json_str)

    json_str = &
         '{"mesh": {"filename": "data/mesh/col100.exo", ' // &
         '          "zones": {"all": {"-": null}},' // &
         '          "minc": {"rock": {"zones": ["all"]}, ' // &
         '                   "geometry": {"fracture": {"volume": 0.1}, ' // &
         '                                "matrix": {"volume": [0.3, 0.6]}}}},' // &
         ' "eos": {"name": "we"}, ' // &
         ' "initial": {"filename": "data/initial/fluid.h5", "minc": false}}'
    call initial_hdf5_test('MINC', json_str)

  contains

!........................................................................

    subroutine primary_variables(z, g, primary)
      !! Return primary variables as function of elevation z.

      PetscReal, intent(in) :: z, g
      PetscReal, intent(out) :: primary(:)
      ! Locals:
      PetscReal, parameter :: P0 = 1.e5_dp, rho = 997._dp
      PetscReal, parameter :: T0 = 20._dp, dtdz = 25._dp / 1.e3_dp

      associate(pressure => primary(1), temperature => primary(2))
        pressure = P0 + rho * g * z
        temperature = T0 - dtdz * z
      end associate

    end subroutine primary_variables

!........................................................................

    subroutine initial_hdf5_test(name, json_str)

      character(*), intent(in) :: name
      character(*), intent(in) :: json_str
      ! Locals:
      type(fson_value), pointer :: json
      type(IAPWS_type) :: thermo
      type(eos_we_type) :: eos
      type(mesh_type) :: mesh
      PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
      PetscViewer :: viewer
      Vec :: fluid_vector, y
      PetscInt :: y_range_start, fluid_range_start, start_cell, end_cell
      PetscInt :: c, ghost, fluid_offset, cell_geom_offset, y_offset
      PetscSection :: fluid_section, cell_geom_section, y_section
      PetscReal, pointer, contiguous :: fluid_array(:), y_array(:)
      DMLabel :: ghost_label
      type(fluid_type) :: fluid
      PetscReal, allocatable :: expected_primary(:)
      PetscReal, pointer, contiguous :: primary(:)
      type(relative_permeability_corey_type) :: relative_permeability
      type(capillary_pressure_zero_type) :: capillary_pressure
      PetscReal, pointer, contiguous :: cell_geom_array(:)
      type(cell_type) :: cell
      PetscReal :: t
      type(logfile_type) :: logfile
      PetscErrorCode :: err
      PetscReal, parameter :: tol = 1.e-3
      PetscInt, parameter :: expected_region = 1

      json => fson_parse_mpi(str = json_str)
      viewer = PETSC_NULL_VIEWER
      call logfile%init('', echo = PETSC_FALSE)

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
      call global_vec_section(y, y_section)
      call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)

      call fson_destroy_mpi(json)

      call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
      CHKERRQ(ierr)
      call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr)
      CHKERRQ(ierr)
      call fluid%init(eos%num_components, eos%num_phases)
      allocate(expected_primary(eos%num_primary_variables))

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
              call primary_variables(z, gravity(3), expected_primary)
            end associate
            call global_section_offset(y_section, c, &
                 y_range_start, y_offset, ierr); CHKERRQ(ierr)
            primary => y_array(y_offset: y_offset + eos%num_primary_variables - 1)
            call assert_equals(expected_primary, primary, &
                 eos%num_primary_variables, tol, name // ': primary')
            call assert_equals(expected_region, nint(fluid%region), &
                 name // ': region')
         end if
      end do

      call VecRestoreArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
      CHKERRQ(ierr)
      call fluid%destroy()
      call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
      call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)

      call VecDestroy(fluid_vector, ierr); CHKERRQ(ierr)
      call VecDestroy(y, ierr); CHKERRQ(ierr)
      call mesh%destroy()
      call eos%destroy()
      call thermo%destroy()
      call relative_permeability%destroy()
      call capillary_pressure%destroy()
      call cell%destroy()
      call logfile%destroy()
      deallocate(expected_primary)

    end subroutine initial_hdf5_test

  end subroutine test_initial_hdf5

!------------------------------------------------------------------------

end module initial_test
