module initial_test

  ! Tests for initial module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use fson
  use initial_module
  use fson_mpi_module

  implicit none
  private

  character(len = 512) :: data_path

  public :: setup, teardown
  public :: test_initial, test_initial_region
  ! public :: test_setup_initial_hdf5

contains

!------------------------------------------------------------------------

  subroutine setup()

    use profiling_module, only: init_profiling

    ! Locals:
    PetscErrorCode :: ierr
    PetscInt :: ios

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call init_profiling()

    call get_environment_variable('WAIWERA_TEST_DATA_PATH', &
         data_path, status = ios)
    if (ios /= 0) data_path = ''

  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown()

    PetscErrorCode :: ierr

    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

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

!------------------------------------------------------------------------

  subroutine test_initial(test)
    ! initial conditions

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    character(:), allocatable :: json_str
    PetscMPIInt :: rank
    character(:), allocatable :: initial_primary_json
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)

    ! HDF5 initial tests:

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/col100.exo"},' // &
         ' "eos": {"name": "we"}, ' // &
         ' "initial": {"filename": "' // trim(adjustl(data_path)) // 'initial/fluid.h5"}}'
    call initial_test_case('single porosity', json_str)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/col100.exo"},' // &
         ' "eos": {"name": "we"}, ' // &
         ' "initial": {"filename": "' // trim(adjustl(data_path)) // 'initial/fluid_minimal.h5"}}'
    call initial_test_case('single porosity minimal', json_str)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/col100.exo", ' // &
         '          "zones": {"all": {"-": null}},' // &
         '          "minc": {"rock": {"zones": ["all"]}, ' // &
         '                   "geometry": {"fracture": {"volume": 0.1}, ' // &
         '                                "matrix": {"volume": [0.3, 0.6]}}}},' // &
         ' "eos": {"name": "we"}, ' // &
         ' "initial": {"filename": "' // trim(adjustl(data_path)) // 'initial/fluid.h5", "minc": false}}'
    call initial_test_case('MINC', json_str)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/col100.exo", ' // &
         '          "zones": {"all": {"-": null}},' // &
         '          "minc": {"rock": {"zones": ["all"]}, ' // &
         '                   "geometry": {"fracture": {"volume": 0.1}, ' // &
         '                                "matrix": {"volume": [0.3, 0.6]}}}},' // &
         ' "eos": {"name": "we"}, ' // &
         ' "initial": {"filename": "' // trim(adjustl(data_path)) // 'initial/fluid_minimal.h5",' // &
         ' "minc": false}}'
    call initial_test_case('MINC minimal false', json_str)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/col100.exo", ' // &
         '          "zones": {"all": {"-": null}},' // &
         '          "minc": {"rock": {"zones": ["all"]}, ' // &
         '                   "geometry": {"fracture": {"volume": 0.1}, ' // &
         '                                "matrix": {"volume": [0.3, 0.6]}}}},' // &
         ' "eos": {"name": "we"}, ' // &
         ' "initial": {"filename": "' // trim(adjustl(data_path)) // 'initial/fluid_minimal_minc.h5",' // &
         ' "minc": true}}'
    call initial_test_case('MINC minimal true', json_str)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/col100.exo", ' // &
         '          "boundaries": [{"faces": {"cells": [0], ' // &
         '                          "normal": [0, 1, 0]}}]},' // &
         ' "eos": {"name": "we"}, ' // &
         ' "initial": {"filename": "' // trim(adjustl(data_path)) // 'initial/fluid_minimal.h5"}}'
    call initial_test_case('single porosity minimal boundary', json_str)

    ! JSON initial tests:

    initial_primary_json = &
         '[ 588530.0, 21.25], ' // &
         '[1565590.0, 23.75], ' // &
         '[2542650.0, 26.25], ' // &
         '[3519710.0, 28.75], ' // &
         '[4496770.0, 31.25], ' // &
         '[5473830.0, 33.75], ' // &
         '[6450890.0, 36.25], ' // &
         '[7427950.0, 38.75], ' // &
         '[8405010.0, 41.25], ' // &
         '[9382070.0, 43.75]'

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/col10.exo"},' // &
         ' "eos": {"name": "we"}, ' // &
         '   "initial": {"primary": [' // &
              initial_primary_json // &
         '    ] , "region": 1}}'
    call initial_test_case('JSON initial', json_str)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/col10.exo"},' // &
         ' "boundaries": [{"faces": {"cells": [0], "normal": [0, 1, 0]}}], ' // &
         ' "eos": {"name": "we"}, ' // &
         '   "initial": {"primary": [' // &
              initial_primary_json // &
         '    ], "region": 1}}'
    call initial_test_case('JSON initial boundary', json_str)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/col10.exo", ' // &
         '          "zones": {"all": {"-": null}},' // &
         '          "minc": {"rock": {"zones": ["all"]}, ' // &
         '                   "geometry": {"fracture": {"volume": 0.1}, ' // &
         '                                "matrix": {"volume": [0.3, 0.6]}}}},' // &
         ' "eos": {"name": "we"}, ' // &
         ' "initial": {"primary": [' // &
            initial_primary_json // &
         '  ], "minc": false}}'
    call initial_test_case('MINC initial JSON false', json_str)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/col10.exo", ' // &
         '          "zones": {"all": {"-": null}},' // &
         '          "minc": {"rock": {"zones": ["all"]}, ' // &
         '                   "geometry": {"fracture": {"volume": 0.1}, ' // &
         '                                "matrix": {"volume": [0.3, 0.6]}}}},' // &
         ' "boundaries": [{"faces": {"cells": [0], "normal": [0, 1, 0]}}], ' // &
         ' "eos": {"name": "we"}, ' // &
         ' "initial": {"primary": [' // &
            initial_primary_json // &
         '  ], "minc": false}}'
    call initial_test_case('MINC initial boundary JSON false', json_str)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/col10.exo", ' // &
         '          "zones": {"all": {"-": null}},' // &
         '          "minc": {"rock": {"zones": ["all"]}, ' // &
         '                   "geometry": {"fracture": {"volume": 0.1}, ' // &
         '                                "matrix": {"volume": 0.9}}}},' // &
         ' "eos": {"name": "we"}, ' // &
         ' "initial": {"primary": [' // &
            initial_primary_json // ', ' // &
            initial_primary_json // &
         '  ], "minc": true}}'
    call initial_test_case('MINC initial JSON true', json_str)

  contains

!........................................................................

    subroutine initial_test_case(name, json_str)

      use IAPWS_module
      use eos_we_module
      use mesh_module
      use eos_module, only: max_component_name_length, &
           max_phase_name_length
      use fluid_module, only: fluid_type, create_fluid_vector
      use dm_utils_module
      use relative_permeability_module, only: relative_permeability_corey_type
      use capillary_pressure_module, only: capillary_pressure_zero_type
      use cell_module, only: cell_type
      use logfile_module

      character(*), intent(in) :: name
      character(*), intent(in) :: json_str
      ! Locals:
      type(fson_value), pointer :: json
      type(IAPWS_type) :: thermo
      type(eos_we_type) :: eos
      type(mesh_type) :: mesh
      PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
      Vec :: fluid_vector, y
      PetscInt :: y_range_start, fluid_range_start
      PetscInt :: start_cell, end_cell, end_interior_cell
      PetscInt :: c, ghost, fluid_offset, cell_geom_offset, y_offset
      PetscSection :: fluid_section, cell_geom_section, y_section
      PetscReal, pointer, contiguous :: fluid_array(:), y_array(:)
      DMLabel :: ghost_label
      type(fluid_type) :: fluid
      PetscReal, allocatable :: expected_primary(:)
      PetscReal, pointer, contiguous :: primary(:)
      PetscReal, pointer, contiguous :: cell_geom_array(:)
      type(cell_type) :: cell
      PetscReal :: t
      type(logfile_type) :: logfile
      PetscErrorCode :: err
      PetscInt, parameter :: expected_region = 1

      json => fson_parse_mpi(str = json_str)
      call logfile%init('', echo = PETSC_FALSE)

      call thermo%init()
      call eos%init(json, thermo)
      call mesh%init(eos, json)
      call mesh%configure(gravity, json, err = err)
      call DMCreateGlobalVector(mesh%dm, y, ierr); CHKERRQ(ierr)
      call PetscObjectSetName(y, "primary", ierr); CHKERRQ(ierr)
      call global_vec_range_start(y, y_range_start)

      call create_fluid_vector(mesh%dm, max_component_name_length, &
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
      call mesh%destroy_distribution_data()

      call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
      CHKERRQ(ierr)
      end_interior_cell = dm_get_end_interior_cell(mesh%dm, end_cell)
      call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr)
      CHKERRQ(ierr)
      call fluid%init(eos%num_components, eos%num_phases)
      allocate(expected_primary(eos%num_primary_variables))

      call local_vec_section(mesh%cell_geom, cell_geom_section)
      call VecGetArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
      CHKERRQ(ierr)
      call cell%init(eos%num_components, eos%num_phases)

      do c = start_cell, end_interior_cell - 1
         call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
         if (ghost < 0) then
            fluid_offset = global_section_offset(fluid_section, c, &
                 fluid_range_start)
            call fluid%assign(fluid_array, fluid_offset)
            cell_geom_offset = section_offset(cell_geom_section, c)
            call cell%assign_geometry(cell_geom_array, cell_geom_offset)
            associate(z => cell%centroid(3))
              call primary_variables(z, gravity(3), expected_primary)
            end associate
            y_offset = global_section_offset(y_section, c, y_range_start)
            primary => y_array(y_offset: y_offset + eos%num_primary_variables - 1)
            call test%assert(expected_primary, primary, &
                 name // ': primary')
            call test%assert(expected_region, nint(fluid%region), &
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
      call cell%destroy()
      call logfile%destroy()
      deallocate(expected_primary)

    end subroutine initial_test_case

  end subroutine test_initial

!------------------------------------------------------------------------

  subroutine write_initial_hdf5(test)
    !! Writes HDF5 reference results for initial conditions unit tests.
    !! Rename to test_setup_initial_hdf5() to execute.

    use IAPWS_module
    use eos_we_module
    use mesh_module
    use eos_module, only: max_component_name_length, &
         max_phase_name_length
    use fluid_module, only: fluid_type, create_fluid_vector
    use rock_module, only: rock_type
    use dm_utils_module, only: global_vec_section, global_section_offset, &
         local_vec_section, section_offset, global_vec_range_start, &
         section_get_field_names
    use relative_permeability_module, only: relative_permeability_corey_type
    use capillary_pressure_module, only: capillary_pressure_zero_type
    use cell_module, only: cell_type
    use unit_test_utils_module, only: vec_write
    use logfile_module
    use hdf5io_module, only: max_field_name_length
    use utils_module, only: str_array_index, str_to_lower

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    character(:), allocatable :: json_str
    PetscMPIInt :: rank
    PetscErrorCode :: ierr
    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    type(mesh_type) :: mesh
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    Vec :: fluid_vector
    PetscInt :: fluid_range_start, start_cell, end_cell, i
    PetscInt :: c, ghost, fluid_offset, cell_geom_offset
    DM :: fluid_dm
    PetscSection :: fluid_section, cell_geom_section, local_fluid_section
    PetscReal, pointer, contiguous :: fluid_array(:)
    DMLabel :: ghost_label
    type(rock_type) :: rock
    type(fluid_type) :: fluid
    PetscReal, allocatable :: primary(:)
    type(relative_permeability_corey_type) :: relative_permeability
    type(capillary_pressure_zero_type) :: capillary_pressure
    PetscReal, pointer, contiguous :: cell_geom_array(:)
    type(cell_type) :: cell
    type(logfile_type) :: logfile
    character(max_field_name_length), allocatable :: fields(:)
    PetscInt, allocatable :: output_field_indices(:)
    PetscErrorCode :: err
    PetscInt, parameter :: region = 1

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/col100.exo",' // &
         '          "zones": {"all": {"-": null}},' // &
         '          "minc": {"rock": {"zones": ["all"]}, ' // &
         '                   "geometry": {"fracture": {"volume": 0.1}, ' // &
         '                                "matrix": {"volume": [0.3, 0.6]}}}},' // &
         ' "eos": "we"}'

    json => fson_parse_mpi(str = json_str)
    call logfile%init('', echo = PETSC_FALSE)

    call thermo%init()
    call eos%init(json, thermo)
    call mesh%init(eos, json)
    call mesh%configure(gravity, json, err = err)

    call create_fluid_vector(mesh%dm, max_component_name_length, &
         eos%component_names, max_phase_name_length, &
         eos%phase_names, fluid_vector, fluid_range_start)
    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call VecGetDM(fluid_vector, fluid_dm, ierr); CHKERRQ(ierr)
    call DMGetSection(fluid_dm, local_fluid_section, ierr); CHKERRQ(ierr)
    call section_get_field_names(local_fluid_section, PETSC_TRUE, fields)
    ! Required fields:
    associate(num_fields => size(eos%required_output_fluid_fields))
      allocate(output_field_indices(num_fields))
      do i = 1, num_fields
         output_field_indices(i) = str_array_index( &
              str_to_lower(eos%required_output_fluid_fields(i)), fields) - 1
         call test%assert(output_field_indices(i) >= 0, "setup field")
      end do
    end associate

    ! All fields:
    ! associate(num_fields => size(fields))
    !   allocate(output_field_indices(num_fields))
    !   output_field_indices = [(i, i = 0, num_fields - 1)]
    ! end associate

    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)
    call fluid%init(eos%num_components, eos%num_phases)
    call local_vec_section(mesh%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)
    call cell%init(eos%num_components, eos%num_phases)
    call rock%init()
    call relative_permeability%init(json, err = err)
    call capillary_pressure%init(json, err = err)
    call rock%assign_relative_permeability(relative_permeability)
    call rock%assign_capillary_pressure(capillary_pressure)
    allocate(primary(eos%num_primary_variables))
    call fson_destroy_mpi(json)
    call mesh%destroy_distribution_data()

    do c = start_cell, end_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          fluid_offset = global_section_offset(fluid_section, c, fluid_range_start)
          call fluid%assign(fluid_array, fluid_offset)
          cell_geom_offset = section_offset(cell_geom_section, c)
          call cell%assign_geometry(cell_geom_array, cell_geom_offset)
          associate(z => cell%centroid(3))
            call primary_variables(z, gravity(3), primary)
          end associate
          fluid%region = region
          call eos%bulk_properties(primary, fluid, err)
          call eos%phase_properties(primary, rock, fluid, err)
       end if
    end do

    call vec_write(fluid_vector, "fluid_minimal_minc", trim(adjustl(data_path)) // "initial/", &
         mesh%cell_index, output_field_indices, "/cell_fields")

    call VecRestoreArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)
    call fluid%destroy()
    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call VecDestroy(fluid_vector, ierr); CHKERRQ(ierr)
    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()
    call relative_permeability%destroy()
    call capillary_pressure%destroy()
    call cell%destroy()
    call logfile%destroy()
    deallocate(primary)
    call rock%destroy()

  end subroutine write_initial_hdf5

!------------------------------------------------------------------------

  subroutine test_initial_region(test)
    ! initial region array

    use IAPWS_module
    use eos_we_module
    use mesh_module
    use dm_utils_module
    use logfile_module
    use eos_module, only: max_component_name_length, &
         max_phase_name_length
    use fluid_module, only: fluid_type, create_fluid_vector
    use cell_module, only: cell_type

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    character(:), allocatable :: json_str
    PetscMPIInt :: rank
    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    type(mesh_type) :: mesh
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    type(logfile_type) :: logfile
    Vec :: y, fluid_vector
    PetscInt :: y_range_start, fluid_range_start, ghost
    PetscInt :: start_cell, end_cell, end_interior_cell, c
    PetscInt :: cell_geom_offset, fluid_offset, expected_region
    PetscReal :: t
    type(fluid_type) :: fluid
    type(cell_type) :: cell
    DMLabel :: ghost_label
    PetscSection :: cell_geom_section, fluid_section
    PetscReal, pointer, contiguous :: fluid_array(:), cell_geom_array(:)
    PetscErrorCode :: ierr, err

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)

    json_str = &
         '{"mesh": {"filename": "' // trim(adjustl(data_path)) // 'mesh/col10.exo"},' // &
         ' "eos": {"name": "we"}, ' // &
         '   "initial": {"primary": [' // &
         '                          [ 1e5, 20], [ 1e5, 20], [ 1e5, 20], ' // &
         '                          [ 1e5, 120], [ 1e5, 120], [ 1e5, 120], ' // &
         '                          [ 1e5, 0.4], [ 1e5, 0.4], [ 1e5, 0.4], [ 1e5, 0.4]], ' // &
         '               "region":  [1, 1, 1, 2, 2, 2, 4, 4, 4, 4]}}'

    json => fson_parse_mpi(str = json_str)
    call logfile%init('', echo = PETSC_FALSE)

    call thermo%init()
    call eos%init(json, thermo)
    call mesh%init(eos, json)
    call mesh%configure(gravity, json, err = err)

    call DMCreateGlobalVector(mesh%dm, y, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(y, "primary", ierr); CHKERRQ(ierr)
    call global_vec_range_start(y, y_range_start)

    call create_fluid_vector(mesh%dm, max_component_name_length, &
         eos%component_names, max_phase_name_length, &
         eos%phase_names, fluid_vector, fluid_range_start)

    t = 0._dp
    call setup_initial(json, mesh, eos, t, y, fluid_vector, &
         y_range_start, fluid_range_start, logfile)

    call fson_destroy_mpi(json)
    call mesh%destroy_distribution_data()

    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    end_interior_cell = dm_get_end_interior_cell(mesh%dm, end_cell)
    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)
    call fluid%init(eos%num_components, eos%num_phases)

    call local_vec_section(mesh%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)
    call cell%init(eos%num_components, eos%num_phases)
    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    do c = start_cell, end_interior_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          fluid_offset = global_section_offset(fluid_section, c, &
               fluid_range_start)
          call fluid%assign(fluid_array, fluid_offset)
          cell_geom_offset = section_offset(cell_geom_section, c)
          call cell%assign_geometry(cell_geom_array, cell_geom_offset)
          associate(z => cell%centroid(3))
            if (z > -300._dp) then
               expected_region = 1
            else if (z > -600._dp) then
               expected_region = 2
            else
               expected_region = 4
            end if
          end associate
        call test%assert(expected_region, nint(fluid%region), &
             'region')
     end if
  end do

  call fluid%destroy()
  call cell%destroy()
  call VecRestoreArrayReadF90(mesh%cell_geom, cell_geom_array, ierr)
  CHKERRQ(ierr)
  call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
  call VecDestroy(fluid_vector, ierr); CHKERRQ(ierr)
  call VecDestroy(y, ierr); CHKERRQ(ierr)
  call mesh%destroy()
  call eos%destroy()
  call thermo%destroy()
  call logfile%destroy()

  end subroutine test_initial_region

!------------------------------------------------------------------------

end module initial_test
