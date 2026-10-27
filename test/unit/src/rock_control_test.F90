module rock_control_test

  ! Test for rock control module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use rock_module
  use control_module
  use rock_control_module
  use rock_setup_module

  implicit none
  private

  character(len = 512) :: data_path

  public :: setup, teardown
  public :: test_rock_control_table
  public :: test_rock_control_temperature_dependent_permeability
  public :: test_rock_control_temperature_dependent_specific_heat

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

  subroutine test_rock_control_table(test)
    ! Table rock control

    use fson
    use fson_mpi_module
    use mesh_module
    use list_module
    use IAPWS_module
    use eos_wge_module
    use dm_utils_module, only: global_vec_section, global_section_offset, &
         global_to_local_vec_section, restore_dm_local_vec

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(IAPWS_type) :: thermo
    type(eos_wge_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(list_type) :: rock_controls
    type(rock_type) :: rock
    Vec :: rock_vector
    PetscReal, pointer, contiguous :: rock_array(:)
    PetscSection :: rock_section
    PetscInt :: rock_range_start, num_rock_types
    PetscErrorCode :: ierr, err
    PetscMPIInt :: rank
    DMLabel :: ghost_label
    PetscReal :: time
    PetscInt, parameter :: expected_num_rocks = 5
    character(12), parameter :: rock_names(expected_num_rocks) = &
         ["constant    ", "scalar      ", "array       ", &
         "value scalar", "value array "]
    PetscReal, parameter :: start_time = 0._dp, test_time = 4500._dp
    PetscReal, parameter :: expected_initial_permeability(expected_num_rocks, 3) = reshape([ &
         1.e-13_dp, 1.e-13_dp, 1.e-14_dp, 1.e-13_dp, 1.e-14_dp, &
         1.e-13_dp, 1.e-13_dp, 2.e-14_dp, 1.e-13_dp, 2.e-14_dp, &
         1.e-13_dp, 1.e-13_dp, 3.e-14_dp, 1.e-13_dp, 3.e-14_dp], &
         [expected_num_rocks, 3])
    PetscReal, parameter :: expected_permeability(expected_num_rocks, 3) = reshape([ &
         1.e-13_dp, 7.e-14_dp, 4.e-15_dp, 7.e-14_dp, 4.e-15_dp, &
         1.e-13_dp, 7.e-14_dp, 5.e-15_dp, 7.e-14_dp, 5.e-15_dp, &
         1.e-13_dp, 7.e-14_dp, 6.e-15_dp, 7.e-14_dp, 6.e-15_dp], &
         [expected_num_rocks, 3])
    PetscReal, parameter :: expected_initial_porosity(expected_num_rocks) = &
         [0.1_dp, 0.1_dp, 0.2_dp, 0.1_dp, 0.2_dp]
    PetscReal, parameter :: expected_porosity(expected_num_rocks) = &
         [0.1_dp, 0.05_dp, 0.2_dp, 0.05_dp, 0.2_dp]
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscReal, parameter :: ktol = 1.e-20_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    json => fson_parse_mpi(trim(adjustl(data_path)) // "rock/test_rock_controls_table.json")

    call thermo%init()
    call eos%init(json, thermo)
    call mesh%init(eos, json)
    call DMCreateLabel(mesh%serial_dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(gravity, json, err = err)

    call setup_rocks(json, mesh%dm, rock_vector, mesh%rock_types, &
         rock_controls, rock_range_start, mesh%ghost_cell, eos, err = err)
    call test%assert(0, err, "rock setup error")
    num_rock_types = mesh%rock_types%count()
    if (rank == 0) then
       call test%assert(expected_num_rocks, num_rock_types, "number of rocks")
    end if

    call VecGetArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
    call global_vec_section(rock_vector, rock_section)
    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    call rock%init()

    time = start_time
    call rock_controls%traverse(rock_control_iterator)
    call rock_test("t = 0", expected_initial_permeability, expected_initial_porosity)

    time = test_time
    call rock_controls%traverse(rock_control_iterator)
    call rock_test("t > 0", expected_permeability, expected_porosity)

    call rock%destroy()
    call rock_controls%destroy(rock_control_list_node_data_destroy, &
         reverse = PETSC_TRUE)
    call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
    call VecDestroy(rock_vector, ierr); CHKERRQ(ierr)
    call mesh%destroy_distribution_data()
    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  contains

    subroutine rock_control_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      select type (rock_control => node%data)
      class is (table_vector_control_type)
         call rock_control%update(time, rock_array, rock_section, rock_range_start)
      end select
      stopped = PETSC_FALSE
    end subroutine rock_control_iterator

     subroutine rock_control_list_node_data_destroy(node)
      type(list_node_type), pointer, intent(in out) :: node
      select type (rock_control => node%data)
      class is (table_vector_control_type)
         call rock_control%destroy()
      end select
    end subroutine rock_control_list_node_data_destroy

    subroutine rock_test(title, permeability, porosity)

      character(*), intent(in) :: title
      PetscReal, intent(in) :: permeability(expected_num_rocks, 3)
      PetscReal, intent(in) :: porosity(expected_num_rocks)
      ! Locals:
      PetscInt :: ir, num_cells, c, ic
      PetscInt :: rock_offset, ghost
      type(list_node_type), pointer :: node
      IS :: cell_IS
      PetscInt, pointer, contiguous :: cells(:)

      do ir = 1, num_rock_types
         node => mesh%rock_types%get(rock_names(ir))
         select type (item => node%data)
         type is (rock_dict_item_type)
            call DMGetStratumSize(mesh%dm, rock_type_label_name, item%label_value, &
                 num_cells, ierr); CHKERRQ(ierr)
            if (num_cells > 0) then
               call DMGetStratumIS(mesh%dm, rock_type_label_name, item%label_value, &
                    cell_IS, ierr); CHKERRQ(ierr)
               call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
               do ic = 1, num_cells
                  c = cells(ic)
                  call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
                  if (ghost < 0) then
                     rock_offset = global_section_offset(rock_section, c, &
                          rock_range_start)
                     call rock%assign(rock_array, rock_offset)
                     call test%assert(permeability(ir, :), rock%permeability, &
                          trim(rock_names(ir)) // " permeability " // &
                          trim(title), tol = ktol)
                     call test%assert(porosity(ir), rock%porosity, &
                          trim(rock_names(ir)) // " porosity " // trim(title))
                  end if
               end do
               call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
            end if
         end select
      end do

    end subroutine rock_test

  end subroutine test_rock_control_table

!------------------------------------------------------------------------

  subroutine test_rock_control_temperature_dependent_permeability(test)
    ! Temperature-dependent permeability rock control

    use fson
    use fson_mpi_module
    use mesh_module
    use list_module
    use IAPWS_module
    use eos_we_module
    use dm_utils_module, only: global_vec_section, global_section_offset, &
         global_to_local_vec_section, restore_dm_local_vec
    use fluid_module, only: create_fluid_vector, fluid_type
    use eos_module, only: max_component_name_length, max_phase_name_length

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(list_type) :: rock_controls
    type(rock_type) :: rock
    type(fluid_type) :: fluid
    Vec :: rock_vector, fluid_vector
    PetscReal, pointer, contiguous :: rock_array(:), fluid_array(:)
    PetscSection :: rock_section, fluid_section
    PetscInt :: rock_range_start, num_rock_types, fluid_range_start, ir, it
    PetscInt :: start_cell, end_cell
    PetscErrorCode :: ierr, err
    PetscMPIInt :: rank
    DMLabel :: ghost_label
    PetscInt, parameter :: dim = 2, num_rocks = 4, num_temperatures = 5
    character(12), parameter :: rock_names(num_rocks) = &
         ["HI scalar   ", "HI array    ", "table scalar", "table array "]
    character(16) :: title(num_temperatures) = &
         ["100 deg C", "370 deg C", "400 deg C", "450 deg C", "600 deg C"]
    PetscReal, parameter :: temperature(num_temperatures) = &
         [100._dp, 370._dp, 400._dp, 450._dp, 600._dp]
    PetscReal, parameter :: expected_log_permeability(dim, num_temperatures, &
         num_rocks) = reshape([ &
         -15._dp, -15._dp, & ! HI scalar
         -15.5_dp, -15.5_dp, &
         -17._dp, -17._dp, &
         -19.5_dp, -19.5_dp, &
         -22._dp, -22._dp, &
         -14._dp, -15._dp, & ! HI array
         -14.75_dp, -15.5_dp, &
         -17._dp, -17._dp, &
         -19.5_dp, -19.5_dp, &
         -22._dp, -22._dp, &
         -15._dp, -15._dp, & ! table scalar
         -16.4_dp, -16.4_dp, &
         -18.5_dp, -18.5_dp, &
         -22._dp, -22._dp, &
         -22._dp, -22._dp, &
         -15._dp, -16._dp, & ! table array
         -16.2_dp, -17.2_dp, &
         -18._dp, -19._dp, &
         -21._dp, -22._dp, &
         -21._dp, -22._dp], &
         [dim, num_temperatures, num_rocks])
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    json => fson_parse_mpi(trim(adjustl(data_path)) // &
         "rock/test_rock_controls_temperature.json")

    call thermo%init()
    call eos%init(json, thermo)
    call mesh%init(eos, json)
    call DMCreateLabel(mesh%serial_dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(gravity, json, err = err)

    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, &
         ierr); CHKERRQ(ierr)

    call setup_rocks(json, mesh%dm, rock_vector, mesh%rock_types, &
         rock_controls, rock_range_start, mesh%ghost_cell, eos, err = err)
    call test%assert(0, err, "rock setup error")
    num_rock_types = mesh%rock_types%count()
    if (rank == 0) then
       call test%assert(num_rocks, num_rock_types, "number of rocks")
    end if

    call VecGetArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
    call global_vec_section(rock_vector, rock_section)
    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    call rock%init()

    call create_fluid_vector(mesh%dm, max_component_name_length, &
         eos%component_names, max_phase_name_length, &
         eos%phase_names, fluid_vector, fluid_range_start)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call global_vec_section(fluid_vector, fluid_section)
    call fluid%init(eos%num_components, eos%num_phases)

    do ir = 1, num_rocks
       do it = 1, num_temperatures
          call rock_test(title(it), ir, temperature(it), &
               expected_log_permeability(:, it, ir))
       end do
    end do

    call fluid%destroy()
    call rock%destroy()
    call rock_controls%destroy(rock_control_list_node_data_destroy, &
         reverse = PETSC_TRUE)
    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
    call VecDestroy(fluid_vector, ierr); CHKERRQ(ierr)
    call VecDestroy(rock_vector, ierr); CHKERRQ(ierr)
    call mesh%destroy_distribution_data()
    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  contains

    subroutine rock_control_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      select type (rock_control => node%data)
      class is (fluid_rock_control_type)
         call rock_control%update(fluid_array, fluid_section, fluid_range_start, &
              rock_array, rock_section, rock_range_start)
      end select
      stopped = PETSC_FALSE
    end subroutine rock_control_iterator

     subroutine rock_control_list_node_data_destroy(node)
      type(list_node_type), pointer, intent(in out) :: node
      select type (rock_control => node%data)
      class is (fluid_rock_control_type)
         call rock_control%destroy()
      end select
    end subroutine rock_control_list_node_data_destroy

    subroutine rock_test(title, ir, temperature, log_permeability)

      character(*), intent(in) :: title
      PetscInt, intent(in) :: ir
      PetscReal, intent(in) :: temperature
      PetscReal, intent(in) :: log_permeability(dim)
      ! Locals:
      PetscInt :: num_cells, c, ic
      PetscInt :: rock_offset, ghost, fluid_offset
      type(list_node_type), pointer :: node
      IS :: cell_IS
      PetscInt, pointer, contiguous :: cells(:)

      do c = start_cell, end_cell - 1
         if (mesh%ghost_cell(c) < 0) then
            fluid_offset = global_section_offset(fluid_section, c, &
                 fluid_range_start)
            call fluid%assign(fluid_array, fluid_offset)
            fluid%temperature = temperature
         end if
      end do

      call rock_controls%traverse(rock_control_iterator)

      node => mesh%rock_types%get(rock_names(ir))
      select type (item => node%data)
      type is (rock_dict_item_type)
         call DMGetStratumSize(mesh%dm, rock_type_label_name, item%label_value, &
              num_cells, ierr); CHKERRQ(ierr)
         if (num_cells > 0) then
            call DMGetStratumIS(mesh%dm, rock_type_label_name, item%label_value, &
                 cell_IS, ierr); CHKERRQ(ierr)
            call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
            do ic = 1, num_cells
               c = cells(ic)
               call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
               if (ghost < 0) then
                  rock_offset = global_section_offset(rock_section, c, &
                       rock_range_start)
                  call rock%assign(rock_array, rock_offset)
                  call test%assert(log_permeability, log10(rock%permeability(1:dim)), &
                       trim(rock_names(ir)) // " log permeability " // &
                       trim(title))
               end if
            end do
            call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
         end if
      end select

    end subroutine rock_test

  end subroutine test_rock_control_temperature_dependent_permeability

!------------------------------------------------------------------------

  subroutine test_rock_control_temperature_dependent_specific_heat(test)
    ! Temperature-dependent specific heat rock control

    use fson
    use fson_mpi_module
    use mesh_module
    use list_module
    use IAPWS_module
    use eos_we_module
    use dm_utils_module, only: global_vec_section, global_section_offset, &
         global_to_local_vec_section, restore_dm_local_vec
    use fluid_module, only: create_fluid_vector, fluid_type
    use eos_module, only: max_component_name_length, max_phase_name_length

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(list_type) :: rock_controls
    type(rock_type) :: rock
    type(fluid_type) :: fluid
    Vec :: rock_vector, fluid_vector
    PetscReal, pointer, contiguous :: rock_array(:), fluid_array(:)
    PetscSection :: rock_section, fluid_section
    PetscInt :: rock_range_start, fluid_range_start, it
    PetscInt :: start_cell, end_cell
    PetscErrorCode :: ierr, err
    PetscMPIInt :: rank
    DMLabel :: ghost_label
    PetscInt, parameter :: num_temperatures = 3
    character(13), parameter :: rock_name = 'specific_heat'
    character(16) :: title(num_temperatures) = &
         ["100 deg C", "500 deg C", "790 deg C"]
    PetscReal, parameter :: temperature(num_temperatures) = &
         [100._dp, 500._dp, 790._dp]
    PetscReal, parameter :: expected_specific_heat(num_temperatures) = [ &
         1.e3_dp, 1.e3_dp, 2.e3_dp]
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    json => fson_parse_mpi(trim(adjustl(data_path)) // &
         "rock/test_rock_controls_specific_heat.json")

    call thermo%init()
    call eos%init(json, thermo)
    call mesh%init(eos, json)
    call DMCreateLabel(mesh%serial_dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(gravity, json, err = err)

    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, &
         ierr); CHKERRQ(ierr)

    call setup_rocks(json, mesh%dm, rock_vector, mesh%rock_types, &
         rock_controls, rock_range_start, mesh%ghost_cell, eos, err = err)
    call test%assert(0, err, "rock setup error")

    call VecGetArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
    call global_vec_section(rock_vector, rock_section)
    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    call rock%init()

    call create_fluid_vector(mesh%dm, max_component_name_length, &
         eos%component_names, max_phase_name_length, &
         eos%phase_names, fluid_vector, fluid_range_start)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call global_vec_section(fluid_vector, fluid_section)
    call fluid%init(eos%num_components, eos%num_phases)

    do it = 1, num_temperatures
       call rock_test(title(it), temperature(it), expected_specific_heat(it))
    end do

    call fluid%destroy()
    call rock%destroy()
    call rock_controls%destroy(rock_control_list_node_data_destroy, &
         reverse = PETSC_TRUE)
    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
    call VecDestroy(fluid_vector, ierr); CHKERRQ(ierr)
    call VecDestroy(rock_vector, ierr); CHKERRQ(ierr)
    call mesh%destroy_distribution_data()
    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  contains

    subroutine rock_control_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      select type (rock_control => node%data)
      class is (fluid_rock_control_type)
         call rock_control%update(fluid_array, fluid_section, fluid_range_start, &
              rock_array, rock_section, rock_range_start)
      end select
      stopped = PETSC_FALSE
    end subroutine rock_control_iterator

     subroutine rock_control_list_node_data_destroy(node)
      type(list_node_type), pointer, intent(in out) :: node
      select type (rock_control => node%data)
      class is (fluid_rock_control_type)
         call rock_control%destroy()
      end select
    end subroutine rock_control_list_node_data_destroy

    subroutine rock_test(title, temperature, specific_heat)

      character(*), intent(in) :: title
      PetscReal, intent(in) :: temperature
      PetscReal, intent(in) :: specific_heat
      ! Locals:
      PetscInt :: num_cells, c, ic
      PetscInt :: rock_offset, ghost, fluid_offset
      type(list_node_type), pointer :: node
      IS :: cell_IS
      PetscInt, pointer, contiguous :: cells(:)

      do c = start_cell, end_cell - 1
         if (mesh%ghost_cell(c) < 0) then
            fluid_offset = global_section_offset(fluid_section, c, &
                 fluid_range_start)
            call fluid%assign(fluid_array, fluid_offset)
            fluid%temperature = temperature
         end if
      end do

      call rock_controls%traverse(rock_control_iterator)

      node => mesh%rock_types%get(rock_name)
      select type (item => node%data)
      type is (rock_dict_item_type)
         call DMGetStratumSize(mesh%dm, rock_type_label_name, item%label_value, &
              num_cells, ierr); CHKERRQ(ierr)
         if (num_cells > 0) then
            call DMGetStratumIS(mesh%dm, rock_type_label_name, item%label_value, &
                 cell_IS, ierr); CHKERRQ(ierr)
            call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
            do ic = 1, num_cells
               c = cells(ic)
               call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
               if (ghost < 0) then
                  rock_offset = global_section_offset(rock_section, c, &
                       rock_range_start)
                  call rock%assign(rock_array, rock_offset)
                  call test%assert(specific_heat, rock%specific_heat, &
                       " specific heat " // trim(title))
               end if
            end do
            call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
         end if
      end select

    end subroutine rock_test

  end subroutine test_rock_control_temperature_dependent_specific_heat

end module rock_control_test
