module source_control_test

  ! Test for source control module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use source_module
  use source_control_module
  use source_setup_module

  implicit none
  private

  character(len = 512) :: data_path

  public :: setup, teardown
  public :: test_source_control_table, test_source_control_pressure_reference

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

  subroutine test_source_control_table(test)
    ! Table source control

    use fson
    use fson_mpi_module
    use mesh_module
    use list_module
    use IAPWS_module
    use eos_module, only: max_component_name_length, max_phase_name_length
    use eos_wge_module
    use fluid_module, only: create_fluid_vector
    use dm_utils_module, only: global_vec_section, global_section_offset, &
         global_to_local_vec_section, restore_dm_local_vec

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(IAPWS_type) :: thermo
    type(eos_wge_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(list_type) :: source_controls
    type(source_type) :: source
    Vec :: fluid_vector, local_fluid_vector, source_vector
    PetscReal, pointer, contiguous :: fluid_array(:), local_fluid_array(:)
    PetscReal, pointer, contiguous :: source_array(:)
    PetscSection :: fluid_section, local_fluid_section, source_section
    PetscInt :: num_sources, total_num_sources, source_vector_size
    PetscInt :: fluid_range_start, source_range_start
    PetscInt :: s, source_offset, source_index
    PetscReal :: t, interval(2)
    character(len = 8) :: srcstr
    PetscErrorCode :: ierr, err
    PetscReal, parameter :: start_time = 0._dp
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscInt, parameter :: expected_num_sources = 8
    PetscMPIInt :: rank
    IS :: source_is

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    json => fson_parse_mpi(trim(adjustl(data_path)) // "source/test_source_controls_table.json")

    call thermo%init()
    call eos%init(json, thermo)
    call mesh%init(eos, json)
    call DMCreateLabel(mesh%serial_dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(gravity, json, err = err)

    call create_fluid_vector(mesh%dm, max_component_name_length, &
         eos%component_names, max_phase_name_length, eos%phase_names, &
         fluid_vector, fluid_range_start)
    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call setup_sources(json, mesh%dm, mesh%cell_natural_global, eos, thermo, start_time, &
         fluid_vector, fluid_range_start, source_vector, source_range_start, &
         num_sources, total_num_sources, source_controls, source_is, err = err)
    call test%assert(0, err, "source setup error")
    call source%init(eos)
    call test%assert(13, source%dof, "source dof")

    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call VecGetSize(source_vector, source_vector_size, ierr); CHKERRQ(ierr)
    if (rank == 0) then
       call test%assert(expected_num_sources, total_num_sources, &
            "number of sources")
       call test%assert(expected_num_sources * source%dof, &
            source_vector_size, "source vector size")
    end if

    call global_to_local_vec_section(fluid_vector, local_fluid_vector, &
         local_fluid_section)
    call VecGetArrayReadF90(local_fluid_vector, local_fluid_array, ierr)
    CHKERRQ(ierr)
    call global_vec_section(source_vector, source_section)
    call VecGetArrayReadF90(source_vector, source_array, ierr); CHKERRQ(ierr)

    t = 120._dp
    interval = [30._dp, t]
    call source_controls%traverse(source_control_iterator)
    do s = 0, num_sources - 1
       source_offset = global_section_offset(source_section, s, &
            source_range_start)
       call source%assign(source_array, source_offset)
       source_index = nint(source%source_index)
       write(srcstr, '(a, i1)') 'source ', source_index
       select case (source_index)
       case (0)
          call test%assert(-2.25_dp, source%rate, trim(srcstr))
       case (1)
          call test%assert(2.25_dp, source%rate, trim(srcstr))
       case (2)
          call test%assert(104.5e3_dp, source%injection_enthalpy, &
               trim(srcstr))
       case (3)
          call test%assert(-2.5_dp * 0.75_dp, source%rate, trim(srcstr))
       case (4)
          call test%assert(-2.5_dp * 0.5_dp, source%rate, trim(srcstr))
       end select
    end do

    call VecRestoreArrayReadF90(local_fluid_vector, local_fluid_array, ierr)
    CHKERRQ(ierr)
    call restore_dm_local_vec(local_fluid_vector)

    call ISDestroy(source_is, ierr); CHKERRQ(ierr)
    call source%destroy()
    call source_controls%destroy(source_control_list_node_data_destroy, &
         reverse = PETSC_TRUE)
    call VecRestoreArrayReadF90(source_vector, source_array, ierr); CHKERRQ(ierr)
    call VecDestroy(source_vector, ierr); CHKERRQ(ierr)
    call VecDestroy(fluid_vector, ierr); CHKERRQ(ierr)
    call mesh%destroy_distribution_data()
    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  contains

    subroutine source_control_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      select type (source_control => node%data)
      class is (source_control_type)
         call source_control%update(t, interval, source_array, &
              source_section, source_range_start, fluid_array, &
              fluid_section, eos)
      end select
      stopped = PETSC_FALSE
    end subroutine source_control_iterator

     subroutine source_control_list_node_data_destroy(node)
      type(list_node_type), pointer, intent(in out) :: node
      select type (source_control => node%data)
      class is (source_control_type)
         call source_control%destroy()
      end select
    end subroutine source_control_list_node_data_destroy

  end subroutine test_source_control_table

!------------------------------------------------------------------------

  subroutine test_source_control_pressure_reference(test)
    ! Pressure reference source controls

    use fson
    use fson_mpi_module
    use mesh_module
    use list_module
    use IAPWS_module
    use eos_module, only: max_component_name_length, max_phase_name_length
    use eos_wge_module
    use fluid_module, only: fluid_type, create_fluid_vector
    use dm_utils_module, only: global_vec_section, global_section_offset, &
         global_to_local_vec_section, restore_dm_local_vec, section_offset
    use interpolation_module, only: INTERP_STEP, INTERP_AVERAGING_ENDPOINT

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(IAPWS_type) :: thermo
    type(eos_wge_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(list_type) :: source_controls
    type(source_type) :: source
    Vec :: source_vector
    Vec :: fluid_vector, local_fluid_vector
    PetscReal, pointer, contiguous :: fluid_array(:), local_fluid_array(:)
    PetscReal, pointer, contiguous :: source_array(:)
    PetscSection :: source_section, fluid_section, local_fluid_section
    type(fluid_type) :: fluid
    PetscInt :: num_sources, total_num_sources, num_source_controls, source_index
    PetscInt :: start_cell, end_cell, c, s, s12, source_range_start, fluid_range_start
    PetscInt :: fluid_offset, source_offset, cell_phase_composition
    PetscReal :: t, interval(2), props(2)
    PetscReal :: cell_temperature, cell_liquid_density, cell_liquid_internal_energy
    PetscReal :: cell_vapour_density, cell_vapour_internal_energy
    PetscReal :: cell_liquid_viscosity, cell_vapour_viscosity
    IS :: source_is
    PetscErrorCode :: ierr, err
    PetscReal, parameter :: cell_pressure = 50.e5_dp, cell_vapour_saturation = 0.8_dp
    PetscInt, parameter :: cell_region = 4
    PetscReal, parameter :: start_time = 0._dp
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscReal, parameter :: expected_rates(0: 14) = [ &
         -12.8728519749_dp, -10._dp, -9.3081349399_dp, -11._dp, -10.1910078135_dp, &
         0.0_dp, 130.0_dp, -13.0_dp, 0.0_dp, -12.9264888581701_dp, &
         -10.3366086953508_dp, 0._dp, -2.25_dp, -12.8728519749_dp * 0.75_dp, &
         -12.8728519749_dp * 0.375_dp]
    PetscMPIInt :: rank
    character(len = 16) :: srcstr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    json => fson_parse_mpi(trim(adjustl(data_path)) // &
         "source/test_source_controls_pressure_reference.json")

    call thermo%init()
    call eos%init(json, thermo)
    call thermo%saturation%temperature(cell_pressure, cell_temperature, ierr)
    cell_phase_composition = thermo%phase_composition(cell_region, &
         cell_pressure, cell_temperature)
    call thermo%water%properties([cell_pressure, cell_temperature], &
         props, ierr)
    cell_liquid_density = props(1)
    cell_liquid_internal_energy = props(2)
    call thermo%water%viscosity(cell_temperature, cell_pressure, cell_liquid_density, &
         cell_liquid_viscosity)
    call thermo%steam%properties([cell_pressure, cell_temperature], &
         props, ierr)
    cell_vapour_density = props(1)
    cell_vapour_internal_energy = props(2)
    call thermo%steam%viscosity(cell_temperature, cell_pressure, cell_vapour_density, &
         cell_vapour_viscosity)

    call mesh%init(eos, json)
    call fluid%init(eos%num_components, eos%num_phases)
    call DMCreateLabel(mesh%serial_dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(gravity, json, err = err)

    call create_fluid_vector(mesh%dm, max_component_name_length, &
         eos%component_names, max_phase_name_length, eos%phase_names, &
         fluid_vector, fluid_range_start)
    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    do c = start_cell, end_cell - 1
       if (mesh%ghost_cell(c) < 0) then
          fluid_offset = global_section_offset(fluid_section, c, fluid_range_start)
          call fluid%assign(fluid_array, fluid_offset)
          fluid%pressure = cell_pressure
          fluid%temperature = cell_temperature
          fluid%region = dble(cell_region)
          fluid%phase_composition = cell_phase_composition
          fluid%phase(1)%density = cell_liquid_density
          fluid%phase(1)%viscosity = cell_liquid_viscosity
          fluid%phase(1)%saturation = 1._dp - cell_vapour_saturation
          fluid%phase(1)%relative_permeability = 1._dp - cell_vapour_saturation
          fluid%phase(1)%specific_enthalpy = cell_liquid_internal_energy + &
               cell_pressure / cell_liquid_density
          fluid%phase(1)%internal_energy = cell_liquid_internal_energy
          fluid%phase(1)%mass_fraction = [0.75_dp, 0.25_dp]
          fluid%phase(2)%density = cell_vapour_density
          fluid%phase(2)%viscosity = cell_vapour_viscosity
          fluid%phase(2)%saturation = cell_vapour_saturation
          fluid%phase(2)%relative_permeability = cell_vapour_saturation
          fluid%phase(2)%specific_enthalpy = cell_vapour_internal_energy + &
               cell_pressure / cell_vapour_density
          fluid%phase(2)%internal_energy = cell_vapour_internal_energy
          fluid%phase(2)%mass_fraction = [0.9_dp, 0.1_dp]
       end if
    end do
    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call setup_sources(json, mesh%dm, mesh%cell_natural_global, eos, thermo, start_time, &
         fluid_vector, fluid_range_start, source_vector, source_range_start, &
         num_sources, total_num_sources, source_controls, source_is, err = err)
    call test%assert(0, err, "source setup error")
    call source%init(eos)

    if (rank == 0) then
      call test%assert(15, total_num_sources, "number of sources")
    end if

    call MPI_reduce(source_controls%count, num_source_controls, 1, &
         MPI_INTEGER, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
      call test%assert(28, num_source_controls, "number of source controls")
    end if

    call global_to_local_vec_section(fluid_vector, local_fluid_vector, &
         local_fluid_section)
    call VecGetArrayF90(local_fluid_vector, local_fluid_array, ierr)
    CHKERRQ(ierr)

    call global_vec_section(source_vector, source_section)
    call VecGetArrayF90(source_vector, source_array, ierr); CHKERRQ(ierr)

    t = 120._dp
    interval = [30._dp, t]

    call source_controls%traverse(source_control_update_iterator)
    call source_controls%traverse(source_control_test_iterator)

    s12 = -1
    do s = 0, num_sources - 1
       source_offset = global_section_offset(source_section, s, source_range_start)
       call source%assign(source_array, source_offset)
       source_index = nint(source%source_index)
       write(srcstr, '(a, i1, a)') 'source ', source_index + 1, ' rate'
       call test%assert(expected_rates(source_index), source%rate, &
            trim(srcstr))
       if (source_index == 12) s12 = s
    end do

    ! Test deliverability threshold control- reduce fluid pressure:
    call reset_fluid_pressures(6.e5_dp)
    call source_controls%traverse(source_control_update_iterator)
    if (s12 >= 0) then
       source_offset = global_section_offset(source_section, s12, source_range_start)
       call source%assign(source_array, source_offset)
       call test%assert(-2.25_dp, source%rate, "source 13 rate P = 6 bar")
    end if

    call reset_fluid_pressures(4.e5_dp)
    call source_controls%traverse(source_control_update_iterator)
    if (s12 >= 0) then
       call test%assert(-1.125_dp, source%rate, "source 13 rate P = 4 bar")
    end if
    call reset_fluid_pressures(3.e5_dp)
    call source_controls%traverse(source_control_update_iterator)
    if (s12 >= 0) then
       call test%assert(-0.5625_dp, source%rate, "source 13 rate P = 3 bar")
    end if

    call VecRestoreArrayF90(local_fluid_vector, local_fluid_array, ierr)
    CHKERRQ(ierr)
    call restore_dm_local_vec(local_fluid_vector)

    call ISDestroy(source_is, ierr); CHKERRQ(ierr)
    call source%destroy()
    call source_controls%destroy(source_control_list_node_data_destroy, &
      reverse = PETSC_TRUE)
    call VecRestoreArrayF90(source_vector, source_array, ierr); CHKERRQ(ierr)
    call VecDestroy(source_vector, ierr); CHKERRQ(ierr)
    call fluid%destroy()
    call VecDestroy(fluid_vector, ierr); CHKERRQ(ierr)
    call mesh%destroy_distribution_data()
    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  contains

    subroutine source_control_update_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      select type (source_control => node%data)
      class is (source_control_type)
         call source_control%update(t, interval, source_array, &
              source_section, source_range_start, local_fluid_array, &
              local_fluid_section, eos)
      end select
      stopped = PETSC_FALSE
    end subroutine source_control_update_iterator

    subroutine source_control_test_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      PetscInt :: s, source_index

      select type (source_control => node%data)

      type is (source_control_deliverability_type)
         s = source_control%local_source_index
         source_offset = global_section_offset(source_section, s, &
              source_range_start)
         call source%assign(source_array, source_offset)
         source_index = nint(source%source_index)
         select case (source_index)
         case (1)
            call test%assert(1.e-12_dp, &
                 source_control%productivity%val(1,1), &
                 "source 2 productivity")
            call test%assert(2.e5_dp, &
                 source_control%reference_pressure%val(1,1), &
                 "source 2 reference pressure")
         case (3)
            call test%assert(8.54511496085953E-13_dp, &
                 source_control%productivity%val(1,1), &
                 "source 4 productivity")
         case (9)
            call test%assert(SRC_PRESSURE_TABLE_COORD_TIME, &
                 source_control%pressure_table_coordinate, &
                 "source 10 pressure table coordinate")
            call test%assert(4, &
                 source_control%reference_pressure%coord%size, &
                 "source 10 pressure table size")
            call test%assert(INTERP_STEP, &
                 source_control%reference_pressure%interpolation_type, &
                 "source 10 interpolation type")
            call test%assert(INTERP_AVERAGING_ENDPOINT, &
                 source_control%reference_pressure%averaging_type, &
                 "source 10 averaging type")
            call test%assert(1.8e5_dp, &
                 source_control%reference_pressure%average(interval, 1), &
                 "source 10 reference pressure")
         case (10)
            call test%assert(SRC_PRESSURE_TABLE_COORD_ENTHALPY, &
                 source_control%pressure_table_coordinate, &
                 "source 11 pressure table coordinate")
         case (12)
            call test%assert(5.e5_dp, source_control%threshold, &
                 "source 13 deliverability threshold")
         end select

      type is (source_control_recharge_type)
         s = source_control%local_source_index
         source_offset = global_section_offset(source_section, s, source_range_start)
         call source%assign(source_array, source_offset)
         source_index = nint(source%source_index)
         select case (source_index)
         case (6)
            call test%assert(1.3e-2_dp, &
                 source_control%coefficient%val(1, 1), &
                 "source 7 recharge coefficient")
            call test%assert(50.1e5_dp, &
                 source_control%reference_pressure%val(1, 1), &
                 "source 7 reference pressure")
         case (11)
            call test%assert(cell_pressure, &
                 source_control%reference_pressure%val(1, 1), &
                 "source 12 reference pressure")
         end select

      type is (source_control_limiter_type)
         select case (source_control%type)
         case (SRC_CONTROL_LIMITER_TYPE_TOTAL)
            call test%assert(10._dp, source_control%limit, "total limiter limit")
         case (SRC_CONTROL_LIMITER_TYPE_STEAM)
            call test%assert(5._dp, source_control%limit, "steam limiter limit")
         end select

      type is (source_control_separator_type)
         call test%assert(10.e5_dp, &
              source_control%separator_pressure, "separator pressure")
         call test%assert(0.5371645375_dp, &
              source_control%steam_fraction, "separator steam fraction")
         call test%assert(-5.9580123975_dp, &
              source_control%water_flow_rate, "separator water flow rate")
         call test%assert(-6.9148395774_dp, &
              source_control%steam_flow_rate, "separator steam flow rate")
      end select
      stopped = PETSC_FALSE

    end subroutine source_control_test_iterator

    subroutine reset_fluid_pressures(P)
      PetscReal, intent(in) :: P
      do c = start_cell, end_cell - 1
         if (mesh%ghost_cell(c) < 0) then
            fluid_offset = section_offset(local_fluid_section, c)
            call fluid%assign(local_fluid_array, fluid_offset)
            fluid%pressure = P
         end if
      end do

    end subroutine reset_fluid_pressures

    subroutine source_list_node_data_destroy(node)
      type(list_node_type), pointer, intent(in out) :: node
      select type (source => node%data)
      type is (source_type)
         call source%destroy()
      end select
    end subroutine source_list_node_data_destroy

     subroutine source_control_list_node_data_destroy(node)
      type(list_node_type), pointer, intent(in out) :: node
      select type (source_control => node%data)
      class is (source_control_type)
         call source_control%destroy()
      end select
    end subroutine source_control_list_node_data_destroy

  end subroutine test_source_control_pressure_reference

!------------------------------------------------------------------------

end module source_control_test
