module source_control_test

  ! Test for source control module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use source_module
  use control_module
  use source_control_module
  use source_network_control_module
  use source_network_group_module
  use source_network_module
  use source_setup_module
  use tracer_module

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
    use tracer_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(IAPWS_type) :: thermo
    type(eos_wge_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(source_type) :: source
    Vec :: fluid_vector, local_fluid_vector
    PetscReal, pointer, contiguous :: fluid_array(:), local_fluid_array(:)
    PetscReal, pointer, contiguous :: source_array(:)
    PetscSection :: fluid_section, local_fluid_section, source_section
    type(source_network_type) :: source_network
    PetscInt :: source_vector_size
    PetscInt :: fluid_range_start
    PetscReal :: t, interval(2)
    PetscErrorCode :: ierr, err
    PetscReal, parameter :: start_time = 0._dp
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    type(tracer_type), allocatable :: tracers(:)
    PetscBool, parameter :: rate_specified = PETSC_FALSE, &
         enthalpy_specified = PETSC_FALSE
    PetscReal, parameter :: specified_rate = -1._dp, specified_enthalpy = 0._dp
    PetscInt, parameter :: expected_num_sources = 11
    PetscMPIInt :: rank

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    json => fson_parse_mpi(trim(adjustl(data_path)) // "source/test_source_controls_table.json")

    call thermo%init()
    call eos%init(json, thermo)
    call setup_tracers(json, eos, tracers, err = err)
    call mesh%init(eos, json)
    call DMCreateLabel(mesh%serial_dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(gravity, json, err = err)

    call create_fluid_vector(mesh%dm, max_component_name_length, &
         eos%component_names, max_phase_name_length, eos%phase_names, &
         fluid_vector, fluid_range_start)
    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call setup_source_network(json, mesh%dm, mesh%cell_natural_global, eos, tracers, &
         thermo, start_time, fluid_vector, fluid_range_start, source_network, &
         err = err)
    call test%assert(0, err, "source network setup error")
    call source%init("", eos, 0, 0, 0._dp, 0, 0, rate_specified, &
         specified_rate, enthalpy_specified, specified_enthalpy, size(tracers))
    call test%assert(13 + size(tracers) * 2, source%dof, "source dof")
    call source%destroy()

    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call VecGetSize(source_network%source, source_vector_size, ierr); CHKERRQ(ierr)
    if (rank == 0) then
       call test%assert(expected_num_sources, source_network%num_sources, &
            "number of sources")
       call test%assert(0, source_network%num_groups, "number of source groups")
       call test%assert(expected_num_sources * source%dof, &
            source_vector_size, "source vector size")
    end if

    call global_to_local_vec_section(fluid_vector, local_fluid_vector, &
         local_fluid_section)
    call VecGetArrayReadF90(local_fluid_vector, local_fluid_array, ierr)
    CHKERRQ(ierr)
    call global_vec_section(source_network%source, source_section)
    call VecGetArrayReadF90(source_network%source, source_array, ierr); CHKERRQ(ierr)

    t = 120._dp
    interval = [30._dp, t]
    call source_network%source_controls%traverse(source_control_iterator)
    call source_network%sources%traverse(source_test_iterator)

    call VecRestoreArrayReadF90(local_fluid_vector, local_fluid_array, ierr)
    CHKERRQ(ierr)
    call restore_dm_local_vec(local_fluid_vector)

    call VecRestoreArrayReadF90(source_network%source, source_array, ierr); CHKERRQ(ierr)
    call source_network%destroy()
    call VecDestroy(fluid_vector, ierr); CHKERRQ(ierr)
    call mesh%destroy_distribution_data()
    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()
    deallocate(tracers)
    call fson_destroy_mpi(json)

  contains

    subroutine source_control_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      stopped = PETSC_FALSE
      select type (source_control => node%data)
      class is (integer_object_control_type)
         call source_control%update()
      class is (interval_update_object_control_type)
         call source_control%update(interval)
      class is (pressure_reference_source_control_type)
         call source_control%update(t, interval, local_fluid_array, &
              local_fluid_section)
      end select
    end subroutine source_control_iterator

    subroutine source_test_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: s, source_offset, source_index
      character(len = 8) :: srcstr

      stopped = PETSC_FALSE
      select type (source => node%data)
      type is (source_type)
         s = source%local_source_index
         source_offset = global_section_offset(source_section, s, &
              source_network%source_range_start)
         call source%assign(source_array, source_offset)
         source_index = nint(source%source_index)
         write(srcstr, '(a, i1)') 'source ', source_index
         select case (source_index)
         case (0)
            call test%assert(-2.25_dp, source%rate, trim(srcstr))
         case (1: 3)
            call test%assert(2.25_dp, source%rate, trim(srcstr))
         case (4: 5)
            call test%assert(104.5e3_dp, source%injection_enthalpy, &
                 trim(srcstr))
         case (6)
            call test%assert(-2.5_dp * 0.75_dp, source%rate, trim(srcstr))
         case (7)
            call test%assert(-2.5_dp * 0.5_dp, source%rate, trim(srcstr))
         case (8)
            call test%assert([0.001_dp / 3._dp, 0.001_dp / 3._dp], &
                 source%tracer_injection_rate, trim(srcstr))
         case (9)
            call test%assert([0.001_dp / 3._dp, 0.001_dp * 7._dp / 9._dp], &
                 source%tracer_injection_rate, trim(srcstr))
         case (10)
            call test%assert([0._dp, 0.001_dp * 7._dp / 9._dp], &
                 source%tracer_injection_rate, trim(srcstr))
         end select
      end select

    end subroutine source_test_iterator

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
    use separator_module, only: SEPARATED_FLOW_TYPE_TOTAL, &
         SEPARATED_FLOW_TYPE_STEAM

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(IAPWS_type) :: thermo
    type(eos_wge_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(source_type) :: source
    Vec :: fluid_vector, local_fluid_vector
    PetscReal, pointer, contiguous :: fluid_array(:), local_fluid_array(:)
    PetscReal, pointer, contiguous :: source_array(:)
    PetscSection :: source_section, fluid_section, local_fluid_section
    type(source_network_type) :: source_network
    type(fluid_type) :: fluid
    PetscInt :: num_separators, num_source_controls
    PetscInt :: start_cell, end_cell, c, s12, fluid_range_start
    PetscInt :: fluid_offset, source_offset, cell_phase_composition
    PetscReal :: t, interval(2), props(2)
    PetscReal :: cell_temperature, cell_liquid_density, cell_liquid_internal_energy
    PetscReal :: cell_vapour_density, cell_vapour_internal_energy
    PetscReal :: cell_liquid_viscosity, cell_vapour_viscosity
    PetscErrorCode :: ierr, err
    PetscReal, parameter :: cell_pressure = 50.e5_dp, cell_vapour_saturation = 0.8_dp
    PetscInt, parameter :: cell_region = 4
    PetscReal, parameter :: start_time = 0._dp
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    type(tracer_type) :: tracers(0)
    PetscBool, parameter :: rate_specified = PETSC_FALSE, &
         enthalpy_specified = PETSC_FALSE
    PetscReal, parameter :: specified_rate = -1._dp, specified_enthalpy = 0._dp
    PetscInt, parameter :: expected_num_sources = 17
    PetscReal, parameter :: expected_rates(0: expected_num_sources - 1) = [ &
         -12.8728519749_dp, -10._dp, -9.3081349399_dp, -11._dp, -10.1910078135_dp, &
         0.0_dp, 130.0_dp, -13.0_dp, 0.0_dp, -12.9264888581701_dp, &
         -10.3366086953508_dp, 0._dp, -2.25_dp, -12.8728519749_dp * 0.75_dp, &
         -12.8728519749_dp * 0.375_dp, -10._dp, -9.3081349399_dp]
    PetscReal, parameter :: expected_steam_rates(0: expected_num_sources - 1) = [ &
         0._dp, 0._dp, -5._dp, 0._dp, 0._dp, &
         0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
         0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
         -5.60996474758954_dp, -5._dp]
    PetscMPIInt :: rank

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
          fluid%permeability_factor = 1._dp
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

    call setup_source_network(json, mesh%dm, mesh%cell_natural_global, eos, tracers, &
         thermo, start_time, fluid_vector, fluid_range_start, source_network, &
         err = err)
    call test%assert(0, err, "source setup error")

    if (rank == 0) then
      call test%assert(expected_num_sources, source_network%num_sources, "number of sources")
      call test%assert(0, source_network%num_groups, "number of source groups")
    end if

    call MPI_reduce(source_network%source_controls%count, num_source_controls, 1, &
         MPI_INTEGER, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
       call test%assert(31, num_source_controls, "number of source controls")
    end if

    call MPI_reduce(source_network%separated_sources%count, num_separators, 1, &
         MPI_INTEGER, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
       call test%assert(3, num_separators, "number of separators")
    end if

    call global_to_local_vec_section(fluid_vector, local_fluid_vector, &
         local_fluid_section)
    call VecGetArrayF90(local_fluid_vector, local_fluid_array, ierr)
    CHKERRQ(ierr)

    call global_vec_section(source_network%source, source_section)
    call VecGetArrayF90(source_network%source, source_array, ierr); CHKERRQ(ierr)

    t = 120._dp
    interval = [30._dp, t]

    call source_network%separated_sources%traverse(source_separator_iterator)
    call source_network%source_controls%traverse(source_control_update_iterator)
    call source_network%source_controls%traverse(source_control_test_iterator)

    s12 = -1
    call source_network%sources%traverse(source_rate_test_iterator)

    ! Test deliverability threshold control- reduce fluid pressure:
    call source%init("source 12", eos, 0, 0, 0._dp, 0, 0, rate_specified, &
         specified_rate, enthalpy_specified, specified_enthalpy, size(tracers))
    call reset_fluid_pressures(6.e5_dp)
    call source_network%source_controls%traverse(source_control_update_iterator)
    if (s12 >= 0) then
       source_offset = global_section_offset(source_section, s12, &
            source_network%source_range_start)
       call source%assign(source_array, source_offset)
       call test%assert(-2.25_dp, source%rate, "source 13 rate P = 6 bar")
    end if

    call reset_fluid_pressures(4.e5_dp)
    call source_network%source_controls%traverse(source_control_update_iterator)
    if (s12 >= 0) then
       call test%assert(-1.125_dp, source%rate, "source 13 rate P = 4 bar")
    end if
    call reset_fluid_pressures(3.e5_dp)
    call source_network%source_controls%traverse(source_control_update_iterator)
    if (s12 >= 0) then
       call test%assert(-0.5625_dp, source%rate, "source 13 rate P = 3 bar")
    end if

    t = t + 90._dp
    interval = [150._dp, t]
    call source_network%source_controls%traverse(source_control_update_iterator)
    if (s12 >= 0) then
       call test%assert(-0.0291666666667_dp, source%rate, &
            "source 13 rate P = 3 bar case 2")
    end if

    call VecRestoreArrayF90(local_fluid_vector, local_fluid_array, ierr)
    CHKERRQ(ierr)
    call restore_dm_local_vec(local_fluid_vector)

    call source%destroy()
    call VecRestoreArrayF90(source_network%source, source_array, ierr); CHKERRQ(ierr)
    call source_network%destroy()
    call fluid%destroy()
    call VecDestroy(fluid_vector, ierr); CHKERRQ(ierr)
    call mesh%destroy_distribution_data()
    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  contains

    subroutine source_separator_iterator(node, stopped)
      !! Updates enthalpy and separated outputs from separators.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: s, source_offset
      PetscReal, allocatable :: phase_flow_fractions(:)

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)

         s = source%local_source_index
         source_offset = global_section_offset(source_section, &
              s, source_network%source_range_start)
         call source%assign(source_array, source_offset)

         call source%assign_fluid_local(local_fluid_array, local_fluid_section)
         allocate(phase_flow_fractions(source%fluid%num_phases))
         phase_flow_fractions = source%fluid%phase_flow_fractions()
         source%enthalpy = source%fluid%specific_enthalpy(phase_flow_fractions)
         deallocate(phase_flow_fractions)

         call source%get_separated_flows()

      end select

    end subroutine source_separator_iterator

    subroutine source_control_update_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      select type (source_control => node%data)
      class is (integer_object_control_type)
         call source_control%update()
      class is (interval_update_object_control_type)
         call source_control%update(interval)
      class is (pressure_reference_source_control_type)
         call source_control%update(t, interval, local_fluid_array, &
              local_fluid_section)
      end select
      stopped = PETSC_FALSE
    end subroutine source_control_update_iterator

    subroutine source_control_test_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      PetscInt :: s, source_index

      select type (source_control => node%data)

      type is (deliverability_source_control_type)
         select type (source => source_control%objects%head%data)
         type is (source_type)
            s = source%local_source_index
            source_offset = global_section_offset(source_section, s, &
              source_network%source_range_start)
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
            case (2)
               call test%assert(source%separator%on, &
                    "source 3 separator on")
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
         end select

      type is (recharge_source_control_type)
         select type (source => source_control%objects%head%data)
         type is (source_type)
            s = source%local_source_index
            source_offset = global_section_offset(source_section, s, &
                 source_network%source_range_start)
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
         end select

      type is (limiter_table_source_network_control_type)
         select case (source_control%flow_type(1))
            case (SEPARATED_FLOW_TYPE_TOTAL)
               call test%assert(10._dp, source_control%table(1)%ptr%average(interval, 1), &
                    "total limiter limit")
            case (SEPARATED_FLOW_TYPE_STEAM)
               call test%assert(5._dp, source_control%table(1)%ptr%average(interval, 1), &
                    "steam limiter limit")
            end select
      end select
      stopped = PETSC_FALSE

    end subroutine source_control_test_iterator

    subroutine source_rate_test_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      PetscInt :: s, source_index
      character(len = 16) :: srcstr

      select type (source => node%data)
      type is (source_type)
         s = source%local_source_index
         source_offset = global_section_offset(source_section, s, &
              source_network%source_range_start)
         call source%assign(source_array, source_offset)
         source_index = nint(source%source_index)
         write(srcstr, '(a, i2)') 'source ', source_index + 1
         call test%assert(expected_rates(source_index), source%rate, &
              trim(srcstr) // ' rate')
         call test%assert(expected_steam_rates(source_index), source%steam_rate, &
              trim(srcstr) // ' steam rate')
         if (source_index == 12) s12 = s
      end select

    end subroutine source_rate_test_iterator

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

  end subroutine test_source_control_pressure_reference

!------------------------------------------------------------------------

end module source_control_test
