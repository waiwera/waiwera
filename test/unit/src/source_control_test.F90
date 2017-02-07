module source_control_test

  ! Test for source control module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use fruit
  use source_module
  use source_control_module
  use source_setup_module

  implicit none
  private

  public :: test_source_control_table, test_source_control_pressure_reference

contains

!------------------------------------------------------------------------

  subroutine test_source_control_table
    ! Table source control

    use fson
    use fson_mpi_module
    use mesh_module
    use list_module
    use IAPWS_module
    use eos_module, only: max_component_name_length, max_phase_name_length
    use eos_test
    use fluid_module, only: setup_fluid_vector
    use dm_utils_module, only: global_vec_section, &
         global_to_local_vec_section, restore_dm_local_vec

    character(16), parameter :: path = "data/source/"
    type(IAPWS_type) :: thermo
    type(eos_test_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(list_type) :: sources, source_controls
    Vec :: fluid_vector, local_fluid_vector
    PetscReal, pointer, contiguous :: fluid_array(:), local_fluid_array(:)
    PetscSection :: fluid_section, local_fluid_section
    PetscInt :: num_sources, fluid_range_start
    PetscReal :: t, interval(2)
    PetscErrorCode :: ierr
    PetscReal, parameter :: start_time = 0._dp
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscMPIInt :: rank

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    json => fson_parse_mpi(trim(path) // "test_source_controls_table.json")

    call thermo%init()
    call eos%init(json, thermo)
    call mesh%init(json)
    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(eos%primary_variable_names, gravity)
    call setup_fluid_vector(mesh%dm, max_component_name_length, &
         eos%component_names, max_phase_name_length, eos%phase_names, &
         fluid_vector, fluid_range_start)
    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call setup_sources(json, mesh%dm, eos, thermo, start_time, fluid_vector, &
         fluid_range_start, sources, source_controls)

    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    call MPI_reduce(sources%count, num_sources, 1, MPI_INTEGER, MPI_SUM, &
         0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
      call assert_equals(6, num_sources, "number of sources")
    end if

    call global_to_local_vec_section(fluid_vector, local_fluid_vector, &
         local_fluid_section)
    call VecGetArrayReadF90(local_fluid_vector, local_fluid_array, ierr)
    CHKERRQ(ierr)

    t = 120._dp
    interval = [30._dp, t]
    call source_controls%traverse(source_control_iterator)
    call sources%traverse(source_test_iterator)

    call VecRestoreArrayReadF90(local_fluid_vector, local_fluid_array, ierr)
    CHKERRQ(ierr)
    call restore_dm_local_vec(local_fluid_vector)

    call sources%destroy(source_list_node_data_destroy)
    call source_controls%destroy(source_control_list_node_data_destroy, &
         reverse = PETSC_TRUE)

    call VecDestroy(fluid_vector, ierr); CHKERRQ(ierr)
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
         call source_control%update(t, interval, fluid_array, fluid_section)
      end select
      stopped = PETSC_FALSE
    end subroutine source_control_iterator

    subroutine source_test_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscReal, parameter :: tol = 1.e-6_dp
      select type (source => node%data)
      class is (source_type)
         select case (node%tag)
         case ("rate table single source")
            call assert_equals(-2.25_dp, source%rate, tol, node%tag)
         case ("rate table 3 sources")
            call assert_equals(2.25_dp, source%rate, tol, node%tag)
         case ("enthalpy table")
            call assert_equals(104.5e3_dp, source%injection_enthalpy, &
                 tol, node%tag)
         end select
      end select
      stopped = PETSC_FALSE
    end subroutine source_test_iterator

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

  end subroutine test_source_control_table

!------------------------------------------------------------------------

  subroutine test_source_control_pressure_reference
    ! Pressure reference source controls

    use fson
    use fson_mpi_module
    use mesh_module
    use list_module
    use IAPWS_module
    use eos_module, only: max_component_name_length, max_phase_name_length
    use eos_test
    use fluid_module, only: fluid_type, setup_fluid_vector
    use dm_utils_module, only: global_vec_section, global_section_offset, &
         global_to_local_vec_section, restore_dm_local_vec
    use interpolation_module, only: INTERP_STEP, INTERP_AVERAGING_ENDPOINT

    character(16), parameter :: path = "data/source/"
    type(IAPWS_type) :: thermo
    type(eos_test_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(list_type) :: sources, source_controls
    Vec :: fluid_vector, local_fluid_vector
    PetscReal, pointer, contiguous :: fluid_array(:), local_fluid_array(:)
    PetscSection :: fluid_section, local_fluid_section
    type(fluid_type) :: fluid
    PetscInt :: num_sources, num_source_controls, fluid_range_start, c
    PetscInt :: fluid_offset, cell_phase_composition
    PetscReal :: t, interval(2), props(2)
    PetscReal :: cell_temperature, cell_liquid_density, cell_liquid_internal_energy
    PetscReal :: cell_vapour_density, cell_vapour_internal_energy
    PetscReal :: cell_liquid_viscosity, cell_vapour_viscosity
    PetscErrorCode :: ierr
    PetscReal, parameter :: cell_pressure = 50.e5_dp, cell_vapour_saturation = 0.8_dp
    PetscInt, parameter :: cell_region = 4
    PetscReal, parameter :: start_time = 0._dp
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscMPIInt :: rank

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    json => fson_parse_mpi(trim(path) // "test_source_controls_pressure_reference.json")

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

    call mesh%init(json)
    call fluid%init(eos%num_components, eos%num_phases)
    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(eos%primary_variable_names, gravity)
    call setup_fluid_vector(mesh%dm, max_component_name_length, &
         eos%component_names, max_phase_name_length, eos%phase_names, &
         fluid_vector, fluid_range_start)
    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    do c = mesh%start_cell, mesh%end_cell - 1
       if (mesh%ghost_cell(c) < 0) then
          call global_section_offset(fluid_section, c, fluid_range_start, &
               fluid_offset, ierr)
          CHKERRQ(ierr)
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

    call setup_sources(json, mesh%dm, eos, thermo, start_time, fluid_vector, &
         fluid_range_start, sources, source_controls)

    call MPI_reduce(sources%count, num_sources, 1, MPI_INTEGER, MPI_SUM, &
         0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
      call assert_equals(12, num_sources, "number of sources")
    end if

    call MPI_reduce(source_controls%count, num_source_controls, 1, &
         MPI_INTEGER, MPI_SUM, 0, PETSC_COMM_WORLD, ierr)
    if (rank == 0) then
      call assert_equals(20, num_source_controls, "number of source controls")
    end if

    call global_to_local_vec_section(fluid_vector, local_fluid_vector, &
         local_fluid_section)
    call VecGetArrayReadF90(local_fluid_vector, local_fluid_array, ierr)
    CHKERRQ(ierr)

    t = 120._dp
    interval = [30._dp, t]

    call source_controls%traverse(source_control_update_iterator)
    call source_controls%traverse(source_control_test_iterator)
    call sources%traverse(source_test_iterator)

    call VecRestoreArrayReadF90(local_fluid_vector, local_fluid_array, ierr)
    CHKERRQ(ierr)
    call restore_dm_local_vec(local_fluid_vector)

    call sources%destroy(source_list_node_data_destroy)
    call source_controls%destroy(source_control_list_node_data_destroy, &
      reverse = PETSC_TRUE)

    call fluid%destroy()
    call VecDestroy(fluid_vector, ierr); CHKERRQ(ierr)
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
         call source_control%update(t, interval, local_fluid_array, &
              local_fluid_section)
      end select
      stopped = PETSC_FALSE
    end subroutine source_control_update_iterator

    subroutine source_control_test_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      PetscReal, parameter :: PI_tol = 1.e-16_dp, tol = 1.e-6_dp
      select type (source_control => node%data)

      type is (source_control_deliverability_type)
         select case (source_control%sources%head%tag)
         case ("source 2")
            call assert_equals(1.e-12_dp, &
                 source_control%productivity%val(1), PI_tol, &
                 "source 2 productivity")
            call assert_equals(2.e5_dp, &
                 source_control%reference_pressure%val(1), tol, &
                 "source 2 reference pressure")
         case ("source 4")
            call assert_equals(8.54511496085953E-13_dp, &
                 source_control%productivity%val(1), PI_tol, &
                 "source 4 productivity")
         case ("source 10")
            call assert_equals(SRC_PRESSURE_TABLE_COORD_TIME, &
                 source_control%pressure_table_coordinate, &
                 "source 10 pressure table coordinate")
            call assert_equals(4, &
                 source_control%reference_pressure%coord%size, &
                 "source 10 pressure table size")
            call assert_equals(INTERP_STEP, &
                 source_control%reference_pressure%interpolation_type, &
                 "source 10 interpolation type")
            call assert_equals(INTERP_AVERAGING_ENDPOINT, &
                 source_control%reference_pressure%averaging_type, &
                 "source 10 averaging type")
            call assert_equals(1.8e5_dp, &
                 source_control%reference_pressure%average(interval), &
                 tol, "source 10 reference pressure")
         case ("source 11")
            call assert_equals(SRC_PRESSURE_TABLE_COORD_ENTHALPY, &
                 source_control%pressure_table_coordinate, &
                 "source 11 pressure table coordinate")
         end select

      type is (source_control_recharge_type)
         select case (source_control%sources%head%tag)
         case ("source 7")
            call assert_equals(1.3e-2_dp, &
                 source_control%coefficient%val(1), tol, &
                 "source 7 recharge coefficient")
            call assert_equals(50.1e5_dp, &
                 source_control%reference_pressure%val(1), &
                 tol, "source 7 reference pressure")
         case ("source 12")
            call assert_equals(cell_pressure, &
                 source_control%reference_pressure%val(1), &
                 tol, "source 12 reference pressure")
         end select

      type is (source_control_limiter_type)
         select case (source_control%type)
         case (SRC_CONTROL_LIMITER_TYPE_TOTAL)
            call assert_equals(10._dp, source_control%limit, tol, "total limiter limit")
         case (SRC_CONTROL_LIMITER_TYPE_STEAM)
            call assert_equals(5._dp, source_control%limit, tol, "steam limiter limit")
         end select

      type is (source_control_separator_type)
         call assert_equals(10.e5_dp, &
              source_control%separator_pressure, tol, "separator pressure")
         call assert_equals(0.5371645375_dp, &
              source_control%steam_fraction, tol, "separator steam fraction")
         call assert_equals(-5.9580123975_dp, &
              source_control%water_flow_rate, tol, "separator water flow rate")
         call assert_equals(-6.9148395774_dp, &
              source_control%steam_flow_rate, tol, "separator steam flow rate")
      end select
      stopped = PETSC_FALSE

    end subroutine source_control_test_iterator

    subroutine source_test_iterator(node, stopped)
      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscReal, parameter :: tol = 1.e-6_dp
      select type (source => node%data)
      class is (source_type)
         select case (node%tag)
         case ("source 1")
            call assert_equals(-12.8728519749_dp, source%rate, tol, "source 1 rate")
         case ("source 2")
            call assert_equals(-10._dp, source%rate, tol, "source 2 rate")
         case ("source 3")
            call assert_equals(-9.3081349399_dp, source%rate, tol, "source 3 rate")
         case ("source 4")
            call assert_equals(-11._dp, source%rate, tol, "source 4 rate")
         case ("source 5")
            call assert_equals(-10.1910078135_dp, source%rate, tol, "source 5 rate")
         case ("source 6")
            call assert_equals(0.0_dp, source%rate, tol, "source 6 rate")
         case ("source 7")
            call assert_equals(130.0_dp, source%rate, tol, "source 7 rate")
         case ("source 8")
            call assert_equals(-13.0_dp, source%rate, tol, "source 8 rate")
         case ("source 9")
            call assert_equals(0.0_dp, source%rate, tol, "source 9 rate")
         case ("source 10")
            call assert_equals(-12.9264888581701_dp, source%rate, tol, "source 10 rate")
         case ("source 11")
            call assert_equals(-10.3366086953508_dp, source%rate, tol, "source 11 rate")
         case ("source 12")
            call assert_equals(0._dp, source%rate, tol, "source 12 rate")
         end select
      end select
      stopped = PETSC_FALSE
    end subroutine source_test_iterator

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
