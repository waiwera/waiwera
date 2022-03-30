module source_network_control_test

  ! Test for source network control module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use source_module
  use control_module
  use source_control_module
  use source_network_group_module
  use source_network_control_module
  use source_setup_module

  implicit none
  private

  character(len = 512) :: data_path

  public :: setup, teardown
  public :: test_source_network_limiter

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

  subroutine test_source_network_limiter(test)
    ! Source network limiter

    use fson
    use fson_mpi_module
    use IAPWS_module
    use eos_we_module
    use mesh_module
    use tracer_module
    use dm_utils_module, only: global_vec_section, global_section_offset

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(IAPWS_type) :: thermo
    type(eos_we_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(tracer_type), allocatable :: tracers(:)
    Vec :: fluid_vector, source_vector, group_vector
    PetscInt :: fluid_range_start, source_range_start, group_range_start
    PetscReal, pointer, contiguous :: source_array(:), group_array(:)
    PetscSection :: source_section, group_section
    type(list_type) :: sources, source_controls, source_network_groups, &
         separated_sources, source_network_controls
    IS :: source_index, source_network_group_index
    PetscInt :: total_num_sources, total_num_source_network_groups
    PetscMPIInt :: rank
    PetscReal, parameter :: start_time = 0._dp, end_time = 100._dp
    PetscReal, parameter :: interval(2) = [start_time, end_time]
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscErrorCode :: err, ierr
    PetscInt, parameter :: expected_num_sources = 22
    PetscInt, parameter :: expected_num_groups = 13

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    json => fson_parse_mpi(trim(adjustl(data_path)) // "source/test_source_network_limiter.json")

    call thermo%init()
    call eos%init(json, thermo)
    call setup_tracers(json, eos, tracers, err = err)
    call mesh%init(eos, json)
    call DMCreateLabel(mesh%serial_dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(gravity, json, err = err)
    call DMGetGlobalVector(mesh%dm, fluid_vector, ierr); CHKERRQ(ierr) ! dummy- not used

    call setup_sources(json, mesh%dm, mesh%cell_natural_global, eos, tracers%name, &
         thermo, start_time, fluid_vector, fluid_range_start, source_vector, &
         source_range_start, group_vector, group_range_start, &
         sources, total_num_sources, total_num_source_network_groups, source_controls, &
         source_index, source_network_group_index, separated_sources, source_network_groups, &
         source_network_controls, err = err)

    call test%assert(0, err, "error")
    if (rank == 0) then
       call test%assert(expected_num_sources, total_num_sources, "number of sources")
       call test%assert(expected_num_groups, total_num_source_network_groups, "number of groups")
    end if

    if (err == 0) then

       call global_vec_section(source_vector, source_section)
       call global_vec_section(group_vector, group_section)
       call VecGetArrayF90(source_vector, source_array, ierr); CHKERRQ(ierr)
       call VecGetArrayF90(group_vector, group_array, ierr); CHKERRQ(ierr)

       call sources%traverse(source_setup_iterator)
       call source_network_groups%traverse(group_assign_iterator)
       call separated_sources%traverse(source_separator_iterator)
       call source_network_groups%traverse(group_sum_iterator)
       call source_network_controls%traverse(network_control_iterator)
       call sources%traverse(source_test_iterator)
       call source_network_groups%traverse(group_test_iterator)

       call VecRestoreArrayF90(group_vector, group_array, ierr); CHKERRQ(ierr)
       call VecRestoreArrayF90(source_vector, source_array, ierr); CHKERRQ(ierr)

    end if

    call ISDestroy(source_index, ierr); CHKERRQ(ierr)
    call ISDestroy(source_network_group_index, ierr); CHKERRQ(ierr)
    call VecDestroy(source_vector, ierr); CHKERRQ(ierr)
    call VecDestroy(group_vector, ierr); CHKERRQ(ierr)
    call separated_sources%destroy()
    call source_network_groups%destroy(source_network_group_list_node_data_destroy)
    call source_controls%destroy(source_control_list_node_data_destroy)
    call source_network_controls%destroy(source_control_list_node_data_destroy)
    call sources%destroy(source_list_node_data_destroy)
    call DMRestoreGlobalVector(mesh%dm, fluid_vector, ierr); CHKERRQ(ierr)
    call mesh%destroy_distribution_data()
    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  contains

    subroutine source_setup_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: s, source_offset

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         s = source%local_source_index
         source_offset = global_section_offset(source_section, &
              s, source_range_start)
         call source%assign(source_array, source_offset)
         select case (source%name)
         case ("s1")
            source%enthalpy = 500.e3_dp
         case ("s2")
            source%enthalpy = 800.e3_dp
         case ("s3")
            source%enthalpy = 500.e3_dp
         case ("s4")
            source%enthalpy = 800.e3_dp
         case ("s5")
            source%enthalpy = 500.e3_dp
         case ("s6")
            source%enthalpy = 800.e3_dp
         case ("s7")
            source%enthalpy = 1200.e3_dp
         case ("s8")
            source%enthalpy = 1750.e3_dp
         case ("s9")
            source%enthalpy = 500.e3_dp
         case ("s10")
            source%enthalpy = 600.e3_dp
         case ("s11")
            source%enthalpy = 450.e3_dp
         case ("s12")
            source%enthalpy = 500.e3_dp
         case ("s13")
            source%enthalpy = 600.e3_dp
         case ("s14")
            source%enthalpy = 475.e3_dp
         case ("s15")
            source%enthalpy = 1200.e3_dp
         case ("s16")
            source%enthalpy = 1750.e3_dp
         case ("s17")
            source%enthalpy = 500.e3_dp
         case ("s18")
            source%enthalpy = 800.e3_dp
         case ("s19")
            source%enthalpy = 1500.e3_dp
         case ("s20")
            source%enthalpy = 1800.e3_dp
         case ("s21")
            source%enthalpy = 1600.e3_dp
         case ("s22")
            source%enthalpy = 1400.e3_dp
         end select
      end select

    end subroutine source_setup_iterator

!........................................................................

    subroutine group_assign_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: g, group_offset

      stopped = PETSC_FALSE
      select type(group => node%data)
      class is (source_network_group_type)
         if (group%rank == 0) then
            g = group%local_group_index
            group_offset = global_section_offset(group_section, &
              g, group_range_start)
            call group%assign(group_array, group_offset)
         end if
      end select

    end subroutine group_assign_iterator

!........................................................................

    subroutine source_separator_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         call source%get_separated_flows()
      end select

    end subroutine source_separator_iterator

!........................................................................

    subroutine group_sum_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (group => node%data)
      class is (source_network_group_type)
         call group%sum()
      end select

    end subroutine group_sum_iterator

!........................................................................

    subroutine network_control_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type(control => node%data)
      class is (interval_update_object_control_type)
         call control%update(interval)
      end select

    end subroutine network_control_iterator

!........................................................................

    subroutine flow_test(node, rate, water_rate, steam_rate)

      class(source_network_node_type), intent(in) :: node
      PetscReal, intent(in) :: rate, water_rate, steam_rate

      call test%assert(rate, node%rate, trim(node%name) // " rate")
      call test%assert(water_rate, node%water_rate, &
           trim(node%name) // " water rate")
      call test%assert(steam_rate, node%steam_rate, &
           trim(node%name) //" steam rate")

    end subroutine flow_test

!........................................................................

    subroutine source_test_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type(source => node%data)
      type is (source_type)
         select case (source%name)
         case ("s1")
            call flow_test(source, -3.2_dp, 0.0_dp, 0.0_dp)
         case ("s2")
            call flow_test(source, -4.8_dp, 0.0_dp, 0.0_dp)
         case ("s3")
            call flow_test(source, -1.5_dp, 0.0_dp, 0.0_dp)
         case ("s4")
            call flow_test(source, -3.5_dp, 0.0_dp, 0.0_dp)
         case ("s5")
            call flow_test(source, -2.5_dp, 0.0_dp, 0.0_dp)
         case ("s6")
            call flow_test(source, -4.0_dp, 0.0_dp, 0.0_dp)
         case ("s7")
            call flow_test(source, -1.2556341098337822_dp, &
                 -0.922167171793673_dp, -0.333466938040109_dp)
         case ("s8")
            call flow_test(source, -1.6741788131117095_dp, &
                 -0.8076457511518187_dp, -0.8665330619598908_dp)
         case ("s9")
            call flow_test(source, -1.2_dp, 0.0_dp, 0.0_dp)
         case ("s10")
            call flow_test(source, -1.6_dp, 0.0_dp, 0.0_dp)
         case ("s11")
            call flow_test(source, -5.0_dp, 0.0_dp, 0.0_dp)
         case ("s12")
            call flow_test(source, -0.3_dp, 0.0_dp, 0.0_dp)
         case ("s13")
            call flow_test(source, -0.6_dp, 0.0_dp, 0.0_dp)
         case ("s14")
            call flow_test(source, -3.6_dp, 0.0_dp, 0.0_dp)
         case ("s15")
            call flow_test(source, -2.768782992805067_dp, 0.0_dp, 0.0_dp)
         case ("s16")
            call flow_test(source, -4.153174489207601_dp, 0.0_dp, 0.0_dp)
         case ("s17")
            call flow_test(source, -6._dp, 0.0_dp, 0.0_dp)
         case ("s18")
            call flow_test(source, -1._dp, 0.0_dp, 0.0_dp)
         case ("s19")
            call flow_test(source, -4.0_dp, -2.3684129664887807_dp, &
                 -1.631587033511219_dp)
         case ("s20")
            call flow_test(source, -3.007165178154069_dp, 0._dp, 0._dp)
         case ("s21")
            call flow_test(source, -3.6085982137848824_dp, 0._dp, 0._dp)
         case ("s22")
            call flow_test(source, 0._dp, 0._dp, 0._dp)
         end select
      end select

    end subroutine source_test_iterator

!........................................................................

    subroutine group_test_iterator(node, stopped)

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type(group => node%data)
      class is (source_network_group_type)
         if (group%rank == 0) then
            select case (group%name)
            case ("group1")
               call flow_test(group, -8.0_dp, 0.0_dp, 0.0_dp)
            case ("group2a")
               call flow_test(group, -5.0_dp, 0.0_dp, 0.0_dp)
            case ("group2b")
               call flow_test(group, -6.5_dp, 0.0_dp, 0.0_dp)
            case ("group2")
               call flow_test(group, -11.5_dp, 0.0_dp, 0.0_dp)
            case ("group3")
               call flow_test(group, -2.9298129229454917_dp, &
                    -1.7298129229454917_dp, -1.2_dp)
            case ("group4a")
               call flow_test(group, -2.8_dp, 0.0_dp, 0.0_dp)
            case ("group4")
               call flow_test(group, -7.8_dp, 0.0_dp, 0.0_dp)
            case ("group5a")
               call flow_test(group, -0.9_dp, 0.0_dp, 0.0_dp)
            case ("group5")
               call flow_test(group, -4.5_dp, 0.0_dp, 0.0_dp)
            case ("group6")
               call flow_test(group, -6.921957482012667_dp, -4.0_dp, &
                    -2.921957482012667_dp)
            case ("group7")
               call flow_test(group, -7.0_dp, 0.0_dp, 0.0_dp)
            case ("group8a")
               call flow_test(group, -6.615763391938951_dp, &
                    -3.2473504254501706_dp, -3.3684129664887807_dp)
            case ("group8")
               call flow_test(group, -10.615763391938952_dp, &
                    -5.615763391938952_dp, -5.0_dp)
            end select
         end if
      end select

    end subroutine group_test_iterator

  end subroutine test_source_network_limiter

!------------------------------------------------------------------------

  subroutine source_control_list_node_data_destroy(node)
    type(list_node_type), pointer, intent(in out) :: node
    select type (source_control => node%data)
    class is (object_control_type)
       call source_control%destroy()
    end select
  end subroutine source_control_list_node_data_destroy

  subroutine source_network_group_list_node_data_destroy(node)
    type(list_node_type), pointer, intent(in out) :: node
    select type (source_network_group => node%data)
    class is (source_network_group_type)
       call source_network_group%destroy()
    end select
  end subroutine source_network_group_list_node_data_destroy

  subroutine source_list_node_data_destroy(node)
    type(list_node_type), pointer, intent(in out) :: node
    select type (source => node%data)
    class is (source_type)
       call source%destroy()
    end select
  end subroutine source_list_node_data_destroy

!------------------------------------------------------------------------

end module source_network_control_test
