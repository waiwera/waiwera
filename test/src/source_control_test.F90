module source_control_test

  ! Test for source control module

  use kinds_module
  use mpi_module
  use fruit
  use source_module
  use source_control_module
  use source_setup_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  public :: test_source_control_table

contains

!------------------------------------------------------------------------

  subroutine test_source_control_table
    ! table source control tests

    use fson
    use fson_mpi_module
    use mesh_module
    use list_module
    use IAPWS_module
    use eos_module, only: max_component_name_length, max_phase_name_length
    use eos_test
    use fluid_module, only: setup_fluid_vector
    use dm_utils_module, only: global_to_local_vec_section, restore_dm_local_vec

    character(16), parameter :: path = "data/source/"
    type(IAPWS_type) :: thermo
    type(eos_test_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(list_type) :: sources, source_controls
    Vec :: global_fluid, local_fluid
    PetscReal, pointer, contiguous :: fluid_array(:)
    PetscSection :: fluid_section
    PetscInt :: num_sources, range_start
    PetscReal :: t, interval(2)
    PetscErrorCode :: ierr

    json => fson_parse_mpi(trim(path) // "test_source_controls_table.json")

    call thermo%init()
    call eos%init(json, thermo)
    call mesh%init(json)
    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(eos%primary_variable_names)
    call setup_fluid_vector(mesh%dm, max_component_name_length, &
         eos%component_names, max_phase_name_length, eos%phase_names, &
         global_fluid, range_start)
    call global_to_local_vec_section(global_fluid, local_fluid, fluid_section)
    call VecGetArrayReadF90(local_fluid, fluid_array, ierr); CHKERRQ(ierr)

    call setup_sources(json, mesh%dm, eos, sources, source_controls)

    call MPI_reduce(sources%count, num_sources, 1, MPI_INTEGER, MPI_SUM, &
         mpi%input_rank, mpi%comm, ierr)
    if (mpi%rank == mpi%input_rank) then
      call assert_equals(6, num_sources, "number of sources")
    end if

    t = 120._dp
    interval = [30._dp, t]
    call source_controls%traverse(source_control_iterator)
    call sources%traverse(source_test_iterator)

    call sources%destroy(source_list_node_data_destroy)
    call source_controls%destroy(source_control_list_node_data_destroy)

    call VecRestoreArrayReadF90(local_fluid, fluid_array, ierr); CHKERRQ(ierr)
    call restore_dm_local_vec(local_fluid)
    call VecDestroy(global_fluid, ierr); CHKERRQ(ierr)
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

end module source_control_test
