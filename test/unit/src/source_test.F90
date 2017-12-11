module source_test

  ! Test for source module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use fruit
  use source_module

  implicit none
  private

public :: test_source_update_flow

contains

!------------------------------------------------------------------------

  subroutine test_source_update_flow
    ! update_flow() test

    use fson
    use fson_mpi_module
    use IAPWS_module, only: IAPWS_type
    use eos_module, only: max_component_name_length, max_phase_name_length
    use IAPWS_module
    use eos_test
    use mesh_module
    use fluid_module, only: fluid_type, setup_fluid_vector
    use dm_utils_module, only: global_vec_section, global_section_offset, &
         global_to_local_vec_section, restore_dm_local_vec

    type(source_type) :: source
    type(fluid_type) :: fluid
    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_test_type) :: eos
    type(mesh_type) :: mesh
    Vec :: fluid_vector, local_fluid_vector
    PetscSection :: fluid_section, local_fluid_section
    PetscInt :: start_cell, end_cell, c, fluid_offset, fluid_range_start
    PetscReal, allocatable :: fluid_cell_data(:)
    PetscReal, pointer, contiguous :: fluid_array(:), local_fluid_array(:)
    PetscErrorCode :: ierr, err
    PetscInt, parameter :: offset = 1
    PetscReal, parameter :: tol = 1.e-6_dp
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscMPIInt :: rank
    PetscViewer :: viewer

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    json => fson_parse_mpi(str = '{"mesh": "data/flow_simulation/mesh/3x3_2d.exo"}')
    call thermo%init()
    call eos%init(json, thermo)
    viewer = PETSC_NULL_VIEWER

    call mesh%init(json)
    call fluid%init(eos%num_components, eos%num_phases)
    call mesh%configure(eos, gravity, json, viewer = viewer, err = err)
    call setup_fluid_vector(mesh%dm, max_component_name_length, &
         eos%component_names, max_phase_name_length, eos%phase_names, &
         fluid_vector, fluid_range_start)
    call DMPlexGetHeightStratum(mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    call global_vec_section(fluid_vector, fluid_section)
    call VecGetArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)

    fluid_cell_data = [2.7e5_dp, 130._dp, 4._dp, 3._dp, &
         935._dp, 1.e-6_dp, 0.8_dp, 0.7_dp, 0._dp, 83.9e3_dp, 5.461e5_dp, 0.7_dp, 0.3_dp, &
         1.5_dp,  2.e-7_dp, 0.2_dp, 0.3_dp, 0._dp, 800.e3_dp, 2.540e6_dp, 0.4_dp, 0.6_dp]

    do c = start_cell, end_cell - 1
       if (mesh%ghost_cell(c) < 0) then
          call global_section_offset(fluid_section, c, fluid_range_start, &
               fluid_offset, ierr)
          CHKERRQ(ierr)
          fluid_array(fluid_offset: fluid_offset + fluid%dof - 1) = fluid_cell_data
       end if
    end do
    call VecRestoreArrayF90(fluid_vector, fluid_array, ierr); CHKERRQ(ierr)
    deallocate(fluid_cell_data)

    call global_to_local_vec_section(fluid_vector, local_fluid_vector, &
         local_fluid_section)
    call VecGetArrayReadF90(local_fluid_vector, local_fluid_array, ierr)
    CHKERRQ(ierr)

    if (rank == 0) then

       call source_flow_test("inject 1", 10._dp, 200.e3_dp, 1, 0, &
            [10._dp, 0._dp, 2.e6_dp], 1)
       call source_flow_test("inject 2", 5._dp, 200.e3_dp, 2, 0, &
            [0._dp, 5._dp, 1.e6_dp], 2)
       call source_flow_test("inject 3", 5._dp, 200.e3_dp, 2, 2, &
            [0._dp, 5._dp, 1.e6_dp], 2)
       call source_flow_test("inject heat", 1000._dp, 0._dp, 3, 0, &
            [0._dp, 0._dp, 1000._dp], 3)

       call source_flow_test("produce all", -5._dp, 0._dp, 0, 0, &
            [-3.4948610582_dp, -1.5051389418_dp, -431766.653977922_dp], 0)
       call source_flow_test("produce 1", -5._dp, 0._dp, 0, 1, &
            [-5._dp, 0._dp, -431766.653977922_dp], 1)
       call source_flow_test("produce heat", -5000._dp, 0._dp, 0, 3, &
            [0._dp, 0._dp, -5000._dp], 3)

       call source_flow_test("no flow 1", 0._dp, 100.e3_dp, 1, 0, &
            [0._dp, 0._dp, 0._dp], 1)
       call source_flow_test("no flow all", 0._dp, 100.e3_dp, 0, 0, &
            [0._dp, 0._dp, 0._dp], 1)

    end if

    call fluid%destroy()
    call VecRestoreArrayReadF90(local_fluid_vector, local_fluid_array, ierr)
    CHKERRQ(ierr)
    call restore_dm_local_vec(local_fluid_vector)
    call VecDestroy(fluid_vector, ierr); CHKERRQ(ierr)
    call mesh%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)

  contains

    subroutine source_flow_test(tag, rate, enthalpy, injection_component, &
         production_component, flow, component)
      !! Runs asserts for single flow update_source() test.

      character(*), intent(in) :: tag
      PetscReal, intent(in) :: rate, enthalpy
      PetscInt, intent(in) :: injection_component, production_component
      PetscReal, intent(in) :: flow(:)
      PetscInt, intent(in) :: component

      call source%init(0, 0, eos, rate, enthalpy, &
           injection_component, production_component)
      call source%update_flow(local_fluid_array, local_fluid_section)
      call assert_equals(flow, source%flow, &
           eos%num_primary_variables, tol, trim(tag) // " flow")
      call assert_equals(component, source%component, &
           trim(tag) // " component")
      call source%destroy()

    end subroutine source_flow_test

  end subroutine test_source_update_flow

!------------------------------------------------------------------------

end module source_test
