module rock_control_test

  ! Test for rock control module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use rock_module
  use rock_control_module
  use rock_setup_module

  implicit none
  private

  character(len = 512) :: data_path

  public :: setup, teardown
  public :: test_rock_control_table

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
    PetscInt, parameter :: expected_num_rocks = 3
    character(10), parameter :: rock_names(expected_num_rocks) = &
         ["constant", "scalar  ", "array   "]
    PetscReal, parameter :: start_time = 0._dp, t = 4500._dp
    PetscReal, parameter :: expected_initial_permeability(expected_num_rocks, 3) = reshape([ &
         1.e-13_dp, 1.e-13_dp, 1.e-14_dp, &
         1.e-13_dp, 1.e-13_dp, 2.e-14_dp, &
         1.e-13_dp, 1.e-13_dp, 3.e-14_dp], &
         [expected_num_rocks, 3])
    PetscReal, parameter :: expected_permeability(expected_num_rocks, 3) = reshape([ &
         1.e-13_dp, 7.e-14_dp, 4.e-15_dp, &
         1.e-13_dp, 7.e-14_dp, 5.e-15_dp, &
         1.e-13_dp, 7.e-14_dp, 6.e-15_dp], &
         [expected_num_rocks, 3])
    PetscReal, parameter :: expected_initial_porosity(expected_num_rocks) = &
         [0.1_dp, 0.1_dp, 0.2_dp]
    PetscReal, parameter :: expected_porosity(expected_num_rocks) = &
         [0.1_dp, 0.05_dp, 0.2_dp]
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscReal, parameter :: ktol = 1.e-20_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    json => fson_parse_mpi(trim(adjustl(data_path)) // "rock/test_rock_controls_table.json")

    call thermo%init()
    call eos%init(json, thermo)
    call mesh%init(eos, json)
    call DMCreateLabel(mesh%serial_dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(gravity, json, err = err)

    call setup_rocks(json, mesh%dm, start_time, rock_vector, mesh%rock_types, &
         rock_controls, rock_range_start, mesh%ghost_cell, err = err)
    call test%assert(0, err, "rock setup error")
    num_rock_types = mesh%rock_types%count()
    if (rank == 0) then
       call test%assert(expected_num_rocks, num_rock_types, "number of rocks")
    end if

    call VecGetArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
    call global_vec_section(rock_vector, rock_section)
    call DMGetLabel(mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    call rock%init()

    call rock_test("t = 0", expected_initial_permeability, expected_initial_porosity)

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
      class is (rock_control_type)
         call rock_control%update(t, rock_array, rock_section, rock_range_start)
      end select
      stopped = PETSC_FALSE
    end subroutine rock_control_iterator

     subroutine rock_control_list_node_data_destroy(node)
      type(list_node_type), pointer, intent(in out) :: node
      select type (rock_control => node%data)
      class is (rock_control_type)
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

end module rock_control_test
