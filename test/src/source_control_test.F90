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

  public :: test_source_control_rate

contains

!------------------------------------------------------------------------

  subroutine test_source_control_rate
    ! rate source control test

    use fson
    use fson_mpi_module
    use mesh_module
    use list_module
    use IAPWS_module
    use eos_test

    character(16), parameter :: path = "data/source/"
    type(IAPWS_type) :: thermo
    type(eos_test_type) :: eos
    type(fson_value), pointer :: json
    type(mesh_type) :: mesh
    type(list_type) :: sources, source_controls
    PetscReal :: t, interval(2)
    PetscErrorCode :: ierr

    json => fson_parse_mpi(trim(path) // "test_source_controls.json")

    call thermo%init()
    call eos%init(json, thermo)
    call mesh%init(json)
    call DMCreateLabel(mesh%dm, open_boundary_label_name, ierr); CHKERRQ(ierr)
    call mesh%configure(eos%primary_variable_names)

    call setup_sources(json, mesh%dm, eos, sources, source_controls)

    t = 0.5_dp
    interval = [0._dp, t]
    call source_controls%traverse(source_control_iterator)
    ! call sources%traverse(source_test_iterator)

    call sources%destroy(source_list_node_data_destroy)
    call source_controls%destroy()
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
         call source_control%update(t, interval)
      end select
      stopped = PETSC_FALSE
    end subroutine source_control_iterator

    subroutine source_list_node_data_destroy(node)
      type(list_node_type), pointer, intent(in out) :: node
      select type (source => node%data)
      type is (source_type)
         call source%destroy()
      end select
    end subroutine source_list_node_data_destroy

  end subroutine test_source_control_rate

!------------------------------------------------------------------------

end module source_control_test
