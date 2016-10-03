module source_test

  ! Test for source module

  use kinds_module
  use mpi_module
  use fruit
  use source_module
  use eos_test, only: eos_test_type

  implicit none
  private

#include <petsc/finclude/petsc.h90>

public :: test_source_update_flow

contains

!------------------------------------------------------------------------

  subroutine test_source_update_flow
    ! update_flow() test

    use fson
    use fluid_module, only: fluid_type
    use IAPWS_module, only: IAPWS_type

    type(source_type) :: source
    type(fluid_type) :: fluid
    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_test_type) :: eos
    PetscInt, parameter :: offset = 1
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal, parameter :: tol = 1.e-6_dp

    if (mpi%rank == mpi%output_rank) then

       json => fson_parse(str = "{}")
       call thermo%init()
       call eos%init(json, thermo)
       call fluid%init(eos%num_components, eos%num_phases)

       allocate(fluid_data(offset - 1 + fluid%dof))
       fluid_data = [2.7e5_dp, 130._dp, 4._dp, 3._dp, &
            935._dp, 1.e-6_dp, 0.8_dp, 0.7_dp, 83.9e3_dp, 5.461e5_dp, 0.7_dp, 0.3_dp, &
            1.5_dp,  2.e-7_dp, 0.2_dp, 0.3_dp, 800.e3_dp, 2.540e6_dp, 0.4_dp, 0.6_dp]

       call fluid%assign(fluid_data, offset)

       call source_flow_test("inject 1", 10._dp, 200.e3_dp, 1, 0, &
            [10._dp, 0._dp, 2.e6_dp])
       call source_flow_test("inject 2", 5._dp, 200.e3_dp, 2, 0, &
            [0._dp, 5._dp, 1.e6_dp])
       call source_flow_test("inject 3", 5._dp, 200.e3_dp, 2, 2, &
            [0._dp, 5._dp, 1.e6_dp])
       call source_flow_test("inject heat", 1000._dp, 0._dp, 3, 0, &
            [0._dp, 0._dp, 1000._dp])

       call source_flow_test("produce all", -5._dp, 0._dp, 0, 0, &
            [-3.4948610582_dp, -1.5051389418_dp, -431766.653977922_dp])
       call source_flow_test("produce 1", -5._dp, 0._dp, 0, 1, &
            [-5._dp, 0._dp, -431766.653977922_dp])
       call source_flow_test("produce heat", -5000._dp, 0._dp, 0, 3, &
            [0._dp, 0._dp, -5000._dp])

       call source_flow_test("no flow 1", 0._dp, 100.e3_dp, 1, 0, &
            [0._dp, 0._dp, 0._dp])
       call source_flow_test("no flow all", 0._dp, 100.e3_dp, 0, 0, &
            [0._dp, 0._dp, 0._dp])

       call fluid%destroy()
       deallocate(fluid_data)
       call eos%destroy()
       call thermo%destroy()
       call fson_destroy(json)

    end if

  contains

    subroutine source_flow_test(tag, rate, enthalpy, injection_component, &
         production_component, flow)
      !! Runs asserts for single flow update_source() test.

      character(*), intent(in) :: tag
      PetscReal, intent(in) :: rate, enthalpy
      PetscInt, intent(in) :: injection_component, production_component
      PetscReal, intent(in) :: flow(:)

      call source%init(0, 0, eos, rate, enthalpy, &
           injection_component, production_component)
      call source%update_flow(fluid_data, offset)
      call assert_equals(flow, source%flow, &
           eos%num_primary_variables, tol, "Source update_flow() " // trim(tag))
      call source%destroy()

    end subroutine source_flow_test

  end subroutine test_source_update_flow

!------------------------------------------------------------------------

end module source_test
