module eos_test
  !! Tests for eos module.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use fson
  use fson_mpi_module
  use eos_module, only: max_component_name_length
  use eos_wge_module

  implicit none
  private

  public :: setup, teardown
  public :: test_eos_component_index, test_scaling

contains

!------------------------------------------------------------------------

  subroutine setup()

    use profiling_module, only: init_profiling

    ! Locals:
    PetscErrorCode :: ierr

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr); CHKERRQ(ierr)
    call init_profiling()

  end subroutine setup

!------------------------------------------------------------------------

  subroutine teardown()

    PetscErrorCode :: ierr

    call PetscFinalize(ierr); CHKERRQ(ierr)

  end subroutine teardown

!------------------------------------------------------------------------

  subroutine test_eos_component_index(test)
    ! eos component_index() test

    use IAPWS_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(IAPWS_type) :: thermo
    type(eos_wge_type) :: eos
    type(fson_value), pointer :: json
    character(max_component_name_length) :: name
    PetscMPIInt :: rank
    PetscInt :: ierr

    json => fson_parse_mpi(str = "{}")
    call thermo%init()
    call eos%init(json, thermo)

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       name = "water"
       call test%assert(1, eos%component_index(name), name)

       name = "Water"
       call test%assert(1, eos%component_index(name), name)

       name = "gas"
       call test%assert(2, eos%component_index(name), name)

       name = "energy"
       call test%assert(eos%num_primary_variables, &
            eos%component_index(name), name)

       name = "fred"
       call test%assert(-1, eos%component_index(name), name)

    end if

    call eos%destroy()
    call thermo%destroy()

  end subroutine test_eos_component_index

!------------------------------------------------------------------------

  subroutine test_scaling(test)
    ! eos scale/unscale

    use eos_module
    use eos_setup_module
    use IAPWS_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '{"eos": {"name": "w"}}')
    call scale_test(json, [3.e5_dp], 1, [0.3_dp], "w region 1")
    call scale_test(json, [4.e5_dp], 2, [0.4_dp], "w region 2")
    call fson_destroy(json)

    json => fson_parse_mpi(str = '{"eos": {"name": "w", "primary": {"scale": {"pressure": 1.e5}}}}')
    call scale_test(json, [3.e5_dp], 1, [3.0_dp], "w region 1 scale")
    call fson_destroy(json)

    json => fson_parse_mpi(str = '{"eos": {"name": "we"}}')
    call scale_test(json, [1.e5_dp, 20._dp], 1, [0.1_dp, 0.2_dp], "we region 1")
    call scale_test(json, [0.9e5_dp, 100._dp], 2, [0.09_dp, 1._dp], "we region 2")
    call scale_test(json, [13.e5_dp, 0.4_dp], 4, [1.3_dp, 0.4_dp], "we region 4")
    call fson_destroy(json)

    json => fson_parse_mpi(str = '{"eos": {"name": "we", ' // &
         '"primary": {"scale": {"pressure": 1.e7, "temperature": 200}}}}')
    call scale_test(json, [1.e5_dp, 20._dp], 1, [0.01_dp, 0.1_dp], "we region 1 scale")
    call scale_test(json, [0.9e5_dp, 100._dp], 2, [0.009_dp, 0.5_dp], "we region 2 scale")
    call scale_test(json, [13.e5_dp, 0.4_dp], 4, [0.13_dp, 0.4_dp], "we region 4 scale")
    call fson_destroy(json)

    json => fson_parse_mpi(str = '{"eos": {"name": "wce"}}')
    call scale_test(json, [15.e5_dp, 40._dp, 3.e5_dp], 1, [1.5_dp, 0.4_dp, 0.2_dp], "wce region 1")
    call scale_test(json, [0.8e5_dp, 110._dp, 0.6e5_dp], 2, [0.08_dp, 1.1_dp, 0.75_dp], "wce region 2")
    call scale_test(json, [20.e5_dp, 0.4_dp, 10.e5_dp], 4, [2.0_dp, 0.4_dp, 0.5_dp], "wce region 4")
    call fson_destroy(json)

    json => fson_parse_mpi(str = '{"eos": {"name": "wce", ' // &
         '"primary": {"scale": {"pressure": 1.e7, "temperature": 200, ' // &
         '"partial_pressure": 1.e5}}}}')
    call scale_test(json, [15.e5_dp, 40._dp, 2.e5_dp], 1, [0.15_dp, 0.2_dp, 2._dp], &
         "wce region 1 scale")
    call scale_test(json, [0.7e5_dp, 110._dp, 0.6e5_dp], 2, [0.007_dp, 0.55_dp, 0.6_dp], &
         "wce region 2 scale")
    call scale_test(json, [13.e5_dp, 0.4_dp, 10.e5_dp], 4, [0.13_dp, 0.4_dp, 10.0_dp], &
         "wce region 4 scale")
    call fson_destroy(json)

    json => fson_parse_mpi(str = '{"eos": {"name": "wae", ' // &
         '"primary": {"scale": {"pressure": 1.e7, "temperature": 200, ' // &
         '"partial_pressure": 1.e5}}}}')
    call scale_test(json, [15.e5_dp, 40._dp, 2.e5_dp], 1, [0.15_dp, 0.2_dp, 2._dp], &
         "wae region 1 scale")
    call scale_test(json, [0.7e5_dp, 110._dp, 0.6e5_dp], 2, [0.007_dp, 0.55_dp, 0.6_dp], &
         "wae region 2 scale")
    call scale_test(json, [13.e5_dp, 0.4_dp, 10.e5_dp], 4, [0.13_dp, 0.4_dp, 10.0_dp], &
         "wae region 4 scale")
    call fson_destroy(json)

    json => fson_parse_mpi(str = '{"eos": {"name": "wce", ' // &
         '"primary": {"scale": {"partial_pressure": "pressure"}}}}')
    call scale_test(json, [15.e5_dp, 40._dp, 3.e5_dp], 1, [1.5_dp, 0.4_dp, 0.2_dp], &
         "wce region 1 adaptive scale")
    call scale_test(json, [0.8e5_dp, 100._dp, 0.6e5_dp], 2, [0.08_dp, 1.0_dp, 0.75_dp], &
         "wce region 2 adaptive scale")
    call scale_test(json, [100.e5_dp, 0.5_dp, 50.e5_dp], 4, [10._dp, 0.5_dp, 0.5_dp], &
         "wce region 4 adaptive scale")

  contains

    subroutine scale_test(json, primary, region, expected_scaled_primary, title)

      type(fson_value), pointer, intent(in) :: json
      PetscReal, intent(in) :: primary(:)
      PetscInt, intent(in) :: region
      PetscReal, intent(in) :: expected_scaled_primary(:)
      character(*), intent(in) :: title
      ! Locals:
      type(IAPWS_type) :: thermo
      class(eos_type), allocatable :: eos
      PetscReal, allocatable :: scaled_primary(:), unscaled_primary(:)

      call thermo%init()
      call setup_eos(json, thermo, eos)
      associate (n => size(primary))
        allocate(scaled_primary(n), unscaled_primary(n))
      end associate

      scaled_primary = eos%scale(primary, region)
      call test%assert(expected_scaled_primary, scaled_primary, title // ": scaled")

      unscaled_primary = eos%unscale(scaled_primary, region)
      call test%assert(primary, unscaled_primary, title // ": unscaled")

      call eos%destroy()
      deallocate(eos, scaled_primary, unscaled_primary)
      call thermo%destroy()

    end subroutine scale_test

  end subroutine test_scaling

!------------------------------------------------------------------------

end module eos_test
