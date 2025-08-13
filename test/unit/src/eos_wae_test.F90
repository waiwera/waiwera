module eos_wae_test_module

  ! Tests for eos_wae module (non-isothermal water and air NCG
  ! equation of state)

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use zofu
  use fluid_module
  use rock_module
  use relative_permeability_module
  use capillary_pressure_module
  use IAPWS_module
  use fson
  use fson_mpi_module
  use eos_wae_module
  use unit_test_utils_module, only: fluid_compare

  implicit none
  private

  public :: setup, teardown, setup_test
  public :: test_eos_wae_fluid_properties, test_eos_wae_bdy_consistency

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

  subroutine setup_test(test)

    class(unit_test_type), intent(in out) :: test

    test%tolerance = 1.e-9

  end subroutine setup_test

!------------------------------------------------------------------------

  subroutine test_eos_wae_fluid_properties(test)

    ! eos_wae fluid properties test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: fluid
    type(rock_type) :: rock
    PetscInt,  parameter :: offset = 1, num_primary_variables = 3
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal, allocatable :: expected(:)
    PetscReal :: primary(num_primary_variables)
    PetscInt :: region, expected_err
    type(eos_wae_type) :: eos
    type(IAPWS_type) :: thermo
    class(relative_permeability_type), allocatable :: rp
    class(capillary_pressure_type), allocatable :: cp
    type(fson_value), pointer :: json
    character(2) :: json_str = '{}'
    character(60) :: title
    PetscErrorCode :: err
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, parameter :: small = 1.e-6_dp

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = json_str)
    call thermo%init()
    call eos%init(json, thermo)
    call setup_relative_permeabilities(json, rp)
    call setup_capillary_pressures(json, cp)
    call rock%assign_relative_permeability(rp)
    call rock%assign_capillary_pressure(cp)
    call fluid%init(eos%num_components, eos%num_phases)
    allocate(fluid_data(fluid%dof), expected(fluid%dof))
    fluid_data = 0._dp
    call fluid%assign(fluid_data, offset)
    call rock%init()

    if (rank == 0) then
       associate (bulk => expected(1: 10), phase1 => expected(11: 19), &
            phase2 => expected(20: 28))

         title = 'cold water, Pg = 0'
         primary = [1.e5_dp, 20._dp, 0._dp]
         region = 1
         expected = 0._dp
         bulk = [1.e5_dp, 20._dp, 1._dp, 0._dp, 1._dp, 1._dp, 1._dp, 0._dp, &
              1.e5_dp, 0._dp]
         phase1 = [998.2054863776967_dp, 0.0010015972622270245_dp, 1._dp, &
              1._dp, 0._dp, 84.01181116713627e3_dp, 83.9116313931672e3_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'steam, Pg = 0'
         primary = [1.e5_dp, 150._dp, 0._dp]
         region = 2
         expected = 0._dp
         bulk = [1.e5_dp, 150._dp, 2._dp, 0._dp, 2._dp, 1._dp, 0._dp, 0._dp, &
              1.e5_dp, 0._dp]
         phase2 = [0.5163351360139934_dp, 1.419241230472252e-5_dp, 1._dp, &
              1._dp, 0._dp, 2776.591815449922e3_dp, 2582.919153984385e3_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = '2-phase, Pg = 0'
         primary = [10.e5_dp, 0.4_dp, 0._dp]
         region = 4
         expected = 0._dp
         bulk = [10.e5_dp, 179.88563239146663_dp, 4._dp, 0._dp, 3._dp, &
              1._dp, 0.6_dp, 0._dp, 10.e5_dp, 0._dp]
         phase1 = [887.1274516747791_dp, 0.00015048492650911237_dp, 0.6_dp, &
              0.6_dp, 0._dp, 762.6828443354106e3_dp, 761.5556105900089e3_dp, 1._dp, 0._dp]
         phase2 = [5.145385853182684_dp, 1.4981316222701134e-5_dp, 0.4_dp, &
              0.4_dp, 0._dp, 2777.1195376846623e3_dp, 2582.77065335727e3_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'liquid water, Pg > 0'
         primary = [10.e5_dp, 90._dp, 1.e5_dp]
         region = 1
         expected = 0._dp
         bulk = [10.e5_dp, 90._dp, 1._dp, 0._dp, 1._dp, 1._dp, 1._dp, 0._dp, &
              9.e5_dp, 1.e5_dp]
         phase1 = [965.7286048999848_dp, 0.00031442392084784024_dp, 1._dp, &
              1._dp, 0._dp, 377.68387768277636e3_dp, 376.6483900779502e3_dp, &
              0.9999858156956914_dp, 1.4184304308575404e-5_dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'steam, Pg > 0'
         primary = [2.e5_dp, 150._dp, 1.e5_dp]
         region = 2
         expected = 0._dp
         bulk = [2.e5_dp, 150._dp, 2._dp, 0._dp, 2._dp, 1._dp, 0._dp, 0._dp, &
              1.e5_dp, 1.e5_dp]
         phase2 = [1.3394684416068174_dp, 2.054627278536341e-05_dp, 1._dp, &
              1._dp, 0._dp, 1162.7795459599963e3_dp, 1013.4665843495822e3_dp, &
              0.3854776417088268_dp, 0.6145223582911732_dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = '2-phase, Pg > 0'
         primary = [10.e5_dp, 0.4_dp, 2.e5_dp]
         region = 4
         expected = 0._dp
         bulk = [10.e5_dp, 170.41351081360017_dp, 4._dp, 0._dp, 3._dp, &
              1._dp, 0.6_dp, 0._dp, 8.e5_dp, 2.e5_dp]
         phase1 = [897.1582523817934_dp, 0.00015941397404013526_dp, 0.6_dp, &
              0.6_dp, 0._dp, 721.1227402064757e3_dp, 720.0081096523397e3_dp, &
              0.999959224263865_dp, 4.0775736135034536e-05_dp]
         phase2 = [5.73149097766101_dp, 1.7406141056041412e-05_dp, 0.4_dp, &
              0.4_dp, 0._dp, 2056.6975008786342e3_dp, 1882.2228303439359e3_dp, &
              0.725987048943491_dp, 0.27401295105650897_dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

       end associate
    end if

    call rock%destroy()
    call fluid%destroy()
    deallocate(fluid_data, expected)
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)
    deallocate(rp, cp)

  contains

    subroutine properties_test(title, primary, region, expected, expected_err)

      character(*), intent(in) :: title
      PetscReal, intent(in) :: primary(:)
      PetscInt, intent(in) :: region
      PetscReal, intent(in) :: expected(:)
      PetscInt, intent(in) :: expected_err

      fluid%region = dble(region)
      call eos%fluid_properties(primary, rock, fluid, err)

      call test%assert(expected_err, err, trim(title) // ' err')
      if (err == 0) then
         call test%assert(expected, fluid_data(offset: offset + fluid%dof - 1), &
              trim(title) // ' data')
      end if

    end subroutine properties_test

  end subroutine test_eos_wae_fluid_properties

!------------------------------------------------------------------------

  subroutine test_eos_wae_bdy_consistency(test)

    ! Test eos_wae boundary consistency

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_wae_type) :: eos
    PetscReal, pointer, contiguous :: fluid_data(:), rock_data(:)
    type(fluid_type) :: fluid1, fluid2
    type(rock_type) :: rock
    class(relative_permeability_type), allocatable :: rp
    class(capillary_pressure_type), allocatable :: cp
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, allocatable :: primary1(:), primary2(:)
    PetscReal :: Pw
    PetscErrorCode :: err

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '{}')
    call thermo%init()
    call eos%init(json, thermo)
    allocate(primary1(eos%num_primary_variables))
    allocate(primary2(eos%num_primary_variables))
    call fluid1%init(eos%num_components, eos%num_phases)
    call fluid2%init(eos%num_components, eos%num_phases)
    call rock%init()
    allocate(fluid_data(fluid1%dof + fluid2%dof))
    allocate(rock_data(rock%dof))
    call setup_relative_permeabilities(json, rp)
    call setup_capillary_pressures(json, cp)

    fluid_data = 0._dp
    rock_data = 0._dp

    call fluid1%assign(fluid_data, 1)
    call fluid2%assign(fluid_data, fluid1%dof + 1)
    call rock%assign(rock_data, 1)
    call rock%assign_relative_permeability(rp)
    call rock%assign_capillary_pressure(cp)

    if (rank == 0) then

       ! region 1 / 4
       associate (P1 => primary1(1), T1 => primary1(2), &
            Pa1 => primary1(3), P2 => primary2(1), &
            Sv2 => primary2(2), Pa2 => primary2(3))
         T1 = 200._dp
         Pa1 = 0.1e5_dp
         call thermo%saturation%pressure(T1, Pw, err)
         P1 = Pw + Pa1
         P2 = P1
         Sv2 = 0._dp
         Pa2 = Pa1
       end associate
       fluid1%region = dble(1)
       fluid2%region = dble(4)
       call eos%fluid_properties(primary1, rock, fluid1, err)
       call eos%fluid_properties(primary2, rock, fluid2, err)
       call fluid_compare(test, fluid1, fluid2, "region 1/4")

       ! region 2 / 4
       associate (P1 => primary1(1), T1 => primary1(2), &
            Pa1 => primary1(3), P2 => primary2(1), &
            Sv2 => primary2(2), Pa2 => primary2(3))
         T1 = 300._dp
         Pa1 = 0.2e5_dp
         call thermo%saturation%pressure(T1, Pw, err)
         P1 = Pw + Pa1
         P2 = P1
         Sv2 = 1._dp
         Pa2 = Pa1
       end associate
       fluid1%region = dble(2)
       fluid2%region = dble(4)
       call eos%fluid_properties(primary1, rock, fluid1, err)
       call eos%fluid_properties(primary2, rock, fluid2, err)
       call fluid_compare(test, fluid1, fluid2, "region 2/4")

    end if

    call fluid1%destroy()
    call fluid2%destroy()
    call rock%destroy()
    call eos%destroy()
    call thermo%destroy()
    call fson_destroy_mpi(json)
    deallocate(fluid_data, primary1, primary2)
    call rp%destroy()
    deallocate(rp)
    call cp%destroy()
    deallocate(cp)

  end subroutine test_eos_wae_bdy_consistency

!------------------------------------------------------------------------

end module eos_wae_test_module
