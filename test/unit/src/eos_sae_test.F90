module eos_sae_test_module

  ! Tests for eos_sae module (supercritical water and air NCG equation
  ! of state)

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
  use eos_sae_module

  implicit none
  private

  public :: setup, teardown, setup_test
  public :: test_eos_sae_fluid_properties

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

  subroutine test_eos_sae_fluid_properties(test)

    ! eos_sae fluid properties test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fluid_type) :: fluid
    type(rock_type) :: rock
    PetscInt,  parameter :: offset = 1, num_primary_variables = 3
    PetscReal, pointer, contiguous :: fluid_data(:)
    PetscReal, allocatable :: expected(:)
    PetscReal :: primary(num_primary_variables)
    PetscInt :: region, expected_err
    type(eos_sae_type) :: eos
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
            phase2 => expected(20: 28), phase3 => expected(29: 37))

         title = 'cold water, Pg = 0'
         primary = [1.e5_dp, 20._dp, 0._dp]
         region = 1
         expected = 0._dp
         bulk = [1.e5_dp, 20._dp, 1._dp, 0._dp, 1._dp, 1._dp, 1._dp, 0._dp, &
              1.e5_dp, 0._dp]
         phase1 = [998.2054863776967_dp, 0.0010015972622270245_dp, 1._dp, &
              1._dp, 0._dp, 84011.81116713627_dp, 83911.6313931672_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'steam, Pg = 0'
         primary = [1.e5_dp, 150._dp, 0._dp]
         region = 2
         expected = 0._dp
         bulk = [1.e5_dp, 150._dp, 2._dp, 0._dp, 2._dp, 1._dp, 0._dp, 0._dp, &
              1.e5_dp, 0._dp]
         phase2 = [0.5163351360139934_dp, 1.419241230472252e-5_dp, 1._dp, &
              1._dp, 0._dp, 2776591.815449922_dp, 2582919.153984385_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = '2-phase, Pg = 0'
         primary = [10.e5_dp, 0.4_dp, 0._dp]
         region = 4
         expected = 0._dp
         bulk = [10.e5_dp, 179.88563239146663_dp, 4._dp, 0._dp, 3._dp, &
              1._dp, 0.6_dp, 0._dp, 10.e5_dp, 0._dp]
         phase1 = [887.1274516747791_dp, 0.00015048492650911237_dp, 0.6_dp, &
              0.6_dp, 0._dp, 762682.8443354106_dp, 761555.6105900089_dp, 1._dp, 0._dp]
         phase2 = [5.145385853182684_dp, 1.4981316222701134e-5_dp, 0.4_dp, &
              0.4_dp, 0._dp, 2777119.5376846623_dp, 2582770.65335727_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 subcritical liquid, Pg = 0'
         primary = [650._dp, 360._dp, 0._dp]
         region = 3
         expected = 0._dp
         bulk = [40.45712306840381e6_dp, 360._dp, 3._dp, 0._dp, 1._dp, &
              1._dp, 1._dp, 0._dp, 40.45712306840381e6_dp, 0._dp]
         phase1 = [650._dp, 7.663779619075069e-5_dp, 1._dp, &
              1._dp, 0._dp, 1646669.5333441563_dp, 1584427.805546612_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 subcritical vapour, Pg = 0'
         primary = [125._dp, 360._dp, 0._dp]
         region = 3
         expected = 0._dp
         bulk = [18.06931379044633e6_dp, 360._dp, 3._dp, 0._dp, 2._dp, &
              1._dp, 0._dp, 0._dp, 18.06931379044633e6_dp, 0._dp]
         phase2 = [125._dp, 2.4792120096153134e-5_dp, 1._dp, &
              1._dp, 0._dp, 2558815.868019601_dp, 2414261.3576960303_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 liquidlike supercritical, Pg = 0'
         primary = [500._dp, 400._dp, 0._dp]
         region = 3
         expected = 0._dp
         bulk = [37.23729144851866e6_dp, 400._dp, 3._dp, 0._dp, 4._dp, &
              1._dp, 1._dp, 1._dp, 37.23729144851866e6_dp, 0._dp]
         phase3 = [500._dp, 5.8837922820019035e-5_dp, 1._dp, &
              1._dp, 0._dp, 1958054.759408955_dp, 1883580.1765119177_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 vapourlike supercritical, Pg = 0'
         primary = [320._dp, 500._dp, 0._dp]
         region = 3
         expected = 0._dp
         bulk = [57.58529377353408e6_dp, 500._dp, 3._dp, 0._dp, 4._dp, &
              1._dp, 0._dp, 2._dp, 57.58529377353408e6_dp, 0._dp]
         phase3 = [320._dp, 4.5929888420621406e-5_dp, 1._dp, &
              1._dp, 0._dp, 2602672.3916362077_dp, 2422718.3485939135_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 2 supercritical steam, Pg = 0'
         primary = [40.e6_dp, 600._dp, 0._dp]
         region = 2
         expected = 0._dp
         bulk = [40.e6_dp, 600._dp, 2._dp, 0._dp, 4._dp, 1._dp, 0._dp, 2._dp, &
              40.e6_dp, 0._dp]
         phase3 = [123.62382847002988_dp, 3.696154052412127e-5_dp, 1._dp, &
              1._dp, 0._dp, 3350432.7456605723_dp, 3026870.528772125_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 Widom delta, Pg = 0'
         primary = [420._dp, 460._dp, 0._dp]
         region = 3
         expected = 0._dp
         bulk = [55.956158699122846e6_dp, 460._dp, 3._dp, 0._dp, 4._dp, &
              1._dp, 0.75747923886852431_dp, 3._dp, 55.956158699122846e6_dp, 0._dp]
         phase3 = [420._dp, 5.314389040730653e-05_dp, 1._dp, &
              1._dp, 0._dp, 2296858.7930659675_dp, 2163629.8437823416_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'liquid water, Pg > 0'
         primary = [10.e5_dp, 90._dp, 1.e5_dp]
         region = 1
         expected = 0._dp
         bulk = [10.e5_dp, 90._dp, 1._dp, 0._dp, 1._dp, 1._dp, 1._dp, 0._dp, &
              9.e5_dp, 1.e5_dp]
         phase1 = [965.7286048999848_dp, 0.00031442392084784024_dp, 1._dp, &
              1._dp, 0._dp, 377683.87768277636_dp, 376648.3900779502_dp, &
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
              1._dp, 0._dp, 1162779.5459599963_dp, 1013466.5843495822_dp, &
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
              0.6_dp, 0._dp, 721122.7402064757_dp, 720008.1096523397_dp, &
              0.999959224263865_dp, 4.0775736135034536e-05_dp]
         phase2 = [5.73149097766101_dp, 1.7406141056041412e-05_dp, 0.4_dp, &
              0.4_dp, 0._dp, 2056697.5008786342_dp, 1882222.8303439359_dp, &
              0.725987048943491_dp, 0.27401295105650897_dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 subcritical liquid, Pg > 0'
         primary = [650._dp, 360._dp, 10.e5_dp]
         region = 3
         expected = 0._dp
         bulk = [41.45712306840381e6_dp, 360._dp, 3._dp, 0._dp, 1._dp, 1._dp, &
              1._dp, 0._dp, 40.45712306840381e6_dp, 10.e5_dp]
         phase1 = [650._dp, 7.663779619075069e-5_dp, 1._dp, 1._dp, 0._dp, &
              1670793.4565923912_dp, 1607013.2672563854_dp, &
              0.995776602363838_dp, 0.004223397636161937_dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 subcritical vapour, Pg > 0'
         primary = [125._dp, 360._dp, 1.e6_dp]
         region = 3
         expected = 0._dp
         bulk = [19.06931379044633e6_dp, 360._dp, 3._dp, 0._dp, 2._dp, 1._dp, &
              0._dp, 0._dp, 18.06931379044633e6_dp, 1.e6_dp]
         phase2 = [130.5012060058691_dp, 2.517178462683078e-5_dp, 1._dp, 1._dp, &
              0._dp, 2466519.6061534435_dp, 2320395.9466557014_dp, &
              0.9578455542731021_dp, 0.042154445726897825_dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 liquidlike supercritical, Pg > 0'
         primary = [500._dp, 400._dp, 2.e6_dp]
         region = 3
         expected = 0._dp
         bulk = [39.23729144851866e6_dp, 400._dp, 3._dp, 0._dp, 4._dp, 1._dp, &
              1._dp, 1._dp, 37.23729144851866e6_dp, 2.e6_dp]
         phase3 = [500._dp, 5.8837922820019035e-5_dp, 1._dp, 1._dp, 0._dp, &
              2157696.3443172574_dp, 2079221.7614202201_dp, &
              0.646883470490572_dp, 0.35311652950942796_dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 vapourlike supercritical, Pg > 0'
         primary = [320._dp, 500._dp, 2.e6_dp]
         region = 3
         expected = 0._dp
         bulk = [59.58529377353408e6_dp, 500._dp, 3._dp, 0._dp, 4._dp, 1._dp, &
              0._dp, 2._dp, 57.58529377353408e6_dp, 2.e6_dp]
         phase3 = [329.01012373437504_dp, 4.5686238006410807e-5_dp, 1._dp, 1._dp, &
              0._dp, 2545648.992343294_dp, 2364544.249737744_dp, &
              0.9726144483576763_dp, 0.027385551642323733_dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 Widom delta, Pg > 0'
         primary = [420._dp, 460._dp, 8.e6_dp]
         region = 3
         expected = 0._dp
         bulk = [63.956158699122846e6_dp, 460._dp, 3._dp, 0._dp, 4._dp, 1._dp, &
              0.7574792388685243_dp, 3._dp, 55.956158699122846e6_dp, 8.e6_dp]
         phase3 = [429.21744575210414_dp, 5.261854000243448e-5_dp, 1._dp, 1._dp, &
              0._dp, 2262231.824241309_dp, 2113225.395606857_dp, &
              0.9783738482549783_dp, 0.021626151745021704_dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 2 supercritical steam, Pg > 0'
         primary = [40.e6_dp, 600._dp, 5.e6_dp]
         region = 2
         expected = 0._dp
         bulk = [40.e6_dp, 600._dp, 2._dp, 0._dp, 4._dp, 1._dp, 0._dp, 2._dp, &
              35.e6_dp, 5.e6_dp]
         phase3 = [124.95373736753187_dp, 3.764443853154568e-5_dp, 1._dp, 1._dp, &
              0._dp, 2957148.420310664_dp, 2637029.944123383_dp, &
              0.8403766521370382_dp, 0.15962334786296184_dp]
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

  end subroutine test_eos_sae_fluid_properties

!------------------------------------------------------------------------

end module eos_sae_test_module
