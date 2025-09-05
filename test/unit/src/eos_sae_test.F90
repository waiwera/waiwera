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
  use unit_test_utils_module, only: fluid_compare

  implicit none
  private

  public :: setup, teardown, setup_test
  public :: test_eos_sae_fluid_properties, test_eos_sae_bdy_consistency

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

         title = 'region 3 subcritical liquid, Pg = 0'
         primary = [650._dp, 360._dp, 0._dp]
         region = 3
         expected = 0._dp
         bulk = [40.45712306840381e6_dp, 360._dp, 3._dp, 0._dp, 1._dp, &
              1._dp, 1._dp, 0._dp, 40.45712306840381e6_dp, 0._dp]
         phase1 = [650._dp, 7.663779619075069e-5_dp, 1._dp, &
              1._dp, 0._dp, 1646.6695333441563e3_dp, 1584.427805546612e3_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 subcritical vapour, Pg = 0'
         primary = [125._dp, 360._dp, 0._dp]
         region = 3
         expected = 0._dp
         bulk = [18.06931379044633e6_dp, 360._dp, 3._dp, 0._dp, 2._dp, &
              1._dp, 0._dp, 0._dp, 18.06931379044633e6_dp, 0._dp]
         phase2 = [125._dp, 2.4792120096153134e-5_dp, 1._dp, &
              1._dp, 0._dp, 2558.815868019601e3_dp, 2414.2613576960303e3_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 liquidlike supercritical, Pg = 0'
         primary = [500._dp, 400._dp, 0._dp]
         region = 3
         expected = 0._dp
         bulk = [37.23729144851866e6_dp, 400._dp, 3._dp, 0._dp, 4._dp, &
              1._dp, 1._dp, 1._dp, 37.23729144851866e6_dp, 0._dp]
         phase3 = [500._dp, 5.8837922820019035e-5_dp, 1._dp, &
              1._dp, 0._dp, 1958.054759408955e3_dp, 1883.5801765119177e3_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 vapourlike supercritical, Pg = 0'
         primary = [320._dp, 500._dp, 0._dp]
         region = 3
         expected = 0._dp
         bulk = [57.58529377353408e6_dp, 500._dp, 3._dp, 0._dp, 4._dp, &
              1._dp, 0._dp, 2._dp, 57.58529377353408e6_dp, 0._dp]
         phase3 = [320._dp, 4.5929888420621406e-5_dp, 1._dp, &
              1._dp, 0._dp, 2602.6723916362077e3_dp, 2422.7183485939135e3_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 2 supercritical steam, Pg = 0'
         primary = [40.e6_dp, 600._dp, 0._dp]
         region = 2
         expected = 0._dp
         bulk = [40.e6_dp, 600._dp, 2._dp, 0._dp, 4._dp, 1._dp, 0._dp, 2._dp, &
              40.e6_dp, 0._dp]
         phase3 = [123.62382847002988_dp, 3.696154052412127e-5_dp, 1._dp, &
              1._dp, 0._dp, 3350.4327456605723e3_dp, 3026.870528772125e3_dp, 1._dp, 0._dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 Widom delta, Pg = 0'
         primary = [420._dp, 460._dp, 0._dp]
         region = 3
         expected = 0._dp
         bulk = [55.956158699122846e6_dp, 460._dp, 3._dp, 0._dp, 4._dp, &
              1._dp, 0.75747923886852431_dp, 3._dp, 55.956158699122846e6_dp, 0._dp]
         phase3 = [420._dp, 5.314389040730653e-05_dp, 1._dp, &
              1._dp, 0._dp, 2296.8587930659675e3_dp, 2163.6298437823416e3_dp, 1._dp, 0._dp]
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

         title = 'region 3 subcritical liquid, Pg > 0'
         primary = [650._dp, 360._dp, 10.e5_dp]
         region = 3
         expected = 0._dp
         bulk = [41.45712306840381e6_dp, 360._dp, 3._dp, 0._dp, 1._dp, 1._dp, &
              1._dp, 0._dp, 40.45712306840381e6_dp, 10.e5_dp]
         phase1 = [653.8020093587036_dp, 7.671989656222635e-5_dp, 1._dp, 1._dp, 0._dp, &
              1656.3103248353535e3_dp, 1592.9010319990285e3_dp, &
              0.9961718022439375_dp, 0.003828197756062468_dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 subcritical vapour, Pg > 0'
         primary = [125._dp, 360._dp, 1.e6_dp]
         region = 3
         expected = 0._dp
         bulk = [19.06931379044633e6_dp, 360._dp, 3._dp, 0._dp, 2._dp, 1._dp, &
              0._dp, 0._dp, 18.06931379044633e6_dp, 1.e6_dp]
         phase2 = [130.5012060058691_dp, 2.517178462683078e-5_dp, 1._dp, 1._dp, &
              0._dp, 2466.5196061534435e3_dp, 2320.3959466557014e3_dp, &
              0.9578455542731021_dp, 0.042154445726897825_dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 liquidlike supercritical, Pg > 0'
         primary = [500._dp, 400._dp, 2.e6_dp]
         region = 3
         expected = 0._dp
         bulk = [39.23729144851866e6_dp, 400._dp, 3._dp, 0._dp, 4._dp, 1._dp, &
              1._dp, 1._dp, 37.23729144851866e6_dp, 2.e6_dp]
         phase3 = [510.3486253661622_dp, 5.784079879654913e-5_dp, 1._dp, 1._dp, 0._dp, &
              1926.7063450531405e3_dp, 1849.8230352952923e3_dp, &
              0.9797224390312851_dp, 0.020277560968714865_dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 vapourlike supercritical, Pg > 0'
         primary = [320._dp, 500._dp, 2.e6_dp]
         region = 3
         expected = 0._dp
         bulk = [59.58529377353408e6_dp, 500._dp, 3._dp, 0._dp, 4._dp, 1._dp, &
              0._dp, 2._dp, 57.58529377353408e6_dp, 2.e6_dp]
         phase3 = [329.01012373437504_dp, 4.5686238006410807e-5_dp, 1._dp, 1._dp, &
              0._dp, 2545.648992343294e3_dp, 2364.544249737744e3_dp, &
              0.9726144483576763_dp, 0.027385551642323733_dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 3 Widom delta, Pg > 0'
         primary = [420._dp, 460._dp, 8.e6_dp]
         region = 3
         expected = 0._dp
         bulk = [63.956158699122846e6_dp, 460._dp, 3._dp, 0._dp, 4._dp, 1._dp, &
              0.7574792388685243_dp, 3._dp, 55.956158699122846e6_dp, 8.e6_dp]
         phase3 = [458.00683170009995_dp, 5.0977682440528604e-5_dp, 1._dp, 1._dp, &
              0._dp, 2145.828286414576e3_dp, 2006.1881014380322e3_dp, &
              0.9170168891170895_dp, 0.08298311088291055_dp]
         expected_err = 0
         call properties_test(title, primary, region, expected, expected_err)

         title = 'region 2 supercritical steam, Pg > 0'
         primary = [40.e6_dp, 600._dp, 5.e6_dp]
         region = 2
         expected = 0._dp
         bulk = [40.e6_dp, 600._dp, 2._dp, 0._dp, 4._dp, 1._dp, 0._dp, 2._dp, &
              35.e6_dp, 5.e6_dp]
         phase3 = [124.95373736753187_dp, 3.764443853154568e-5_dp, 1._dp, 1._dp, &
              0._dp, 2957.148420310664e3_dp, 2637.029944123383e3_dp, &
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

  subroutine test_eos_sae_bdy_consistency(test)

    ! Test eos_sae boundary consistency

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(fson_value), pointer :: json
    type(IAPWS_type) :: thermo
    type(eos_sae_type) :: eos
    PetscReal, pointer, contiguous :: fluid_data(:), rock_data(:)
    type(fluid_type) :: fluid1, fluid2
    type(rock_type) :: rock
    class(relative_permeability_type), allocatable :: rp
    class(capillary_pressure_type), allocatable :: cp
    PetscMPIInt :: rank
    PetscInt :: ierr
    PetscReal, allocatable :: primary1(:), primary2(:)
    PetscReal :: props(2), t, d, Pw1

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

       ! region 1 / 3
       associate (P1 => primary1(1), T1 => primary1(2), &
            Pa1 => primary1(3), d2 => primary2(1), &
            T2 => primary2(2), Pa2 => primary2(3), &
            density => props(1))
         P1 = 40.e6_dp
         T1 = thermo%temperature_bdy_1_3
         Pa1 = 1.e5_dp
         call thermo%water%properties([P1 - Pa1, T1], props, err)
         d2 = density
         T2 = T1
         Pa2 = Pa1
       end associate
       fluid1%region = dble(1)
       fluid2%region = dble(3)
       call eos%fluid_properties(primary1, rock, fluid1, err)
       call eos%fluid_properties(primary2, rock, fluid2, err)
       call fluid_compare(test, fluid1, fluid2, "region 1/3")

       ! region 2 / 3
       associate (P1 => primary1(1), T1 => primary1(2), &
            Pa1 => primary1(3), d2 => primary2(1), &
            T2 => primary2(2), Pa2 => primary2(3), &
            density => props(1))
         T1 = 400._dp
         Pa1 = 2.e5_dp
         call thermo%boundary23%pressure(T1, Pw1)
         P1 = Pw1 + Pa1
         call thermo%steam%properties([Pw1, T1], props, err)
         d2 = density
         T2 = T1
         Pa2 = Pa1
       end associate
       fluid1%region = dble(2)
       fluid2%region = dble(3)
       call eos%fluid_properties(primary1, rock, fluid1, err)
       call eos%fluid_properties(primary2, rock, fluid2, err)
       call fluid_compare(test, fluid1, fluid2, "region 2/3")

       ! region 4 at 350 deg C
       associate (P1 => primary1(1), Sv1 => primary1(2), &
            Pa1 => primary1(3), P2 => primary2(1), &
            Sv2 => primary2(2), Pa2 => primary2(3))
         Pw1 = thermo%saturation_pressure_bdy_1_3
         Pa1 = 5.e5_dp
         P1 = Pw1 + Pa1
         Sv1 = 0.5_dp
         P2 = P1 + 1.e-9_dp
         Sv2 = Sv1
         Pa2 = Pa1
       end associate
       fluid1%region = dble(4)
       fluid2%region = dble(4)
       call eos%fluid_properties(primary1, rock, fluid1, err)
       call eos%fluid_properties(primary2, rock, fluid2, err)
       call fluid_compare(test, fluid1, fluid2, "region 4 350 C")

       ! region 4 / 3 liquid
       associate (P1 => primary1(1), Sv1 => primary1(2), &
            Pa1 => primary1(3), d2 => primary2(1), &
            T2 => primary2(2), Pa2 => primary2(3))
         Pw1 = 18.e6_dp
         call thermo%saturation%temperature(Pw1, t, err)
         Pa1 = 7.e5_dp
         P1 = Pw1 + Pa1
         Sv1 = 0._dp
         select type (region => thermo%supercritical)
         type is (IAPWS_region3_type)
            call region%saturation_density([Pw1, t], &
                 PETSC_TRUE, d, err, PETSC_TRUE)
         end select
         d2 = d + 1.e-9_dp
         T2 = t
         Pa2 = Pa1
       end associate
       fluid1%region = dble(4)
       fluid2%region = dble(3)
       call eos%fluid_properties(primary1, rock, fluid1, err)
       call eos%fluid_properties(primary2, rock, fluid2, err)
       call fluid_compare(test, fluid1, fluid2, "region 4/3 L")

       ! region 4 / 3 liquid, near-critical
       associate (P1 => primary1(1), Sv1 => primary1(2), &
            Pa1 => primary1(3), d2 => primary2(1), &
            T2 => primary2(2), Pa2 => primary2(3))
         Pw1 = 21.99e6_dp
         call thermo%saturation%temperature(Pw1, t, err)
         Pa1 = 2.e5_dp
         P1 = Pw1 + Pa1
         Sv1 = 0._dp
         select type (region => thermo%supercritical)
         type is (IAPWS_region3_type)
            call region%saturation_density([Pw1, t], &
                 PETSC_TRUE, d, err, PETSC_TRUE)
         end select
         d2 = d + 1.e-6_dp
         T2 = t
         Pa2 = Pa1
       end associate
       fluid1%region = dble(4)
       fluid2%region = dble(3)
       call eos%fluid_properties(primary1, rock, fluid1, err)
       call eos%fluid_properties(primary2, rock, fluid2, err)
       call fluid_compare(test, fluid1, fluid2, &
            "region 4/3 near-critical L", tol = 1.e-8_dp)

       ! region 4 / 3 vapour
       associate (P1 => primary1(1), Sv1 => primary1(2), &
            Pa1 => primary1(3), d2 => primary2(1), &
            T2 => primary2(2), Pa2 => primary2(3))
         Pw1 = 18.e6_dp
         Sv1 = 1._dp
         Pa1 = 1.5e5_dp
         P1 = Pw1 + Pa1
         call thermo%saturation%temperature(Pw1, t, err)
         select type (region => thermo%supercritical)
         type is (IAPWS_region3_type)
            call region%saturation_density([Pw1, t], &
                 PETSC_FALSE, d, err, PETSC_TRUE)
         end select
         d2 = d - 1.e-9_dp
         T2 = t
         Pa2 = Pa1
       end associate
       fluid1%region = dble(4)
       fluid2%region = dble(3)
       call eos%fluid_properties(primary1, rock, fluid1, err)
       call eos%fluid_properties(primary2, rock, fluid2, err)
       call fluid_compare(test, fluid1, fluid2, "region 4/3 V")

       ! region 4 / 3 vapour, near-critical
       associate (P1 => primary1(1), Sv1 => primary1(2), &
            Pa1 => primary1(3), d2 => primary2(1), &
            T2 => primary2(2), Pa2 => primary2(3))
         Pw1 = 21.99e6_dp
         Sv1 = 1._dp
         Pa1 = 1.3e5_dp
         P1 = Pw1 + Pa1
         call thermo%saturation%temperature(Pw1, t, err)
         select type (region => thermo%supercritical)
         type is (IAPWS_region3_type)
            call region%saturation_density([Pw1, t], &
                 PETSC_FALSE, d, err, PETSC_TRUE)
         end select
         d2 = d
         T2 = t + 1.e-9_dp
         Pa2 = Pa1
       end associate
       fluid1%region = dble(4)
       fluid2%region = dble(3)
       call eos%fluid_properties(primary1, rock, fluid1, err)
       call eos%fluid_properties(primary2, rock, fluid2, err)
       call fluid_compare(test, fluid1, fluid2, "region 4/3 near-critical V")

       ! region 4 liquid / vapour at critical point
       associate (P1 => primary1(1), Sv1 => primary1(2), &
            Pa1 => primary1(3), P2 => primary2(1), &
            Sv2 => primary2(2), Pa2 => primary2(3))
         Pw1 = thermo%critical%pressure - 1.e-9_dp
         Sv1 = 0._dp
         Pa1 = 1.e5_dp
         P1 = Pw1 + Pa1
         P2 = P1
         Sv2 = 1._dp
         Pa2 = Pa1
       end associate
       fluid1%region = dble(4)
       fluid2%region = dble(4)
       call eos%fluid_properties(primary1, rock, fluid1, err)
       call eos%fluid_properties(primary2, rock, fluid2, err)
       call test%assert(fluid1%phase(1)%density, fluid2%phase(2)%density, &
            "region 4 L/V CP density", tol = 1.e-10_dp)
       call test%assert(fluid1%phase(1)%viscosity, fluid2%phase(2)%viscosity, &
            "region 4 L/V CP viscosity", tol = 1.e-10_dp)
       call test%assert(fluid1%phase(1)%specific_enthalpy, &
            fluid2%phase(2)%specific_enthalpy, &
            "region 4 L/V CP enthalpy", tol = 1.e-10_dp)

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

  end subroutine test_eos_sae_bdy_consistency

!------------------------------------------------------------------------

end module eos_sae_test_module
