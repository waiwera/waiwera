module face_test

  ! Test for face module

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use zofu
  use face_module
  use IAPWS_module
  use fson

  implicit none
  private

  public :: setup, teardown
  public :: test_face_assign, test_face_distances, &
       test_face_permeability_direction, test_face_normal_gradient, &
       test_face_harmonic_average, test_face_flux_zero_horizontal, &
       test_face_flux_vertical_gravity, test_face_flux_hydrostatic, &
       test_face_flux_two_phase_vertical

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

  subroutine test_face_assign(test)

    ! Face assign() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(face_type) :: face
    PetscReal, parameter :: area = 300._dp
    PetscReal, parameter :: distance(2) = [20._dp, 30._dp]
    PetscReal, parameter :: normal(3) = [0.5_dp, -0.25_dp, 0.75_dp]
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscReal, parameter :: gravity_normal = dot_product(gravity, normal)
    PetscReal, parameter :: centroid(3) = [-1250._dp, 3560._dp, -2530._dp]
    PetscReal, parameter :: permeability_direction = dble(2)
    PetscInt, parameter :: offset = 6
    PetscReal :: offset_padding(offset-1) = 0._dp
    PetscReal, pointer, contiguous :: face_data(:)
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call face%init()

       allocate(face_data(offset - 1 + face%dof))
       face_data = [offset_padding, area, distance, sum(distance), normal, &
            gravity_normal, centroid, permeability_direction]

       call test%assert(face%dof, size(face_data) - (offset-1), "face dof")

       call face%assign_geometry(face_data, offset)

       call test%assert(area, face%area, "area")
       call test%assert(distance, face%distance, "distances")
       call test%assert(normal, face%normal, "normal")
       call test%assert(-7.35_dp, face%gravity_normal, "gravity normal")
       call test%assert(centroid, face%centroid, "centroid")
       call test%assert(permeability_direction, face%permeability_direction, &
            "permeability direction")

       call face%destroy()
       deallocate(face_data)

    end if

  end subroutine test_face_assign

!------------------------------------------------------------------------

  subroutine test_face_distances(test)

    ! Face distance test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(face_type) :: face
    PetscInt, parameter :: num_tests = 2
    PetscReal, parameter :: centroid(3, num_tests) = reshape([ &
         0._dp, 200._dp, 50._dp, &
         0._dp, 200._dp, 50._dp &
         ], [3, num_tests])
    PetscReal, parameter :: normal(3, num_tests) = reshape([ &
         1._dp, 0._dp, 0._dp, &
         -1._dp, 0._dp, 0._dp &
         ], [3, num_tests])
    PetscReal, parameter :: expected_distances(2, num_tests) = reshape([ &
         80._dp, 100._dp, &
         -80._dp, -100._dp &
         ], [2, num_tests])
    PetscReal, pointer, contiguous :: cell_data(:), face_data(:)
    PetscInt :: i
    PetscInt, parameter :: offset = 1
    character(len = 32) :: msg
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       allocate(cell_data(8))
       cell_data = [-80._dp, 200._dp, 50._dp, 0._dp, &
            100._dp, 200._dp, 50._dp, 0._dp]

       do i = 1, num_tests

          call face%init()
          allocate(face_data(offset - 1 + face%dof))
          face_data = 0._dp
          face_data(offset + 4: offset + 6) = normal(:, i)
          face_data(offset + 8: offset + 10) = centroid(:, i)
          call face%assign_geometry(face_data, offset)
          call face%assign_cell_geometry(cell_data, [1, 5])

          call face%calculate_distances()
          write(msg, '(a, i2)') "Face distances", i
          call test%assert(expected_distances(:, i), face%distance, msg)
          write(msg, '(a, i2)') "Face distance12", i
          call test%assert(sum(expected_distances(:, i)), face%distance12, msg)

          call face%destroy()
          deallocate(face_data)

       end do

       deallocate(cell_data)

    end if

  end subroutine test_face_distances

!------------------------------------------------------------------------

  subroutine test_face_permeability_direction(test)

    ! Face permeability_direction() test

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(face_type) :: face
    PetscReal, parameter :: area = 10._dp
    PetscReal, parameter :: distance(2) = [10._dp, 10._dp]
    PetscReal, parameter :: centroid(3) = [0._dp, 0._dp, 0._dp]
    PetscReal, parameter :: initial_permeability_direction = dble(0)
    PetscInt,  parameter :: num_tests = 3
    PetscReal, parameter :: normal(3, num_tests) = reshape( &
         [  1._dp, -0.25_dp,  0.3_dp, &
          -0.1_dp,   2.1_dp,  0.5_dp,&
            1._dp,  -1.4_dp, -1.6_dp], [3, num_tests])
    PetscReal, parameter :: expected_permeability_direction(num_tests) = &
         [dble(1), dble(2), dble(3)]
    PetscReal, parameter :: gravity(3) = [0._dp, 0._dp, -9.8_dp]
    PetscReal :: gravity_normal
    PetscReal, pointer, contiguous :: face_data(:)
    PetscInt :: i
    PetscInt :: offset = 1
    PetscReal :: rotation(3, 3) = reshape([ &
         1._dp, 0._dp, 0._dp, &
         0._dp, 1._dp, 0._dp, &
         0._dp, 0._dp, 1._dp], [3, 3])

    character(len = 32) :: msg
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       do i = 1, num_tests

          call face%init()
          gravity_normal = dot_product(gravity, normal(:,i))
          allocate(face_data(offset - 1 + face%dof))
          face_data = [area, distance, sum(distance), normal(:,i), &
               gravity_normal, centroid, initial_permeability_direction]
          call face%assign_geometry(face_data, offset)
          call face%calculate_permeability_direction(rotation)

          write(msg, '(a, i2)') "Permeability direction test ", i
          call test%assert(expected_permeability_direction(i), &
               face%permeability_direction, trim(msg))

          call face%destroy()
          deallocate(face_data)

       end do

    end if

  end subroutine test_face_permeability_direction

!------------------------------------------------------------------------

  subroutine test_face_normal_gradient(test)

    ! Face normal_gradient() test

    use cell_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(face_type) :: face
    type(cell_type) :: cell
    PetscReal, parameter :: distance(2) = [25._dp, 32._dp]
    PetscReal, parameter :: x(2) = [240._dp, 170._dp]
    PetscReal, parameter :: expected_d12 = 57._dp
    PetscReal, parameter :: expected_g = -1.22807017544_dp
    PetscReal, pointer, contiguous :: face_data(:)
    PetscReal, pointer, contiguous :: cell_data(:)
    PetscInt :: face_offset, cell_offsets(2)
    PetscReal :: g
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call cell%init(1,1) ! dummy argument values
       call face%init()

       face_offset = 1
       cell_offsets = [1, 1 + cell%dof]
       allocate(face_data(face_offset - 1 + face%dof), &
            cell_data(cell%dof * 2))
       face_data = 0._dp
       face_data(2:3) = distance
       face_data(4) = sum(distance)
       cell_data = 0._dp

       call face%assign_geometry(face_data, face_offset)
       call face%assign_cell_geometry(cell_data, cell_offsets)

       g = face%normal_gradient(x)

       call test%assert(expected_d12, face%distance12, "face distance12")
       call test%assert(expected_g, g, "face normal gradient")

       call cell%destroy()
       call face%destroy()
       deallocate(face_data, cell_data)

    end if

  end subroutine test_face_normal_gradient

!------------------------------------------------------------------------

  subroutine test_face_harmonic_average(test)

    ! Face harmonic_average() test

    use cell_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    type(face_type) :: face
    type(cell_type) :: cell
    PetscReal, pointer, contiguous :: face_data(:)
    PetscReal, pointer, contiguous :: cell_data(:)
    PetscInt :: face_offset, cell_offsets(2)
    PetscReal :: xh
    PetscInt,  parameter :: num_tests = 6
    PetscReal, parameter :: x(2, num_tests) = reshape( &
         [240._dp, 170._dp, &
          240._dp, 170._dp, &
          240._dp, 170._dp, &
            0._dp, 170._dp, &
          240._dp,   0._dp, &
            0._dp,   0._dp], [2, num_tests])
    PetscReal, parameter :: distance(2, num_tests) = reshape( &
         [25._dp, 32._dp, &
           0._dp, 10._dp, &
          22._dp,  0._dp, &
          25._dp, 32._dp, &
          25._dp, 32._dp, &
          25._dp, 32._dp], [2, num_tests])
    PetscReal, parameter :: expected_xh(num_tests) = &
         [194.937133277_dp, 170._dp, 240._dp, 0._dp, 0._dp, 0._dp]
    PetscInt :: i
    character(len = 32) :: msg
    PetscMPIInt :: rank
    PetscInt :: ierr

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call cell%init(1,1) ! dummy argument values
       call face%init()
       face_offset = 1
       cell_offsets = [1, 1 + cell%dof]
       allocate(face_data(face_offset - 1 + face%dof), &
            cell_data(cell%dof * 2))
       face_data = 0._dp
       cell_data = 0._dp

       do i = 1, num_tests
          face_data(2:3) = distance(:, i)
          face_data(4) = sum(distance(:, i))
          call face%assign_geometry(face_data, face_offset)
          call face%assign_cell_geometry(cell_data, cell_offsets)
          xh = face%harmonic_average(x(:, i))
          write(msg, '(a, i2)') "Face harmonic average test ", i
          call test%assert(expected_xh(i), xh, msg)
       end do

       call cell%destroy()
       call face%destroy()
       deallocate(face_data, cell_data)

    end if

  end subroutine test_face_harmonic_average

!------------------------------------------------------------------------

  subroutine test_face_flux_zero_horizontal(test)

    ! Face flux() test, 1-phase horizontal
    ! Fluid properties in the two cells are identical, so the fluxes
    ! should be zero.

    use cell_module
    use rock_module
    use fluid_module
    use eos_we_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: nc = 1, num_phases = 1, num_primary = 2
    type(face_type) :: face
    type(cell_type) :: cell
    type(rock_type) :: rock
    type(fluid_type) :: fluid
    type(eos_we_type) :: eos
    type(IAPWS_type) :: thermo
    type(fson_value), pointer :: json
    PetscReal, pointer, contiguous :: face_data(:), cell_data(:)
    PetscReal, pointer, contiguous :: rock_data(:), fluid_data(:)
    PetscReal, allocatable :: flux(:)
    PetscInt :: face_offset, cell_offsets(2)
    PetscInt :: rock_offsets(2), fluid_offsets(2)
    PetscReal, parameter :: expected_mass_flux = 0._dp
    PetscReal, parameter :: expected_heat_flux = 0._dp
    PetscMPIInt :: rank
    PetscInt :: ierr

    call thermo%init()
    json => fson_parse(str = '{}')
    call eos%init(json, thermo)

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call cell%init(nc, num_phases)
       call face%init(nc, num_phases)
       call fluid%init(nc, num_phases)
       call rock%init()
       allocate(face_data(face%dof), cell_data(cell%dof))
       allocate(rock_data(rock%dof), fluid_data(fluid%dof))
       allocate(flux(num_primary))
       face_offset = 1
       cell_offsets = [1, 1]
       rock_offsets = [1, 1]
       fluid_offsets = [1, 1]
       face_data = [0._dp,  25._dp, 35._dp,  60._dp, 1._dp, 0._dp, 0._dp, &
            0._dp, 0._dp, 0._dp, 0._dp, 1._dp]
       cell_data = 0._dp ! not needed
       rock_data = [ &
            1.e-14_dp, 2.e-14_dp, 3.e-15_dp,  2.5_dp,  2.5_dp, 0.1_dp, &
            2200._dp, 1000._dp]
       fluid_data = [ &
            1.e5_dp, 20._dp, 1._dp, 1._dp, 0._dp, &
            998.2_dp, 1.e-3_dp, 1._dp, 1._dp, 0._dp, &
            84011.8_dp, 83911.6_dp, 1._dp]

       call face%assign_geometry(face_data, face_offset)
       call face%assign_cell_geometry(cell_data, cell_offsets)
       call face%assign_cell_rock(rock_data, rock_offsets)
       call face%assign_cell_fluid(fluid_data, fluid_offsets)

       flux = face%flux(eos)

       call test%assert(num_primary, size(flux), "Flux array size")
       call test%assert(expected_mass_flux, flux(1), "Mass flux")
       call test%assert(expected_heat_flux, flux(2), "Heat flux")

       call cell%destroy()
       call face%destroy()
       call fluid%destroy()
       call rock%destroy()
       deallocate(face_data, cell_data, rock_data, fluid_data, flux)

    end if

    call eos%destroy()
    call fson_destroy(json)
    call thermo%destroy()

  end subroutine test_face_flux_zero_horizontal

!------------------------------------------------------------------------

  subroutine test_face_flux_vertical_gravity(test)

    ! Face flux() test, 1-phase vertical, gravity only
    ! Fluid properties in both cells are identical, so the only flow
    ! is from gravity.

    use cell_module
    use rock_module
    use fluid_module
    use eos_we_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: nc = 1, num_phases = 1, num_primary = 2
    type(face_type) :: face
    type(cell_type) :: cell
    type(rock_type) :: rock
    type(fluid_type) :: fluid
    type(eos_we_type) :: eos
    type(IAPWS_type) :: thermo
    type(fson_value), pointer :: json
    PetscReal, pointer, contiguous :: face_data(:), cell_data(:)
    PetscReal, pointer, contiguous :: rock_data(:), fluid_data(:)
    PetscReal, allocatable :: flux(:)
    PetscInt :: face_offset, cell_offsets(2)
    PetscInt :: rock_offsets(2), fluid_offsets(2)
    PetscReal, parameter :: expected_mass_flux = 2.9294255256e-5_dp
    PetscReal, parameter :: expected_heat_flux = 2.4610631137_dp
    PetscMPIInt :: rank
    PetscInt :: ierr

    call thermo%init()
    json => fson_parse(str = '{}')
    call eos%init(json, thermo)

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call cell%init(nc, num_phases)
       call face%init(nc, num_phases)
       call fluid%init(nc, num_phases)
       call rock%init()
       allocate(face_data(face%dof), cell_data(cell%dof))
       allocate(rock_data(rock%dof), fluid_data(fluid%dof))
       allocate(flux(num_primary))
       face_offset = 1
       cell_offsets = [1, 1]
       rock_offsets = [1, 1]
       fluid_offsets = [1, 1]
       face_data = [0._dp,  25._dp, 35._dp,  60._dp, 0._dp, 0._dp, -1._dp, &
            9.8_dp, 0._dp, 0._dp, 0._dp, 3._dp]
       cell_data = 0._dp ! not needed
       rock_data = [ &
            1.e-14_dp, 2.e-14_dp, 3.e-15_dp,  2.5_dp,  2.5_dp, 0.1_dp, &
            2200._dp, 1000._dp]
       fluid_data = [ &
            1.e5_dp, 20._dp, 1._dp, 1._dp, 0._dp, &
            998.2_dp, 1.e-3_dp, 1._dp, 1._dp, 0._dp, &
            84011.8_dp, 83911.6_dp, 1._dp]

       call face%assign_geometry(face_data, face_offset)
       call face%assign_cell_geometry(cell_data, cell_offsets)
       call face%assign_cell_rock(rock_data, rock_offsets)
       call face%assign_cell_fluid(fluid_data, fluid_offsets)

       flux = face%flux(eos)

       call test%assert(expected_mass_flux, flux(1), "Mass flux")
       call test%assert(expected_heat_flux, flux(2), "Heat flux")

       call cell%destroy()
       call face%destroy()
       call fluid%destroy()
       call rock%destroy()
       deallocate(face_data, cell_data, rock_data, fluid_data, flux)

    end if

    call eos%destroy()
    call fson_destroy(json)
    call thermo%destroy()

  end subroutine test_face_flux_vertical_gravity

!------------------------------------------------------------------------

  subroutine test_face_flux_hydrostatic(test)

    ! Face flux() test, vertical hydrostatic
    ! Pressure in cell 2 is chosen (by solving a nonlinear equation)
    ! to make the pressure gradient balance the gravity term exactly.

    use cell_module
    use rock_module
    use fluid_module
    use eos_we_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: nc = 1, num_phases = 1, num_primary = 2
    type(face_type) :: face
    type(cell_type) :: cell
    type(rock_type) :: rock
    type(fluid_type) :: fluid
    type(eos_we_type) :: eos
    type(IAPWS_type) :: thermo
    type(fson_value), pointer :: json
    PetscReal, pointer, contiguous :: face_data(:), cell_data(:)
    PetscReal, pointer, contiguous :: rock_data(:), fluid_data(:)
    PetscReal, allocatable :: flux(:)
    PetscInt :: face_offset, cell_offsets(2)
    PetscInt :: rock_offsets(2), fluid_offsets(2)
    PetscReal, parameter :: expected_mass_flux = 0._dp
    PetscReal, parameter :: expected_heat_flux = 0._dp
    PetscMPIInt :: rank
    PetscInt :: ierr

    call thermo%init()
    json => fson_parse(str = '{}')
    call eos%init(json, thermo)

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call cell%init(nc, num_phases)
       call face%init(nc, num_phases)
       call fluid%init(nc, num_phases)
       call rock%init()
       allocate(face_data(face%dof), cell_data(cell%dof))
       allocate(rock_data(rock%dof), fluid_data(fluid%dof*2))
       allocate(flux(num_primary))
       face_offset = 1
       cell_offsets = [1, 1]
       rock_offsets = [1, 1]
       fluid_offsets = [1, 1 + fluid%dof]
       face_data = [0._dp,  25._dp, 35._dp,  60._dp, 0._dp, 0._dp, -1._dp, &
            9.8_dp, 0._dp, 0._dp, 0._dp, 3._dp]
       cell_data = 0._dp ! not needed
       rock_data = [ &
            1.e-14_dp, 2.e-14_dp, 3.e-15_dp,  2.5_dp,  2.5_dp, 0.1_dp, &
            2200._dp, 1000._dp]
       fluid_data = [ &
            2.e5_dp, 20._dp, 1._dp, 1._dp, 0._dp, &                ! cell 1
            998.2512244888_dp, 0.00100156652270771_dp, 1._dp, 1._dp, 0._dp, &
            84105.9189422008_dp, 83905.5685743839_dp, 1._dp, &
            7.87050606076185e5_dp, 20._dp, 1._dp, 1._dp, 0._dp, &  ! cell 2
            998.5195444779_dp, 0.00100138700807062_dp, 1._dp, 1._dp, 0._dp, &
            84658.2021844106_dp, 83869.9846573438_dp, 1._dp]

       call face%assign_geometry(face_data, face_offset)
       call face%assign_cell_geometry(cell_data, cell_offsets)
       call face%assign_cell_rock(rock_data, rock_offsets)
       call face%assign_cell_fluid(fluid_data, fluid_offsets)

       flux = face%flux(eos)

       call test%assert(expected_mass_flux, flux(1), "Mass flux")
       call test%assert(expected_heat_flux, flux(2), "Heat flux")

       call cell%destroy()
       call face%destroy()
       call fluid%destroy()
       call rock%destroy()
       deallocate(face_data, cell_data, rock_data, fluid_data, flux)

    end if

    call eos%destroy()
    call fson_destroy(json)
    call thermo%destroy()

  end subroutine test_face_flux_hydrostatic

!------------------------------------------------------------------------

  subroutine test_face_flux_two_phase_vertical(test)

    ! Face flux() test, 2-phase vertical
    ! The cells have different rock and two-phase fluid properties.

    use cell_module
    use rock_module
    use fluid_module
    use eos_we_module

    class(unit_test_type), intent(in out) :: test
    ! Locals:
    PetscInt, parameter :: nc = 1, num_phases = 2, num_primary = 2
    type(face_type) :: face
    type(cell_type) :: cell
    type(rock_type) :: rock
    type(fluid_type) :: fluid
    type(eos_we_type) :: eos
    type(IAPWS_type) :: thermo
    type(fson_value), pointer :: json
    PetscReal, pointer, contiguous :: face_data(:), cell_data(:)
    PetscReal, pointer, contiguous :: rock_data(:), fluid_data(:)
    PetscReal, allocatable :: flux(:)
    PetscInt :: face_offset, cell_offsets(2)
    PetscInt :: rock_offsets(2), fluid_offsets(2)
    PetscReal :: density
    PetscReal, parameter :: expected_liquid_density = 900.3915384615_dp
    PetscReal, parameter :: expected_vapour_density = 3.7044444444_dp
    PetscReal, parameter :: expected_mass_flux = 9.14772841429594e-5_dp
    PetscReal, parameter :: expected_heat_flux = 57.9124776818_dp
    PetscMPIInt :: rank
    PetscInt :: ierr

    call thermo%init()
    json => fson_parse(str = '{}')
    call eos%init(json, thermo)

    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    if (rank == 0) then

       call cell%init(nc, num_phases)
       call face%init(nc, num_phases)
       call fluid%init(nc, num_phases)
       call rock%init()
       allocate(face_data(face%dof), cell_data(cell%dof))
       allocate(rock_data(rock%dof * 2), fluid_data(fluid%dof * 2))
       allocate(flux(num_primary))
       face_offset = 1
       cell_offsets = [1, 1]
       rock_offsets = [1, 1 + rock%dof]
       fluid_offsets = [1, 1 + fluid%dof]
       face_data = [0._dp,  25._dp, 35._dp,  60._dp, 0._dp, 0._dp, -1._dp, 9.8_dp, &
            0._dp, 0._dp, 0._dp, 3._dp]
       cell_data = 0._dp ! not needed
       rock_data = [ &
            1.e-14_dp, 2.e-14_dp, 3.e-15_dp,  2.5_dp, 2.5_dp, 0.1_dp, &  ! cell 1
            2200._dp, 1000._dp, &
            2.e-14_dp, 3.e-14_dp, 6.e-15_dp,  2.7_dp, 2.7_dp, 0.05_dp, & ! cell 2
            2300._dp, 995._dp]
       fluid_data = [ &
            6.2e5_dp, 160._dp, 4._dp, 3._dp, 0._dp, &        ! cell 1
            907.45_dp, 1.7e-4_dp, 0.25_dp, 0.75_dp, 0._dp, & ! liquid
            675574.7_dp, 674893.5_dp, 1._dp, &
            3.26_dp, 1.43e-5_dp, 0.75_dp, 0.25_dp, 0._dp,  & ! vapour
            2757430.53_dp, 2567774.0_dp, 1._dp, &
            8.2e5_dp, 171.44_dp, 4._dp, 3._dp, 0._dp, &      ! cell 2
            895.98_dp, 1.58e-4_dp, 0.4_dp, 0.6_dp, 0._dp,  & ! liquid
            725517.1_dp, 724601.9_dp, 1._dp, &
            4.26_dp, 1.47e-5_dp, 0.6_dp, 0.4_dp, 0._dp,    & ! vapour
            2769308.8_dp, 2576807.25_dp, 1._dp]

       call face%assign_geometry(face_data, face_offset)
       call face%assign_cell_geometry(cell_data, cell_offsets)
       call face%assign_cell_rock(rock_data, rock_offsets)
       call face%assign_cell_fluid(fluid_data, fluid_offsets)

       density = face%phase_density(1)
       call test%assert(expected_liquid_density, density, "Liquid density")
       density = face%phase_density(2)
       call test%assert(expected_vapour_density, density, "Vapour density")
       flux = face%flux(eos)

       call test%assert(expected_mass_flux, flux(1), "Mass flux")
       call test%assert(expected_heat_flux, flux(2), "Heat flux")

       call cell%destroy()
       call face%destroy()
       call fluid%destroy()
       call rock%destroy()
       deallocate(face_data, cell_data, rock_data, fluid_data, flux)

    end if

    call eos%destroy()
    call fson_destroy(json)
    call thermo%destroy()

  end subroutine test_face_flux_two_phase_vertical

!------------------------------------------------------------------------

end module face_test
