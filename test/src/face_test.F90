module face_test

  ! Test for face module

  use kinds_module
  use mpi_module
  use fruit
  use face_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>

public :: test_face_assign, test_face_permeability_direction, &
     test_face_normal_gradient, test_face_harmonic_average, &
     test_face_flux_zero_horizontal, test_face_flux_vertical_gravity, &
     test_face_flux_hydrostatic, test_face_flux_two_phase_vertical

PetscReal, parameter :: tol = 1.e-6_dp
PetscReal, parameter :: mass_tol = 1.e-10_dp, heat_tol = 1.e-6

contains
  
!------------------------------------------------------------------------

  subroutine test_face_assign

    ! Face assign() test

    type(face_type) :: face
    PetscReal, parameter :: area = 300._dp
    PetscReal, parameter :: distance(2) = [20._dp, 30._dp]
    PetscReal, parameter :: normal(3) = [0.5_dp, -0.25_dp, 0.75_dp]
    PetscReal, parameter :: centroid(3) = [-1250._dp, 3560._dp, -2530._dp]
    PetscReal, parameter :: permeability_direction = dble(2)
    PetscInt, parameter :: offset = 6
    PetscReal :: offset_padding(offset-1) = 0._dp
    PetscReal, allocatable :: face_data(:)

    if (mpi%rank == mpi%output_rank) then

       call face%init()

       face_data = [offset_padding, area, distance, normal, centroid, &
            permeability_direction]

       call assert_equals(face%dof(), size(face_data) - (offset-1), "face dof")

       call face%assign(face_data, offset)

       call assert_equals(area, face%area, tol, "area")
       call assert_equals(0._dp, norm2(face%distance - distance), tol, "distances")
       call assert_equals(0._dp, norm2(face%normal - normal), tol, "normal")
       call assert_equals(0._dp, norm2(face%centroid - centroid), tol, "centroid")
       call assert_equals(permeability_direction, face%permeability_direction, &
            tol, "permeability direction")

       call face%destroy()
       deallocate(face_data)

    end if

  end subroutine test_face_assign

!------------------------------------------------------------------------

  subroutine test_face_permeability_direction

    ! Face permeability_direction() test

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
    PetscReal, allocatable :: face_data(:)
    PetscInt :: i
    PetscInt :: offset = 1
    character(len = 32) :: msg

    if (mpi%rank == mpi%output_rank) then

       call face%init()

       do i = 1, num_tests

          face_data = [area, distance, normal(:,i), centroid, &
               initial_permeability_direction]
          call face%assign(face_data, offset)
          call face%calculate_permeability_direction()

          write(msg, '(a, i2)') "Permeability direction test ", i
          call assert_equals(expected_permeability_direction(i), &
               face%permeability_direction, tol, trim(msg))

          call face%destroy()
          deallocate(face_data)

       end do

    end if

  end subroutine test_face_permeability_direction

!------------------------------------------------------------------------

  subroutine test_face_normal_gradient

    ! Face normal_gradient() test

    use cell_module

    type(face_type) :: face
    type(cell_type) :: cell
    PetscInt :: face_dof, cell_dof
    PetscReal, parameter :: distance(2) = [25._dp, 32._dp]
    PetscReal, parameter :: x(2) = [240._dp, 170._dp]
    PetscReal, parameter :: expected_d12 = 57._dp
    PetscReal, parameter :: expected_g = -1.22807017544_dp
    PetscReal, allocatable :: face_data(:)
    PetscReal, allocatable :: cell_data(:)
    PetscInt :: face_offset, cell_offsets(2)
    PetscReal :: g

    if (mpi%rank == mpi%output_rank) then

       call face%init()
       face_dof = face%dof()
       cell_dof = cell%dof()

       allocate(face_data(face_dof), cell_data(cell_dof * 2))
       face_offset = 1
       cell_offsets = [1, 1 + cell_dof]
       face_data = 0._dp
       face_data(2:3) = distance
       cell_data = 0._dp

       call face%assign(face_data, face_offset, cell_data, cell_offsets)

       g = face%normal_gradient(x)

       call assert_equals(expected_d12, face%distance12, tol, "face distance12")
       call assert_equals(expected_g, g, tol, "face normal gradient")

       deallocate(face_data, cell_data)

    end if

  end subroutine test_face_normal_gradient

!------------------------------------------------------------------------

  subroutine test_face_harmonic_average

    ! Face harmonic_average() test

    use cell_module

    type(face_type) :: face
    type(cell_type) :: cell
    PetscInt :: face_dof, cell_dof
    PetscReal, allocatable :: face_data(:)
    PetscReal, allocatable :: cell_data(:)
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

    if (mpi%rank == mpi%output_rank) then

       call face%init()
       face_dof = face%dof()
       cell_dof = cell%dof()
       allocate(face_data(face_dof), cell_data(cell_dof * 2))
       face_offset = 1
       cell_offsets = [1, 1 + cell_dof]
       face_data = 0._dp
       cell_data = 0._dp

       do i = 1, num_tests
          face_data(2:3) = distance(:, i)
          call face%assign(face_data, face_offset, cell_data, cell_offsets)
          xh = face%harmonic_average(x(:, i))
          write(msg, '(a, i2)') "Face harmonic average test ", i
          call assert_equals(expected_xh(i), xh, tol, msg)
       end do

       deallocate(face_data, cell_data)

    end if

  end subroutine test_face_harmonic_average

!------------------------------------------------------------------------

  subroutine test_face_flux_zero_horizontal

    ! Face flux() test, 1-phase horizontal
    ! Fluid properties in the two cells are identical, so the fluxes
    ! should be zero.

    use cell_module
    use rock_module
    use fluid_module

    PetscInt, parameter :: nc = 1, np = 1, num_primary = 2
    PetscReal, parameter :: gravity = 9.8_dp
    type(face_type) :: face
    type(cell_type) :: cell
    type(rock_type) :: rock
    type(fluid_type) :: fluid
    PetscReal, allocatable :: face_data(:), cell_data(:)
    PetscReal, allocatable :: rock_data(:), fluid_data(:)
    PetscReal, allocatable :: flux(:)
    PetscInt :: face_offset, cell_offsets(2)
    PetscInt :: rock_offsets(2), fluid_offsets(2)
    PetscReal, parameter :: expected_mass_flux = 0._dp
    PetscReal, parameter :: expected_heat_flux = 0._dp

    if (mpi%rank == mpi%output_rank) then

       call face%init(nc, np)
       call fluid%init(nc, np)
       allocate(face_data(face%dof()), cell_data(cell%dof()))
       allocate(rock_data(rock%dof()), fluid_data(fluid%dof()))
       allocate(flux(num_primary))
       face_offset = 1
       cell_offsets = [1, 1]
       rock_offsets = [1, 1]
       fluid_offsets = [1, 1]
       face_data = [0._dp,  25._dp, 35._dp,  1._dp, 0._dp, 0._dp, &
            0._dp, 0._dp, 0._dp, 1._dp]
       cell_data = 0._dp ! not needed
       rock_data = [ &
            1.e-14_dp, 2.e-14_dp, 3.e-15_dp,  2.5_dp,  0.1_dp, &
            2200._dp, 1000._dp]
       fluid_data = [ &
            1.e5_dp, 20._dp, 1._dp, & 
            998.2_dp, 1.e-3_dp, 1._dp, 1._dp, &
            84011.8_dp, 83911.6_dp, 1._dp]

       call face%assign(face_data, face_offset, cell_data, cell_offsets, &
            rock_data, rock_offsets, fluid_data, fluid_offsets)

       flux = face%flux(num_primary, gravity)

       call assert_equals(num_primary, size(flux), "Flux array size")
       call assert_equals(expected_mass_flux, flux(1), tol, "Mass flux")
       call assert_equals(expected_heat_flux, flux(2), tol, "Heat flux")

       call face%destroy()
       deallocate(face_data, cell_data, rock_data, fluid_data)
       deallocate(flux)
       call fluid%destroy()

    end if

  end subroutine test_face_flux_zero_horizontal

!------------------------------------------------------------------------

  subroutine test_face_flux_vertical_gravity

    ! Face flux() test, 1-phase vertical, gravity only
    ! Fluid properties in both cells are identical, so the only flow
    ! is from gravity.

    use cell_module
    use rock_module
    use fluid_module

    PetscInt, parameter :: nc = 1, np = 1, num_primary = 2
    PetscReal, parameter :: gravity = 9.8_dp
    type(face_type) :: face
    type(cell_type) :: cell
    type(rock_type) :: rock
    type(fluid_type) :: fluid
    PetscReal, allocatable :: face_data(:), cell_data(:)
    PetscReal, allocatable :: rock_data(:), fluid_data(:)
    PetscReal, allocatable :: flux(:)
    PetscInt :: face_offset, cell_offsets(2)
    PetscInt :: rock_offsets(2), fluid_offsets(2)
    PetscReal, parameter :: expected_mass_flux = 2.9294255256e-5_dp
    PetscReal, parameter :: expected_heat_flux = 2.4610631137_dp

    if (mpi%rank == mpi%output_rank) then

       call face%init(nc, np)
       call fluid%init(nc, np)
       allocate(face_data(face%dof()), cell_data(cell%dof()))
       allocate(rock_data(rock%dof()), fluid_data(fluid%dof()))
       allocate(flux(num_primary))
       face_offset = 1
       cell_offsets = [1, 1]
       rock_offsets = [1, 1]
       fluid_offsets = [1, 1]
       face_data = [0._dp,  25._dp, 35._dp,  0._dp, 0._dp, -1._dp, &
            0._dp, 0._dp, 0._dp, 3._dp]
       cell_data = 0._dp ! not needed
       rock_data = [ &
            1.e-14_dp, 2.e-14_dp, 3.e-15_dp,  2.5_dp,  0.1_dp, &
            2200._dp, 1000._dp]
       fluid_data = [ &
            1.e5_dp, 20._dp, 1._dp, & 
            998.2_dp, 1.e-3_dp, 1._dp, 1._dp, &
            84011.8_dp, 83911.6_dp, 1._dp]

       call face%assign(face_data, face_offset, cell_data, cell_offsets, &
            rock_data, rock_offsets, fluid_data, fluid_offsets)

       flux = face%flux(num_primary, gravity)

       call assert_equals(expected_mass_flux, flux(1), mass_tol, "Mass flux")
       call assert_equals(expected_heat_flux, flux(2), heat_tol, "Heat flux")

       call face%destroy()
       deallocate(face_data, cell_data, rock_data, fluid_data)
       deallocate(flux)
       call fluid%destroy()

    end if

  end subroutine test_face_flux_vertical_gravity

!------------------------------------------------------------------------

  subroutine test_face_flux_hydrostatic

    ! Face flux() test, vertical hydrostatic
    ! Pressure in cell 2 is chosen (by solving a nonlinear equation)
    ! to make the pressure gradient balance the gravity term exactly.

    use cell_module
    use rock_module
    use fluid_module

    PetscInt, parameter :: nc = 1, np = 1, num_primary = 2
    PetscReal, parameter :: gravity = 9.8_dp
    type(face_type) :: face
    type(cell_type) :: cell
    type(rock_type) :: rock
    type(fluid_type) :: fluid
    PetscReal, allocatable :: face_data(:), cell_data(:)
    PetscReal, allocatable :: rock_data(:), fluid_data(:)
    PetscReal, allocatable :: flux(:)
    PetscInt :: face_offset, cell_offsets(2)
    PetscInt :: rock_offsets(2), fluid_offsets(2)
    PetscReal, parameter :: expected_mass_flux = 0._dp
    PetscReal, parameter :: expected_heat_flux = 0._dp

    if (mpi%rank == mpi%output_rank) then

       call face%init(nc, np)
       call fluid%init(nc, np)
       allocate(face_data(face%dof()), cell_data(cell%dof()))
       allocate(rock_data(rock%dof()), fluid_data(fluid%dof()*2))
       allocate(flux(num_primary))
       face_offset = 1
       cell_offsets = [1, 1]
       rock_offsets = [1, 1]
       fluid_offsets = [1, 1 + fluid%dof()]
       face_data = [0._dp,  25._dp, 35._dp,  0._dp, 0._dp, -1._dp, &
            0._dp, 0._dp, 0._dp, 3._dp]
       cell_data = 0._dp ! not needed
       rock_data = [ &
            1.e-14_dp, 2.e-14_dp, 3.e-15_dp,  2.5_dp,  0.1_dp, &
            2200._dp, 1000._dp]
       fluid_data = [ &
            2.e5_dp, 20._dp, 1._dp, &                ! cell 1
            998.2512244888_dp, 0.00100156652270771_dp, 1._dp, 1._dp, &
            84105.9189422008_dp, 83905.5685743839_dp, 1._dp, &
            7.87050606076185e5_dp, 20._dp, 1._dp, &  ! cell 2
            998.5195444779_dp, 0.00100138700807062_dp, 1._dp, 1._dp, &
            84658.2021844106_dp, 83869.9846573438_dp, 1._dp]

       call face%assign(face_data, face_offset, cell_data, cell_offsets, &
            rock_data, rock_offsets, fluid_data, fluid_offsets)

       flux = face%flux(num_primary, gravity)

       call assert_equals(expected_mass_flux, flux(1), mass_tol, "Mass flux")
       call assert_equals(expected_heat_flux, flux(2), heat_tol, "Heat flux")

       call face%destroy()
       deallocate(face_data, cell_data, rock_data, fluid_data)
       deallocate(flux)
       call fluid%destroy()

    end if

  end subroutine test_face_flux_hydrostatic

!------------------------------------------------------------------------

  subroutine test_face_flux_two_phase_vertical

    ! Face flux() test, 2-phase vertical
    ! The cells have different rock and two-phase fluid properties.

    use cell_module
    use rock_module
    use fluid_module

    PetscInt, parameter :: nc = 1, np = 2, num_primary = 2
    PetscReal, parameter :: gravity = 9.8_dp
    type(face_type) :: face
    type(cell_type) :: cell
    type(rock_type) :: rock
    type(fluid_type) :: fluid
    PetscReal, allocatable :: face_data(:), cell_data(:)
    PetscReal, allocatable :: rock_data(:), fluid_data(:)
    PetscReal, allocatable :: flux(:)
    PetscInt :: face_offset, cell_offsets(2)
    PetscInt :: rock_offsets(2), fluid_offsets(2)
    PetscReal, parameter :: expected_mass_flux = 9.16974670293235e-5_dp
    PetscReal, parameter :: expected_heat_flux = 58.061787312_dp

    if (mpi%rank == mpi%output_rank) then

       call face%init(nc, np)
       call fluid%init(nc, np)
       allocate(face_data(face%dof()), cell_data(cell%dof()))
       allocate(rock_data(rock%dof() * 2), fluid_data(fluid%dof() * 2))
       allocate(flux(num_primary))
       face_offset = 1
       cell_offsets = [1, 1]
       rock_offsets = [1, 1 + rock%dof()]
       fluid_offsets = [1, 1 + fluid%dof()]
       face_data = [0._dp,  25._dp, 35._dp,  0._dp, 0._dp, -1._dp, 0._dp, &
            0._dp, 0._dp, 3._dp]
       cell_data = 0._dp ! not needed
       rock_data = [ &
            1.e-14_dp, 2.e-14_dp, 3.e-15_dp,  2.5_dp,  0.1_dp, &  ! cell 1
            2200._dp, 1000._dp, &
            2.e-14_dp, 3.e-14_dp, 6.e-15_dp,  2.7_dp,  0.05_dp, & ! cell 2
            2300._dp, 995._dp]
       fluid_data = [ &
            6.2e5_dp, 160._dp, 4._dp, &               ! cell 1
            907.45_dp, 1.7e-4_dp, 0.25_dp, 0.75_dp, & ! liquid
            675574.7_dp, 674893.5_dp, 1._dp, &
            3.26_dp, 1.43e-5_dp, 0.75_dp, 0.25_dp, &  ! vapour
            2757430.53_dp, 2567774.0_dp, 1._dp, &
            8.2e5_dp, 171.44_dp, 4._dp, &             ! cell 2
            895.98_dp, 1.58e-4_dp, 0.4_dp, 0.6_dp, &  ! liquid
            725517.1_dp, 724601.9_dp, 1._dp, &
            4.26_dp, 1.47e-5_dp, 0.6_dp, 0.4_dp, &    ! vapour
            2769308.8_dp, 2576807.25_dp, 1._dp]

       call face%assign(face_data, face_offset, cell_data, cell_offsets, &
            rock_data, rock_offsets, fluid_data, fluid_offsets)

       flux = face%flux(num_primary, gravity)

       call assert_equals(expected_mass_flux, flux(1), mass_tol, "Mass flux")
       call assert_equals(expected_heat_flux, flux(2), heat_tol, "Heat flux")

       call face%destroy()
       deallocate(face_data, cell_data, rock_data, fluid_data)
       deallocate(flux)
       call fluid%destroy()

    end if

  end subroutine test_face_flux_two_phase_vertical

!------------------------------------------------------------------------

end module face_test
