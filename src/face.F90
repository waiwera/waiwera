!   Copyright 2016 University of Auckland.

!   This file is part of Waiwera.

!   Waiwera is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.

!   Waiwera is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.

!   You should have received a copy of the GNU Lesser General Public License
!   along with Waiwera.  If not, see <http://www.gnu.org/licenses/>.

module face_module
  !! Defines type for accessing local quantities defined on a mesh face.

#include <petsc/finclude/petscsys.h>

  use petscsys
  use kinds_module
  use cell_module

  implicit none
  private

  type face_type
     !! Type for accessing local face properties.
     private
     PetscReal, pointer, public :: area !! face area
     PetscReal, pointer, contiguous, public :: distance(:) !! cell centroid normal distances on either side of the face
     PetscReal, pointer, public :: distance12 !! normal distance between cell centroids
     PetscReal, pointer, contiguous, public :: normal(:) !! normal vector to face
     PetscReal, pointer, public :: gravity_normal !! dot product of normal with gravity vector
     PetscReal, pointer, contiguous, public :: centroid(:) !! centroid of face
     PetscReal, pointer, public :: permeability_direction !! direction of permeability (1.. 3)
     type(cell_type), allocatable, public :: cell(:) !! cells on either side of face
     PetscInt, public :: dof !! Number of degrees of freedom
   contains
     private
     procedure, public :: init => face_init
     procedure, public :: assign_geometry => face_assign_geometry
     procedure, public :: assign_cell_geometry => face_assign_cell_geometry
     procedure, public :: assign_cell_rock => face_assign_cell_rock
     procedure, public :: assign_cell_fluid => face_assign_cell_fluid
     procedure, public :: destroy => face_destroy
     procedure, public :: calculate_permeability_direction => &
          face_calculate_permeability_direction
     procedure, public :: calculate_distances => face_calculate_distances
     procedure, public :: reversed_orientation => face_reversed_orientation
     procedure, public :: reverse_geometry => face_reverse_geometry
     procedure, public :: normal_gradient => face_normal_gradient
     procedure, public :: pressure_gradient => face_pressure_gradient
     procedure, public :: temperature_gradient => face_temperature_gradient
     procedure, public :: phase_density => face_phase_density
     procedure, public :: harmonic_average => face_harmonic_average
     procedure, public :: permeability => face_permeability
     procedure, public :: heat_conductivity => face_heat_conductivity
     procedure, public :: upstream_index => face_upstream_index
     procedure, public :: flux => face_flux
  end type face_type

  PetscInt, parameter, public :: num_face_variables = 7
  PetscInt, parameter, public :: &
       face_variable_num_components(num_face_variables) = &
       [1, 2, 1, 3, 1, 3, 1]
  PetscInt, parameter, public :: max_face_variable_name_length = 24
  character(max_face_variable_name_length), parameter, public :: &
       face_variable_names(num_face_variables) = &
       [character(max_face_variable_name_length):: &
       "area", "distance", "distance12", "normal", "gravity_normal", "centroid", &
       "permeability_direction"]

  type petsc_face_type
     !! Type for accessing face geometry parameters calculated by
     !! PETSc DMPlexTSGetGeometryFVM().
     private
     PetscReal, pointer, contiguous, public :: area_normal(:) !! normal vector multiplied by area
     PetscReal, pointer, contiguous, public :: centroid(:) !! centroid of face
   contains
     private
     procedure, public :: assign_geometry => petsc_face_assign_geometry
     procedure, public :: destroy => petsc_face_destroy
  end type petsc_face_type

  public :: face_type, petsc_face_type

contains

!------------------------------------------------------------------------

  subroutine face_init(self, num_components, num_phases)
    !! Initialises a face.

    class(face_type), intent(in out) :: self
    PetscInt, intent(in), optional :: num_components !! Number of fluid components
    PetscInt, intent(in), optional :: num_phases !! Number of fluid phases
    ! Locals:
    PetscInt, parameter :: num_cells = 2
    PetscInt :: i

    allocate(self%cell(num_cells))
    if ((present(num_components)) .and. (present(num_phases))) then
       do i = 1, num_cells
          call self%cell(i)%init(num_components, num_phases)
       end do
    end if

    self%dof = sum(face_variable_num_components)

  end subroutine face_init

!------------------------------------------------------------------------

  subroutine face_assign_geometry(self, data, offset)
    !! Assigns geometry pointers in a face to elements of the specified data array,
    !! starting from the given offset.

    class(face_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! array with face geometry data
    PetscInt, intent(in) :: offset  !! face geometry array offset for this face

    self%area => data(offset)
    self%distance => data(offset + 1: offset + 2)
    self%distance12 => data(offset + 3)
    self%normal => data(offset + 4: offset + 6)
    self%gravity_normal => data(offset + 7)
    self%centroid => data(offset + 8: offset + 10)
    self%permeability_direction => data(offset + 11)

  end subroutine face_assign_geometry

!------------------------------------------------------------------------

  subroutine face_assign_cell_geometry(self, data, offsets)
    !! Assigns cell geometry pointers for both cells on the face.

    class(face_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! array with cell geometry data
    PetscInt, intent(in) :: offsets(:)  !! cell geometry array offsets for the face cells
    ! Locals:
    PetscInt :: i

    do i = 1, 2
       call self%cell(i)%assign_geometry(data, offsets(i))
    end do

  end subroutine face_assign_cell_geometry

!------------------------------------------------------------------------

  subroutine face_assign_cell_rock(self, data, offsets)
    !! Assigns rock pointers for both cells on the face.

    class(face_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! array with rock data
    PetscInt, intent(in) :: offsets(:)  !! rock array offsets for the face cells
    ! Locals:
    PetscInt :: i

    do i = 1, 2
       call self%cell(i)%rock%assign(data, offsets(i))
    end do

  end subroutine face_assign_cell_rock

!------------------------------------------------------------------------

  subroutine face_assign_cell_fluid(self, data, offsets)
    !! Assigns fluid pointers for both cells on the face.

    class(face_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! array with fluid data
    PetscInt, intent(in) :: offsets(:)  !! fluid array offsets for the face cells
    ! Locals:
    PetscInt :: i

    do i = 1, 2
       call self%cell(i)%fluid%assign(data, offsets(i))
    end do

  end subroutine face_assign_cell_fluid

!------------------------------------------------------------------------

  subroutine face_destroy(self)
    !! Destroys a face (nullifies all pointer components).

    class(face_type), intent(in out) :: self

    nullify(self%area)
    nullify(self%distance)
    nullify(self%distance12)
    nullify(self%normal)
    nullify(self%gravity_normal)
    nullify(self%centroid)
    nullify(self%permeability_direction)
    if (allocated(self%cell)) then
       deallocate(self%cell)
    end if

  end subroutine face_destroy

!------------------------------------------------------------------------

  subroutine face_calculate_permeability_direction(self, rotation)
    !! Calculates permeability direction for the face, being the
    !! coordinate axis most closely aligned with the face normal
    !! vector. The rotation matrix corresponds to the rotation
    !! transformation of the first horizontal permeability direction.

    class(face_type), intent(in out) :: self
    PetscReal, intent(in) :: rotation(3, 3)
    ! Locals:
    PetscReal :: d(3)
    PetscInt :: index

    d = matmul(rotation, self%normal)
    index = maxloc(abs(d), 1)
    self%permeability_direction = dble(index)

  end subroutine face_calculate_permeability_direction

!------------------------------------------------------------------------

  subroutine face_calculate_distances(self)
    !! Calculates normal distances from cells to face, and normal
    !! distance between cells.

    class(face_type), intent(in out) :: self
    ! Locals:
    PetscReal :: correction

    associate(centroid1 => self%cell(1)%centroid, &
         centroid2 => self%cell(2)%centroid)
      self%distance(1) = dot_product(self%centroid - centroid1, self%normal)
      self%distance(2) = dot_product(centroid2 - self%centroid, self%normal)
      self%distance12 = dot_product(centroid2 - centroid1, self%normal)
    end associate

    ! Correction for non-orthogonal meshes:
    correction = self%distance12 / sum(self%distance)
    self%distance = self%distance * correction

  end subroutine face_calculate_distances

!------------------------------------------------------------------------

  PetscBool function face_reversed_orientation(self) result(reversed)
    !! Returns true if the face does not have the canonical
    !! orientation (i.e. order of cells in its support) but has been
    !! reversed.

    class(face_type), intent(in) :: self

    associate(d12 => self%cell(2)%centroid - self%cell(1)%centroid)
      reversed = (dot_product(d12, self%normal) < 0._dp)
    end associate

  end function face_reversed_orientation

!------------------------------------------------------------------------

  subroutine face_reverse_geometry(self)
    !! Reverses face geometry: flips normal vector by 180 degrees and
    !! reverses centroid distance array.

    class(face_type), intent(in out) :: self

    self%normal = -self%normal
    self%gravity_normal = -self%gravity_normal
    self%distance = self%distance(2:1:-1)

  end subroutine face_reverse_geometry

!------------------------------------------------------------------------

  PetscReal function face_normal_gradient(self, x) result(grad)
    !! Returns gradient along normal vector of the two cell values
    !! x. This assumes that the line joining the two cell centroids is
    !! orthogonal to the face.

    class(face_type), intent(in) :: self
    PetscReal, intent(in) :: x(2)

    grad = (x(2) - x(1)) / self%distance12

  end function face_normal_gradient

!------------------------------------------------------------------------

  PetscReal function face_pressure_gradient(self, p) result(dpdn)
    !! Returns effective pressure gradient for phase p across the
    !! face, including capillary pressure effects.

    class(face_type), intent(in) :: self
    PetscInt, intent(in) :: p !! Phase index
    ! Locals:
    PetscReal :: pressure(2)
    PetscInt :: i
    
    do i = 1, 2
       associate (fluid => self%cell(i)%fluid)
         pressure(i) = fluid%pressure + fluid%phase(p)%capillary_pressure
       end associate
    end do
    dpdn = self%normal_gradient(pressure)

  end function face_pressure_gradient

!------------------------------------------------------------------------

  PetscReal function face_temperature_gradient(self) result(dtdn)
    !! Returns temperature gradient across the face.

    class(face_type), intent(in) :: self
    ! Locals:
    PetscReal :: t(2)
    PetscInt :: i

    do i = 1, 2
       t(i) = self%cell(i)%fluid%temperature
    end do
    dtdn = self%normal_gradient(t)

  end function face_temperature_gradient

!------------------------------------------------------------------------

  PetscReal function face_phase_density(self, p) result(rho)
    !! Returns phase density on the face for a given phase p. It is
    !! assumed that the phase is present in at least one of the cells.

    class(face_type), intent(in) :: self
    PetscInt, intent(in) :: p !! Phase index
    ! Locals:
    PetscInt :: i
    PetscReal :: weight

    rho = 0._dp
    weight = 0._dp
    do i = 1, 2
       associate(phase => self%cell(i)%fluid%phase(p))
         rho = rho + phase%saturation * phase%density
         weight = weight + phase%saturation
       end associate
    end do
    rho = rho / weight

  end function face_phase_density

!------------------------------------------------------------------------

  PetscReal function face_harmonic_average(self, x) result(xh)
    !! Returns harmonic average of the two cell values x, based on the
    !! distances of the cell centroids from the face.

    class(face_type), intent(in) :: self
    PetscReal, intent(in) :: x(2)
    ! Locals:
    PetscReal :: wx
    PetscReal, parameter :: tol = 1.e-30_dp

    wx = (self%distance(1) * x(2) + self%distance(2) * x(1)) / &
         self%distance12

    if (abs(wx) > tol) then
       xh = x(1) * x(2) / wx
    else
       xh = 0._dp
    end if

  end function face_harmonic_average

!------------------------------------------------------------------------

  PetscReal function face_permeability(self) result(k)
    !! Returns effective permeability on the face, harmonic weighted
    !! between the two cells.

    class(face_type), intent(in) :: self
    ! Locals:
    PetscInt :: direction, i
    PetscReal :: perm(2)

    direction = nint(self%permeability_direction)
    do i = 1, 2
       perm(i) = self%cell(i)%rock%permeability(direction)
    end do
    k = self%harmonic_average(perm)

  end function face_permeability

!------------------------------------------------------------------------

  PetscReal function face_heat_conductivity(self, eos) result(K)
    !! Returns effective heat conductivity on the face, harmonically
    !! averaged between the two cells.

    use eos_module, only: eos_type

    class(face_type), intent(in) :: self
    class(eos_type), intent(in) :: eos

    ! Locals:
    PetscReal :: kcell(2)
    PetscInt :: i

    do i = 1, 2
       associate(cell => self%cell(i))
         kcell(i) = eos%conductivity(cell%rock, cell%fluid)
       end associate
    end do
    K = self%harmonic_average(kcell)

  end function face_heat_conductivity

!------------------------------------------------------------------------

  PetscInt function face_upstream_index(self, gradient) result(index)
    !! Returns index of upstream cell for a given effective pressure
    !! gradient (including gravity term).

    class(face_type), intent(in) :: self
    PetscReal, intent(in) :: gradient

    if (gradient <= 0._dp) then
       index = 1
    else
       index = 2
    end if

  end function face_upstream_index

!------------------------------------------------------------------------

  function face_flux(self, eos) result(flux)
    !! Returns array containing the mass fluxes for each component
    !! through the face, from cell(1) to cell(2), and energy flux
    !! for non-isothermal simulations.

    use eos_module, only: eos_type

    class(face_type), intent(in) :: self
    class(eos_type), intent(in) :: eos
    PetscReal :: flux(eos%num_primary_variables)
    ! Locals:
    PetscInt :: nc, np
    PetscInt :: i, p, up
    PetscReal :: dpdn, dtdn, G, face_density, F
    PetscReal :: phase_flux(self%cell(1)%fluid%num_components)
    PetscReal :: k, h, cond
    PetscInt :: phases(2), phase_present

    nc = eos%num_components
    np = eos%num_primary_variables

    if (.not. eos%isothermal) then
       ! Heat conduction:
       cond = self%heat_conductivity(eos)
       dtdn = self%temperature_gradient()
       flux(np) = -cond * dtdn
    end if
    flux(1: nc) = 0._dp

    do i = 1, 2
       phases(i) = nint(self%cell(i)%fluid%phase_composition)
    end do
    phase_present = ior(phases(1), phases(2))

    do p = 1, eos%num_phases

       if (btest(phase_present, p - 1)) then

          face_density = self%phase_density(p)
          dpdn = self%pressure_gradient(p)
          G = dpdn - face_density * self%gravity_normal

          up = self%upstream_index(G)

          if (btest(phases(up), p - 1)) then

             k = self%permeability()
             associate(upstream => self%cell(up)%fluid%phase(p))
               ! Mass flows:
               F = -k * upstream%mobility() * G
               phase_flux = F * upstream%mass_fraction
               flux(1:nc) = flux(1:nc) + phase_flux
               if (.not. eos%isothermal) then
                  ! Heat convection:
                  h = upstream%specific_enthalpy
                  flux(np) = flux(np) + h * F
               end if
             end associate

          end if

       end if

    end do

  end function face_flux

!------------------------------------------------------------------------
! petsc_face_type routines
!------------------------------------------------------------------------

  subroutine petsc_face_assign_geometry(self, data, offset)
    !! Assigns geometry pointers in a petsc_face to elements of the
    !! specified data array, starting from the given offset.

    class(petsc_face_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! array with face geometry data
    PetscInt, intent(in)  :: offset  !! face geometry array offset for this face

    self%area_normal => data(offset: offset + 2)
    self%centroid => data(offset + 3: offset + 5)

  end subroutine petsc_face_assign_geometry

!------------------------------------------------------------------------

  subroutine petsc_face_destroy(self)
    !! Destroys a petsc_face (nullifies all pointer components).

    class(petsc_face_type), intent(in out) :: self

    nullify(self%area_normal)
    nullify(self%centroid)

  end subroutine petsc_face_destroy

!------------------------------------------------------------------------

end module face_module
