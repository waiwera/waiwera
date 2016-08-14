module face_module
  !! Defines type for accessing local quantities defined on a mesh face.

  use kinds_module
  use cell_module

  implicit none
  private

#include <petsc/finclude/petscsys.h>

  type face_type
     !! Type for accessing local face properties.
     private
     PetscReal, pointer, public :: area !! face area
     PetscReal, pointer, contiguous, public :: distance(:) !! cell centroid distances on either side of the face
     PetscReal, pointer, contiguous, public :: normal(:) !! normal vector to face
     PetscReal, pointer, contiguous, public :: centroid(:) !! centroid of face
     PetscReal, pointer, public :: permeability_direction !! direction of permeability (1.. 3)
     type(cell_type), allocatable, public :: cell(:) !! cells on either side of face
     PetscReal, public :: distance12 !! distance between cell centroids
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
     procedure, public :: normal_gradient => face_normal_gradient
     procedure, public :: pressure_gradient => face_pressure_gradient
     procedure, public :: temperature_gradient => face_temperature_gradient
     procedure, public :: phase_density => face_phase_density
     procedure, public :: harmonic_average => face_harmonic_average
     procedure, public :: permeability => face_permeability
     procedure, public :: heat_conductivity => face_heat_conductivity
     procedure, public :: flow_indices => face_flow_indices
     procedure, public :: upstream_weight => face_upstream_weight
     procedure, public :: transport => face_transport
     procedure, public :: flux => face_flux
  end type face_type

  PetscInt, parameter :: num_face_variables = 5
  PetscInt, parameter, public :: &
       face_variable_num_components(num_face_variables) = &
       [1, 2, 3, 3, 1]
  PetscInt, parameter, public :: &
       face_variable_dim(num_face_variables) = &
       [2, 2, 2, 2, 2]
  PetscInt, parameter :: max_face_variable_name_length = 24
  character(max_face_variable_name_length), parameter, public :: &
       face_variable_names(num_face_variables) = &
       [character(max_face_variable_name_length):: &
       "area", "distance", "normal", "centroid", &
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
    self%normal => data(offset + 3: offset + 5)
    self%centroid => data(offset + 6: offset + 8)
    self%permeability_direction => data(offset + 9)
    self%distance12 = sum(self%distance)

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
    nullify(self%normal)
    nullify(self%centroid)
    nullify(self%permeability_direction)
    if (allocated(self%cell)) then
       deallocate(self%cell)
    end if

  end subroutine face_destroy

!------------------------------------------------------------------------

  subroutine face_calculate_permeability_direction(self)
    !! Calculates permeability direction for the face, being the
    !! coordinate axis most closely aligned with the face normal
    !! vector.

    class(face_type), intent(in out) :: self
    ! Locals:
    PetscInt :: index

    index = maxloc(abs(self%normal), 1)
    self%permeability_direction = dble(index)

  end subroutine face_calculate_permeability_direction

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

  PetscReal function face_pressure_gradient(self) result(dpdn)
    !! Returns pressure gradient across the face.

    class(face_type), intent(in) :: self
    ! Locals:
    PetscReal :: p(2)
    PetscInt :: i
    
    do i = 1, 2
       p(i) = self%cell(i)%fluid%pressure
    end do
    dpdn = self%normal_gradient(p)

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
    PetscInt, intent(in) :: p
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

  subroutine face_flow_indices(self, gradient, up, down)
    !! Returns indices of upstream and downstream cells for a given
    !! effective pressure gradient (including gravity term).

    class(face_type), intent(in) :: self
    PetscReal, intent(in) :: gradient
    PetscInt, intent(out) :: up, down

    if (gradient <= 0._dp) then
       up = 1; down = 2
    else
       up = 2; down = 1
    end if

  end subroutine face_flow_indices

!------------------------------------------------------------------------

  PetscReal function face_upstream_weight(self, gradient) result (w)
    !! Returns upstream weight for computing effective mobility etc.
    !! at the face.

    class(face_type), intent(in) :: self
    PetscReal, intent(in) :: gradient
    ! Locals:
    PetscReal, parameter :: upstream_threshold = 1.e-6_dp

    if (abs(gradient) > upstream_threshold) then
       ! Upstream weighting:
       w = 1._dp
    else
       ! Smoothed upstream weighting to avoid discontinous
       ! mobilities etc. when gradient changes sign:
       w = cubic(gradient / upstream_threshold)
    end if

  contains

    PetscReal function cubic(x) result(f)
      !! Cubic interpolant with the properties:
      !! f(-1) = 1; f(1) = 0; f'(-1) = f'(1) = 0
      PetscReal, intent(in) :: x
      f = 0.5_dp - 0.25_dp * x * (3._dp - x * x)
    end function cubic

  end function face_upstream_weight

!------------------------------------------------------------------------

  subroutine face_transport(self, p, w_up, up, down, &
       mobility, mass_fraction, h)
    !! Returns effective transport quantities (mobility, mass_fractions
    !! and enthalpy) on the face, for phase p, according to the given
    !! upstream weight w_up.

    class(face_type), intent(in) :: self
    PetscInt, intent(in) :: p !! Phase index
    PetscReal, intent(in) :: w_up !! Upstream weight
    PetscInt, intent(in) :: up, down
    PetscReal, intent(out) :: mobility !! Mobility
    PetscReal, dimension(self%cell(1)%fluid%num_components), intent(out) &
         :: mass_fraction
    PetscReal, intent(out) :: h !! Enthalpy
    ! Locals:
    PetscInt :: i
    PetscReal :: w(2)
    PetscReal, parameter :: tol = 1.e-12_dp

    if (w_up >= 1._dp - tol) then
       associate(upstream => self%cell(up)%fluid%phase(p))
         mobility = upstream%mobility()
         mass_fraction = upstream%mass_fraction
         h = upstream%specific_enthalpy
       end associate
    else
       mobility = 0._dp
       mass_fraction = 0._dp
       h = 0._dp
       w(up) = w_up; w(down) = 1._dp - w_up
       do i = 1, 2
          mobility = mobility + w(i) * self%cell(i)%fluid%phase(p)%mobility()
          mass_fraction = mass_fraction + w(i) * &
               self%cell(i)%fluid%phase(p)%mass_fraction
          h = h + w(i) * self%cell(i)%fluid%phase(p)%specific_enthalpy
       end do
    end if

  end subroutine face_transport

!------------------------------------------------------------------------

  function face_flux(self, eos, gravity) result(flux)
    !! Returns array containing the mass fluxes for each component
    !! through the face, from cell(1) to cell(2), and energy flux
    !! for non-isothermal simulations.

    use eos_module, only: eos_type

    class(face_type), intent(in) :: self
    class(eos_type), intent(in) :: eos
    PetscReal, intent(in) :: gravity
    PetscReal :: flux(eos%num_primary_variables)
    ! Locals:
    PetscInt :: nc, np
    PetscInt :: i, p, up, down
    PetscReal :: dpdn, dtdn, gn, G, face_density, F
    PetscReal :: k, h, cond, w_up, mobility
    PetscReal :: mass_fraction(self%cell(1)%fluid%num_components)
    PetscInt :: phases(2), phase_present

    nc = eos%num_components
    np = eos%num_primary_variables
    dpdn = self%pressure_gradient()
    gn = gravity * self%normal(3)

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
          G = dpdn + face_density * gn
          call self%flow_indices(G, up, down)

          if (btest(phases(up), p - 1)) then
             k = self%permeability()
             if (btest(phases(down), p - 1)) then
                w_up = self%upstream_weight(G)
             else
                w_up = 1.0_dp
             end if
             call self%transport(p, w_up, up, down, &
                  mobility, mass_fraction, h)
             ! Mass flows:
             F = -k * mobility * G
             flux(1:nc) = flux(1:nc) + F * mass_fraction
             if (.not. eos%isothermal) then
                ! Heat convection:
                flux(np) = flux(np) + h * F
             end if

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
