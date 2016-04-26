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
     PetscReal, pointer, public :: distance(:) !! cell centroid distances on either side of the face
     PetscReal, pointer, public :: normal(:) !! normal vector to face
     PetscReal, pointer, public :: centroid(:) !! centroid of face
     PetscReal, pointer, public :: permeability_direction !! direction of permeability (1.. 3)
     type(cell_type), allocatable, public :: cell(:) !! cells on either side of face
     PetscReal, public :: distance12 !! distance between cell centroids
   contains
     private
     procedure, public :: init => face_init
     procedure, public :: assign_geometry => face_assign_geometry
     procedure, public :: assign_cell_geometry => face_assign_cell_geometry
     procedure, public :: assign_cell_rock => face_assign_cell_rock
     procedure, public :: assign_cell_fluid => face_assign_cell_fluid
     procedure, public :: dof => face_dof
     procedure, public :: destroy => face_destroy
     procedure, public :: calculate_permeability_direction => &
          face_calculate_permeability_direction
     procedure, public :: normal_gradient => face_normal_gradient
     procedure, public :: pressure_gradient => face_pressure_gradient
     procedure, public :: temperature_gradient => face_temperature_gradient
     procedure, public :: phase_density => face_phase_density
     procedure, public :: harmonic_average => face_harmonic_average
     procedure, public :: permeability => face_permeability
     procedure, public :: mobility => face_mobility
     procedure, public :: heat_conductivity => face_heat_conductivity
     procedure, public :: upstream_index => face_upstream_index
     procedure, public :: flux => face_flux
  end type face_type

  PetscInt, parameter :: num_face_variables = 5
  PetscInt, parameter, public :: &
       face_variable_num_components(num_face_variables) = &
       [1, 2, 3, 3, 1]

  type petsc_face_type
     !! Type for accessing face geometry parameters calculated by
     !! PETSc DMPlexTSGetGeometryFVM().
     private
     PetscReal, pointer, public :: area_normal(:) !! normal vector multiplied by area
     PetscReal, pointer, public :: centroid(:) !! centroid of face
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

  end subroutine face_init

!------------------------------------------------------------------------

  subroutine face_assign_geometry(self, face_geom_data, face_geom_offset)
    !! Assigns geometry pointers in a face to elements of the specified data array,
    !! starting from the given offset.

    use profiling_module, only: assign_pointers_event

    class(face_type), intent(in out) :: self
    PetscReal, target, intent(in), optional :: face_geom_data(:)  !! array with face geometry data
    PetscInt, intent(in), optional  :: face_geom_offset  !! face geometry array offset for this face
    ! Locals:
    PetscErrorCode :: ierr

    call PetscLogEventBegin(assign_pointers_event, ierr); CHKERRQ(ierr)
    
    self%area => face_geom_data(face_geom_offset)
    self%distance => face_geom_data(face_geom_offset + 1: face_geom_offset + 2)
    self%normal => face_geom_data(face_geom_offset + 3: face_geom_offset + 5)
    self%centroid => face_geom_data(face_geom_offset + 6: face_geom_offset + 8)
    self%permeability_direction => face_geom_data(face_geom_offset + 9)
    self%distance12 = sum(self%distance)

    call PetscLogEventEnd(assign_pointers_event, ierr); CHKERRQ(ierr)

  end subroutine face_assign_geometry

!------------------------------------------------------------------------

  subroutine face_assign_cell_geometry(self, cell_geom_data, &
       cell_geom_offsets)
    !! Assigns cell geometry pointers for both cells on the face.

    class(face_type), intent(in out) :: self
    PetscReal, target, intent(in), optional :: cell_geom_data(:)  !! array with cell geometry data
    PetscInt, intent(in), optional  :: cell_geom_offsets(:)  !! cell geometry array offsets for the face cells
    ! Locals:
    PetscInt :: i

    do i = 1, 2
       call self%cell(i)%assign_geometry(cell_geom_data, cell_geom_offsets(i))
    end do

  end subroutine face_assign_cell_geometry

!------------------------------------------------------------------------

  subroutine face_assign_cell_rock(self, rock_data, rock_offsets)
    !! Assigns rock pointers for both cells on the face.

    class(face_type), intent(in out) :: self
    PetscReal, target, intent(in), optional :: rock_data(:)  !! array with rock data
    PetscInt, intent(in), optional  :: rock_offsets(:)  !! rock array offsets for the face cells
    ! Locals:
    PetscInt :: i

    do i = 1, 2
       call self%cell(i)%rock%assign(rock_data, rock_offsets(i))
    end do

  end subroutine face_assign_cell_rock

!------------------------------------------------------------------------

  subroutine face_assign_cell_fluid(self, fluid_data, fluid_offsets)
    !! Assigns fluid pointers for both cells on the face.

    class(face_type), intent(in out) :: self
    PetscReal, target, intent(in), optional :: fluid_data(:)  !! array with fluid data
    PetscInt, intent(in), optional  :: fluid_offsets(:)  !! fluid array offsets for the face cells
    ! Locals:
    PetscInt :: i

    do i = 1, 2
       call self%cell(i)%fluid%assign(fluid_data, fluid_offsets(i))
    end do

  end subroutine face_assign_cell_fluid

!------------------------------------------------------------------------

  PetscInt function face_dof(self)
    !! Returns number of degrees of freedom in a face object.

    class(face_type), intent(in) :: self

    face_dof = sum(face_variable_num_components)

  end function face_dof

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
       rho = rho + self%cell(i)%fluid%phase(p)%saturation * &
            self%cell(i)%fluid%phase(p)%density
       weight = weight + self%cell(i)%fluid%phase(p)%saturation
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

  PetscReal function face_mobility(self, p, up) result(mobility)
    !! Returns effective mobility of phase p on the face, upstream
    !! weighted between the two cells.

    class(face_type), intent(in) :: self
    PetscInt, intent(in) :: p  !! Phase index
    PetscInt, intent(in) :: up !! Upstream cell index

    mobility = self%cell(up)%fluid%phase(p)%mobility()

  end function face_mobility

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
       kcell(i) = eos%conductivity(self%cell(i)%rock, self%cell(i)%fluid)
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

  function face_flux(self, eos, gravity) result(flux)
    !! Returns array containing the mass fluxes for each component
    !! through the face, from cell(1) to cell(2), and energy flux
    !! for non-isothermal simulations.

    use eos_module, only: eos_type
    use profiling_module, only: face_flux_event

    class(face_type), intent(in) :: self
    class(eos_type), intent(in) :: eos
    PetscReal, intent(in) :: gravity
    PetscReal :: flux(eos%num_primary_variables)
    ! Locals:
    PetscInt :: nc, np
    PetscInt :: i, p, up
    PetscReal :: dpdn, dtdn, gn, G, face_density, F
    PetscReal :: phase_flux(self%cell(1)%fluid%num_components)
    PetscReal :: k, h, cond, mobility
    PetscInt :: phases(2), phase_present
    PetscErrorCode :: ierr

    call PetscLogEventBegin(face_flux_event, ierr); CHKERRQ(ierr)

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

    do p = 1, self%cell(1)%fluid%num_phases

       if (btest(phase_present, p - 1)) then

          face_density = self%phase_density(p)
          G = dpdn + face_density * gn

          up = self%upstream_index(G)

          if (btest(phases(up), p - 1)) then

             k = self%permeability()
             mobility = self%mobility(p, up)

             ! Mass flows:
             F = -k * mobility * G
             phase_flux = F * self%cell(up)%fluid%phase(p)%mass_fraction
             flux(1:nc) = flux(1:nc) + phase_flux

             if (.not. eos%isothermal) then
                ! Heat convection:
                h = self%cell(up)%fluid%phase(p)%specific_enthalpy
                flux(np) = flux(np) + h * F
             end if

          end if

       end if

    end do

    call PetscLogEventEnd(face_flux_event, ierr); CHKERRQ(ierr)

  end function face_flux

!------------------------------------------------------------------------
! petsc_face_type routines
!------------------------------------------------------------------------

  subroutine petsc_face_assign_geometry(self, face_geom_data, face_geom_offset)
    !! Assigns geometry pointers in a petsc_face to elements of the
    !! specified data array, starting from the given offset.

    class(petsc_face_type), intent(in out) :: self
    PetscReal, target, intent(in) :: face_geom_data(:)  !! array with face geometry data
    PetscInt, intent(in)  :: face_geom_offset  !! face geometry array offset for this face

    self%area_normal => face_geom_data(face_geom_offset: face_geom_offset + 2)
    self%centroid => face_geom_data(face_geom_offset + 3: face_geom_offset + 5)

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
