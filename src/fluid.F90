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

module fluid_module
  !! Module for accessing fluid properties.
  !!
  !! A [[fluid_type(type)]] object has two kinds of properties:
  !! **bulk** properties and **phase** properties. The bulk properties
  !! relate to the fluid as a whole. The phase properties describe the
  !! state of each fluid phase and are stored in an array component of
  !! [[phase_type]] objects.

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module

  implicit none
  private

  PetscInt, parameter, public :: num_fluid_variables = 5
  PetscInt, parameter, public :: num_phase_variables = 8
  PetscInt, parameter, public :: max_fluid_variable_name_length = 16
  character(max_fluid_variable_name_length), public :: &
       fluid_variable_names(num_fluid_variables) = [ &
       "pressure        ", "temperature     ", &
       "region          ", "phases          ", &
       "partial_pressure"]
  PetscInt, parameter, public :: max_phase_variable_name_length = 21
  character(max_phase_variable_name_length), public :: &
       phase_variable_names(num_phase_variables) = [ &
       "density              ", "viscosity            ", &
       "saturation           ", "relative_permeability", &
       "capillary_pressure   ", &
       "specific_enthalpy    ", "internal_energy      ", &
       "mass_fraction        "]

  type phase_type
     !! Type for accessing local fluid properties for a particular phase.
     PetscReal, pointer :: density    !! Phase density
     PetscReal, pointer :: viscosity  !! Viscosity
     PetscReal, pointer :: saturation !! Phase saturation
     PetscReal, pointer :: relative_permeability !! Relative permeability
     PetscReal, pointer :: capillary_pressure !! Capillary pressure
     PetscReal, pointer :: specific_enthalpy !! Specific enthalpy
     PetscReal, pointer :: internal_energy !! Internal energy
     PetscReal, pointer, contiguous :: mass_fraction(:) !! Component mass fractions
   contains
     private
     procedure, public :: destroy => phase_destroy
     procedure, public :: mobility => phase_mobility
  end type phase_type

  type fluid_type
     !! Type for accessing local fluid properties.
     !!
     !! The main bulk fluid properties are [[fluid_type:pressure]] and
     !! [[fluid_type:temperature]]. The [[fluid_type:phase]] component
     !! is an array of [[phase_type]] objects containing the fluid
     !! properties of individual phases.  The thermodynamic
     !! [[fluid_type:region]] and [[fluid_type:phase_composition]]
     !! contain integer data (though stored as real in the underlying
     !! data vector).
     private
     PetscReal, pointer, public :: pressure    !! Pressure
     PetscReal, pointer, public :: temperature !! Temperature
     PetscReal, pointer, public :: region      !! Thermodynamic region
     PetscReal, pointer, public :: phase_composition   !! Phase composition
     PetscReal, pointer, contiguous, public :: partial_pressure(:) !! Component partial pressures
     type(phase_type), allocatable, public :: phase(:) !! Phase variables
     PetscInt, public :: num_phases !! Number of phases
     PetscInt, public :: num_components !! Number of mass components
     PetscInt, public :: dof !! Number of degrees of freedom
     PetscInt, public :: phase_dof !! Number of degrees of freedom per phase
   contains
     private
     procedure, public :: init => fluid_init
     procedure, public :: assign => fluid_assign
     procedure, public :: destroy => fluid_destroy
     procedure, public :: component_density => fluid_component_density
     procedure, public :: component_mass_fraction => fluid_component_mass_fraction
     procedure, public :: energy => fluid_energy
     procedure, public :: phase_mobilities => fluid_phase_mobilities
     procedure, public :: phase_flow_fractions => fluid_phase_flow_fractions
     procedure, public :: component_flow_fractions => fluid_component_flow_fractions
     procedure, public :: specific_enthalpy => fluid_specific_enthalpy
     procedure, public :: update_phase_composition => &
          fluid_update_phase_composition
  end type fluid_type

  public :: fluid_type, phase_type, create_fluid_vector

contains

!------------------------------------------------------------------------
! Phase procedures
!------------------------------------------------------------------------

  subroutine phase_destroy(self)
    !! Destroys a phase object.

    class(phase_type), intent(in out) :: self

    self%density => null()
    self%viscosity => null()
    self%saturation => null()
    self%relative_permeability => null()
    self%capillary_pressure => null()
    self%specific_enthalpy => null()
    self%internal_energy => null()
    self%mass_fraction => null()

  end subroutine phase_destroy

!------------------------------------------------------------------------

  PetscReal function phase_mobility(self)
    !! Returns mobility of the phase:
    !! mobility = \(k_r \rho / \mu \)

    class(phase_type), intent(in) :: self

    phase_mobility = self%relative_permeability * self%density / &
         self%viscosity

  end function phase_mobility

!------------------------------------------------------------------------
! Fluid procedures
!------------------------------------------------------------------------

  subroutine fluid_init(self, num_components, num_phases)
    !! Initialises a fluid object.

    class(fluid_type), intent(in out) :: self
    PetscInt, intent(in) :: num_components !! Number of fluid components
    PetscInt, intent(in) :: num_phases !! Number of fluid phases

    self%num_phases = num_phases
    self%num_components = num_components
    allocate(self%phase(num_phases))

    self%phase_dof = num_phase_variables + self%num_components - 1
    associate(bulk_dof => num_fluid_variables + self%num_components - 1)
      self%dof = bulk_dof + self%num_phases * self%phase_dof
    end associate

  end subroutine fluid_init
    
!------------------------------------------------------------------------

  subroutine fluid_assign(self, data, offset)
    !! Assigns pointers in a fluid object to elements in the data array,
    !! starting from the specified offset.

    class(fluid_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: data(:)  !! fluid data array
    PetscInt, intent(in) :: offset  !! fluid array offset
    ! Locals:
    PetscInt :: i, p

    self%pressure => data(offset)
    self%temperature => data(offset + 1)
    self%region => data(offset + 2)
    self%phase_composition => data(offset + 3)
    self%partial_pressure => data(offset + 4: offset + 4 + self%num_components - 1)
    
    associate(bulk_dof => num_fluid_variables + self%num_components - 1)
      i = offset + bulk_dof
    end associate
    do p = 1, self%num_phases
       associate(phase => self%phase(p))
         phase%density => data(i)
         phase%viscosity => data(i + 1)
         phase%saturation => data(i + 2)
         phase%relative_permeability => data(i + 3)
         phase%capillary_pressure => data(i + 4)
         phase%specific_enthalpy => data(i + 5)
         phase%internal_energy => data(i + 6)
         phase%mass_fraction => data(i + 7: i + 7 + self%num_components - 1)
         i = i + self%phase_dof
       end associate
    end do

  end subroutine fluid_assign

!------------------------------------------------------------------------

  subroutine fluid_destroy(self)
    !! Destroys a fluid object.
    
    class(fluid_type), intent(in out) :: self
    ! Locals:
    PetscInt :: p

    nullify(self%pressure)
    nullify(self%temperature)
    nullify(self%region)
    nullify(self%phase_composition)
    nullify(self%partial_pressure)

    do p = 1, self%num_phases
       call self%phase(p)%destroy()
    end do
    deallocate(self%phase)

  end subroutine fluid_destroy

!------------------------------------------------------------------------

  function fluid_component_density(self) result(d)
    !! Returns total fluid density for each mass component, over all
    !! phases.

    class(fluid_type), intent(in) :: self
    PetscReal :: d(self%num_components)
    ! Locals:
    PetscInt :: p, c, phases
    PetscReal :: ds

    phases = nint(self%phase_composition)
    d = 0._dp
    do p = 1, self%num_phases
       if (btest(phases, p - 1)) then
          associate(phase => self%phase(p))
            ds = phase%density * phase%saturation
            do c = 1, self%num_components
               d(c) = d(c) + ds * phase%mass_fraction(c)
            end do
          end associate
       end if
    end do

  end function fluid_component_density

!------------------------------------------------------------------------

  PetscReal function fluid_component_mass_fraction(self, c) result(xc)
    !! Returns total mass fraction for specified component, with
    !! contributions from all phases.

    class(fluid_type), intent(in) :: self
    PetscInt, intent(in) :: c !! Mass component
    ! Locals:
    PetscInt :: p, phases
    PetscReal :: ds, total
    PetscReal, parameter :: small = 1.e-30_dp

    xc = 0._dp
    total = 0._dp

    phases = nint(self%phase_composition)
    do p = 1, self%num_phases
       if (btest(phases, p - 1)) then
          associate(phase => self%phase(p))
            ds = phase%density * phase%saturation
            xc = xc + ds * phase%mass_fraction(c)
            total = total + ds
          end associate
       end if
    end do

    if (total > small) then
       xc = xc / total
    else
       xc = 0._dp
    end if

  end function fluid_component_mass_fraction

!------------------------------------------------------------------------

  PetscReal function fluid_energy(self) result(ef)
    !! Returns total fluid energy density.

    class(fluid_type), intent(in) :: self
    ! Locals:
    PetscInt :: p, phases
    PetscReal :: ds

    phases = nint(self%phase_composition)
    ef = 0._dp
    do p = 1, self%num_phases
       if (btest(phases, p - 1)) then
          associate(phase => self%phase(p))
            ds = phase%density * phase%saturation
            ef = ef + ds * phase%internal_energy
          end associate
       end if
    end do

  end function fluid_energy

!------------------------------------------------------------------------

  function fluid_phase_mobilities(self) result(mobilities)
    !! Returns array containing the mobility for each phase.

    class(fluid_type), intent(in) :: self
    PetscReal :: mobilities(self%num_phases)
    ! Locals:
    PetscInt :: p, phases

    phases = nint(self%phase_composition)

    mobilities = 0._dp
    do p = 1, self%num_phases
       if (btest(phases, p - 1)) then
          mobilities(p) = self%phase(p)%mobility()
       end if
    end do

  end function fluid_phase_mobilities

!------------------------------------------------------------------------

  function fluid_phase_flow_fractions(self) result(f)
    !! Returns array containing the flow fractions for each
    !! phase. There are in proportion to the mobility of each phase,
    !! scaled to sum to 1.

    class(fluid_type), intent(in) :: self
    PetscReal :: f(self%num_phases)

    f = self%phase_mobilities()
    f = f / sum(f)

  end function fluid_phase_flow_fractions

!------------------------------------------------------------------------

  function fluid_component_flow_fractions(self, phase_flow_fractions) result(f)
    !! Returns array containing the flow fractions for each
    !! component, given the array of phase flow fractions (calculated
    !! using the flow_fractions() method).

    class(fluid_type), intent(in) :: self
    PetscReal, intent(in) :: phase_flow_fractions(self%num_phases)
    PetscReal :: f(self%num_components)
    ! Locals:
    PetscInt :: c, p, phases

    phases = nint(self%phase_composition)
    do c = 1, self%num_components
       f(c) = 0._dp
       do p = 1, self%num_phases
          if (btest(phases, p - 1)) then
             f(c) = f(c) + phase_flow_fractions(p) * self%phase(p)%mass_fraction(c)
          end if
       end do
    end do
    f = f / sum(f)

  end function fluid_component_flow_fractions

!------------------------------------------------------------------------

  PetscReal function fluid_specific_enthalpy(self, phase_flow_fractions) result(h)
    !! Returns total specific enthalpy, with contributions from all
    !! phases, given the array of phase flow fractions (calculated
    !! using the flow_fractions() method).

    class(fluid_type), intent(in) :: self
    PetscReal, intent(in) :: phase_flow_fractions(self%num_phases)
    ! Locals:
    PetscInt :: p, phases

    h = 0._dp

    phases = nint(self%phase_composition)
    do p = 1, self%num_phases
       if (btest(phases, p - 1)) then
          h = h + phase_flow_fractions(p) * &
               self%phase(p)%specific_enthalpy
       end if
    end do

  end function fluid_specific_enthalpy

!------------------------------------------------------------------------

  subroutine fluid_update_phase_composition(self, thermo)
    !! Updates fluid phase composition from thermodynamic region,
    !! pressure and temperature, according to specified thermodynamic
    !! formulation.

    use thermodynamics_module

    class(fluid_type), intent(in out) :: self
    class(thermodynamics_type), intent(in) :: thermo
    ! Locals:
    PetscInt :: region, phases

    region = nint(self%region)
    phases = thermo%phase_composition(region, self%pressure, &
         self%temperature)
    self%phase_composition = dble(phases)

  end subroutine fluid_update_phase_composition

!------------------------------------------------------------------------
! Fluid vector setup routine
!------------------------------------------------------------------------

  subroutine create_fluid_vector(dm, max_component_name_length, &
       component_names, max_phase_name_length, phase_names, &
       fluid_vec, range_start)
    !! Creates and returns global vector for fluid properties, with
    !! specified numbers of components and phases.

    use dm_utils_module, only: dm_set_data_layout, global_vec_range_start

    DM, intent(in) :: dm
    PetscInt, intent(in) :: max_component_name_length
    character(max_component_name_length), intent(in) :: component_names(:)
    PetscInt, intent(in) :: max_phase_name_length
    character(max_phase_name_length), intent(in) :: phase_names(:)
    Vec, intent(in out) :: fluid_vec
    PetscInt, intent(out) :: range_start
    ! Locals:
    PetscInt :: num_components, num_phases
    PetscInt :: p, i, j, dim
    PetscInt, allocatable :: num_field_components(:), field_dim(:)
    PetscErrorCode :: ierr
    PetscInt, parameter :: max_field_name_length = 40
    character(max_field_name_length), allocatable :: field_names(:)
    DM :: fluid_dm
    type(fluid_type) :: fluid

    num_components = size(component_names)
    num_phases = size(phase_names)
    call fluid%init(num_components, num_phases)

    allocate(num_field_components(fluid%dof), field_dim(fluid%dof), &
         field_names(fluid%dof))

    num_field_components = 1

    ! Assemble field names:
    field_names(1: num_fluid_variables - 1) = &
         fluid_variable_names(1: num_fluid_variables - 1) ! scalar bulk properties
    i = num_fluid_variables
    ! array bulk properties (partial pressures):
    do j = 1, num_components
       field_names(i) = trim(component_names(j)) // '_' // &
            trim(fluid_variable_names(num_fluid_variables))
       i = i + 1
    end do
    do p = 1, num_phases
       ! scalar phase properties:
       do j = 1, num_phase_variables - 1
          field_names(i) = trim(phase_names(p)) &
               // '_' // trim(phase_variable_names(j))
          i = i + 1
       end do
       ! array phase properties (mass fractions):
       do j = 1, num_components
          field_names(i) = trim(phase_names(p)) &
               // '_' // trim(component_names(j)) // '_' // &
               trim(phase_variable_names(num_phase_variables))
          i = i + 1
       end do
    end do

    call DMClone(dm, fluid_dm, ierr); CHKERRQ(ierr)

    call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)
    field_dim = dim
    call dm_set_data_layout(fluid_dm, num_field_components, field_dim, &
         field_names)

    call DMCreateGlobalVector(fluid_dm, fluid_vec, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(fluid_vec, "fluid", ierr); CHKERRQ(ierr)
    call global_vec_range_start(fluid_vec, range_start)

    deallocate(num_field_components, field_dim, field_names)
    call DMDestroy(fluid_dm, ierr); CHKERRQ(ierr)
    call fluid%destroy()

  end subroutine create_fluid_vector

!------------------------------------------------------------------------

end module fluid_module
