module fluid_module

  use kinds_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  PetscInt, parameter, public :: num_phase_variables = 6
  PetscInt, parameter, public :: num_fluid_variables = 4

  type phase_type
     !! Type for accessing local fluid properties for a particular phase.
     PetscReal, pointer :: density    !! Phase density
     PetscReal, pointer :: viscosity  !! Viscosity
     PetscReal, pointer :: saturation !! Phase saturation
     PetscReal, pointer :: relative_permeability !! Relative permeability
     PetscReal, pointer :: specific_enthalpy !! Specific enthalpy
     PetscReal, pointer :: internal_energy !! Internal energy
     PetscReal, pointer :: mass_fraction(:)
   contains
     private
     procedure, public :: init => phase_init
     procedure, public :: destroy => phase_destroy
     procedure, public :: dof => phase_dof
     procedure, public :: mobility => phase_mobility
  end type phase_type

  type fluid_type
     !! Type for accessing local fluid properties.
     private
     PetscReal, pointer, public :: pressure    !! Bulk pressure
     PetscReal, pointer, public :: temperature !! Temperature
     PetscReal, pointer, public :: region      !! Thermodynamic region
     PetscReal, pointer, public :: phase_composition   !! Phase composition
     type(phase_type), allocatable, public :: phase(:) !! Phase variables
     PetscInt, public :: num_phases !! Number of phases
     PetscInt, public :: num_concurrent_phases !! Number of phases that can be present at once
     PetscInt, public :: num_components !! Number of mass components
   contains
     private
     procedure, public :: init => fluid_init
     procedure, public :: assign => fluid_assign
     procedure, public :: destroy => fluid_destroy
     procedure, public :: dof => fluid_dof
     procedure, public :: component_density => fluid_component_density
     procedure, public :: energy => fluid_energy
     procedure, public :: flow_fractions => fluid_flow_fractions
     procedure, public :: energy_production => fluid_energy_production
     procedure, public :: update_phase_composition => &
          fluid_update_phase_composition
  end type fluid_type

  public :: fluid_type, setup_fluid_vector, initialise_fluid_regions

contains

!------------------------------------------------------------------------
! Phase procedures
!------------------------------------------------------------------------

  subroutine phase_init(self, num_components)
    !! Initialises a phase object.
    
    class(phase_type), intent(in out) :: self
    PetscInt, intent(in) :: num_components !! Number of fluid components

    allocate(self%mass_fraction(num_components))
    
  end subroutine phase_init

!------------------------------------------------------------------------

  subroutine phase_destroy(self)
    !! Destroys a phase object.

    class(phase_type), intent(in out) :: self

    nullify(self%density)
    nullify(self%viscosity)
    nullify(self%saturation)
    nullify(self%relative_permeability)
    nullify(self%specific_enthalpy)
    nullify(self%internal_energy)
    nullify(self%mass_fraction)

  end subroutine phase_destroy

!------------------------------------------------------------------------

  PetscInt function phase_dof(self)
    !! Returns number of degrees of freedom in the phase object.

    class(phase_type), intent(in) :: self

    phase_dof = num_phase_variables + size(self%mass_fraction)

  end function phase_dof

!------------------------------------------------------------------------

  PetscReal function phase_mobility(self)
    !! Returns mobility of the phase.

    class(phase_type), intent(in) :: self

    phase_mobility = self%relative_permeability * self%density / &
         self%viscosity

  end function phase_mobility

!------------------------------------------------------------------------
! Fluid procedures
!------------------------------------------------------------------------

  subroutine fluid_init(self, num_components, num_phases, &
       num_concurrent_phases)
    !! Initialises a fluid object.
    
    class(fluid_type), intent(in out) :: self
    PetscInt, intent(in) :: num_components !! Number of fluid components
    PetscInt, intent(in) :: num_phases     !! Number of fluid phases
    PetscInt, intent(in) :: num_concurrent_phases !! Number of phases that can be present at once
    ! Locals:
    PetscInt :: ip

    self%num_phases = num_phases
    self%num_concurrent_phases = num_concurrent_phases
    self%num_components = num_components
    allocate(self%phase(num_concurrent_phases))
    do ip = 1, num_concurrent_phases
       call self%phase(ip)%init(num_components)
    end do

  end subroutine fluid_init
    
!------------------------------------------------------------------------

  subroutine fluid_assign(self, data, offset)
    !! Assigns pointers in a fluid object to elements in the data array,
    !! starting from the specified offset.

    class(fluid_type), intent(in out) :: self
    PetscReal, target, intent(in) :: data(:)  !! fluid data array
    PetscInt, intent(in) :: offset  !! fluid array offset
    ! Locals:
    PetscInt :: i, ip

    self%pressure => data(offset)
    self%temperature => data(offset + 1)
    self%region => data(offset + 2)
    self%phase_composition => data(offset + 3)
    
    i = offset + num_fluid_variables
    do ip = 1, self%num_concurrent_phases
       self%phase(ip)%density => data(i)
       self%phase(ip)%viscosity => data(i+1)
       self%phase(ip)%saturation => data(i+2)
       self%phase(ip)%relative_permeability => data(i+3)
       self%phase(ip)%specific_enthalpy => data(i+4)
       self%phase(ip)%internal_energy => data(i+5)
       self%phase(ip)%mass_fraction => data(i+6: i+6 + &
            self%num_components-1)
       i = i + self%phase(ip)%dof()
    end do

  end subroutine fluid_assign

!------------------------------------------------------------------------

  subroutine fluid_destroy(self)
    !! Destroys a fluid object.
    
    class(fluid_type), intent(in out) :: self
    ! Locals:
    PetscInt :: ip

    nullify(self%pressure)
    nullify(self%temperature)
    nullify(self%region)
    nullify(self%phase_composition)

    do ip = 1, self%num_concurrent_phases
       call self%phase(ip)%destroy()
    end do
    deallocate(self%phase)

  end subroutine fluid_destroy

!------------------------------------------------------------------------

  PetscInt function fluid_dof(self)
    !! Returns number of degrees of freedom in a fluid object.

    class(fluid_type), intent(in) :: self
    ! Locals:
    PetscInt :: ip

    fluid_dof = num_fluid_variables
    do ip = 1, self%num_concurrent_phases
       fluid_dof = fluid_dof + self%phase(ip)%dof()
    end do

  end function fluid_dof

!------------------------------------------------------------------------

  function fluid_component_density(self, phase_index) result(d)
    !! Returns total fluid density for each mass component, over all
    !! phases.

    class(fluid_type), intent(in) :: self
    PetscInt, intent(in) :: phase_index(:)
    PetscReal :: d(self%num_components)
    ! Locals:
    PetscInt :: p, ip, c, phases
    PetscReal :: ds

    phases = nint(self%phase_composition)
    d = 0._dp
    do p = 1, self%num_phases
       if (btest(phases, p - 1)) then
          ip = phase_index(p)
          ds = self%phase(ip)%density * self%phase(ip)%saturation
          do c = 1, self%num_components
             d(c) = d(c) + ds * self%phase(ip)%mass_fraction(c)
          end do
       end if
    end do

  end function fluid_component_density

!------------------------------------------------------------------------

  PetscReal function fluid_energy(self, phase_index) result(ef)
    !! Returns total fluid energy density.

    class(fluid_type), intent(in) :: self
    PetscInt, intent(in) :: phase_index(:)
    ! Locals:
    PetscInt :: p, ip, phases
    PetscReal :: ds

    phases = nint(self%phase_composition)
    ef = 0._dp
    do p = 1, self%num_phases
       if (btest(phases, p - 1)) then
          ip = phase_index(p)
          ds = self%phase(ip)%density * self%phase(ip)%saturation
          ef = ef + ds * self%phase(ip)%internal_energy
       end if
    end do

  end function fluid_energy

!------------------------------------------------------------------------

  function fluid_flow_fractions(self, phase_index) result(f)
    !! Returns array containing the flow fractions for each
    !! phase. There are in proportion to the mobility of each phase,
    !! scaled to sum to 1.

    class(fluid_type), intent(in) :: self
    PetscInt, intent(in) :: phase_index(:)
    PetscReal :: f(self%num_concurrent_phases)
    ! Locals:
    PetscInt :: p, ip, phases

    phases = nint(self%phase_composition)

    do p = 1, self%num_phases
       if (btest(phases, p - 1)) then
          ip = phase_index(p)
          f(ip) = self%phase(ip)%mobility()
       else
          f(ip) = 0._dp
       end if
    end do
    f = f / sum(f)

  end function fluid_flow_fractions

!------------------------------------------------------------------------

  subroutine fluid_energy_production(self, source, phase_index, isothermal)
    !! If source array contains production, and EOS is
    !! non-isothermal, calculate associated energy production.

    class(fluid_type), intent(in) :: self
    PetscReal, target, intent(in out) :: source(:)
    PetscInt, intent(in) :: phase_index(:)
    PetscBool, intent(in) :: isothermal
    ! Locals:
    PetscInt :: p, ip, np, phases, c
    PetscReal :: flow_fractions(self%num_concurrent_phases), hc
    PetscReal, pointer :: q, qenergy

    if (.not. isothermal) then

       np = size(source)
       qenergy => source(np)
       phases = nint(self%phase_composition)
       flow_fractions = self%flow_fractions(phase_index)

       do c = 1, self%num_components
          q => source(c)
          if (q < 0._dp) then
             hc = 0._dp
             do p = 1, self%num_phases
                if (btest(phases, p - 1)) then
                   ip = phase_index(p)
                   hc = hc + flow_fractions(ip) * &
                        self%phase(ip)%specific_enthalpy * &
                        self%phase(ip)%mass_fraction(c)
                end if
             end do
             qenergy = qenergy + q * hc
          end if
       end do

    end if

  end subroutine fluid_energy_production

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

  subroutine setup_fluid_vector(dm, num_concurrent_phases, num_components, &
       fluid, range_start)
    !! Sets up global vector for fluid properties, with specified
    !! numbers of components and phases.

    use dm_utils_module, only: set_dm_data_layout, global_vec_range_start

    DM, intent(in) :: dm
    PetscInt, intent(in) :: num_concurrent_phases, num_components
    Vec, intent(out) :: fluid
    PetscInt, intent(out) :: range_start
    ! Locals:
    PetscInt :: num_vars
    PetscInt, allocatable :: num_field_components(:), field_dim(:)
    DM :: dm_fluid
    PetscErrorCode :: ierr

    num_vars = num_fluid_variables + num_concurrent_phases * &
         (num_phase_variables + num_components)

    allocate(num_field_components(num_vars), field_dim(num_vars))

    ! All fluid variables are scalars defined on cells:
    num_field_components = 1
    field_dim = 3

    call DMClone(dm, dm_fluid, ierr); CHKERRQ(ierr)

    call set_dm_data_layout(dm_fluid, num_field_components, field_dim)

    call DMCreateGlobalVector(dm_fluid, fluid, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(fluid, "fluid", ierr); CHKERRQ(ierr)
    call global_vec_range_start(fluid, range_start)

    deallocate(num_field_components, field_dim)
    call DMDestroy(dm_fluid, ierr); CHKERRQ(ierr)

  end subroutine setup_fluid_vector

!------------------------------------------------------------------------

  subroutine initialise_fluid_regions(dm, fluid, start_cell, end_cell, &
       range_start, num_components, num_phases, num_concurrent_phases)
    !! Initialise fluid regions in each cell. For now, just assume all
    !! cells are initially region 1 (liquid).

    use dm_utils_module, only: global_vec_section, global_section_offset

    DM, intent(in) :: dm
    Vec, intent(in out) :: fluid
    PetscInt, intent(in) :: start_cell, end_cell
    PetscInt, intent(in) :: range_start
    PetscInt, intent(in) :: num_components
    PetscInt, intent(in) :: num_phases, num_concurrent_phases
    ! Locals:
    PetscSection :: fluid_section
    PetscReal, pointer :: fluid_array(:)
    type(fluid_type) :: f
    DMLabel :: ghost_label
    PetscInt :: ghost, fluid_offset, c
    PetscErrorCode :: ierr

    call global_vec_section(fluid, fluid_section)
    call VecGetArrayF90(fluid, fluid_array, ierr); CHKERRQ(ierr)
    call f%init(num_components, num_phases, num_concurrent_phases)

    call DMPlexGetLabel(dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    do c = start_cell, end_cell - 1
       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then
          call global_section_offset(fluid_section, c, &
               range_start, fluid_offset, ierr); CHKERRQ(ierr)
          call f%assign(fluid_array, fluid_offset)
          f%region = 1
       end if
    end do

    call f%destroy()
    call VecRestoreArrayF90(fluid, fluid_array, ierr); CHKERRQ(ierr)

  end subroutine initialise_fluid_regions

!------------------------------------------------------------------------

end module fluid_module
