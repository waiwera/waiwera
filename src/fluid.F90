module fluid_module
  !! Module for accessing fluid properties.

  use kinds_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  PetscInt, parameter, public :: num_fluid_variables = 4
  PetscInt, parameter, public :: num_phase_variables = 6  ! (excluding mass fractions)
  PetscInt, parameter, public :: max_fluid_variable_name_length = 11
  character(max_fluid_variable_name_length), public :: &
       fluid_variable_names(num_fluid_variables) = [ &
       "pressure   ", "temperature", &
       "region     ", "phases     "]
  PetscInt, parameter, public :: max_phase_variable_name_length = 21
  character(max_phase_variable_name_length), public :: &
       phase_variable_names(num_phase_variables) = [ &
       "density              ", "viscosity            ", &
       "saturation           ", "relative_permeability", &
       "specific_enthalpy    ", "internal_energy      "]

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
     procedure, public :: destroy => phase_destroy
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
     PetscInt, public :: num_components !! Number of mass components
     PetscInt, public :: dof !! Number of degrees of freedom
     PetscInt, public :: phase_dof !! Number of degrees of freedom per phase
   contains
     private
     procedure, public :: init => fluid_init
     procedure, public :: assign => fluid_assign
     procedure, public :: destroy => fluid_destroy
     procedure, public :: component_density => fluid_component_density
     procedure, public :: energy => fluid_energy
     procedure, public :: flow_fractions => fluid_flow_fractions
     procedure, public :: energy_production => fluid_energy_production
     procedure, public :: update_phase_composition => &
          fluid_update_phase_composition
  end type fluid_type

  public :: fluid_type, setup_fluid_vector

contains

!------------------------------------------------------------------------
! Phase procedures
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

  PetscReal function phase_mobility(self)
    !! Returns mobility of the phase.

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

    self%phase_dof = num_phase_variables + self%num_components
    self%dof = num_fluid_variables + self%num_phases * self%phase_dof

  end subroutine fluid_init
    
!------------------------------------------------------------------------

  subroutine fluid_assign(self, data, offset)
    !! Assigns pointers in a fluid object to elements in the data array,
    !! starting from the specified offset.

    use profiling_module, only: assign_pointers_event

    class(fluid_type), intent(in out) :: self
    PetscReal, target, intent(in) :: data(:)  !! fluid data array
    PetscInt, intent(in) :: offset  !! fluid array offset
    ! Locals:
    PetscInt :: i, p
    PetscErrorCode :: ierr

    call PetscLogEventBegin(assign_pointers_event, ierr); CHKERRQ(ierr)

    self%pressure => data(offset)
    self%temperature => data(offset + 1)
    self%region => data(offset + 2)
    self%phase_composition => data(offset + 3)
    
    i = offset + num_fluid_variables
    do p = 1, self%num_phases
       associate(phase => self%phase(p))
         phase%density => data(i)
         phase%viscosity => data(i+1)
         phase%saturation => data(i+2)
         phase%relative_permeability => data(i+3)
         phase%specific_enthalpy => data(i+4)
         phase%internal_energy => data(i+5)
         phase%mass_fraction => data(i+6: i+6 + self%num_components-1)
         i = i + self%phase_dof
       end associate
    end do

    call PetscLogEventEnd(assign_pointers_event, ierr); CHKERRQ(ierr)

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

  function fluid_flow_fractions(self) result(f)
    !! Returns array containing the flow fractions for each
    !! phase. There are in proportion to the mobility of each phase,
    !! scaled to sum to 1.

    class(fluid_type), intent(in) :: self
    PetscReal :: f(self%num_phases)
    ! Locals:
    PetscInt :: p, phases

    phases = nint(self%phase_composition)

    f = 0._dp
    do p = 1, self%num_phases
       if (btest(phases, p - 1)) then
          f(p) = self%phase(p)%mobility()
       end if
    end do
    f = f / sum(f)

  end function fluid_flow_fractions

!------------------------------------------------------------------------

  subroutine fluid_energy_production(self, source, isothermal)
    !! If source array contains production, and EOS is
    !! non-isothermal, calculate associated energy production.

    class(fluid_type), intent(in) :: self
    PetscReal, intent(in out) :: source(:)
    PetscBool, intent(in) :: isothermal
    ! Locals:
    PetscInt :: p, np, phases, c
    PetscReal :: flow_fractions(self%num_phases), hc

    if (.not. isothermal) then
       np = size(source)
       associate (qenergy => source(np))

         phases = nint(self%phase_composition)
         flow_fractions = self%flow_fractions()

         do c = 1, self%num_components
            associate(q => source(c))
              if (q < 0._dp) then
                 hc = 0._dp
                 do p = 1, self%num_phases
                    if (btest(phases, p - 1)) then
                       associate(phase => self%phase(p))
                         hc = hc + flow_fractions(p) * &
                              phase%specific_enthalpy * &
                              phase%mass_fraction(c)
                       end associate
                    end if
                 end do
                 qenergy = qenergy + q * hc
              end if
            end associate
         end do

       end associate
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

  subroutine setup_fluid_vector(dm, max_component_name_length, &
       component_names, max_phase_name_length, phase_names, &
       fluid_vec, range_start)
    !! Sets up global vector and DM for fluid properties, with specified
    !! numbers of components and phases.

    use dm_utils_module, only: set_dm_data_layout, global_vec_range_start

    DM, intent(in) :: dm
    PetscInt, intent(in) :: max_component_name_length
    character(max_component_name_length), intent(in) :: component_names(:)
    PetscInt, intent(in) :: max_phase_name_length
    character(max_phase_name_length), intent(in) :: phase_names(:)
    Vec, intent(out) :: fluid_vec
    PetscInt, intent(out) :: range_start
    ! Locals:
    PetscInt :: num_components, num_phases
    PetscInt :: p, i, j
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
    field_dim = 3

    ! Assemble field names:
    field_names(1: num_fluid_variables) = fluid_variable_names
    i = num_fluid_variables + 1
    do p = 1, num_phases
       do j = 1, num_phase_variables
          field_names(i) = trim(phase_names(p)) &
               // '_' // trim(phase_variable_names(j))
          i = i + 1
       end do
       do j = 1, num_components
          field_names(i) = trim(phase_names(p)) &
               // '_' // trim(component_names(j)) // '_mass_fraction'
          i = i + 1
       end do
    end do

    call DMClone(dm, fluid_dm, ierr); CHKERRQ(ierr)

    call set_dm_data_layout(fluid_dm, num_field_components, field_dim, &
         field_names)

    call DMCreateGlobalVector(fluid_dm, fluid_vec, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(fluid_vec, "fluid", ierr); CHKERRQ(ierr)
    call global_vec_range_start(fluid_vec, range_start)

    deallocate(num_field_components, field_dim, field_names)
    call DMDestroy(fluid_dm, ierr); CHKERRQ(ierr)
    call fluid%destroy()

  end subroutine setup_fluid_vector

!------------------------------------------------------------------------

end module fluid_module
