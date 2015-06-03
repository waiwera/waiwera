module flow_simulation_module
  !! Module for high-level representation of a flow simulation ODE.

  use ode_module
  use mesh_module
  use thermodynamics_module
  use eos_module

  implicit none

  private

#include <petsc-finclude/petsc.h90>

  PetscInt, parameter, public :: max_title_length = 120

  type, public, extends(ode_type) :: flow_simulation_type
     !! Simulation type.
     private
     character(max_title_length), public :: title
     Vec, public :: rock
     Vec, public :: fluid
     class(thermodynamics_type), allocatable, public :: thermo
     class(eos_type), allocatable, public :: eos
   contains
     private
     procedure, public :: init => flow_simulation_init
     procedure, public :: destroy => flow_simulation_destroy
     procedure, public :: lhs => flow_simulation_cell_balances
     procedure, public :: rhs => flow_simulation_cell_inflows
     procedure, public :: pre_eval => flow_simulation_fluid_properties
     procedure, public :: output => flow_simulation_output
  end type flow_simulation_type

contains

!------------------------------------------------------------------------

  subroutine flow_simulation_init(self, json)
    !! Initializes a flow simulation using data from the specified JSON object.

    use mpi_module
    use fson
    use fson_mpi_module
    use thermodynamics_setup_module, only: setup_thermodynamics
    use eos_setup_module, only: setup_eos
    use initial_module, only: setup_initial
    use fluid_module, only: setup_fluid_vector
    use rock_module, only: setup_rock_vector

    class(flow_simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    character(len = max_title_length), parameter :: default_title = ""
    PetscErrorCode :: ierr

    call fson_get_mpi(json, "title", default_title, self%title)
    call setup_thermodynamics(json, self%thermo)
    call setup_eos(json, self%thermo, self%eos)
    call self%mesh%init(json)
    call setup_labels(json, self%mesh%dm)
    call self%mesh%configure(self%eos%primary_variable_names)
    call DMCreateGlobalVector(self%mesh%dm, self%solution, ierr)
    CHKERRQ(ierr)
    call PetscObjectSetName(self%solution, "primary", ierr); CHKERRQ(ierr)
    call setup_initial(json, self%time, self%solution)
    call setup_fluid_vector(self%mesh%dm, self%eos%num_phases, &
         self%eos%num_components, self%fluid)
    call setup_rock_vector(json, self%mesh%dm, self%rock)

  end subroutine flow_simulation_init

!------------------------------------------------------------------------

  subroutine flow_simulation_destroy(self)
    !! Destroys the simulation.

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call VecDestroy(self%solution, ierr); CHKERRQ(ierr)
    call VecDestroy(self%fluid, ierr); CHKERRQ(ierr)
    call VecDestroy(self%rock, ierr); CHKERRQ(ierr)
    call self%mesh%destroy()
    call self%thermo%destroy()
    deallocate(self%thermo)
    deallocate(self%eos)

  end subroutine flow_simulation_destroy

!------------------------------------------------------------------------

  subroutine flow_simulation_cell_balances(self, t, y, lhs)
    !! Computes mass and energy balance for each cell, for the given
    !! primary thermodynamic variables and time.

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time (s)
    Vec, intent(in) :: y !! global primary variables vector
    Vec, intent(out) :: lhs
    ! Locals:
    PetscErrorCode :: ierr

    ! Dummy identity function:
    call VecCopy(y, lhs, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_cell_balances

!------------------------------------------------------------------------

  subroutine flow_simulation_cell_inflows(self, t, y, rhs)
    !! Computes net inflow into each cell, from flows through faces and
    !! source terms, for the given primary thermodynamic variables and
    !! time.

    use kinds_module

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time
    Vec, intent(in) :: y !! global primary variables vector
    Vec, intent(out) :: rhs
    ! Locals:
    PetscErrorCode :: ierr

    ! Dummy zero inflows:
    call VecSet(rhs, 0.0_dp, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_cell_inflows

!------------------------------------------------------------------------

  subroutine flow_simulation_fluid_properties(self, t, y)
    !! Computes fluid properties in all cells, based on the current time
    !! and primary thermodynamic variables.

    use dm_utils_module, only: section_offset
    use fluid_module, only: fluid_type

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time
    Vec, intent(in) :: y !! global primary variables vector
    ! Locals:
    PetscInt :: c
    DM :: dm_fluid
    PetscSection :: y_section, fluid_section
    PetscInt :: y_offset, fluid_offset
    PetscReal, pointer :: y_array(:), cell_primary(:)
    PetscReal, pointer :: fluid_array(:)
    type(fluid_type) :: fluid
    PetscInt :: region
    DMLabel :: ghost_label
    PetscInt :: ghost
    PetscErrorCode :: ierr

    ! Need read-only access to primary as it is locked by the SNES:
    call VecGetArrayReadF90(y, y_array, ierr); CHKERRQ(ierr)
    call DMGetDefaultSection(self%mesh%dm, y_section, ierr)
    CHKERRQ(ierr)

    call VecGetDM(self%fluid, dm_fluid, ierr); CHKERRQ(ierr)
    call DMGetDefaultSection(dm_fluid, fluid_section, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)

    call fluid%init(self%eos%num_components, self%eos%num_phases)

    call DMPlexGetLabel(self%mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    do c = self%mesh%start_cell, self%mesh%end_interior_cell - 1

       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then

          call section_offset(y_section, c, y_offset, ierr)
          CHKERRQ(ierr)
          cell_primary => y_array(y_offset : &
               y_offset + self%eos%num_primary_variables - 1)

          call section_offset(fluid_section, c, fluid_offset, ierr)
          CHKERRQ(ierr)

          call fluid%assign(fluid_array, fluid_offset)
          region = nint(fluid%region)

          call self%eos%fluid_properties(region, cell_primary, fluid)

       end if

    end do

    call VecRestoreArrayF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(y, y_array, ierr); CHKERRQ(ierr)
    call fluid%destroy()

  end subroutine flow_simulation_fluid_properties

!------------------------------------------------------------------------

  subroutine flow_simulation_output(self)
    !! Output from flow simulation.

    class(flow_simulation_type), intent(in out) :: self

    continue

  end subroutine flow_simulation_output

!------------------------------------------------------------------------

end module flow_simulation_module
