module simulation_module
  !! Module for high-level representation of a simulation.

  use mesh_module
  use timestepping_module
  use thermodynamics_module
  use eos_module

  implicit none

  private

#include <petsc-finclude/petsc.h90>

  PetscInt, parameter, public :: max_filename_length = 200
  PetscInt, parameter, public :: max_title_length = 120

  type, public :: simulation_type
     !! Simulation type.
     private
     type(mesh_type), public :: mesh
     Vec, public :: initial
     Vec, public :: rock
     Vec, public :: fluid
     type(timestepper_type), public :: timestepper
     class(thermodynamics_type), allocatable, public :: thermo
     class(eos_type), allocatable, public :: eos
     character(max_filename_length), public :: input_filename
     character(max_title_length), public :: title
   contains
     private
     procedure, public :: init => simulation_init
     procedure, public :: run => simulation_run
     procedure, public :: destroy => simulation_destroy
  end type simulation_type

  type(simulation_type), public :: sim

contains

!------------------------------------------------------------------------

  subroutine simulation_init(self, filename, json_str)
    !! Initializes a simulation using data from the input file with 
    !! specified name.

    use mpi_module
    use fson
    use fson_mpi_module
    use thermodynamics_setup_module, only: setup_thermodynamics
    use eos_setup_module, only: setup_eos
    use initial_module, only: setup_initial
    use fluid_module, only: setup_fluid_vector
    use rock_module, only: setup_rock_vector

    class(simulation_type), intent(in out) :: self
    character(*), intent(in), optional :: filename !! Input file name
    character(*), intent(in), optional :: json_str !! JSON string for alternative input
    ! Locals:
    type(fson_value), pointer :: json
    character(len = max_title_length), parameter :: default_title = ""

    if (present(filename)) then
       self%input_filename = filename
    end if

    if (mpi%rank == mpi%input_rank) then
       if (present(json_str)) then
          json => fson_parse(str = json_str)
       else
          json => fson_parse(file = filename)
       end if
    end if

    call fson_get_mpi(json, "title", default_title, self%title)
    call setup_thermodynamics(json, self%thermo)
    call setup_eos(json, self%thermo, self%eos)
    call self%mesh%init(json)
    call setup_labels(json, self%mesh%dm)
    call self%mesh%configure(self%eos%primary_variable_names)
    call setup_initial(json, self%mesh%dm, self%initial)
    call setup_fluid_vector(self%mesh%dm, self%eos%num_phases, &
         self%eos%num_components, self%fluid)
    call setup_rock_vector(json, self%mesh%dm, self%rock)
    call setup_timestepper(json, self%mesh%dm, self%initial, &
         self%timestepper)

    self%timestepper%steps%lhs_func => simulation_cell_balances
    self%timestepper%steps%rhs_func => simulation_cell_inflows
    self%timestepper%steps%pre_eval_proc => simulation_fluid_properties

    if (mpi%rank == mpi%input_rank) then
       call fson_destroy(json)
    end if

  end subroutine simulation_init

!------------------------------------------------------------------------

  subroutine simulation_run(self)
    !! Runs the simulation.

    class(simulation_type), intent(in out) :: self

    call self%timestepper%run()

    ! maybe some final output?

  end subroutine simulation_run

!------------------------------------------------------------------------

  subroutine simulation_destroy(self)
    !! Destroys the simulation.

    class(simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call VecDestroy(self%initial, ierr); CHKERRQ(ierr)
    call VecDestroy(self%fluid, ierr); CHKERRQ(ierr)
    call VecDestroy(self%rock, ierr); CHKERRQ(ierr)
    call self%mesh%destroy()
    call self%thermo%destroy()
    deallocate(self%thermo)
    deallocate(self%eos)
    call self%timestepper%destroy()

  end subroutine simulation_destroy

!------------------------------------------------------------------------
! Cell mass and energy balances
!------------------------------------------------------------------------

  subroutine simulation_cell_balances(t, primary, balance)
    !! Computes mass and energy balance for each cell, for the given
    !! primary thermodynamic variables and time.
    PetscReal, intent(in) :: t !! time (s)
    Vec, intent(in) :: primary !! global primary variables vector
    Vec, intent(out) :: balance
    ! Locals:
    PetscErrorCode :: ierr

    ! Dummy identity function:
    call VecCopy(primary, balance, ierr); CHKERRQ(ierr)

  end subroutine simulation_cell_balances

!------------------------------------------------------------------------
! Cell net inflows and sources
!------------------------------------------------------------------------

  subroutine simulation_cell_inflows(t, primary, inflow)
    !! Computes net inflow into each cell, from flows through faces and
    !! source terms, for the given primary thermodynamic variables and
    !! time.

    use kinds_module

    PetscReal, intent(in) :: t !! time
    Vec, intent(in) :: primary !! global primary variables vector
    Vec, intent(out) :: inflow
    ! Locals:
    PetscErrorCode :: ierr

    ! Dummy zero inflows:
    call VecSet(inflow, 0.0_dp, ierr); CHKERRQ(ierr)

  end subroutine simulation_cell_inflows

!------------------------------------------------------------------------
! Fluid properties
!------------------------------------------------------------------------

  subroutine simulation_fluid_properties(t, primary)
    !! Computes fluid properties in all cells, based on the current time
    !! and primary thermodynamic variables.

    use dm_utils_module, only: section_offset
    use fluid_module, only: fluid_type

    PetscReal, intent(in) :: t !! time
    Vec, intent(in) :: primary !! global primary variables vector
    ! Locals:
    PetscInt :: c
    DM :: dm_fluid
    PetscSection :: primary_section, fluid_section
    PetscInt :: primary_offset, fluid_offset
    PetscReal, pointer :: primary_array(:), cell_primary(:)
    PetscReal, pointer :: fluid_array(:)
    type(fluid_type) :: fluid
    PetscInt :: region
    DMLabel :: ghost_label
    PetscInt :: ghost
    PetscErrorCode :: ierr

    ! Need read-only access to primary as it is locked by the SNES:
    call VecGetArrayReadF90(primary, primary_array, ierr); CHKERRQ(ierr)
    call DMGetDefaultSection(sim%mesh%dm, primary_section, ierr)
    CHKERRQ(ierr)

    call VecGetDM(sim%fluid, dm_fluid, ierr); CHKERRQ(ierr)
    call DMGetDefaultSection(dm_fluid, fluid_section, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(sim%fluid, fluid_array, ierr); CHKERRQ(ierr)

    call DMPlexGetLabel(sim%mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    call fluid%init(sim%eos%num_components, sim%eos%num_phases)

    do c = sim%mesh%start_cell, sim%mesh%end_interior_cell - 1

       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then

          call section_offset(primary_section, c, primary_offset, ierr)
          CHKERRQ(ierr)
          cell_primary => primary_array(primary_offset : &
               primary_offset + sim%eos%num_primary_variables -1)

          call section_offset(fluid_section, c, fluid_offset, ierr)
          CHKERRQ(ierr)

          call fluid%assign(fluid_array, fluid_offset)
          region = nint(fluid%region)

          call sim%eos%fluid_properties(region, cell_primary, fluid)

       end if

    end do

    call VecRestoreArrayF90(sim%fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(primary, primary_array, ierr); CHKERRQ(ierr)
    call fluid%destroy()

  end subroutine simulation_fluid_properties

!------------------------------------------------------------------------

end module simulation_module
