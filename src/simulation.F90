module simulation_module
  !! Module for high-level representation of a simulation.

  use mpi_module
  use kinds_module
  use mesh_module
  use timestepping_module
  use thermodynamics_module
  use IAPWS_module
  use IFC67_module
  use eos_w_module
  use eos_module
  use fson
  use fson_mpi_module

  implicit none

  private

#include <petsc-finclude/petsc.h90>

  integer, parameter, public :: max_filename_length = 200
  integer, parameter, public :: max_title_length = 120

  type, public :: simulation_type
     !! Simulation type.
     private
     type(mesh_type), public :: mesh
     Vec, public :: initial
     Vec, public :: rock
     Vec, public :: fluid
     type(timestepper_type), public :: timestepper
     class(thermodynamics_type), pointer, public :: thermo
     class(eos_type), pointer, public :: eos
     character(max_filename_length), public :: input_filename
     character(max_title_length), public :: title
   contains
     private
     procedure :: setup_thermodynamics => simulation_setup_thermodynamics
     procedure :: setup_eos => simulation_setup_eos
     procedure, public :: init => simulation_init
     procedure, public :: run => simulation_run
     procedure, public :: destroy => simulation_destroy
  end type simulation_type

  type(simulation_type), public :: sim

contains

!------------------------------------------------------------------------

  subroutine simulation_setup_thermodynamics(self, json)
    !! Reads simulation thermodynamic formulation from JSON input file.
    !! If not present, a default value is assigned.

    use utils_module, only : str_to_lower

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    integer, parameter :: max_thermo_ID_length = 8
    character(max_thermo_ID_length) :: thermo_ID
    character(max_thermo_ID_length), parameter :: &
         default_thermo_ID = "IAPWS"

    call fson_get_mpi(json, "thermodynamics", default_thermo_ID, thermo_ID)
    thermo_ID = str_to_lower(thermo_ID)

    select case (thermo_ID)
    case ("ifc67")
       self%thermo => IFC67
    case default
       self%thermo => IAPWS
    end select

    call self%thermo%init()

  end subroutine simulation_setup_thermodynamics

!------------------------------------------------------------------------

  subroutine simulation_setup_eos(self, json)
    !! Reads simulation equation of state from JSON input file.
    !! If not present, a default value is assigned.

    use utils_module, only : str_to_lower

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    integer, parameter :: max_eos_ID_length = 12
    character(max_eos_ID_length) :: eos_ID
    character(max_eos_ID_length), parameter :: &
         default_eos_ID = "W"

    call fson_get_mpi(json, "eos", default_eos_ID, eos_ID)
    eos_ID = str_to_lower(eos_ID)

    select case (eos_ID)
    case ("ew")
       self%eos => eos_w  ! change to eos_ew when it's ready
    case default
       self%eos => eos_w
    end select

    call self%eos%init(self%thermo)

  end subroutine simulation_setup_eos

!------------------------------------------------------------------------

  subroutine simulation_init(self, filename, json_str)
    !! Initializes a simulation using data from the input file with 
    !! specified name.

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
    call self%setup_thermodynamics(json)
    call self%setup_eos(json)
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
    nullify(self%thermo)
    nullify(self%eos)
    call self%timestepper%destroy()

  end subroutine simulation_destroy

!------------------------------------------------------------------------

  subroutine simulation_cell_balances(t, primary, balance)
    !! Computes mass and energy balance for each cell, for the given
    !! primary thermodynamic variables and time.
    PetscReal, intent(in) :: t !! time (s)
    Vec, intent(in) :: primary
    Vec, intent(out) :: balance
    ! Locals:
    PetscErrorCode :: ierr

    ! Dummy identity function:
    call VecCopy(primary, balance, ierr); CHKERRQ(ierr)

  end subroutine simulation_cell_balances

!------------------------------------------------------------------------

  subroutine simulation_cell_inflows(t, primary, inflow)
    !! Computes net inflow into each cell, from flows through faces and
    !! source terms, for the given primary thermodynamic variables and
    !! time.

    PetscReal, intent(in) :: t
    Vec, intent(in) :: primary
    Vec, intent(out) :: inflow
    ! Locals:
    PetscErrorCode :: ierr

    ! Dummy zero inflows:
    call VecSet(inflow, 0.0_dp, ierr); CHKERRQ(ierr)

  end subroutine simulation_cell_inflows

!------------------------------------------------------------------------

end module simulation_module
