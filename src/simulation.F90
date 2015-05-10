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
     procedure :: read_title => simulation_read_title
     procedure :: read_thermodynamics => simulation_read_thermodynamics
     procedure :: read_eos => simulation_read_eos
     procedure :: read_initial => simulation_read_initial
     procedure :: setup_fluid => simulation_setup_fluid
     procedure :: read_rock_types => simulation_read_rock_types
     procedure :: read_rock_properties => simulation_read_rock_properties
     procedure :: read_timestepping => simulation_read_timestepping
     procedure :: setup_rocktype_labels => simulation_setup_rocktype_labels
     procedure :: setup_boundary_labels => simulation_setup_boundary_labels
     procedure :: setup_labels => simulation_setup_labels
     procedure, public :: init => simulation_init
     procedure, public :: run => simulation_run
     procedure, public :: destroy => simulation_destroy
  end type simulation_type

contains

!------------------------------------------------------------------------

  subroutine simulation_read_title(self, json)
    !! Reads simulation title from JSON input file.
    !! If not present, a default value is assigned.

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    character(len = max_title_length), parameter :: default_title = ""

    call fson_get_mpi(json, "title", default_title, self%title)

  end subroutine simulation_read_title

!------------------------------------------------------------------------

  subroutine simulation_read_thermodynamics(self, json)
    !! Reads simulation thermodynamic formulation from JSON input file.
    !! If not present, a default value is assigned.

    use utils_module, only : str_to_lower

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscErrorCode :: ierr
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

  end subroutine simulation_read_thermodynamics

!------------------------------------------------------------------------

  subroutine simulation_read_eos(self, json)
    !! Reads simulation equation of state from JSON input file.
    !! If not present, a default value is assigned.

    use utils_module, only : str_to_lower

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscErrorCode :: ierr
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

  end subroutine simulation_read_eos

!------------------------------------------------------------------------

  subroutine simulation_read_initial(self, json)
    !! Reads initial conditions from JSON input. These may be specified
    !! as a constant value or as an array. The array may contain a complete
    !! of initial conditions for all cells, or if a shorter array is 
    !! given, this is repeated over initial conditions vector.

    use fson_value_m, only : TYPE_REAL, TYPE_INTEGER, TYPE_ARRAY

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    type(fson_value), pointer :: initial
    PetscErrorCode :: ierr
    real(dp) :: const_initial_value
    integer :: int_const_initial_value
    integer, allocatable :: indices(:)
    real(dp), allocatable :: initial_input(:), initial_data(:)
    integer :: i, np, count
    logical :: const
    real(dp), parameter :: default_initial_value = 0.0_dp

    call DMCreateGlobalVector(self%mesh%dm, self%initial, ierr); CHKERRQ(ierr)
    call VecGetSize(self%initial, count, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(self%initial, "initial", ierr); CHKERRQ(ierr)

    const = .true.

    if (fson_has_mpi(json, "initial")) then

       select case (fson_type_mpi(json, "initial"))
          case (TYPE_REAL)
             call fson_get_mpi(json, "initial", val = const_initial_value)
          case (TYPE_INTEGER)
             call fson_get_mpi(json, "initial", val = int_const_initial_value)
             const_initial_value = real(int_const_initial_value)
          case (TYPE_ARRAY)
             const = .false.
             call fson_get_mpi(json, "initial", val = initial_input)
             np = size(initial_input)
             if (np >= count) then
                initial_data = initial_input(1:count)
             else ! repeat input over array:
                do i = 1, np
                   initial_data(i:count:np) = initial_input(i)
                end do
             end if
             deallocate(initial_input)
       end select
    else
       const_initial_value = default_initial_value
    end if

    if (const) then
       call VecSet(self%initial, const_initial_value, ierr); CHKERRQ(ierr)
    else
       allocate(indices(count))
       do i = 1, count
          indices(i) = i-1
       end do
       call VecSetValues(self%initial, count, indices, &
            initial_data, INSERT_VALUES, ierr); CHKERRQ(ierr)
       call VecAssemblyBegin(self%initial, ierr); CHKERRQ(ierr)
       call VecAssemblyEnd(self%initial, ierr); CHKERRQ(ierr)
       deallocate(indices, initial_data)
    end if

  end subroutine simulation_read_initial

!------------------------------------------------------------------------

  subroutine simulation_setup_fluid(self)
    !! Sets up global vector for fluid properties.

    use fluid_module, only: num_fluid_variables, num_phase_variables

    class(simulation_type), intent(in out) :: self
    ! Locals:
    PetscInt :: num_vars
    PetscInt, allocatable :: num_components(:), field_dim(:)
    DM :: dm_fluid
    PetscErrorCode :: ierr

    num_vars = num_fluid_variables + self%eos%num_phases * &
         (num_phase_variables + self%eos%num_components)

    allocate(num_components(num_vars), field_dim(num_vars))

    ! All fluid variables are scalars defined on cells:
    num_components = 1
    field_dim = 3

    call DMClone(self%mesh%dm, dm_fluid, ierr); CHKERRQ(ierr)

    call set_dm_data_layout(dm_fluid, num_components, field_dim)

    call DMCreateGlobalVector(dm_fluid, self%fluid, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(self%fluid, "fluid", ierr); CHKERRQ(ierr)

    deallocate(num_components, field_dim)
    call DMDestroy(dm_fluid, ierr); CHKERRQ(ierr)

  end subroutine simulation_setup_fluid

!------------------------------------------------------------------------

  subroutine simulation_read_rock_types(self, json)
    !! Reads rock properties from rock types in JSON input.

    use rock_module

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscInt :: num_rocktypes, ir, ic, c, num_cells, offset, ghost
    DM :: dm_rock
    type(fson_value), pointer :: rocktypes, r
    IS :: rock_IS
    PetscInt, pointer :: rock_cells(:)
    DMLabel :: ghost_label
    type(rock_type) :: rock
    character(max_rockname_length) :: name
    PetscReal :: porosity, density, specific_heat, heat_conductivity
    PetscReal, allocatable :: permeability(:)
    PetscReal, pointer :: rock_array(:)
    PetscSection :: section
    PetscErrorCode :: ierr

    call VecGetDM(self%rock, dm_rock, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(self%rock, rock_array, ierr); CHKERRQ(ierr)
    call DMGetDefaultSection(dm_rock, section, ierr); CHKERRQ(ierr)
    call DMPlexGetLabel(self%mesh%dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)
    
    call fson_get_mpi(json, "rock.types", rocktypes)
    num_rocktypes = fson_value_count_mpi(rocktypes, ".")
    do ir = 1, num_rocktypes
       r => fson_value_get_mpi(rocktypes, ir)
       call fson_get_mpi(r, "name", "", name)
       call fson_get_mpi(r, "permeability", default_permeability, permeability)
       call fson_get_mpi(r, "heat conductivity", default_heat_conductivity, heat_conductivity)
       call fson_get_mpi(r, "porosity", default_porosity, porosity)
       call fson_get_mpi(r, "density", default_density, density)
       call fson_get_mpi(r, "specific heat", default_specific_heat, specific_heat)
       call DMPlexGetStratumIS(self%mesh%dm, rocktype_label_name, ir, rock_IS, ierr); CHKERRQ(ierr)
       if (rock_IS /= 0) then
          call ISGetIndicesF90(rock_IS, rock_cells, ierr); CHKERRQ(ierr)
          num_cells = size(rock_cells)
          do ic = 1, num_cells
             c = rock_cells(ic)
             call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
             if (ghost < 0) then
                call section_offset(section, c, offset, ierr); CHKERRQ(ierr)
                call rock%assign(rock_array, offset)
                rock%permeability = permeability
                rock%heat_conductivity = heat_conductivity
                rock%porosity = porosity
                rock%density = density
                rock%specific_heat = specific_heat
             end if
          end do
          call ISRestoreIndicesF90(rock_IS, rock_cells, ierr); CHKERRQ(ierr)
       end if
    end do
    call rock%destroy()
    call VecRestoreArrayF90(self%rock, rock_array, ierr); CHKERRQ(ierr)

  end subroutine simulation_read_rock_types

!------------------------------------------------------------------------

  subroutine simulation_read_rock_properties(self, json)
    !! Reads rock properties from JSON input.

    use rock_module, only: rock_variable_num_components, &
         rock_variable_dim, rock_variable_names

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    DM :: dm_rock
    PetscErrorCode :: ierr

    call DMClone(self%mesh%dm, dm_rock, ierr); CHKERRQ(ierr)

    call set_dm_data_layout(dm_rock, rock_variable_num_components, &
         rock_variable_dim, rock_variable_names)

    call DMCreateGlobalVector(dm_rock, self%rock, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(self%rock, "rock", ierr); CHKERRQ(ierr)

    ! TODO: set default rock properties everywhere here? in case of cells with
    ! no properties specified

    if (fson_has_mpi(json, "rock")) then

       if (fson_has_mpi(json, "rock.types")) then

          call self%read_rock_types(json)

       else
          ! other types of rock initialization here- TODO
          ! e.g. read from a rock section in initial conditions HDF5 file
       end if

    end if

    call DMDestroy(dm_rock, ierr); CHKERRQ(ierr)

  end subroutine simulation_read_rock_properties

!------------------------------------------------------------------------

  subroutine simulation_read_timestepping(self, json)
    !! Reads time stepping data from JSON input.

    use utils_module, only : str_to_lower

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    type(fson_value), pointer :: time
    PetscReal :: start_time, stop_time, initial_stepsize, max_stepsize
    PetscReal, allocatable :: steps(:)
    integer :: int_initial_stepsize, max_num_steps
    PetscReal, parameter :: default_start_time = 0.0_dp, &
         default_stop_time = 1.0_dp
    PetscReal, parameter :: default_initial_stepsize = 0.1_dp
    PetscReal, parameter :: default_max_stepsize = 1.0_dp
    integer :: method
    integer, parameter :: max_method_str_len = 12
    character(max_method_str_len) :: method_str
    character(max_method_str_len), parameter :: default_method_str = "beuler"
    integer, parameter :: default_max_num_steps = 100

    !! TODO: steady state simulation input

    call fson_get_mpi(json, "time.start", default_start_time, start_time)
    call fson_get_mpi(json, "time.stop", default_stop_time, stop_time)

    call fson_get_mpi(json, "time.step.method", &
         default_method_str, method_str)
    select case (str_to_lower(method_str))
    case ("beuler")
       method = TS_BEULER
    case ("bdf2")
       method = TS_BDF2
    case ("directss")
       method = TS_DIRECTSS
    case default
       method = TS_BEULER
    end select

    call fson_get_mpi(json, "time.step.initial", &
         default_initial_stepsize, initial_stepsize)
    call fson_get_mpi(json, "time.step.maximum.size", &
         default_max_stepsize, max_stepsize)
    call fson_get_mpi(json, "time.step.maximum.number", &
         default_max_num_steps, max_num_steps)

    call self%timestepper%init(method, self%mesh%dm, simulation_cell_balances, &
         simulation_cell_inflows, start_time, self%initial, initial_stepsize, &
         stop_time, max_num_steps, max_stepsize)

    if (allocated(steps)) deallocate(steps)

  end subroutine simulation_read_timestepping

!------------------------------------------------------------------------

  subroutine simulation_setup_rocktype_labels(self, json)
    !! Sets up rocktype label on the mesh. The values of the rock type
    !! label are the indices (1-based) of the rocktypes specified in the 
    !! JSON input file.

    use rock_module, only : rocktype_label_name

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscErrorCode :: ierr
    type(fson_value), pointer :: rocktypes, r
    PetscInt :: num_rocktypes, num_cells, ir, ic, c
    PetscInt, allocatable :: cells(:)
    PetscInt, allocatable :: default_cells(:)

    default_cells = [PetscInt::] ! empty integer array

    if (fson_has_mpi(json, "rock.types")) then
       call fson_get_mpi(json, "rock.types", rocktypes)
       call DMPlexCreateLabel(self%mesh%dm, rocktype_label_name, ierr); CHKERRQ(ierr)
       num_rocktypes = fson_value_count_mpi(rocktypes, ".")
       do ir = 1, num_rocktypes
          r => fson_value_get_mpi(rocktypes, ir)
          call fson_get_mpi(r, "cells", default_cells, cells)
          num_cells = size(cells)
          do ic = 1, num_cells
             c = cells(ic)
             call DMPlexSetLabelValue(self%mesh%dm, rocktype_label_name, &
                  c, ir, ierr); CHKERRQ(ierr)
          end do
       end do
    end if

  end subroutine simulation_setup_rocktype_labels

!------------------------------------------------------------------------

  subroutine simulation_setup_boundary_labels(self)
    !! Sets up labels identifying boundaries of the mesh (for e.g.
    !! applying boundary conditions.

    class(simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    PetscBool :: has_label

    call DMPlexHasLabel(self%mesh%dm, open_boundary_label_name, has_label, &
         ierr); CHKERRQ(ierr)
    if (.not.(has_label)) then
       call DMPlexCreateLabel(self%mesh%dm, open_boundary_label_name, &
            ierr); CHKERRQ(ierr)
       ! could read boundary faces from input here if needed- i.e. if labels
       ! not present in mesh file
    end if

  end subroutine simulation_setup_boundary_labels

!------------------------------------------------------------------------

  subroutine simulation_setup_labels(self, json)
    !! Sets up labels on the mesh for a simulation.

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json

    call self%setup_rocktype_labels(json)
    call self%setup_boundary_labels()

  end subroutine simulation_setup_labels

!------------------------------------------------------------------------

  subroutine simulation_init(self, filename, json_str)
    !! Initializes a simulation using data from the input file with 
    !! specified name.

    class(simulation_type), intent(in out) :: self
    character(*), intent(in), optional :: filename !! Input file name
    character(*), intent(in), optional :: json_str !! JSON string for alternative input
    ! Locals:
    type(fson_value), pointer :: json
    PetscErrorCode :: ierr

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

    call self%read_title(json)

    call self%read_thermodynamics(json)
    call self%read_eos(json)

    call self%mesh%init(json)

    call self%setup_labels(json)

    call self%mesh%configure(self%eos%primary_variable_names)

    call self%read_initial(json)

    call self%setup_fluid()

    call self%read_rock_properties(json)

    call self%read_timestepping(json)

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
