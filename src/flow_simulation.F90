module flow_simulation_module
  !! Module for high-level representation of a flow simulation ODE.

  use ode_module
  use mesh_module
  use thermodynamics_module
  use eos_module
  use relative_permeability_module
  use logfile_module
  use mpi_module

  implicit none

  private

#include <petsc/finclude/petsc.h90>

  PetscInt, parameter, public :: max_title_length = 120
  PetscInt, parameter, public :: max_flow_simulation_filename_length = 200
  PetscInt, parameter :: max_output_filename_length = 200

  type, public, extends(ode_type) :: flow_simulation_type
     !! Simulation type.
     private
     PetscInt :: solution_range_start, rock_range_start, fluid_range_start
     character(max_flow_simulation_filename_length), public :: filename
     character(max_title_length), public :: title
     Vec, public :: rock
     Vec, public :: fluid, last_timestep_fluid, last_iteration_fluid
     Vec, public :: source
     class(thermodynamics_type), allocatable, public :: thermo
     class(eos_type), allocatable, public :: eos
     PetscReal, public :: gravity
     class(relative_permeability_type), allocatable, public :: relative_permeability
     character(max_output_filename_length), public :: output_filename
     PetscViewer :: hdf5_viewer
     PetscLogDouble :: start_wall_time
   contains
     private
     procedure :: setup_solution_vector => flow_simulation_setup_solution_vector
     procedure :: setup_logfile => flow_simulation_setup_logfile
     procedure :: setup_output => flow_simulation_setup_output
     procedure :: destroy_output => flow_simulation_destroy_output
     procedure, public :: input_summary => flow_simulation_input_summary
     procedure, public :: run_info => flow_simulation_run_info
     procedure, public :: init => flow_simulation_init
     procedure, public :: destroy => flow_simulation_destroy
     procedure, public :: lhs => flow_simulation_cell_balances
     procedure, public :: rhs => flow_simulation_cell_inflows
     procedure, public :: pre_solve => flow_simulation_pre_solve
     procedure, public :: pre_timestep => flow_simulation_pre_timestep
     procedure, public :: pre_retry_timestep => flow_simulation_pre_retry_timestep
     procedure, public :: pre_iteration => flow_simulation_pre_iteration
     procedure, public :: pre_eval => flow_simulation_pre_eval
     procedure, public :: post_linesearch => flow_simulation_post_linesearch
     procedure, public :: fluid_init => flow_simulation_fluid_init
     procedure, public :: fluid_transitions => flow_simulation_fluid_transitions
     procedure, public :: fluid_properties => flow_simulation_fluid_properties
     procedure, public :: output_mesh_geometry => flow_simulation_output_mesh_geometry
     procedure, public :: output => flow_simulation_output
  end type flow_simulation_type

contains

!------------------------------------------------------------------------

  subroutine flow_simulation_setup_solution_vector(self)
    !! Sets up solution vector.

    use dm_utils_module, only: global_vec_range_start

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call DMCreateGlobalVector(self%mesh%dm, self%solution, ierr)
    CHKERRQ(ierr)
    call PetscObjectSetName(self%solution, "primary", ierr); CHKERRQ(ierr)
    call global_vec_range_start(self%solution, self%solution_range_start)

  end subroutine flow_simulation_setup_solution_vector

!------------------------------------------------------------------------

  subroutine flow_simulation_setup_logfile(self, json, datetimestr)
    !! Sets up logfile output.

    use fson
    use fson_value_m, only : TYPE_LOGICAL
    use fson_mpi_module
    use utils_module, only: change_filename_extension

    class(flow_simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    character(len = *), intent(in) :: datetimestr
    ! Locals:
    character(max_logfile_name_length) :: logfile_name, &
         assumed_logfile_name
    character(max_logfile_name_length), parameter :: &
         default_logfile_name = "output.yaml"
    PetscInt, parameter :: default_max_num_length = 12
    PetscInt, parameter :: default_num_log_real_digits = 6
    PetscInt :: max_log_num_length, num_log_real_digits
    PetscBool, parameter :: default_logfile_echo = PETSC_TRUE
    PetscBool :: output_log, default_log, no_logfile, echo

    if (self%filename /= "") then
       assumed_logfile_name = &
            change_filename_extension(self%filename, "yaml")
    else
       assumed_logfile_name = default_logfile_name
    end if

    default_log = PETSC_FALSE
    no_logfile = PETSC_FALSE

    if (fson_has_mpi(json, "logfile")) then
       if (fson_type_mpi(json, "logfile") == TYPE_LOGICAL) then
          call fson_get_mpi(json, "logfile", val = output_log)
          if (output_log) then
             logfile_name = assumed_logfile_name
             default_log = PETSC_TRUE
          else
             logfile_name = ""
             no_logfile = PETSC_TRUE
          end if
       else
          if (fson_has_mpi(json, "logfile.filename")) then
             call fson_get_mpi(json, "logfile.filename", &
                  val = logfile_name)
             no_logfile = (logfile_name == "")
          else
             logfile_name = assumed_logfile_name
             default_log = PETSC_TRUE
          end if
       end if
    else
       logfile_name = assumed_logfile_name
       default_log = PETSC_TRUE
    end if

    call fson_get_mpi(json, "logfile.format.max_num_length", &
         default_max_num_length, max_log_num_length)
    call fson_get_mpi(json, "logfile.format.num_real_digits", &
         default_num_log_real_digits, num_log_real_digits)
    call fson_get_mpi(json, "logfile.echo", &
         default_logfile_echo, echo)

    call self%logfile%init(logfile_name, max_log_num_length, &
         num_log_real_digits, echo)

    call self%run_info()

    call self%logfile%write(LOG_LEVEL_INFO, 'simulation', 'init', &
         str_key = 'time', str_value = datetimestr)
    call self%logfile%write_blank()

    if (default_log) then
       call self%logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
            str_key = 'logfile.filename', &
            str_value = logfile_name)
    end if

    if (no_logfile) then
       call self%logfile%write(LOG_LEVEL_WARN, 'input', 'no logfile')
    end if

  end subroutine flow_simulation_setup_logfile

!------------------------------------------------------------------------

  subroutine flow_simulation_setup_output(self, json)
    !! Sets up simulation output to HDF5 file.

    use fson
    use fson_value_m, only : TYPE_LOGICAL
    use fson_mpi_module
    use utils_module, only: change_filename_extension

    class(flow_simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    character(max_output_filename_length), parameter :: &
         default_output_filename = "output.h5"
    character(max_output_filename_length) :: assumed_output_filename
    PetscErrorCode :: ierr
    PetscBool :: output, default_output, no_output

    default_output = PETSC_FALSE
    no_output = PETSC_FALSE

    if (self%filename /= "") then
       assumed_output_filename = &
            change_filename_extension(self%filename, "h5")
    else
       assumed_output_filename = default_output_filename
    end if

    if (fson_has_mpi(json, "output")) then
       if (fson_type_mpi(json, "output") == TYPE_LOGICAL) then
          call fson_get_mpi(json, "output", val = output)
          if (output) then
             self%output_filename = assumed_output_filename
             default_output = PETSC_TRUE
          else
             self%output_filename = ""
             no_output = PETSC_TRUE
          end if
       else
          if (fson_has_mpi(json, "output.filename")) then
             call fson_get_mpi(json, "output.filename", &
                  val = self%output_filename)
             no_output = (self%output_filename == "")
          else
             self%output_filename = assumed_output_filename
             default_output = PETSC_TRUE
          end if
       end if
    else
       self%output_filename = assumed_output_filename
       default_output = PETSC_TRUE
    end if

    if (self%output_filename /= "") then
       call PetscViewerHDF5Open(mpi%comm, self%output_filename, &
            FILE_MODE_WRITE, self%hdf5_viewer, ierr); CHKERRQ(ierr)
       call PetscViewerHDF5PushGroup(self%hdf5_viewer, "/", ierr)
       CHKERRQ(ierr)
    end if

    if (default_output) then
       call self%logfile%write(LOG_LEVEL_INFO, 'input', 'default', &
            str_key = "output.filename", &
            str_value = self%output_filename)
    end if

    if (no_output) then
       call self%logfile%write(LOG_LEVEL_WARN, 'input', 'no output')
    end if

  end subroutine flow_simulation_setup_output

!------------------------------------------------------------------------

  subroutine flow_simulation_destroy_output(self)
    !! Finalizes simulation output.

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    if (self%output_filename /= "") then
       call PetscViewerDestroy(self%hdf5_viewer, ierr); CHKERRQ(ierr)
    end if

  end subroutine flow_simulation_destroy_output

!------------------------------------------------------------------------

  subroutine flow_simulation_input_summary(self)
    !! Writes summary of important inputs to logfile.

    class(flow_simulation_type), intent(in out) :: self

    call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'input.filename', str_value = self%filename)
    call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'logfile.filename', &
         str_value = self%logfile%filename)
    call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'output.filename', &
         str_value = self%output_filename)
    call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'title', str_value = trim(self%title))
    call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'mesh', str_value = self%mesh%filename)
    call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'eos.name', str_value = self%eos%name)
    call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'thermodynamics', str_value = self%thermo%name)

    call self%logfile%write_blank()

  end subroutine flow_simulation_input_summary

!------------------------------------------------------------------------

subroutine flow_simulation_run_info(self)
  !! Writes run information to logfile.

  use iso_fortran_env
  use version_module

    class(flow_simulation_type), intent(in out) :: self

    call self%logfile%write(LOG_LEVEL_INFO, 'run', 'start', &
         str_key = 'software', str_value = 'Supermodel' // &
         ' version ' // trim(supermodel_version))
    call self%logfile%write(LOG_LEVEL_INFO, 'run', 'start', &
         str_key = 'compiler', str_value = compiler_version())
    call self%logfile%write(LOG_LEVEL_INFO, 'run', 'start', &
         int_keys = ['num_processors'], int_values = [mpi%size])

    call self%logfile%write_blank()

end subroutine flow_simulation_run_info

!------------------------------------------------------------------------

  subroutine flow_simulation_init(self, json, filename)
    !! Initializes a flow simulation using data from the specified JSON object.

    use kinds_module
    use fson
    use fson_mpi_module
    use thermodynamics_setup_module, only: setup_thermodynamics
    use eos_module, only: setup_eos, max_component_name_length, &
         max_phase_name_length
    use initial_module, only: setup_initial
    use fluid_module, only: setup_fluid_vector
    use rock_module, only: setup_rock_vector, setup_rocktype_labels
    use source_module, only: setup_source_vector
    use utils_module, only: date_time_str
    use profiling_module, only: flops, simulation_init_event

    class(flow_simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    character(len = *), intent(in), optional :: filename
    ! Locals:
    character(len = max_title_length), parameter :: default_title = ""
    character(25) :: datetimestr
    PetscReal, parameter :: default_gravity = 9.8_dp
    PetscErrorCode :: ierr

    call PetscLogEventBegin(simulation_init_event, ierr); CHKERRQ(ierr)

    call PetscTime(self%start_wall_time, ierr); CHKERRQ(ierr)
    datetimestr = date_time_str()

    if (present(filename)) then
       self%filename = filename
    else
       self%filename = ""
    end if

    call self%setup_logfile(json, datetimestr)
    call self%setup_output(json)

    call fson_get_mpi(json, "title", default_title, self%title, &
         self%logfile)
    call fson_get_mpi(json, "gravity", default_gravity, self%gravity, &
         self%logfile)

    call setup_thermodynamics(json, self%thermo, self%logfile)
    call setup_eos(json, self%thermo, self%eos, self%logfile)

    call self%mesh%init(json, self%logfile)
    call setup_rocktype_labels(json, self%mesh%dm, self%logfile)
    call self%mesh%setup_boundaries(json, self%eos, self%logfile)
    call self%mesh%configure(self%eos%primary_variable_names, self%hdf5_viewer)
    call self%output_mesh_geometry()

    call self%setup_solution_vector()
    call setup_relative_permeabilities(json, &
         self%relative_permeability, self%logfile)
    call setup_rock_vector(json, self%mesh%dm, self%rock, &
         self%rock_range_start, self%logfile)
    call setup_fluid_vector(self%mesh%dm, max_component_name_length, &
         self%eos%component_names, max_phase_name_length, &
         self%eos%phase_names, self%fluid, self%fluid_range_start)
    call VecDuplicate(self%fluid, self%last_timestep_fluid, ierr)
    CHKERRQ(ierr)
    call VecDuplicate(self%fluid, self%last_iteration_fluid, ierr)
    CHKERRQ(ierr)

    call setup_initial(json, self%mesh, self%eos, &
         self%time, self%solution, self%rock, self%fluid, &
         self%solution_range_start,  self%rock_range_start, &
         self%fluid_range_start, self%logfile)
    call self%mesh%set_boundary_values(self%solution, self%fluid, &
         self%rock, self%eos, self%solution_range_start, &
         self%fluid_range_start, self%rock_range_start)
    call setup_source_vector(json, self%mesh%dm, &
         self%eos%num_primary_variables, self%eos%isothermal, &
         self%source, self%solution_range_start, self%logfile)

    call self%logfile%flush()

    call PetscLogFlops(flops, ierr); CHKERRQ(ierr)
    call PetscLogEventEnd(simulation_init_event, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_init

!------------------------------------------------------------------------

  subroutine flow_simulation_destroy(self)
    !! Destroys the simulation.

    use utils_module, only : date_time_str

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    PetscLogDouble :: end_wall_time, elapsed_time

    call self%destroy_output()
    call self%logfile%destroy()

    call VecDestroy(self%solution, ierr); CHKERRQ(ierr)
    call VecDestroy(self%fluid, ierr); CHKERRQ(ierr)
    call VecDestroy(self%last_timestep_fluid, ierr); CHKERRQ(ierr)
    call VecDestroy(self%last_iteration_fluid, ierr); CHKERRQ(ierr)
    call VecDestroy(self%rock, ierr); CHKERRQ(ierr)
    call VecDestroy(self%source, ierr); CHKERRQ(ierr)
    call self%mesh%destroy()
    call self%thermo%destroy()
    call self%eos%destroy()
    deallocate(self%thermo)
    deallocate(self%eos)
    deallocate(self%relative_permeability)

    call PetscTime(end_wall_time, ierr); CHKERRQ(ierr)
    elapsed_time = end_wall_time - self%start_wall_time
    call self%logfile%write(LOG_LEVEL_INFO, 'simulation', 'destroy', &
         real_keys = ['elapsed_seconds'], real_values = [elapsed_time], &
         str_key = 'time', str_value = date_time_str())

  end subroutine flow_simulation_destroy

!------------------------------------------------------------------------

  subroutine flow_simulation_cell_balances(self, t, y, lhs, err)
    !! Computes mass and energy balance for each cell, for the given
    !! primary thermodynamic variables and time.

    use dm_utils_module, only: global_section_offset, global_vec_section
    use cell_module, only: cell_type
    use profiling_module, only: flops, lhs_fn_event

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time (s)
    Vec, intent(in) :: y !! global primary variables vector
    Vec, intent(out) :: lhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: c, ghost, np, nc
    PetscSection :: fluid_section, rock_section, lhs_section
    PetscInt :: fluid_offset, rock_offset, lhs_offset
    PetscReal, pointer :: fluid_array(:), rock_array(:), lhs_array(:)
    PetscReal, pointer :: balance(:)
    type(cell_type) :: cell
    DMLabel :: ghost_label
    PetscErrorCode :: ierr

    call PetscLogEventBegin(lhs_fn_event, ierr); CHKERRQ(ierr)

    err = 0
    np = self%eos%num_primary_variables
    nc = self%eos%num_components

    call global_vec_section(lhs, lhs_section)
    call VecGetArrayF90(lhs, lhs_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%fluid, fluid_section)
    call VecGetArrayReadF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%rock, rock_section)
    call VecGetArrayReadF90(self%rock, rock_array, ierr); CHKERRQ(ierr)

    call cell%init(nc, self%eos%num_phases)

    call DMGetLabel(self%mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    do c = self%mesh%start_cell, self%mesh%end_cell - 1

       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then

          call global_section_offset(lhs_section, c, &
               self%solution_range_start, lhs_offset, ierr); CHKERRQ(ierr)
          balance => lhs_array(lhs_offset : lhs_offset + np - 1)

          call global_section_offset(fluid_section, c, &
               self%fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)
          call global_section_offset(rock_section, c, &
               self%rock_range_start, rock_offset, ierr); CHKERRQ(ierr)

          call cell%assign( &
               rock_data = rock_array, rock_offset = rock_offset, &
               fluid_data = fluid_array, fluid_offset = fluid_offset)

          balance = cell%balance(np)

       end if

    end do

    call cell%destroy()
    nullify(balance)
    call VecRestoreArrayReadF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(lhs, lhs_array, ierr); CHKERRQ(ierr)

    call PetscLogFlops(flops, ierr); CHKERRQ(ierr)
    call PetscLogEventEnd(lhs_fn_event, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_cell_balances

!------------------------------------------------------------------------

  subroutine flow_simulation_cell_inflows(self, t, y, rhs, err)
    !! Computes net inflow (per unit volume) into each cell, from
    !! flows through faces and source terms, for the given primary
    !! thermodynamic variables and time.

    use kinds_module
    use dm_utils_module
    use cell_module, only: cell_type
    use face_module, only: face_type
    use profiling_module, only: flops, rhs_fn_event

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time
    Vec, intent(in) :: y !! global primary variables vector
    Vec, intent(out) :: rhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: f, c, ghost_cell, ghost_face, i, np, nc
    Vec :: local_fluid, local_rock
    PetscReal, pointer :: rhs_array(:)
    PetscReal, pointer :: cell_geom_array(:), face_geom_array(:)
    PetscReal, pointer :: fluid_array(:), rock_array(:)
    PetscReal, pointer :: source_array(:)
    PetscSection :: rhs_section, rock_section, fluid_section
    PetscSection :: source_section
    PetscSection :: cell_geom_section, face_geom_section
    type(cell_type) :: cell
    type(face_type) :: face
    PetscInt :: face_geom_offset, cell_geom_offsets(2)
    PetscInt :: rock_offsets(2), fluid_offsets(2), rhs_offsets(2)
    PetscInt :: cell_geom_offset, rhs_offset, source_offset
    PetscInt :: fluid_offset
    DMLabel :: ghost_label
    PetscInt, pointer :: cells(:)
    PetscReal, pointer :: inflow(:)
    PetscReal, allocatable :: face_flow(:), source(:)
    PetscReal, parameter :: flux_sign(2) = [-1._dp, 1._dp]
    PetscReal, allocatable :: primary(:)
    PetscErrorCode :: ierr

    call PetscLogEventBegin(rhs_fn_event, ierr); CHKERRQ(ierr)

    err = 0
    np = self%eos%num_primary_variables
    allocate(face_flow(np), primary(np))

    call global_vec_section(rhs, rhs_section)
    call VecGetArrayF90(rhs, rhs_array, ierr); CHKERRQ(ierr)
    rhs_array = 0._dp

    call local_vec_section(self%mesh%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(self%mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)
    call local_vec_section(self%mesh%face_geom, face_geom_section)
    call VecGetArrayReadF90(self%mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)

    call global_to_local_vec_section(self%fluid, local_fluid, fluid_section)
    call VecGetArrayReadF90(local_fluid, fluid_array, ierr); CHKERRQ(ierr)

    call global_to_local_vec_section(self%rock, local_rock, rock_section)
    call VecGetArrayReadF90(local_rock, rock_array, ierr); CHKERRQ(ierr)

    call face%init(self%eos%num_components, self%eos%num_phases)

    call DMGetLabel(self%mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    do f = self%mesh%start_face, self%mesh%end_face - 1

       call DMLabelGetValue(ghost_label, f, ghost_face, ierr); CHKERRQ(ierr)
       if (ghost_face < 0) then

          call section_offset(face_geom_section, f, face_geom_offset, &
               ierr); CHKERRQ(ierr)

          call DMPlexGetSupport(self%mesh%dm, f, cells, ierr); CHKERRQ(ierr)
          do i = 1, 2
             call section_offset(cell_geom_section, cells(i), &
                  cell_geom_offsets(i), ierr); CHKERRQ(ierr)
             call section_offset(fluid_section, cells(i), &
                  fluid_offsets(i), ierr); CHKERRQ(ierr)
             call section_offset(rock_section, cells(i), &
                  rock_offsets(i), ierr); CHKERRQ(ierr)
             call global_section_offset(rhs_section, cells(i), &
                  self%solution_range_start, rhs_offsets(i), ierr)
             CHKERRQ(ierr)
          end do

          call face%assign(face_geom_array, face_geom_offset, &
               cell_geom_array, cell_geom_offsets, &
               rock_array, rock_offsets, fluid_array, fluid_offsets)

          face_flow = face%flux(self%eos, self%gravity) * face%area

          do i = 1, 2
             call DMLabelGetValue(ghost_label, cells(i), ghost_cell, &
                  ierr); CHKERRQ(ierr)
             if ((ghost_cell < 0) .and. &
                  (cells(i) <= self%mesh%end_interior_cell - 1)) then
                inflow => rhs_array(rhs_offsets(i) : rhs_offsets(i) + np - 1)
                inflow = inflow + flux_sign(i) * face_flow / &
                     face%cell(i)%volume
             end if
          end do

       end if
    end do

    call face%destroy()
    call VecRestoreArrayReadF90(self%mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)

    ! Source/ sink terms:
    nc = self%eos%num_components
    call cell%init(nc, self%eos%num_phases)
    call VecGetArrayReadF90(self%source, source_array, ierr); CHKERRQ(ierr)
    call global_vec_section(self%source, source_section)
    allocate(source(np))

    do c = self%mesh%start_cell, self%mesh%end_cell - 1

       call DMLabelGetValue(ghost_label, c, ghost_cell, ierr); CHKERRQ(ierr)
       if (ghost_cell < 0) then
          call global_section_offset(rhs_section, c, &
               self%solution_range_start, rhs_offset, ierr)
          CHKERRQ(ierr)
          call section_offset(cell_geom_section, c, &
               cell_geom_offset, ierr); CHKERRQ(ierr)
          call section_offset(fluid_section, c, &
               fluid_offset, ierr); CHKERRQ(ierr)
          call cell%assign(cell_geom_array, cell_geom_offset, &
               fluid_data = fluid_array, fluid_offset = fluid_offset)
          call global_section_offset(source_section, c, &
               self%solution_range_start, source_offset, ierr)
          CHKERRQ(ierr)
          inflow => rhs_array(rhs_offset : rhs_offset + np - 1)
          source = source_array(source_offset : source_offset + np - 1)
          call cell%fluid%energy_production(source, self%eos%isothermal)
          inflow = inflow + source / cell%volume
       end if

    end do

    nullify(inflow)
    deallocate(source)
    call cell%destroy()
    call VecRestoreArrayReadF90(local_rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(local_fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%source, source_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(rhs, rhs_array, ierr); CHKERRQ(ierr)
    call restore_dm_local_vec(local_fluid)
    call restore_dm_local_vec(local_rock)
    deallocate(face_flow, primary)

    call PetscLogFlops(flops, ierr); CHKERRQ(ierr)
    call PetscLogEventEnd(rhs_fn_event, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_cell_inflows

!------------------------------------------------------------------------

  subroutine flow_simulation_pre_solve(self, t, y, err)
    !! Routine to be called before the timestepper starts to run.
    !! Here the initial fluid properties in all cells are computed.

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time
    Vec, intent(in) :: y !! global primary variables vector
    PetscErrorCode, intent(out) :: err

    err = 0
    call self%fluid_init(t, y, err)

  end subroutine flow_simulation_pre_solve

!------------------------------------------------------------------------

  subroutine flow_simulation_pre_timestep(self)
    !! Routine to be called before starting each time step. Here the
    !! last_timestep_fluid vector is initialized from the fluid
    !! vector, to be used to revert to the previous state when
    !! re-trying a timestep with a smaller step size.

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call VecCopy(self%fluid, self%last_timestep_fluid, ierr)
    CHKERRQ(ierr)

  end subroutine flow_simulation_pre_timestep

!------------------------------------------------------------------------

  subroutine flow_simulation_pre_retry_timestep(self)
    !! Routine to be called before re-trying a time step. Here the
    !! fluid vector is reset to the last_timestep_fluid vector.

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call VecCopy(self%last_timestep_fluid, self%fluid, ierr)
    CHKERRQ(ierr)

  end subroutine flow_simulation_pre_retry_timestep

!------------------------------------------------------------------------

  subroutine flow_simulation_pre_iteration(self, y, err)
    !! Routine to be called at the start of each nonlinear solver
    !! iteration during each time step. Here the fluid vector for the
    !! start of the current iteration is saved.

    class(flow_simulation_type), intent(in out) :: self
    Vec, intent(in out) :: y !! Global primary variables vector
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscErrorCode :: ierr

    err = 0
    call VecCopy(self%fluid, self%last_iteration_fluid, ierr)

  end subroutine flow_simulation_pre_iteration

!------------------------------------------------------------------------

  subroutine flow_simulation_pre_eval(self, t, y, err)
    !! Routine to be called before each function evaluation during the
    !! nonlinear solve at each time step. Here the fluid properties
    !! (excluding phase composition) are updated.

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time
    Vec, intent(in) :: y !! global primary variables vector
    PetscErrorCode, intent(out) :: err !! error code

    err = 0
    call self%fluid_properties(t, y, err)

  end subroutine flow_simulation_pre_eval

!------------------------------------------------------------------------

  subroutine flow_simulation_post_linesearch(self, y_old, search, y, &
       changed_search, changed_y, err)
    !! Routine to be called after each nonlinear solve line search.
    !! Here we check primary variables and make and necessary region
    !! transitions.

    class(flow_simulation_type), intent(in out) :: self
    Vec, intent(in) :: y_old
    Vec, intent(in out) :: search, y
    PetscBool, intent(out) :: changed_search, changed_y
    PetscErrorCode, intent(out) :: err

    err = 0
    call self%fluid_transitions(y_old, search, y, changed_search, &
         changed_y, err)

  end subroutine flow_simulation_post_linesearch

!------------------------------------------------------------------------

  subroutine flow_simulation_fluid_init(self, t, y, err)
    !! Computes fluid properties in all cells, including phase
    !! composition, based on the current time and primary
    !! thermodynamic variables. This is called before the timestepper
    !! starts to run.

    use dm_utils_module, only: global_section_offset, global_vec_section
    use fluid_module, only: fluid_type
    use rock_module, only: rock_type
    use profiling_module, only: flops, fluid_init_event

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time
    Vec, intent(in) :: y !! global primary variables vector
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: c, np, nc, ghost, order
    PetscSection :: y_section, fluid_section, rock_section
    PetscInt :: y_offset, fluid_offset, rock_offset
    PetscReal, pointer :: y_array(:), cell_primary(:)
    PetscReal, pointer :: fluid_array(:), rock_array(:)
    type(fluid_type) :: fluid
    type(rock_type) :: rock
    DMLabel :: ghost_label, order_label
    PetscErrorCode :: ierr

    call PetscLogEventBegin(fluid_init_event, ierr); CHKERRQ(ierr)

    err = 0
    np = self%eos%num_primary_variables
    nc = self%eos%num_components

    call global_vec_section(y, y_section)
    call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%fluid, fluid_section)
    call VecGetArrayF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%rock, rock_section)
    call VecGetArrayF90(self%rock, rock_array, ierr); CHKERRQ(ierr)

    call fluid%init(nc, self%eos%num_phases)

    call DMGetLabel(self%mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)
    call DMGetLabel(self%mesh%dm, cell_order_label_name, order_label, ierr)
    CHKERRQ(ierr)

    do c = self%mesh%start_cell, self%mesh%end_cell - 1

       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then

          call global_section_offset(y_section, c, &
               self%solution_range_start, y_offset, ierr); CHKERRQ(ierr)
          cell_primary => y_array(y_offset : y_offset + np - 1)

          call global_section_offset(fluid_section, c, &
               self%fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)

          call fluid%assign(fluid_array, fluid_offset)

          call global_section_offset(rock_section, c, &
               self%rock_range_start, rock_offset, ierr); CHKERRQ(ierr)

          call rock%assign(rock_array, rock_offset, &
               self%relative_permeability)

          call DMLabelGetValue(order_label, c, order, ierr); CHKERRQ(ierr)

          call self%eos%bulk_properties(cell_primary, fluid, err)

          if (err == 0) then
             call self%eos%phase_composition(fluid, err)
             if (err == 0) then
                call self%eos%phase_properties(cell_primary, rock, &
                     fluid, err)
                if (err > 0) then
                   call self%logfile%write(LOG_LEVEL_ERR, 'initialize', &
                        'fluid', ['cell            '], [order], &
                        rank = mpi%rank)
                   exit
                end if
             else
                call self%logfile%write(LOG_LEVEL_ERR, 'initialize', &
                     'fluid', ['cell            '], [order], &
                     rank = mpi%rank)
                exit
             end if
          else
             call self%logfile%write(LOG_LEVEL_ERR, 'initialize', &
                  'fluid', ['cell            '], [order], &
                  rank = mpi%rank)
             exit
          end if

       end if

    end do

    call VecRestoreArrayF90(self%rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call rock%destroy()
    call fluid%destroy()

    call mpi%broadcast_error_flag(err)

    call PetscLogFlops(flops, ierr); CHKERRQ(ierr)
    call PetscLogEventEnd(fluid_init_event, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_fluid_init

!------------------------------------------------------------------------

  subroutine flow_simulation_fluid_properties(self, t, y, err)
    !! Computes fluid properties in all cells, excluding phase
    !! composition, based on the current time and primary
    !! thermodynamic variables.

    use dm_utils_module, only: global_section_offset, global_vec_section
    use fluid_module, only: fluid_type
    use rock_module, only: rock_type
    use profiling_module, only: flops, fluid_properties_event

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time
    Vec, intent(in) :: y !! global primary variables vector
    PetscErrorCode, intent(out) :: err !! error code
    ! Locals:
    PetscInt :: c, np, nc, ghost, order
    PetscSection :: y_section, fluid_section, rock_section
    PetscInt :: y_offset, fluid_offset, rock_offset
    PetscReal, pointer :: y_array(:), cell_primary(:)
    PetscReal, pointer :: fluid_array(:), rock_array(:)
    type(fluid_type) :: fluid
    type(rock_type) :: rock
    DMLabel :: ghost_label, order_label
    PetscErrorCode :: ierr

    call PetscLogEventBegin(fluid_properties_event, ierr); CHKERRQ(ierr)

    err = 0
    np = self%eos%num_primary_variables
    nc = self%eos%num_components

    call global_vec_section(y, y_section)
    call VecGetArrayReadF90(y, y_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%fluid, fluid_section)
    call VecGetArrayF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%rock, rock_section)
    call VecGetArrayF90(self%rock, rock_array, ierr); CHKERRQ(ierr)

    call fluid%init(nc, self%eos%num_phases)

    call DMGetLabel(self%mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)
    call DMGetLabel(self%mesh%dm, cell_order_label_name, order_label, ierr)
    CHKERRQ(ierr)

    do c = self%mesh%start_cell, self%mesh%end_cell - 1

       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then

          call global_section_offset(y_section, c, &
               self%solution_range_start, y_offset, ierr); CHKERRQ(ierr)
          cell_primary => y_array(y_offset : y_offset + np - 1)

          call global_section_offset(fluid_section, c, &
               self%fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)

          call fluid%assign(fluid_array, fluid_offset)

          call global_section_offset(rock_section, c, &
               self%rock_range_start, rock_offset, ierr); CHKERRQ(ierr)

          call rock%assign(rock_array, rock_offset, &
               self%relative_permeability)

          call DMLabelGetValue(order_label, c, order, ierr); CHKERRQ(ierr)

          call self%eos%bulk_properties(cell_primary, fluid, err)

          if (err == 0) then
             call self%eos%phase_properties(cell_primary, rock, fluid, err)
             if (err > 0) then
                call self%logfile%write(LOG_LEVEL_WARN, 'fluid', &
                     'properties_not_found', &
                     ['cell            '], [order], &
                     real_array_key = 'primary         ', &
                     real_array_value = cell_primary, rank = mpi%rank)
                exit
             end if
          else
             call self%logfile%write(LOG_LEVEL_WARN, 'fluid', &
                  'properties_not_found', &
                  ['cell            '], [order], &
                  real_array_key = 'primary         ', &
                  real_array_value = cell_primary, rank = mpi%rank)
             exit
          end if

       end if

    end do

    call VecRestoreArrayF90(self%rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(y, y_array, ierr); CHKERRQ(ierr)
    call rock%destroy()
    call fluid%destroy()

    call mpi%broadcast_error_flag(err)

    call PetscLogFlops(flops, ierr); CHKERRQ(ierr)
    call PetscLogEventEnd(fluid_properties_event, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_fluid_properties

!------------------------------------------------------------------------

  subroutine flow_simulation_fluid_transitions(self, y_old, search, y, &
       changed_search, changed_y, err)
    !! Checks primary variables and thermodynamic regions in all mesh
    !! cells and updates if region transitions have occurred.

    use dm_utils_module, only: global_section_offset, global_vec_section
    use fluid_module, only: fluid_type

    class(flow_simulation_type), intent(in out) :: self
    Vec, intent(in) :: y_old !! Previous global primary variables vector
    Vec, intent(in out) :: y !! Global primary variables vector
    Vec, intent(in out) :: search !! Global nonlinear solver search direction vector
    PetscBool, intent(out) :: changed_search, changed_y
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscInt :: c, np, nc, ghost, order
    PetscSection :: primary_section, fluid_section
    PetscInt :: primary_offset, fluid_offset
    PetscReal, pointer :: primary_array(:), old_primary_array(:), search_array(:)
    PetscReal, pointer :: cell_primary(:), old_cell_primary(:), cell_search(:)
    PetscReal, pointer :: last_iteration_fluid_array(:), fluid_array(:)
    type(fluid_type) :: last_iteration_fluid, fluid
    DMLabel :: ghost_label, order_label
    PetscBool :: transition
    PetscErrorCode :: ierr

    err = 0
    changed_search = PETSC_FALSE
    changed_y = PETSC_FALSE
    np = self%eos%num_primary_variables
    nc = self%eos%num_components

    call global_vec_section(y, primary_section)
    call VecGetArrayF90(y, primary_array, ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(y_old, old_primary_array, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(search, search_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%fluid, fluid_section)
    call VecGetArrayF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(self%last_iteration_fluid, &
         last_iteration_fluid_array, ierr); CHKERRQ(ierr)

    call last_iteration_fluid%init(nc, self%eos%num_phases)
    call fluid%init(nc, self%eos%num_phases)

    call DMGetLabel(self%mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)
    call DMGetLabel(self%mesh%dm, cell_order_label_name, order_label, ierr)
    CHKERRQ(ierr)

    do c = self%mesh%start_cell, self%mesh%end_cell - 1

       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then

          call global_section_offset(primary_section, c, &
               self%solution_range_start, primary_offset, ierr); CHKERRQ(ierr)
          cell_primary => primary_array(primary_offset : primary_offset + np - 1)
          old_cell_primary => old_primary_array(primary_offset : &
               primary_offset + np - 1)
          cell_search => search_array(primary_offset : primary_offset + np - 1)

          call global_section_offset(fluid_section, c, &
               self%fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)

          call last_iteration_fluid%assign(last_iteration_fluid_array, &
               fluid_offset)
          call fluid%assign(fluid_array, fluid_offset)

          call DMLabelGetValue(order_label, c, order, ierr); CHKERRQ(ierr)

          call self%eos%transition(cell_primary, last_iteration_fluid, &
               fluid, transition, err)

          if (err == 0) then
             err = self%eos%check_primary_variables(fluid, cell_primary)
             if (err == 0) then
                if (transition) then
                   cell_search = old_cell_primary - cell_primary
                   changed_y = PETSC_TRUE
                   changed_search = PETSC_TRUE
                   call self%logfile%write(LOG_LEVEL_INFO, 'fluid', &
                        'transition', &
                        ['cell            ', &
                        'old_region      ', 'new_region      '], &
                        [order, &
                        nint(last_iteration_fluid%region), nint(fluid%region)], &
                        real_array_key = 'new_primary     ', &
                        real_array_value = cell_primary, rank = mpi%rank)
                end if
             else
                call self%logfile%write(LOG_LEVEL_WARN, 'fluid', &
                     'out_of_range', &
                     ['cell            ', 'region          '], &
                     [order, nint(fluid%region)], &
                     real_array_key = 'primary         ', &
                     real_array_value = cell_primary, rank = mpi%rank)
                exit
             end if

          else
             call self%logfile%write(LOG_LEVEL_WARN, 'fluid', &
                  'transition_failed', &
                  ['cell            '], [order], rank = mpi%rank)
             exit
          end if

       end if

    end do

    call VecRestoreArrayReadF90(self%last_iteration_fluid, &
         last_iteration_fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(y, primary_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(y_old, old_primary_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(search, search_array, ierr); CHKERRQ(ierr)
    call last_iteration_fluid%destroy()
    call fluid%destroy()

    call mpi%broadcast_error_flag(err)
    call mpi%broadcast_logical(changed_y)
    call mpi%broadcast_logical(changed_search)

  end subroutine flow_simulation_fluid_transitions

!------------------------------------------------------------------------

  subroutine flow_simulation_output_mesh_geometry(self)
    !! Writes mesh geometry data to output.

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    DM :: geom_dm
    Vec :: global_cell_geom

    if (self%output_filename /= "") then

       call VecGetDM(self%mesh%cell_geom, geom_dm, ierr); CHKERRQ(ierr)
       call DMGetGlobalVector(geom_dm, global_cell_geom, ierr); CHKERRQ(ierr)
       call PetscObjectSetName(global_cell_geom, "cell_geometry", ierr)
       CHKERRQ(ierr)

       call DMLocalToGlobalBegin(geom_dm, self%mesh%cell_geom, &
            INSERT_VALUES, global_cell_geom, ierr); CHKERRQ(ierr)
       call DMLocalToGlobalEnd(geom_dm, self%mesh%cell_geom, &
            INSERT_VALUES, global_cell_geom, ierr); CHKERRQ(ierr)

       call VecView(global_cell_geom, self%hdf5_viewer, ierr); CHKERRQ(ierr)

       call DMRestoreGlobalVector(geom_dm, global_cell_geom, ierr)
       CHKERRQ(ierr)

    end if

  end subroutine flow_simulation_output_mesh_geometry

!------------------------------------------------------------------------

  subroutine flow_simulation_output(self, time_index, time)
    !! Output from flow simulation.

    class(flow_simulation_type), intent(in out) :: self
    PetscInt, intent(in) :: time_index
    PetscReal, intent(in) :: time
    ! Locals:
    DM :: fluid_dm
    PetscErrorCode :: ierr

    if (self%output_filename /= "") then
       call VecGetDM(self%fluid, fluid_dm, ierr); CHKERRQ(ierr)
       call DMSetOutputSequenceNumber(fluid_dm, time_index, time, &
            ierr); CHKERRQ(ierr)
       call VecView(self%fluid, self%hdf5_viewer, ierr); CHKERRQ(ierr)
    end if

  end subroutine flow_simulation_output

!------------------------------------------------------------------------

end module flow_simulation_module
