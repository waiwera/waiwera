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

module flow_simulation_module
  !! Module for high-level representation of a flow simulation ODE.

#include <petsc/finclude/petsc.h>

  use petsc
  use ode_module
  use mesh_module
  use thermodynamics_module
  use eos_module
  use relative_permeability_module
  use logfile_module
  use list_module

  implicit none

  private

  PetscInt, parameter, public :: max_title_length = 120
  PetscInt, parameter, public :: max_flow_simulation_filename_length = 200
  PetscInt, parameter :: max_output_filename_length = 200

  type, public, extends(ode_type) :: flow_simulation_type
     !! Type for simulation of fluid mass and energy flows in porous media.
     private
     PetscInt :: solution_range_start, rock_range_start, fluid_range_start
     character(max_flow_simulation_filename_length), public :: filename !! JSON input filename
     character(max_title_length), public :: title !! Descriptive title for the simulation
     Vec, public :: rock !! Rock properties in each cell
     Vec, public :: fluid !! Fluid properties in each cell
     Vec, public :: last_timestep_fluid !! Fluid properties at previous timestep
     Vec, public :: last_iteration_fluid !! Fluid properties at previous nonlinear solver iteration
     type(list_type), public :: sources !! Source/sink terms
     type(list_type), public :: source_controls !! Source/sink controls
     class(thermodynamics_type), allocatable, public :: thermo !! Fluid thermodynamic formulation
     class(eos_type), allocatable, public :: eos !! Fluid equation of state
     PetscReal, public :: gravity(3) !! Acceleration of gravity vector (\(m.s^{-1}\))
     class(relative_permeability_type), allocatable, public :: relative_permeability !! Rock relative permeability function
     character(max_output_filename_length), public :: output_filename !! HDF5 output filename
     PetscViewer :: hdf5_viewer
     PetscLogDouble :: start_wall_time
   contains
     private
     procedure :: setup_solution_vector => flow_simulation_setup_solution_vector
     procedure :: setup_logfile => flow_simulation_setup_logfile
     procedure :: setup_output => flow_simulation_setup_output
     procedure :: destroy_output => flow_simulation_destroy_output
     procedure, public :: setup_gravity => flow_simulation_setup_gravity
     procedure, public :: input_summary => flow_simulation_input_summary
     procedure, public :: run_info => flow_simulation_run_info
     procedure, public :: init => flow_simulation_init
     procedure, public :: destroy => flow_simulation_destroy
     procedure, public :: lhs => flow_simulation_cell_balances
     procedure, public :: rhs => flow_simulation_cell_inflows
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
     procedure, public :: boundary_residuals => flow_simulation_boundary_residuals
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
         str_key = 'wall_time', str_value = datetimestr)
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
       call PetscViewerHDF5Open(PETSC_COMM_WORLD, self%output_filename, &
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
         str_key = 'title', str_value = trim(self%title))
    call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'logfile.filename', &
         str_value = self%logfile%filename)
    call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'output.filename', &
         str_value = self%output_filename)
    call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'mesh.filename', str_value = self%mesh%filename)
    call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'eos.name', str_value = self%eos%name)
    call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'thermodynamics', str_value = self%thermo%name)

  end subroutine flow_simulation_input_summary

!------------------------------------------------------------------------

  subroutine flow_simulation_run_info(self)
    !! Writes run information to logfile, e.g. software name and
    !! version, compiler details, number of processors.

    use iso_fortran_env
    use version_module

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    character(len = 6) :: major_str, minor_str, subminor_str
    PetscMPIInt :: num_procs
    PetscInt :: ierr

    call self%logfile%write(LOG_LEVEL_INFO, 'run', 'start', &
         str_key = 'software', str_value = 'Waiwera' // &
         ' version ' // trim(waiwera_version))

    write(major_str, '(i0)') PETSC_VERSION_MAJOR
    write(minor_str, '(i0)') PETSC_VERSION_MINOR
    write(subminor_str, '(i0)') PETSC_VERSION_SUBMINOR
    call self%logfile%write(LOG_LEVEL_INFO, 'run', 'start', &
         str_key = 'software', str_value = 'PETSc' // &
         ' version ' // trim(major_str) // '.' // trim(minor_str) // &
         '.' // subminor_str)

    call self%logfile%write(LOG_LEVEL_INFO, 'run', 'start', &
         str_key = 'compiler', str_value = compiler_version())
    call MPI_COMM_SIZE(PETSC_COMM_WORLD, num_procs, ierr)
    call self%logfile%write(LOG_LEVEL_INFO, 'run', 'start', &
         int_keys = ['num_processors'], int_values = [num_procs])

    call self%logfile%write_blank()

  end subroutine flow_simulation_run_info

!------------------------------------------------------------------------

  subroutine flow_simulation_setup_gravity(self, json)
    !! Sets up gravity vector from JSON input. Gravity may be
    !! specified as a scalar or array. If an array is specified, the
    !! gravity array is set to the specified value. If a scalar is
    !! specified, it is treated as the gravity magnitude and applied
    !! in the negative direction of the last dimension of the mesh. If
    !! no gravity is specified, the default value used depends on the
    !! mesh dimension. 2D meshes are effectively assumed horizontal by
    !! default, so no gravity is applied. For 3D meshes, a default
    !! gravity magnitude is applied in the third dimension.

    use kinds_module
    use fson
    use fson_mpi_module
    use fson_value_m, only: TYPE_REAL, TYPE_ARRAY, TYPE_NULL

    class(flow_simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscReal :: gravity_magnitude
    PetscReal, allocatable :: gravity(:)
    PetscInt :: gravity_type, ng, dim
    PetscErrorCode :: ierr

    call DMGetDimension(self%mesh%dm, dim, ierr); CHKERRQ(ierr)
    self%gravity = 0._dp
    if (fson_has_mpi(json, "gravity")) then
       gravity_type = fson_type_mpi(json, "gravity")
       select case (gravity_type)
       case (TYPE_REAL)
          call fson_get_mpi(json, "gravity", val = gravity_magnitude)
          self%gravity(dim) = -gravity_magnitude
       case (TYPE_ARRAY)
          call fson_get_mpi(json, "gravity", val = gravity)
          ng = size(gravity)
          self%gravity(1: ng) = gravity
          deallocate(gravity)
       case (TYPE_NULL)
          call set_default_gravity()
       case default
          call self%logfile%write(LOG_LEVEL_ERR, 'simulation', &
               'init', str_key = 'stop', &
               str_value = 'unrecognised gravity type', &
               rank = 0)
          stop
       end select
    else
       call set_default_gravity()
    end if

  contains

    subroutine set_default_gravity()
      !! Sets default gravity array.
      PetscReal :: default_gravity
      PetscReal, parameter :: default_gravity_2D = 0.0_dp
      PetscReal, parameter :: default_gravity_3D = 9.8_dp
      select case (dim)
      case(2)
         default_gravity = default_gravity_2D
      case(3)
         default_gravity = default_gravity_3D
      end select
      call fson_get_mpi(json, "gravity", default_gravity, &
           gravity_magnitude, self%logfile) ! for logging purposes
      self%gravity(dim) = -gravity_magnitude
    end subroutine set_default_gravity

  end subroutine flow_simulation_setup_gravity

!------------------------------------------------------------------------

  subroutine flow_simulation_init(self, json, filename)
    !! Initializes a flow simulation using data from the specified JSON object.

    use kinds_module
    use fson
    use fson_mpi_module
    use thermodynamics_setup_module, only: setup_thermodynamics
    use eos_module, only: max_component_name_length, &
         max_phase_name_length
    use eos_setup_module, only: setup_eos
    use initial_module, only: setup_initial
    use fluid_module, only: setup_fluid_vector
    use rock_module, only: setup_rock_vector, setup_rocktype_labels
    use source_setup_module, only: setup_sources
    use utils_module, only: date_time_str
    use profiling_module, only: simulation_init_event

    class(flow_simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    character(len = *), intent(in), optional :: filename
    ! Locals:
    character(len = max_title_length), parameter :: default_title = ""
    character(25) :: datetimestr
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

    call setup_thermodynamics(json, self%thermo, self%logfile)
    call setup_eos(json, self%thermo, self%eos, self%logfile)

    call self%mesh%init(json, self%logfile)
    call self%setup_gravity(json)
    call setup_rocktype_labels(json, self%mesh%dm, self%logfile)
    call self%mesh%setup_boundaries(json, self%eos, self%logfile)
    if (self%output_filename == '') then
       call self%mesh%configure(self%eos%primary_variable_names, self%gravity)
    else
       call self%mesh%configure(self%eos%primary_variable_names, &
            self%gravity, self%hdf5_viewer)
    end if
    call self%mesh%override_face_properties(json, self%logfile)
    call self%output_mesh_geometry()

    call self%setup_solution_vector()
    call setup_relative_permeabilities(json, &
         self%relative_permeability, self%logfile)
    call setup_rock_vector(json, self%mesh%dm, self%rock, &
         self%rock_range_start, self%mesh%ghost_cell, self%logfile)
    call setup_fluid_vector(self%mesh%dm, max_component_name_length, &
         self%eos%component_names, max_phase_name_length, &
         self%eos%phase_names, self%fluid, self%fluid_range_start)
    call VecDuplicate(self%fluid, self%last_timestep_fluid, ierr)
    CHKERRQ(ierr)
    call VecDuplicate(self%fluid, self%last_iteration_fluid, ierr)
    CHKERRQ(ierr)

    call setup_initial(json, self%mesh, self%eos, &
         self%time, self%solution, self%fluid, &
         self%solution_range_start, self%fluid_range_start, self%logfile)
    call self%mesh%set_boundary_values(self%solution, self%fluid, &
         self%rock, self%eos, self%solution_range_start, &
         self%fluid_range_start, self%rock_range_start)
    call self%fluid_init(self%time, self%solution, ierr)
    call setup_sources(json, self%mesh%dm, self%eos, self%thermo, &
         self%time, self%fluid, self%fluid_range_start, &
         self%sources, self%source_controls, self%logfile)

    call self%logfile%flush()

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

    call VecDestroy(self%solution, ierr); CHKERRQ(ierr)
    call VecDestroy(self%fluid, ierr); CHKERRQ(ierr)
    call VecDestroy(self%last_timestep_fluid, ierr); CHKERRQ(ierr)
    call VecDestroy(self%last_iteration_fluid, ierr); CHKERRQ(ierr)
    call VecDestroy(self%rock, ierr); CHKERRQ(ierr)
    call self%source_controls%destroy(source_control_list_node_data_destroy, &
         reverse = PETSC_TRUE)
    call self%sources%destroy(source_list_node_data_destroy)
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
         str_key = 'wall_time', str_value = date_time_str())

    call self%logfile%destroy()

  contains

    subroutine source_list_node_data_destroy(node)
      ! Destroys source in each list node.

      use source_module, only: source_type

      type(list_node_type), pointer, intent(in out) :: node

      select type (source => node%data)
      type is (source_type)
         call source%destroy()
      end select
    end subroutine source_list_node_data_destroy

    subroutine source_control_list_node_data_destroy(node)
      ! Destroys source control in each list node.

      use source_control_module, only: source_control_type

      type(list_node_type), pointer, intent(in out) :: node

      select type (source_control => node%data)
      class is (source_control_type)
         call source_control%destroy()
      end select
    end subroutine source_control_list_node_data_destroy

  end subroutine flow_simulation_destroy

!------------------------------------------------------------------------

  subroutine flow_simulation_cell_balances(self, t, interval, y, lhs, err)
    !! Computes mass and energy balance for each cell, for the given
    !! primary thermodynamic variables and time.

    use dm_utils_module, only: global_section_offset, global_vec_section
    use cell_module, only: cell_type
    use profiling_module, only: cell_balances_event

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time (s)
    PetscReal, intent(in) :: interval(2) !! time interval bounds
    Vec, intent(in) :: y !! global primary variables vector
    Vec, intent(out) :: lhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: c, np, nc
    PetscSection :: fluid_section, rock_section, lhs_section
    PetscInt :: fluid_offset, rock_offset, lhs_offset
    PetscReal, pointer, contiguous :: fluid_array(:), rock_array(:), lhs_array(:)
    PetscReal, pointer, contiguous :: balance(:)
    type(cell_type) :: cell
    PetscErrorCode :: ierr

    call PetscLogEventBegin(cell_balances_event, ierr); CHKERRQ(ierr)

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

    do c = self%mesh%start_cell, self%mesh%end_cell - 1

       if (self%mesh%ghost_cell(c) < 0) then

          call global_section_offset(lhs_section, c, &
               self%solution_range_start, lhs_offset, ierr); CHKERRQ(ierr)
          balance => lhs_array(lhs_offset : lhs_offset + np - 1)

          call global_section_offset(fluid_section, c, &
               self%fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)
          call global_section_offset(rock_section, c, &
               self%rock_range_start, rock_offset, ierr); CHKERRQ(ierr)

          call cell%rock%assign(rock_array, rock_offset)
          call cell%fluid%assign(fluid_array, fluid_offset)

          balance = cell%balance(np)

       end if

    end do

    call cell%destroy()
    nullify(balance)
    call VecRestoreArrayReadF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(lhs, lhs_array, ierr); CHKERRQ(ierr)

    call PetscLogEventEnd(cell_balances_event, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_cell_balances

!------------------------------------------------------------------------

  subroutine flow_simulation_cell_inflows(self, t, interval, y, rhs, err)
    !! Computes net inflow (per unit volume) into each cell, from
    !! flows through faces and source terms, for the given primary
    !! thermodynamic variables and time.

    use kinds_module
    use dm_utils_module
    use cell_module, only: cell_type
    use face_module, only: face_type
    use source_module, only: source_type
    use source_control_module, only: source_control_type
    use profiling_module, only: cell_inflows_event, sources_event

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time (s)
    PetscReal, intent(in) :: interval(2) !! time interval bounds
    Vec, intent(in) :: y !! global primary variables vector
    Vec, intent(out) :: rhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: f, i, np
    Vec :: local_fluid, local_rock
    PetscReal, pointer, contiguous :: rhs_array(:)
    PetscReal, pointer, contiguous :: cell_geom_array(:), face_geom_array(:)
    PetscReal, pointer, contiguous :: fluid_array(:), rock_array(:)
    PetscSection :: rhs_section, rock_section, fluid_section
    PetscSection :: cell_geom_section, face_geom_section
    type(cell_type) :: cell
    type(face_type) :: face
    PetscInt :: face_geom_offset, cell_geom_offsets(2)
    PetscInt :: rock_offsets(2), fluid_offsets(2), rhs_offsets(2)
    PetscInt :: cell_geom_offset, rhs_offset
    PetscInt, pointer :: cells(:)
    PetscReal, pointer, contiguous :: inflow(:)
    PetscReal, allocatable :: face_flow(:)
    PetscReal, parameter :: flux_sign(2) = [-1._dp, 1._dp]
    PetscReal, allocatable :: primary(:)
    PetscErrorCode :: ierr

    call PetscLogEventBegin(cell_inflows_event, ierr); CHKERRQ(ierr)
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

    do f = self%mesh%start_face, self%mesh%end_face - 1

       if (self%mesh%ghost_face(f) < 0) then

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

          call face%assign_geometry(face_geom_array, face_geom_offset)
          call face%assign_cell_geometry(cell_geom_array, cell_geom_offsets)
          call face%assign_cell_rock(rock_array, rock_offsets)
          call face%assign_cell_fluid(fluid_array, fluid_offsets)

          face_flow = face%flux(self%eos) * face%area

          do i = 1, 2
             if ((self%mesh%ghost_cell(cells(i)) < 0) .and. &
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
    call PetscLogEventEnd(cell_inflows_event, ierr); CHKERRQ(ierr)

    ! Source / sink terms:
    call PetscLogEventBegin(sources_event, ierr); CHKERRQ(ierr)
    call cell%init(self%eos%num_components, self%eos%num_phases)
    call self%source_controls%traverse(source_control_iterator)
    call self%sources%traverse(source_iterator)
    call PetscLogEventEnd(sources_event, ierr); CHKERRQ(ierr)

    nullify(inflow)
    call cell%destroy()
    call VecRestoreArrayReadF90(local_rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(local_fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(rhs, rhs_array, ierr); CHKERRQ(ierr)
    call restore_dm_local_vec(local_fluid)
    call restore_dm_local_vec(local_rock)
    deallocate(face_flow, primary)

  contains

    subroutine source_control_iterator(node, stopped)
      !! Applies source controls.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      select type (source_control => node%data)
      class is (source_control_type)
         call source_control%update(t, interval, fluid_array, fluid_section)
      end select
      stopped = PETSC_FALSE

    end subroutine source_control_iterator

    subroutine source_iterator(node, stopped)
      !! Assembles source contribution from source list node to global
      !! RHS array.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: c

      select type (source => node%data)
      type is (source_type)

         c = source%cell_index
         if (self%mesh%ghost_cell(c) < 0) then

            call global_section_offset(rhs_section, c, &
                 self%solution_range_start, rhs_offset, ierr)
            CHKERRQ(ierr)
            inflow => rhs_array(rhs_offset : rhs_offset + np - 1)

            call section_offset(cell_geom_section, c, &
                 cell_geom_offset, ierr); CHKERRQ(ierr)
            call cell%assign_geometry(cell_geom_array, cell_geom_offset)

            call source%update_flow(fluid_array, fluid_section)
            inflow = inflow + source%flow / cell%volume

         end if

      end select

      stopped = PETSC_FALSE

    end subroutine source_iterator

  end subroutine flow_simulation_cell_inflows

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
    !! Here we check primary variables and make any necessary region
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
    use cell_module, only: cell_type
    use profiling_module, only: fluid_init_event
    use mpi_utils_module, only: mpi_broadcast_error_flag

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time
    Vec, intent(in) :: y !! global primary variables vector
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: c, np, nc, order
    PetscSection :: y_section, fluid_section, rock_section
    PetscInt :: y_offset, fluid_offset, rock_offset
    PetscReal, pointer, contiguous :: y_array(:), cell_primary(:)
    PetscReal, pointer, contiguous :: fluid_array(:), rock_array(:)
    type(cell_type) :: cell
    DMLabel :: order_label
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call PetscLogEventBegin(fluid_init_event, ierr); CHKERRQ(ierr)
    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    err = 0
    np = self%eos%num_primary_variables
    nc = self%eos%num_components

    call global_vec_section(y, y_section)
    call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%fluid, fluid_section)
    call VecGetArrayF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%rock, rock_section)
    call VecGetArrayF90(self%rock, rock_array, ierr); CHKERRQ(ierr)

    call cell%init(nc, self%eos%num_phases)

    call DMGetLabel(self%mesh%dm, cell_order_label_name, order_label, ierr)
    CHKERRQ(ierr)

    do c = self%mesh%start_cell, self%mesh%end_cell - 1

       if (self%mesh%ghost_cell(c) < 0) then

          call global_section_offset(y_section, c, &
               self%solution_range_start, y_offset, ierr); CHKERRQ(ierr)
          cell_primary => y_array(y_offset : y_offset + np - 1)

          call global_section_offset(fluid_section, c, &
               self%fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)

          call global_section_offset(rock_section, c, &
               self%rock_range_start, rock_offset, ierr); CHKERRQ(ierr)

          call cell%rock%assign(rock_array, rock_offset)
          call cell%rock%assign_relative_permeability(self%relative_permeability)
          call cell%fluid%assign(fluid_array, fluid_offset)

          call self%eos%bulk_properties(cell_primary, cell%fluid, err)

          if (err == 0) then
             call self%eos%phase_properties(cell_primary, cell%rock, &
                  cell%fluid, err)
             if (err > 0) then
                call DMLabelGetValue(order_label, c, order, ierr)
                CHKERRQ(ierr)
                call self%logfile%write(LOG_LEVEL_ERR, 'initialize', &
                     'fluid', ['cell  ', 'region'], [order, int(cell%fluid%region)], &
                     real_array_key = 'primary', real_array_value = cell_primary, &
                     rank = rank)
                exit
             end if
          else
             call DMLabelGetValue(order_label, c, order, ierr)
             CHKERRQ(ierr)
             call self%logfile%write(LOG_LEVEL_ERR, 'initialize', &
                  'fluid', ['cell  ', 'region'], [order, int(cell%fluid%region)], &
                  real_array_key = 'primary', real_array_value = cell_primary, &
                  rank = rank)
             exit
          end if

       end if

    end do

    call VecRestoreArrayF90(self%rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call cell%destroy()

    call mpi_broadcast_error_flag(err)

    call PetscLogEventEnd(fluid_init_event, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_fluid_init

!------------------------------------------------------------------------

  subroutine flow_simulation_fluid_properties(self, t, y, err)
    !! Computes fluid properties in all cells, excluding phase
    !! composition, based on the current time and primary
    !! thermodynamic variables.

    use dm_utils_module, only: global_section_offset, global_vec_section
    use cell_module, only: cell_type
    use profiling_module, only: fluid_properties_event
    use mpi_utils_module, only: mpi_broadcast_error_flag

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time
    Vec, intent(in) :: y !! global primary variables vector
    PetscErrorCode, intent(out) :: err !! error code
    ! Locals:
    PetscInt :: c, np, nc, order
    PetscSection :: y_section, fluid_section, rock_section
    PetscInt :: y_offset, fluid_offset, rock_offset
    PetscReal, pointer, contiguous :: y_array(:), cell_primary(:)
    PetscReal, pointer, contiguous :: fluid_array(:), rock_array(:)
    type(cell_type) :: cell
    DMLabel :: order_label
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call PetscLogEventBegin(fluid_properties_event, ierr); CHKERRQ(ierr)
    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    err = 0
    np = self%eos%num_primary_variables
    nc = self%eos%num_components

    call global_vec_section(y, y_section)
    call VecGetArrayReadF90(y, y_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%fluid, fluid_section)
    call VecGetArrayF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%rock, rock_section)
    call VecGetArrayF90(self%rock, rock_array, ierr); CHKERRQ(ierr)

    call cell%init(nc, self%eos%num_phases)

    call DMGetLabel(self%mesh%dm, cell_order_label_name, order_label, ierr)
    CHKERRQ(ierr)

    do c = self%mesh%start_cell, self%mesh%end_cell - 1

       if (self%mesh%ghost_cell(c) < 0) then

          call global_section_offset(y_section, c, &
               self%solution_range_start, y_offset, ierr); CHKERRQ(ierr)
          cell_primary => y_array(y_offset : y_offset + np - 1)

          call global_section_offset(fluid_section, c, &
               self%fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)

          call global_section_offset(rock_section, c, &
               self%rock_range_start, rock_offset, ierr); CHKERRQ(ierr)

          call cell%rock%assign(rock_array, rock_offset)
          call cell%rock%assign_relative_permeability(self%relative_permeability)
          call cell%fluid%assign(fluid_array, fluid_offset)

          call self%eos%bulk_properties(cell_primary, cell%fluid, err)

          if (err == 0) then
             call self%eos%phase_properties(cell_primary, cell%rock, &
                  cell%fluid, err)
             if (err > 0) then
                call DMLabelGetValue(order_label, c, order, ierr)
                CHKERRQ(ierr)
                call self%logfile%write(LOG_LEVEL_WARN, 'fluid', &
                     'properties_not_found', &
                     ['cell            '], [order], &
                     real_array_key = 'primary         ', &
                     real_array_value = cell_primary, rank = rank)
                exit
             end if
          else
             call DMLabelGetValue(order_label, c, order, ierr); CHKERRQ(ierr)
             call self%logfile%write(LOG_LEVEL_WARN, 'fluid', &
                  'properties_not_found', &
                  ['cell            '], [order], &
                  real_array_key = 'primary         ', &
                  real_array_value = cell_primary, rank = rank)
             exit
          end if

       end if

    end do

    call VecRestoreArrayF90(self%rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(y, y_array, ierr); CHKERRQ(ierr)
    call cell%destroy()

    call mpi_broadcast_error_flag(err)

    call PetscLogEventEnd(fluid_properties_event, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_fluid_properties

!------------------------------------------------------------------------

  subroutine flow_simulation_fluid_transitions(self, y_old, search, y, &
       changed_search, changed_y, err)
    !! Checks primary variables and thermodynamic regions in all mesh
    !! cells and updates if region transitions have occurred.

    use dm_utils_module, only: global_section_offset, global_vec_section
    use fluid_module, only: fluid_type
    use profiling_module, only: fluid_transitions_event
    use mpi_utils_module, only: mpi_broadcast_error_flag, mpi_broadcast_logical

    class(flow_simulation_type), intent(in out) :: self
    Vec, intent(in) :: y_old !! Previous global primary variables vector
    Vec, intent(in out) :: y !! Global primary variables vector
    Vec, intent(in out) :: search !! Global nonlinear solver search direction vector
    PetscBool, intent(out) :: changed_search, changed_y
    PetscErrorCode, intent(out) :: err !! Error code
    ! Locals:
    PetscInt :: c, np, nc, order
    PetscSection :: primary_section, fluid_section
    PetscInt :: primary_offset, fluid_offset
    PetscReal, pointer, contiguous :: primary_array(:), old_primary_array(:), search_array(:)
    PetscReal, pointer, contiguous :: cell_primary(:), old_cell_primary(:), cell_search(:)
    PetscReal, pointer, contiguous :: last_iteration_fluid_array(:), fluid_array(:)
    type(fluid_type) :: old_fluid, fluid
    DMLabel :: order_label
    PetscBool :: transition
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call PetscLogEventBegin(fluid_transitions_event, ierr); CHKERRQ(ierr)
    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

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

    call old_fluid%init(nc, self%eos%num_phases)
    call fluid%init(nc, self%eos%num_phases)

    call DMGetLabel(self%mesh%dm, cell_order_label_name, order_label, ierr)
    CHKERRQ(ierr)

    do c = self%mesh%start_cell, self%mesh%end_cell - 1

       if (self%mesh%ghost_cell(c) < 0) then

          call global_section_offset(primary_section, c, &
               self%solution_range_start, primary_offset, ierr); CHKERRQ(ierr)
          cell_primary => primary_array(primary_offset : primary_offset + np - 1)
          old_cell_primary => old_primary_array(primary_offset : &
               primary_offset + np - 1)
          cell_search => search_array(primary_offset : primary_offset + np - 1)

          call global_section_offset(fluid_section, c, &
               self%fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)

          call old_fluid%assign(last_iteration_fluid_array, &
               fluid_offset)
          call fluid%assign(fluid_array, fluid_offset)

          call self%eos%transition(old_cell_primary, cell_primary, &
               old_fluid, fluid, transition, err)

          if (err == 0) then
             err = self%eos%check_primary_variables(fluid, cell_primary)
             if (err == 0) then
                if (transition) then
                   cell_search = old_cell_primary - cell_primary
                   changed_y = PETSC_TRUE
                   changed_search = PETSC_TRUE
                   call DMLabelGetValue(order_label, c, order, ierr)
                   CHKERRQ(ierr)
                   call self%logfile%write(LOG_LEVEL_INFO, 'fluid', &
                        'transition', &
                        ['cell            ', &
                        'old_region      ', 'new_region      '], &
                        [order, &
                        nint(old_fluid%region), nint(fluid%region)], &
                        real_array_key = 'new_primary     ', &
                        real_array_value = cell_primary, rank = rank)
                end if
             else
                call DMLabelGetValue(order_label, c, order, ierr)
                CHKERRQ(ierr)
                call self%logfile%write(LOG_LEVEL_WARN, 'fluid', &
                     'out_of_range', &
                     ['cell            ', 'region          '], &
                     [order, nint(fluid%region)], &
                     real_array_key = 'primary         ', &
                     real_array_value = cell_primary, rank = rank)
                exit
             end if

          else
             call DMLabelGetValue(order_label, c, order, ierr); CHKERRQ(ierr)
             call self%logfile%write(LOG_LEVEL_WARN, 'fluid', &
                  'transition_failed', &
                  ['cell  ', 'region'], [order, nint(fluid%region)], &
                  real_array_key = 'primary', real_array_value = cell_primary, &
                  rank = rank)
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
    call old_fluid%destroy()
    call fluid%destroy()

    call mpi_broadcast_error_flag(err)
    call mpi_broadcast_logical(changed_y)
    call mpi_broadcast_logical(changed_search)

    call PetscLogEventEnd(fluid_transitions_event, ierr); CHKERRQ(ierr)

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
    !! Checkpoint output from flow simulation at a particular time.

    use profiling_module, only: output_event

    class(flow_simulation_type), intent(in out) :: self
    PetscInt, intent(in) :: time_index
    PetscReal, intent(in) :: time
    ! Locals:
    DM :: fluid_dm
    PetscErrorCode :: ierr

    call PetscLogEventBegin(output_event, ierr); CHKERRQ(ierr)

    if (self%output_filename /= "") then
       call VecGetDM(self%fluid, fluid_dm, ierr); CHKERRQ(ierr)
       call DMSetOutputSequenceNumber(fluid_dm, time_index, time, &
            ierr); CHKERRQ(ierr)
       call VecView(self%fluid, self%hdf5_viewer, ierr); CHKERRQ(ierr)
    end if

    call PetscLogEventEnd(output_event, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_output

!------------------------------------------------------------------------

  subroutine flow_simulation_boundary_residuals(self, y, lhs, residual, err)
    !! Computes residual terms for boundary ghost cells.

    use dm_utils_module, only: global_section_offset, global_vec_section
    use cell_module, only: cell_type

    class(flow_simulation_type), intent(in out) :: self
    Vec, intent(in) :: y !! primary variables
    Vec, intent(in) :: lhs !! initial LHS vector
    Vec, intent(in out) :: residual !! residual vector
    PetscErrorCode, intent(out) :: err !! error code
    ! Locals:
    PetscInt :: c, np, nc
    PetscSection :: fluid_section, rock_section, lhs_section
    PetscInt :: fluid_offset, rock_offset, lhs_offset
    PetscReal, pointer, contiguous :: fluid_array(:), rock_array(:)
    PetscReal, pointer, contiguous :: lhs_array(:), residual_array(:)
    PetscReal, pointer, contiguous :: cell_lhs(:), cell_residual(:)
    type(cell_type) :: cell
    PetscErrorCode :: ierr

    err = 0
    np = self%eos%num_primary_variables
    nc = self%eos%num_components

    call global_vec_section(lhs, lhs_section)
    call VecGetArrayReadF90(lhs, lhs_array, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(residual, residual_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%fluid, fluid_section)
    call VecGetArrayReadF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%rock, rock_section)
    call VecGetArrayReadF90(self%rock, rock_array, ierr); CHKERRQ(ierr)

    call cell%init(nc, self%eos%num_phases)

    do c = self%mesh%end_interior_cell, self%mesh%end_cell - 1

       if (self%mesh%ghost_cell(c) < 0) then

          call global_section_offset(lhs_section, c, &
               self%solution_range_start, lhs_offset, ierr); CHKERRQ(ierr)
          cell_lhs => lhs_array(lhs_offset : lhs_offset + np - 1)
          cell_residual => residual_array(lhs_offset : lhs_offset + np - 1)

          call global_section_offset(fluid_section, c, &
               self%fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)
          call global_section_offset(rock_section, c, &
               self%rock_range_start, rock_offset, ierr); CHKERRQ(ierr)

          call cell%rock%assign(rock_array, rock_offset)
          call cell%fluid%assign(fluid_array, fluid_offset)

          cell_residual = cell%balance(np) - cell_lhs

       end if

    end do

    call cell%destroy()
    call VecRestoreArrayReadF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(lhs, lhs_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(residual, residual_array, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_boundary_residuals

!------------------------------------------------------------------------

end module flow_simulation_module
