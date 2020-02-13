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
  use iso_fortran_env, only: int32
  use kinds_module
  use ode_module
  use mesh_module
  use thermodynamics_module
  use eos_module
  use relative_permeability_module
  use capillary_pressure_module
  use logfile_module
  use list_module

  implicit none

  private

  PetscInt, parameter, public :: max_title_length = 120
  PetscInt, parameter :: max_output_filename_length = 200

  type, public, extends(ode_type) :: flow_simulation_type
     !! Type for simulation of fluid mass and energy flows in porous media.
     private
     PetscInt :: solution_range_start, rock_range_start
     PetscInt :: fluid_range_start, update_cell_range_start, source_range_start
     character(:), allocatable, public :: filename !! JSON input filename
     character(max_title_length), public :: title !! Descriptive title for the simulation
     Vec, public :: rock !! Rock properties in each cell
     Vec, public :: fluid !! Fluid properties in each cell, for unperturbed primary variables
     Vec, public :: current_fluid !! Fluid properties in each cell for current primary variables
     Vec, public :: last_timestep_fluid !! Fluid properties at previous timestep
     Vec, public :: last_iteration_fluid !! Fluid properties at previous nonlinear solver iteration
     Vec, public :: balances !! Mass and energy balances for unperturbed primary variables
     Vec, public :: update_cell !! Which cells have primary variables being updated
     Vec, public :: flux !! Mass or energy fluxes through cell faces for each component
     Vec, public :: source !! Source/sink terms
     PetscInt, public :: num_local_sources !! Number of source/sink terms on current process
     PetscInt, public :: num_sources !! Total number of source/sink terms on all processes
     IS, public :: source_index !! Index set defining natural to global source ordering
     type(list_type), public :: source_controls !! Source/sink controls
     class(thermodynamics_type), allocatable, public :: thermo !! Fluid thermodynamic formulation
     class(eos_type), allocatable, public :: eos !! Fluid equation of state
     PetscReal, public :: gravity(3) !! Acceleration of gravity vector (\(m.s^{-1}\))
     class(relative_permeability_type), allocatable, public :: relative_permeability !! Rock relative permeability function
     class(capillary_pressure_type), allocatable, public :: capillary_pressure !! Rock capillary pressure function
     character(max_output_filename_length), public :: output_filename !! HDF5 output filename
     PetscViewer :: hdf5_viewer !! Viewer for HDF5 output
     PetscInt, allocatable :: output_fluid_field_indices(:) !! Field indices for fluid output
     PetscInt, allocatable :: output_source_field_indices(:) !! Field indices for source output
     integer(int32) :: start_clock !! Start wall clock time of simulation
     PetscBool :: unperturbed !! Whether any primary variables are being perturbed for Jacobian calculation
   contains
     private
     procedure :: create_solution_vector => flow_simulation_create_solution_vector
     procedure :: setup_logfile => flow_simulation_setup_logfile
     procedure :: setup_output => flow_simulation_setup_output
     procedure :: setup_output_fields => flow_simulation_setup_output_fields
     procedure :: setup_update_cell => flow_simulation_setup_update_cell
     procedure :: setup_flux_vector => flow_simulation_setup_flux_vector
     procedure :: identify_update_cells => flow_simulation_identify_update_cells
     procedure :: redistribute => flow_simulation_redistribute
     procedure :: add_boundary_ghost_cells => flow_simulation_add_boundary_ghost_cells
     procedure :: destroy_output => flow_simulation_destroy_output
     procedure, public :: setup_gravity => flow_simulation_setup_gravity
     procedure, public :: input_summary => flow_simulation_input_summary
     procedure, public :: log_statistics => flow_simulation_log_statistics
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
     procedure, public :: output_source_indices => flow_simulation_output_source_indices
     procedure, public :: output_source_cell_indices => flow_simulation_output_source_cell_indices
     procedure, public :: output => flow_simulation_output
     procedure, public :: boundary_residuals => flow_simulation_boundary_residuals
     procedure, public :: get_dof => flow_simulation_get_dof
  end type flow_simulation_type

contains

!------------------------------------------------------------------------

  subroutine flow_simulation_create_solution_vector(self, solution, range_start)
    !! Creates and returns solution vector, and corresponding range_start.

    use dm_utils_module, only: global_vec_range_start

    class(flow_simulation_type), intent(in out) :: self
    Vec, intent(out) :: solution !! Solution vector
    PetscInt, intent(out) :: range_start !! Range start, used for computing offsets
    ! Locals:
    PetscErrorCode :: ierr

    call DMCreateGlobalVector(self%mesh%dm, solution, ierr)
    CHKERRQ(ierr)
    call PetscObjectSetName(solution, "primary", ierr); CHKERRQ(ierr)
    call global_vec_range_start(solution, range_start)

  end subroutine flow_simulation_create_solution_vector

!------------------------------------------------------------------------

  subroutine flow_simulation_setup_flux_vector(self)
    !! Sets up flux vector, for storing mass and energy fluxes through
    !! cell faces.

    use dm_utils_module, only: dm_set_data_layout, global_vec_range_start

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    DM :: dm_flux
    PetscInt :: num_variables
    PetscInt, allocatable :: flux_variable_num_components(:), &
         flux_variable_dim(:)
    character(max_primary_variable_name_length), allocatable :: &
         flux_variable_names(:)
    PetscErrorCode :: ierr

    num_variables = self%eos%num_primary_variables
    allocate(flux_variable_num_components(num_variables), &
         flux_variable_dim(num_variables), flux_variable_names(num_variables))
    flux_variable_num_components = 1

    call DMClone(self%mesh%dm, dm_flux, ierr); CHKERRQ(ierr)
    flux_variable_dim = self%mesh%dim - 1

    flux_variable_names(1: self%eos%num_components) = self%eos%component_names
    if (.not. (self%eos%isothermal)) then
       flux_variable_names(self%eos%num_primary_variables) = energy_component_name
    end if

    call dm_set_data_layout(dm_flux, flux_variable_num_components, &
         flux_variable_dim, flux_variable_names)

    call DMCreateLocalVector(dm_flux, self%flux, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(self%flux, "flux", ierr); CHKERRQ(ierr)
    call VecSet(self%flux, 0._dp, ierr); CHKERRQ(ierr)

    call DMDestroy(dm_flux, ierr); CHKERRQ(ierr)
    deallocate(flux_variable_dim, flux_variable_num_components, &
         flux_variable_names)

  end subroutine flow_simulation_setup_flux_vector

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
    else
       self%hdf5_viewer = PETSC_NULL_VIEWER
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
    if (allocated(self%output_fluid_field_indices)) then
       deallocate(self%output_fluid_field_indices)
    end if
    if (allocated(self%output_source_field_indices)) then
       deallocate(self%output_source_field_indices)
    end if

  end subroutine flow_simulation_destroy_output

!------------------------------------------------------------------------

  subroutine flow_simulation_setup_output_fields(self, json)
    !! Sets up output field indices for writing to HDF5 file.

    use fson
    use source_module, only: default_output_source_fields, &
         required_output_source_fields

    class(flow_simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json

    call setup_vector_output_fields("fluid", self%fluid, &
         self%eos%default_output_fluid_fields, &
         self%eos%required_output_fluid_fields, &
         self%output_fluid_field_indices)

    call setup_vector_output_fields("source", self%source, &
         default_output_source_fields, required_output_source_fields, &
         self%output_source_field_indices)

  contains

    subroutine setup_vector_output_fields(name, v, default_fields, &
         required_fields, field_indices)

      use fson_mpi_module
      use hdf5io_module, only: max_field_name_length
      use dm_utils_module, only: section_get_field_names
      use utils_module, only: str_to_lower, str_array_index
      use fson_value_m, only: TYPE_ARRAY, TYPE_STRING

      character(*), intent(in) :: name
      Vec, intent(in) :: v
      character(max_field_name_length), intent(in) :: &
           default_fields(:), required_fields(:)
      PetscInt, allocatable, intent(out) :: field_indices(:)
      ! Locals:
      character(max_field_name_length), allocatable :: &
           output_fields(:), fields(:), lower_required_fields(:)
      DM :: dm
      PetscSection :: section
      PetscInt :: i
      PetscBool, allocatable :: required_missing(:)
      character(2) :: field_str
      PetscInt :: fields_type
      character(3) :: fields_str
      PetscErrorCode :: ierr

      call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
      call DMGetSection(dm, section, ierr); CHKERRQ(ierr)
      call section_get_field_names(section, PETSC_TRUE, fields)

      if (fson_has_mpi(json, "output.fields." // trim(name))) then
         fields_type = fson_type_mpi(json, "output.fields." // trim(name))
         select case (fields_type)
         case (TYPE_ARRAY)
            call fson_get_mpi(json, "output.fields." // trim(name), &
                 default_fields, max_field_name_length, &
                 output_fields, self%logfile)
         case (TYPE_STRING)
            call fson_get_mpi(json, "output.fields." // trim(name), &
                 val = fields_str)
            if (str_to_lower(fields_str) == "all") then
               output_fields = fields
            else
               call self%logfile%write(LOG_LEVEL_WARN, 'input', 'unrecognised', &
                    str_key = "output.fields." // trim(name), &
                    str_value = fields_str)
               output_fields = default_fields
            end if
         end select
      else ! default:
         call fson_get_mpi(json, "output.fields." // trim(name), &
              default_fields, max_field_name_length, &
              output_fields, self%logfile)
      end if
      output_fields = str_to_lower(output_fields)

      ! Check if required fields are present:
      associate(num_required => size(required_fields))
        allocate(required_missing(num_required), &
             lower_required_fields(num_required))
        lower_required_fields = str_to_lower(required_fields)
        do i = 1, num_required
           required_missing(i) = (str_array_index( &
                lower_required_fields(i), output_fields) == -1)
        end do
      end associate
      output_fields = [output_fields, &
           pack(lower_required_fields, required_missing)]

      associate(num_fields => size(output_fields))

        allocate(field_indices(num_fields))
        do i = 1, num_fields
           field_indices(i) = str_array_index( &
                output_fields(i), fields)
           if (field_indices(i) == -1) then
              write(field_str, '(i2)') i - 1
              call self%logfile%write(LOG_LEVEL_WARN, 'input', 'unrecognised', &
                   str_key = "output.fields." // trim(name) // &
                   "[" // trim(field_str) // "]", &
                   str_value = output_fields(i))
           end if
        end do
        field_indices = pack(field_indices, field_indices > -1)
        field_indices = field_indices - 1 ! zero-based indices

      end associate

      deallocate(fields, required_missing, lower_required_fields)

    end subroutine setup_vector_output_fields

  end subroutine flow_simulation_setup_output_fields

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
    if (self%mesh%has_minc) then
       call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
            str_key = 'mesh.filename', str_value = self%mesh%filename, &
            logical_keys = ['mesh.minc'], logical_values = [PETSC_TRUE])
    else
       call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
            str_key = 'mesh.filename', str_value = self%mesh%filename)
    end if
    call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'eos.name', str_value = self%eos%name)
    call self%logfile%write(LOG_LEVEL_INFO, 'input', 'summary', &
         str_key = 'thermodynamics', str_value = self%thermo%name)

  end subroutine flow_simulation_input_summary

!------------------------------------------------------------------------

  subroutine flow_simulation_log_statistics(self)
    !! Writes simulation statistics to logfile.

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    PetscInt :: dof_total, dof_min, dof_max
    PetscReal :: dof_imbalance

    call self%get_dof(dof_total, dof_min, dof_max, dof_imbalance)
    call self%logfile%write(LOG_LEVEL_INFO, 'simulation', 'dof', &
         int_keys = ['total', 'min  ', 'max  '], &
         int_values = [dof_total, dof_min, dof_max], &
         real_keys = ['imbalance'], &
         real_values = [dof_imbalance], rank = 0)
    call self%logfile%write(LOG_LEVEL_INFO, 'simulation', 'sources', &
         int_keys = ['count'], int_values = [self%num_sources])

    call self%logfile%write_blank()

  end subroutine flow_simulation_log_statistics

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

    use fson
    use fson_mpi_module
    use fson_value_m, only: TYPE_REAL, TYPE_INTEGER, TYPE_ARRAY, TYPE_NULL

    class(flow_simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscReal :: gravity_magnitude
    PetscReal, allocatable :: gravity(:)
    PetscInt :: gravity_type, ng

    self%gravity = 0._dp
    if (fson_has_mpi(json, "gravity")) then
       gravity_type = fson_type_mpi(json, "gravity")
       select case (gravity_type)
       case (TYPE_REAL, TYPE_INTEGER)
          call fson_get_mpi(json, "gravity", val = gravity_magnitude)
          self%gravity(self%mesh%dim) = -gravity_magnitude
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
      select case (self%mesh%dim)
      case(2)
         default_gravity = default_gravity_2D
      case(3)
         default_gravity = default_gravity_3D
      end select
      call fson_get_mpi(json, "gravity", default_gravity, &
           gravity_magnitude, self%logfile) ! for logging purposes
      self%gravity(self%mesh%dim) = -gravity_magnitude
    end subroutine set_default_gravity

  end subroutine flow_simulation_setup_gravity

!------------------------------------------------------------------------

  subroutine flow_simulation_setup_update_cell(self)
    !! Sets up update_cell vector, which stores update status of each
    !! cell:- -1 if cell is not being updated, or 1 if it is.  For
    !! function evaluations where the primary variables are not being
    !! perturbed to calculate finite differences for the Jacobian, all
    !! values are 1. During Jacobian calculation, all values are -1
    !! except those for cells in which variables are being perturbed,
    !! which have the value 1.

    use dm_utils_module, only: dm_set_data_layout, global_vec_range_start

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    DM :: dm_update
    PetscErrorCode :: ierr

    call DMClone(self%mesh%dm, dm_update, ierr); CHKERRQ(ierr)
    call dm_set_data_layout(dm_update, [1], [self%mesh%dim], ["update"])
    call DMCreateGlobalVector(dm_update, self%update_cell, ierr); CHKERRQ(ierr)
    call PetscObjectSetName(self%update_cell, "update_cell", ierr); CHKERRQ(ierr)
    call global_vec_range_start(self%update_cell, self%update_cell_range_start)
    call VecSet(self%update_cell, 1._dp, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_setup_update_cell

!------------------------------------------------------------------------

  subroutine flow_simulation_init(self, json, filename, err)
    !! Initializes a flow simulation using data from the specified JSON object.

    use fson
    use fson_mpi_module
    use thermodynamics_setup_module, only: setup_thermodynamics
    use eos_module, only: max_component_name_length, &
         max_phase_name_length
    use eos_setup_module, only: setup_eos
    use initial_module
    use fluid_module, only: create_fluid_vector
    use rock_module, only: setup_rock_vector
    use source_setup_module, only: setup_sources
    use utils_module, only: date_time_str
    use profiling_module, only: simulation_init_event
    use mpi_utils_module, only: mpi_broadcast_error_flag
    use dm_utils_module, only: dm_get_cell_index

    class(flow_simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    character(len = *), intent(in), optional :: filename
    PetscErrorCode, intent(out) :: err
    ! Locals:
    character(len = max_title_length), parameter :: default_title = ""
    character(25) :: datetimestr
    PetscErrorCode :: ierr, redist_err

    err = 0
    call PetscLogEventBegin(simulation_init_event, ierr); CHKERRQ(ierr)
    call system_clock(self%start_clock)
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

    call self%mesh%init(self%eos, json, self%logfile)
    call self%setup_gravity(json)

    call self%mesh%configure(self%gravity, json, self%logfile, err)
    if (err == 0) then
       call self%mesh%override_face_properties()
       call self%create_solution_vector(self%solution, self%solution_range_start)
       call setup_relative_permeabilities(json, &
            self%relative_permeability, self%logfile, err)
       if (err == 0) then
          call setup_capillary_pressures(json, &
               self%capillary_pressure, self%logfile, err)
          if (err == 0) then
             call setup_rock_vector(json, self%mesh%dm, self%rock, &
                  self%mesh%rock_types, self%rock_range_start, self%mesh%ghost_cell, &
                  self%logfile, err)
             if (err == 0) then
                if (self%mesh%has_minc) then
                   call self%mesh%setup_minc_rock_properties(json, self%rock, &
                        self%rock_range_start, self%logfile, err)
                end if
                if (err == 0) then
                   call create_fluid_vector(self%mesh%dm, max_component_name_length, &
                        self%eos%component_names, max_phase_name_length, &
                        self%eos%phase_names, self%fluid, self%fluid_range_start)

                   call setup_initial(json, self%mesh, self%eos, &
                        self%time, self%solution, self%fluid, &
                        self%solution_range_start, self%fluid_range_start, self%logfile)

                   if (self%mesh%rebalance) then
                      call self%redistribute(redist_err)
                      if (redist_err > 0) then
                         call self%logfile%write(LOG_LEVEL_WARN, 'simulation', &
                              'redistribution', logical_keys = ['fail'], &
                              logical_values = [PETSC_TRUE])
                      end if
                   end if

                   if (self%hdf5_viewer /= PETSC_NULL_VIEWER) then
                      call ISView(self%mesh%cell_index, self%hdf5_viewer, &
                           ierr); CHKERRQ(ierr)
                   end if

                   call self%add_boundary_ghost_cells()

                   call VecDuplicate(self%solution, self%balances, ierr); CHKERRQ(ierr)
                   call VecDuplicate(self%fluid, self%current_fluid, ierr); CHKERRQ(ierr)
                   call VecDuplicate(self%fluid, self%last_timestep_fluid, ierr)
                   CHKERRQ(ierr)
                   call VecDuplicate(self%fluid, self%last_iteration_fluid, ierr)
                   CHKERRQ(ierr)
                   call self%setup_flux_vector()

                   call self%setup_update_cell()
                   call self%mesh%set_boundary_conditions(json, self%solution, self%fluid, &
                        self%rock, self%eos, self%solution_range_start, &
                        self%fluid_range_start, self%rock_range_start, self%logfile)
                   call scale_initial_primary(self%mesh, self%eos, self%solution, self%fluid, &
                        self%solution_range_start, self%fluid_range_start)
                   call self%fluid_init(self%time, self%solution, err)
                   if (err == 0) then
                      call setup_sources(json, self%mesh%dm, self%mesh%cell_natural_global, &
                           self%eos, self%thermo, self%time, self%fluid, &
                           self%fluid_range_start, self%source, self%source_range_start, &
                           self%num_local_sources, self%num_sources, self%source_controls, &
                           self%source_index, self%logfile, err)
                      if (err == 0) then
                         call self%output_mesh_geometry()
                         call self%output_source_indices()
                         call self%output_source_cell_indices()
                         call self%setup_output_fields(json)
                      end if
                   end if
                end if
             end if
          end if
       end if
    end if

    call self%mesh%destroy_distribution_data()

    if (self%mesh%has_minc) then
       call DMDestroy(self%mesh%original_dm, ierr); CHKERRQ(ierr)
    end if

    self%unperturbed = PETSC_TRUE
    call mpi_broadcast_error_flag(err)
    call self%logfile%flush()

    call PetscLogEventEnd(simulation_init_event, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_init

!------------------------------------------------------------------------

  subroutine flow_simulation_destroy(self)
    !! Destroys the simulation.

    use utils_module, only : date_time_str, clock_elapsed_time

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    PetscReal :: elapsed_time

    call self%destroy_output()

    call VecDestroy(self%solution, ierr); CHKERRQ(ierr)
    call VecDestroy(self%balances, ierr); CHKERRQ(ierr)
    call VecDestroy(self%fluid, ierr); CHKERRQ(ierr)
    call VecDestroy(self%current_fluid, ierr); CHKERRQ(ierr)
    call VecDestroy(self%last_timestep_fluid, ierr); CHKERRQ(ierr)
    call VecDestroy(self%last_iteration_fluid, ierr); CHKERRQ(ierr)
    call VecDestroy(self%flux, ierr); CHKERRQ(ierr)
    call VecDestroy(self%rock, ierr); CHKERRQ(ierr)
    call VecDestroy(self%update_cell, ierr); CHKERRQ(ierr)
    call VecDestroy(self%source, ierr); CHKERRQ(ierr)
    call ISDestroy(self%source_index, ierr); CHKERRQ(ierr)
    call self%source_controls%destroy(source_control_list_node_data_destroy, &
         reverse = PETSC_TRUE)
    call self%mesh%destroy()
    call self%thermo%destroy()
    call self%eos%destroy()
    deallocate(self%thermo)
    deallocate(self%eos)
    call self%relative_permeability%destroy()
    deallocate(self%relative_permeability)
    call self%capillary_pressure%destroy()
    deallocate(self%capillary_pressure)

    elapsed_time = clock_elapsed_time(self%start_clock)
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

  subroutine flow_simulation_identify_update_cells(self, perturbed_columns)
    !! Identify which cells have primary variables that are currently
    !! being updated. For a straight function evaluation, all cells
    !! are updated. For Jacobian calculation, only cells with
    !! perturbed variables are updated.

    class(flow_simulation_type), intent(in out) :: self
    PetscInt, optional :: perturbed_columns(:)
    ! Locals:
    PetscInt :: local_num_perturbed
    PetscReal, allocatable :: update(:)
    PetscErrorCode :: ierr

    if (present(perturbed_columns)) then
       local_num_perturbed = size(perturbed_columns)
    else
       local_num_perturbed = 0
    end if
    call MPI_allreduce(local_num_perturbed == 0, self%unperturbed, 1, &
         MPI_LOGICAL, MPI_LAND, PETSC_COMM_WORLD, ierr)

    if (self%unperturbed) then ! update all
       call VecSet(self%update_cell, 1._dp, ierr); CHKERRQ(ierr)
    else
       call VecSet(self%update_cell, -1._dp, ierr); CHKERRQ(ierr)
       allocate(update(local_num_perturbed))
       update = 1._dp
       call VecSetValues(self%update_cell, local_num_perturbed, &
            perturbed_columns, update, INSERT_VALUES, ierr)
       CHKERRQ(ierr)
       call VecAssemblyBegin(self%update_cell, ierr); CHKERRQ(ierr)
       call VecAssemblyEnd(self%update_cell, ierr); CHKERRQ(ierr)
       deallocate(update)
    end if

  end subroutine flow_simulation_identify_update_cells

!------------------------------------------------------------------------

  subroutine flow_simulation_redistribute(self, err)
    !! Redistributes simulation mesh and data vectors.

    use dm_utils_module, only: dm_distribute_global_vec, &
         global_vec_range_start

    class(flow_simulation_type), intent(in out) :: self
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscSF :: sf
    PetscMPIInt :: np
    PetscErrorCode :: ierr

    err = 0
    call MPI_comm_size(PETSC_COMM_WORLD, np, ierr)

    if (np > 1) then

       call self%mesh%redistribute(sf)

       if (sf .ne. PETSC_NULL_SF) then

          call dm_distribute_global_vec(self%mesh%dm, sf, self%solution)
          call global_vec_range_start(self%solution, self%solution_range_start)

          call dm_distribute_global_vec(self%mesh%dm, sf, self%rock)
          call global_vec_range_start(self%rock, self%rock_range_start)

          call dm_distribute_global_vec(self%mesh%dm, sf, self%fluid)
          call global_vec_range_start(self%fluid, self%fluid_range_start)

          call PetscSFDestroy(sf, ierr); CHKERRQ(ierr)

       else
          err = 1
       end if

    end if

  end subroutine flow_simulation_redistribute

!------------------------------------------------------------------------

  subroutine flow_simulation_add_boundary_ghost_cells(self)
    !! Adds ghost cells for Dirichlet boundary conditions to the mesh,
    !! and adds space for these cells to already-created simulation
    !! vectors.

    use rock_module, only: create_rock_vector
    use fluid_module, only: create_fluid_vector
    use eos_module, only: max_component_name_length, &
         max_phase_name_length
    use dm_utils_module, only: vec_copy_common_local

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    Vec :: solution, rock, fluid
    PetscInt :: range_start
    PetscErrorCode :: ierr

    call self%mesh%construct_ghost_cells(self%gravity)

    call self%create_solution_vector(solution, range_start)
    call vec_copy_common_local(self%solution, solution)
    call VecDestroy(self%solution, ierr); CHKERRQ(ierr)
    self%solution = solution
    self%solution_range_start = range_start

    call create_rock_vector(self%mesh%dm, rock, range_start)
    call vec_copy_common_local(self%rock, rock)
    call VecDestroy(self%rock, ierr); CHKERRQ(ierr)
    self%rock = rock
    self%rock_range_start = range_start

    call create_fluid_vector(self%mesh%dm, max_component_name_length, &
         self%eos%component_names, max_phase_name_length, &
         self%eos%phase_names, fluid, range_start)
    call vec_copy_common_local(self%fluid, fluid)
    call VecDestroy(self%fluid, ierr); CHKERRQ(ierr)
    self%fluid = fluid
    self%fluid_range_start = range_start

  end subroutine flow_simulation_add_boundary_ghost_cells

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
    Vec, intent(in out) :: lhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: c, np, nc, start_cell, end_cell
    PetscSection :: fluid_section, rock_section, lhs_section, update_section
    PetscInt :: fluid_offset, rock_offset, lhs_offset, update_offset
    PetscReal, pointer, contiguous :: fluid_array(:), rock_array(:), &
         lhs_array(:), update(:)
    PetscReal, pointer, contiguous :: balance(:)
    type(cell_type) :: cell
    PetscErrorCode :: ierr

    call PetscLogEventBegin(cell_balances_event, ierr); CHKERRQ(ierr)

    err = 0
    np = self%eos%num_primary_variables
    nc = self%eos%num_components

    call VecCopy(self%balances, lhs, ierr); CHKERRQ(ierr)
    call global_vec_section(lhs, lhs_section)
    call VecGetArrayF90(lhs, lhs_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%current_fluid, fluid_section)
    call VecGetArrayReadF90(self%current_fluid, fluid_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%update_cell, update_section)
    call VecGetArrayReadF90(self%update_cell, update, ierr); CHKERRQ(ierr)

    call global_vec_section(self%rock, rock_section)
    call VecGetArrayReadF90(self%rock, rock_array, ierr); CHKERRQ(ierr)

    call cell%init(nc, self%eos%num_phases)
    call DMPlexGetHeightStratum(self%mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)

    do c = start_cell, end_cell - 1

       if (self%mesh%ghost_cell(c) < 0) then

          update_offset = global_section_offset(update_section, c, &
               self%update_cell_range_start)
          if (update(update_offset) > 0) then

             lhs_offset = global_section_offset(lhs_section, c, &
                  self%solution_range_start)
             balance => lhs_array(lhs_offset : lhs_offset + np - 1)

             fluid_offset = global_section_offset(fluid_section, c, &
                  self%fluid_range_start)
             rock_offset = global_section_offset(rock_section, c, &
                  self%rock_range_start)

             call cell%rock%assign(rock_array, rock_offset)
             call cell%fluid%assign(fluid_array, fluid_offset)

             balance = cell%balance(np)

          end if
       end if

    end do

    call cell%destroy()
    nullify(balance)
    call VecRestoreArrayReadF90(self%current_fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%update_cell, update, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(lhs, lhs_array, ierr); CHKERRQ(ierr)
    if (self%unperturbed) then
       call VecCopy(lhs, self%balances, ierr); CHKERRQ(ierr)
    end if

    call PetscLogEventEnd(cell_balances_event, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_cell_balances

!------------------------------------------------------------------------

  subroutine flow_simulation_cell_inflows(self, t, interval, y, rhs, err)
    !! Computes net inflow (per unit volume) into each cell, from
    !! flows through faces and source terms, for the given primary
    !! thermodynamic variables and time.

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
    Vec, intent(in out) :: rhs
    PetscErrorCode, intent(out) :: err
    ! Locals:
    PetscInt :: f, i, np
    PetscInt :: start_cell, end_cell, end_interior_cell
    PetscInt :: start_face, end_face
    Vec :: local_fluid, local_rock, local_update
    PetscReal, pointer, contiguous :: rhs_array(:)
    PetscReal, pointer, contiguous :: cell_geom_array(:), face_geom_array(:)
    PetscReal, pointer, contiguous :: fluid_array(:), rock_array(:)
    PetscReal, pointer, contiguous :: update(:), flux_array(:)
    PetscReal, pointer, contiguous :: source_data(:)
    PetscSection :: rhs_section, rock_section, fluid_section, update_section
    PetscSection :: cell_geom_section, face_geom_section, flux_section
    PetscSection :: source_section
    type(face_type) :: face
    PetscInt :: face_geom_offset, cell_geom_offsets(2), update_offsets(2)
    PetscInt :: rock_offsets(2), fluid_offsets(2), rhs_offsets(2)
    PetscInt :: cell_geom_offset, rhs_offset, flux_offset
    PetscInt, pointer :: cells(:)
    PetscReal, pointer, contiguous :: inflow(:)
    PetscReal, allocatable :: face_flux(:), face_flow(:)
    PetscReal, parameter :: flux_sign(2) = [-1._dp, 1._dp]
    PetscErrorCode :: ierr

    call PetscLogEventBegin(cell_inflows_event, ierr); CHKERRQ(ierr)
    err = 0
    np = self%eos%num_primary_variables
    allocate(face_flux(np), face_flow(np))

    call global_vec_section(rhs, rhs_section)
    call VecGetArrayF90(rhs, rhs_array, ierr); CHKERRQ(ierr)
    rhs_array = 0._dp

    call local_vec_section(self%mesh%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(self%mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)
    call local_vec_section(self%mesh%face_geom, face_geom_section)
    call VecGetArrayReadF90(self%mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)

    call global_to_local_vec_section(self%current_fluid, local_fluid, &
         fluid_section)
    call VecGetArrayReadF90(local_fluid, fluid_array, ierr); CHKERRQ(ierr)

    call global_to_local_vec_section(self%update_cell, local_update, &
         update_section)
    call VecGetArrayReadF90(local_update, update, ierr); CHKERRQ(ierr)

    call global_to_local_vec_section(self%rock, local_rock, rock_section)
    call VecGetArrayReadF90(local_rock, rock_array, ierr); CHKERRQ(ierr)

    call local_vec_section(self%flux, flux_section)
    call VecGetArrayF90(self%flux, flux_array, ierr); CHKERRQ(ierr)

    call face%init(self%eos%num_components, self%eos%num_phases)
    call DMPlexGetHeightStratum(self%mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    call DMPlexGetHeightStratum(self%mesh%dm, 1, start_face, end_face, ierr)
    CHKERRQ(ierr)
    end_interior_cell = dm_get_end_interior_cell(self%mesh%dm, end_cell)

    do f = start_face, end_face - 1

       if (self%mesh%ghost_face(f) < 0) then

          call DMPlexGetSupport(self%mesh%dm, f, cells, ierr); CHKERRQ(ierr)
          do i = 1, 2
             update_offsets(i) = section_offset(update_section, cells(i))
             cell_geom_offsets(i) = section_offset(cell_geom_section, cells(i))
             rhs_offsets(i) = global_section_offset(rhs_section, cells(i), &
                  self%solution_range_start)
          end do
          face_geom_offset = section_offset(face_geom_section, f)
          call face%assign_geometry(face_geom_array, face_geom_offset)
          call face%assign_cell_geometry(cell_geom_array, cell_geom_offsets)

          if ((update(update_offsets(1)) > 0) .or. &
               (update(update_offsets(2)) > 0)) then
             do i = 1, 2
                fluid_offsets(i) = section_offset(fluid_section, cells(i))
                rock_offsets(i) = section_offset(rock_section, cells(i))
             end do
             call face%assign_cell_fluid(fluid_array, fluid_offsets)
             call face%assign_cell_rock(rock_array, rock_offsets)
             face_flux = face%flux(self%eos)
             if (self%unperturbed) then
                flux_offset = section_offset(flux_section, f)
                flux_array(flux_offset : flux_offset + np - 1) = face_flux
             end if
          else
             flux_offset = section_offset(flux_section, f)
             face_flux = flux_array(flux_offset : flux_offset + np - 1)
          end if

          face_flow = face_flux * face%area

          do i = 1, 2
             if ((self%mesh%ghost_cell(cells(i)) < 0) .and. &
                  (cells(i) <= end_interior_cell - 1)) then
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
    call VecRestoreArrayReadF90(local_update, update, ierr); CHKERRQ(ierr)
    call restore_dm_local_vec(local_update)
    call PetscLogEventEnd(cell_inflows_event, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(self%flux, flux_array, ierr); CHKERRQ(ierr)

    ! Source / sink terms:
    call PetscLogEventBegin(sources_event, ierr); CHKERRQ(ierr)
    call global_vec_section(self%source, source_section)
    call VecGetArrayF90(self%source, source_data, ierr); CHKERRQ(ierr)
    call self%source_controls%traverse(source_control_iterator)
    call apply_sources()
    call VecRestoreArrayF90(self%source, source_data, ierr); CHKERRQ(ierr)
    call PetscLogEventEnd(sources_event, ierr); CHKERRQ(ierr)

    nullify(inflow)
    call VecRestoreArrayReadF90(local_rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(local_fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(rhs, rhs_array, ierr); CHKERRQ(ierr)
    call restore_dm_local_vec(local_fluid)
    call restore_dm_local_vec(local_rock)
    deallocate(face_flux, face_flow)

  contains

!........................................................................

    subroutine source_control_iterator(node, stopped)
      !! Applies source controls.

      type(list_node_type), pointer, intent(in out) :: node
      PetscBool, intent(out) :: stopped

      stopped = PETSC_FALSE
      select type (source_control => node%data)
      class is (source_control_type)
         call source_control%update(t, interval, source_data, &
              source_section, self%source_range_start, fluid_array, &
              fluid_section, self%eos)
      end select

    end subroutine source_control_iterator

!........................................................................

    subroutine apply_sources()
      !! Assembles contributions from sources to global RHS array.

      ! Locals:
      PetscInt :: s, c
      type(cell_type) :: cell
      type(source_type) :: source
      PetscInt :: source_offset

      call cell%init(self%eos%num_components, self%eos%num_phases)
      call source%init(self%eos)

      do s = 0, self%num_local_sources - 1

         source_offset = global_section_offset(source_section, s, &
              self%source_range_start)
         call source%assign(source_data, source_offset)
         c = nint(source%local_cell_index)

         rhs_offset = global_section_offset(rhs_section, c, &
              self%solution_range_start)
         inflow => rhs_array(rhs_offset : rhs_offset + np - 1)

         cell_geom_offset = section_offset(cell_geom_section, c)
         call cell%assign_geometry(cell_geom_array, cell_geom_offset)

         call source%update_flow(fluid_array, fluid_section)
         inflow = inflow + source%flow / cell%volume

      end do

      call source%destroy()
      call cell%destroy()

    end subroutine apply_sources

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

  subroutine flow_simulation_pre_eval(self, t, y, perturbed_columns, err)
    !! Routine to be called before each function evaluation during the
    !! nonlinear solve at each time step. Here the update_cell vector is
    !! constructed and fluid properties (excluding phase composition)
    !! are updated.

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time
    Vec, intent(in) :: y !! global primary variables vector
    PetscInt, intent(in), optional :: perturbed_columns(:)
    PetscErrorCode, intent(out) :: err !! error code
    ! Locals:
    PetscErrorCode :: ierr

    err = 0
    call self%identify_update_cells(perturbed_columns)
    call self%fluid_properties(t, y, err)
    if (self%unperturbed) then
       call VecCopy(self%current_fluid, self%fluid, ierr); CHKERRQ(ierr)
    end if

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
    PetscInt :: c, np, nc, natural, minc_level
    PetscInt :: start_cell, end_cell
    PetscSection :: y_section, fluid_section, rock_section
    PetscInt :: y_offset, fluid_offset, rock_offset
    PetscReal, pointer, contiguous :: y_array(:), scaled_cell_primary(:)
    PetscReal, pointer, contiguous :: fluid_array(:), rock_array(:)
    PetscReal :: cell_primary(self%eos%num_primary_variables)
    type(cell_type) :: cell
    ISLocalToGlobalMapping :: l2g
    character(len = 6), allocatable :: cell_keys(:)
    PetscInt, allocatable :: cell_values(:)
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call PetscLogEventBegin(fluid_init_event, ierr); CHKERRQ(ierr)
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
    call DMGetLocalToGlobalMapping(self%mesh%dm, l2g, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(self%mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)

    do c = start_cell, end_cell - 1

       if (self%mesh%ghost_cell(c) < 0) then

          y_offset = global_section_offset(y_section, c, &
               self%solution_range_start)
          scaled_cell_primary => y_array(y_offset : y_offset + np - 1)

          fluid_offset = global_section_offset(fluid_section, c, self%fluid_range_start)
          rock_offset = global_section_offset(rock_section, c, self%rock_range_start)

          call cell%rock%assign(rock_array, rock_offset)
          call cell%rock%assign_relative_permeability(self%relative_permeability)
          call cell%rock%assign_capillary_pressure(self%capillary_pressure)
          call cell%fluid%assign(fluid_array, fluid_offset)
          cell_primary = self%eos%unscale(scaled_cell_primary, nint(cell%fluid%region))

          call self%eos%bulk_properties(cell_primary, cell%fluid, err)

          if (err == 0) then
             call self%eos%phase_properties(cell_primary, cell%rock, &
                  cell%fluid, err)
             if (err > 0) then
                natural = self%mesh%local_to_parent_natural(c)
                minc_level = self%mesh%local_cell_minc_level(c)
                call self%mesh%natural_cell_output_arrays( &
                     natural, minc_level, cell_keys, cell_values)
                call self%logfile%write(LOG_LEVEL_ERR, &
                     'initialize', 'fluid_phase_properties', &
                     [cell_keys, ['region']], &
                     [cell_values, [int(cell%fluid%region)]], &
                     real_array_key = 'primary', real_array_value = cell_primary, &
                     rank = rank)
                exit
             end if
          else
             natural = self%mesh%local_to_parent_natural(c)
             minc_level = self%mesh%local_cell_minc_level(c)
             call self%mesh%natural_cell_output_arrays( &
                  natural, minc_level, cell_keys, cell_values)
             call self%logfile%write(LOG_LEVEL_ERR, &
                  'initialize', 'fluid_bulk_properties', &
                  [cell_keys, ['region']], &
                  [cell_values, [int(cell%fluid%region)]], &
                  real_array_key = 'primary', real_array_value = cell_primary, &
                  rank = rank)
             exit
          end if

       end if

    end do

    call VecRestoreArrayF90(self%rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecCopy(self%fluid, self%current_fluid, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(y, y_array, ierr); CHKERRQ(ierr)
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
    PetscInt :: c, np, nc, natural, minc_level
    PetscInt :: start_cell, end_cell
    PetscSection :: y_section, fluid_section, rock_section, update_section
    PetscInt :: y_offset, fluid_offset, rock_offset, update_offset
    PetscReal, pointer, contiguous :: y_array(:), scaled_cell_primary(:)
    PetscReal, pointer, contiguous :: fluid_array(:), rock_array(:), update(:)
    PetscReal :: cell_primary(self%eos%num_primary_variables)
    type(cell_type) :: cell
    ISLocalToGlobalMapping :: l2g
    character(len = 4), allocatable :: cell_keys(:)
    PetscInt, allocatable :: cell_values(:)
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call PetscLogEventBegin(fluid_properties_event, ierr); CHKERRQ(ierr)
    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)

    err = 0
    np = self%eos%num_primary_variables
    nc = self%eos%num_components

    call global_vec_section(y, y_section)
    call VecGetArrayReadF90(y, y_array, ierr); CHKERRQ(ierr)

    call VecCopy(self%fluid, self%current_fluid, ierr); CHKERRQ(ierr)
    call global_vec_section(self%current_fluid, fluid_section)
    call VecGetArrayF90(self%current_fluid, fluid_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%update_cell, update_section)
    call VecGetArrayReadF90(self%update_cell, update, ierr); CHKERRQ(ierr)

    call global_vec_section(self%rock, rock_section)
    call VecGetArrayF90(self%rock, rock_array, ierr); CHKERRQ(ierr)

    call cell%init(nc, self%eos%num_phases)
    call DMGetLocalToGlobalMapping(self%mesh%dm, l2g, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(self%mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)

    do c = start_cell, end_cell - 1

       if ((self%mesh%ghost_cell(c) < 0)) then

          update_offset = global_section_offset(update_section, c, &
               self%update_cell_range_start)
          if (update(update_offset) > 0) then

             y_offset = global_section_offset(y_section, c, &
                  self%solution_range_start)
             scaled_cell_primary => y_array(y_offset : y_offset + np - 1)

             fluid_offset = global_section_offset(fluid_section, c, &
                  self%fluid_range_start)
             rock_offset = global_section_offset(rock_section, c, &
                  self%rock_range_start)

             call cell%rock%assign(rock_array, rock_offset)
             call cell%rock%assign_relative_permeability(self%relative_permeability)
             call cell%rock%assign_capillary_pressure(self%capillary_pressure)
             call cell%fluid%assign(fluid_array, fluid_offset)
             cell_primary = self%eos%unscale(scaled_cell_primary, nint(cell%fluid%region))

             call self%eos%bulk_properties(cell_primary, cell%fluid, err)

             if (err == 0) then
                call self%eos%phase_properties(cell_primary, cell%rock, &
                     cell%fluid, err)
                if (err > 0) then
                   natural = self%mesh%local_to_parent_natural(c)
                   minc_level = self%mesh%local_cell_minc_level(c)
                   call self%mesh%natural_cell_output_arrays( &
                        natural, minc_level, cell_keys, cell_values)
                   call self%logfile%write(LOG_LEVEL_WARN, 'fluid', &
                        'phase_properties_not_found', &
                        cell_keys, cell_values, &
                        real_array_key = 'primary         ', &
                        real_array_value = cell_primary, rank = rank)
                   exit
                end if
             else
                natural = self%mesh%local_to_parent_natural(c)
                minc_level = self%mesh%local_cell_minc_level(c)
                call self%mesh%natural_cell_output_arrays( &
                     natural, minc_level, cell_keys, cell_values)
                call self%logfile%write(LOG_LEVEL_WARN, 'fluid', &
                     'bulk_properties_not_found', &
                     cell_keys, cell_values, &
                     real_array_key = 'primary         ', &
                     real_array_value = cell_primary, rank = rank)
                exit
             end if

          end if
       end if

    end do

    call VecRestoreArrayF90(self%rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(self%current_fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%update_cell, update, ierr); CHKERRQ(ierr)
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
    PetscInt :: c, np, nc, minc_level, natural
    PetscInt :: start_cell, end_cell
    PetscSection :: primary_section, fluid_section
    PetscInt :: primary_offset, fluid_offset
    PetscReal, pointer, contiguous :: primary_array(:), old_primary_array(:), search_array(:)
    PetscReal, pointer, contiguous :: scaled_cell_primary(:), old_scaled_cell_primary(:)
    PetscReal, pointer, contiguous :: scaled_cell_search(:)
    PetscReal, pointer, contiguous :: last_iteration_fluid_array(:), fluid_array(:)
    PetscReal, dimension(self%eos%num_primary_variables) :: cell_primary, old_cell_primary
    type(fluid_type) :: old_fluid, fluid
    PetscBool :: transition
    character(len = 16), allocatable :: cell_keys(:)
    PetscInt, allocatable :: cell_values(:)
    ISLocalToGlobalMapping :: l2g
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

    call DMGetLocalToGlobalMapping(self%mesh%dm, l2g, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(self%mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)

    do c = start_cell, end_cell - 1

       if (self%mesh%ghost_cell(c) < 0) then

          primary_offset = global_section_offset(primary_section, c, &
               self%solution_range_start)
          scaled_cell_primary => primary_array(primary_offset : primary_offset + np - 1)
          old_scaled_cell_primary => old_primary_array(primary_offset : &
               primary_offset + np - 1)
          scaled_cell_search => search_array(primary_offset : primary_offset + np - 1)

          fluid_offset = global_section_offset(fluid_section, c, &
               self%fluid_range_start)

          call old_fluid%assign(last_iteration_fluid_array, &
               fluid_offset)
          call fluid%assign(fluid_array, fluid_offset)
          cell_primary = self%eos%unscale(scaled_cell_primary, nint(fluid%region))
          old_cell_primary = self%eos%unscale(old_scaled_cell_primary, &
               nint(old_fluid%region))

          call self%eos%transition(old_cell_primary, cell_primary, &
               old_fluid, fluid, transition, err)

          if (err == 0) then
             call self%eos%check_primary_variables(fluid, cell_primary, &
                  changed_y, err)
             if (err == 0) then
                if (transition) then
                   changed_y = PETSC_TRUE
                   natural = self%mesh%local_to_parent_natural(c)
                   minc_level = self%mesh%local_cell_minc_level(c)
                   call self%mesh%natural_cell_output_arrays( &
                        natural, minc_level, cell_keys, cell_values)
                   call self%logfile%write(LOG_LEVEL_INFO, 'fluid', &
                        'transition', &
                        [cell_keys, &
                        ['old_region      ', 'new_region      ']], &
                        [cell_values, &
                        [nint(old_fluid%region), nint(fluid%region)]], &
                        real_array_key = 'new_primary     ', &
                        real_array_value = cell_primary, rank = rank)
                end if
             else
                natural = self%mesh%local_to_parent_natural(c)
                minc_level = self%mesh%local_cell_minc_level(c)
                call self%mesh%natural_cell_output_arrays( &
                     natural, minc_level, cell_keys, cell_values)
                call self%logfile%write(LOG_LEVEL_WARN, 'fluid', &
                     'out_of_range', &
                     [cell_keys, ['region          ']], &
                     [cell_values, [nint(fluid%region)]], &
                     real_array_key = 'primary         ', &
                     real_array_value = cell_primary, rank = rank)
                exit
             end if
             if (changed_y) then
                changed_search = PETSC_TRUE
                scaled_cell_primary = self%eos%scale(cell_primary, nint(fluid%region))
                scaled_cell_search = old_scaled_cell_primary - scaled_cell_primary
             end if
          else
             natural = self%mesh%local_to_parent_natural(c)
             minc_level = self%mesh%local_cell_minc_level(c)
             call self%mesh%natural_cell_output_arrays( &
                  natural, minc_level, cell_keys, cell_values)
             call self%logfile%write(LOG_LEVEL_WARN, 'fluid', &
                  'transition_failed', &
                  [cell_keys, ['region          ']], &
                  [cell_values, [nint(fluid%region)]], &
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

    use hdf5io_module, only: vec_view_fields_hdf5

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr
    DM :: geom_dm
    Vec :: global_cell_geom
    PetscInt, parameter :: cell_geom_indices(2) = [0, 1]

    if (self%output_filename /= "") then

       call VecGetDM(self%mesh%cell_geom, geom_dm, ierr); CHKERRQ(ierr)
       call DMGetGlobalVector(geom_dm, global_cell_geom, ierr); CHKERRQ(ierr)
       call PetscObjectSetName(global_cell_geom, "cell_geometry", ierr)
       CHKERRQ(ierr)

       call DMLocalToGlobal(geom_dm, self%mesh%cell_geom, &
            INSERT_VALUES, global_cell_geom, ierr); CHKERRQ(ierr)

       call vec_view_fields_hdf5(global_cell_geom, cell_geom_indices, &
            "/cell_fields", self%hdf5_viewer)

       call DMRestoreGlobalVector(geom_dm, global_cell_geom, ierr)
       CHKERRQ(ierr)

    end if

  end subroutine flow_simulation_output_mesh_geometry

!------------------------------------------------------------------------

  subroutine flow_simulation_output_source_indices(self)
    !! Writes source indices to output.

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    if ((self%output_filename /= "") .and. (self%num_sources > 0)) then
       call ISView(self%source_index, self%hdf5_viewer, ierr); CHKERRQ(ierr)
    end if

  end subroutine flow_simulation_output_source_indices

!------------------------------------------------------------------------

  subroutine flow_simulation_output_source_cell_indices(self)
    !! Writes source cell natural indices to output.

    use source_module, only: source_type
    use dm_utils_module, only: global_vec_section, global_section_offset

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    PetscInt :: i, source_offset
    PetscSection :: source_section
    PetscReal, pointer, contiguous :: source_data(:)
    PetscInt, allocatable :: source_cell_indices(:)
    type(source_type) :: source
    IS :: is_cell_indices
    PetscErrorCode :: ierr

    if ((self%output_filename /= "") .and. (self%num_sources > 0)) then

       call global_vec_section(self%source, source_section)
       call VecGetArrayReadF90(self%source, source_data, ierr); CHKERRQ(ierr)
       call source%init(self%eos)
       allocate(source_cell_indices(self%num_local_sources))
       do i = 1, self%num_local_sources
          source_offset = global_section_offset(source_section, i - 1, &
               self%source_range_start)
          call source%assign(source_data, source_offset)
          source_cell_indices(i) = nint(source%natural_cell_index)
       end do
       call VecRestoreArrayReadF90(self%source, source_data, ierr); CHKERRQ(ierr)
       call source%destroy()
       call ISCreateGeneral(PETSC_COMM_WORLD, self%num_local_sources, &
            source_cell_indices, PETSC_COPY_VALUES, is_cell_indices, ierr)
       CHKERRQ(ierr)
       deallocate(source_cell_indices)
       call PetscObjectSetName(is_cell_indices, "source_natural_cell_index", &
            ierr); CHKERRQ(ierr)
       call PetscViewerHDF5PushGroup(self%hdf5_viewer, "/source_fields", &
            ierr); CHKERRQ(ierr)
       call ISView(is_cell_indices, self%hdf5_viewer, ierr); CHKERRQ(ierr)
       call PetscViewerHDF5PopGroup(self%hdf5_viewer, ierr); CHKERRQ(ierr)
       call ISDestroy(is_cell_indices, ierr); CHKERRQ(ierr)

    end if

  end subroutine flow_simulation_output_source_cell_indices

!------------------------------------------------------------------------

  subroutine flow_simulation_output(self, time_index, time)
    !! Checkpoint output from flow simulation at a particular time.

    use profiling_module, only: output_event
    use hdf5io_module, only: vec_sequence_view_hdf5

    class(flow_simulation_type), intent(in out) :: self
    PetscInt, intent(in) :: time_index
    PetscReal, intent(in) :: time
    ! Locals:
    PetscErrorCode :: ierr

    call PetscLogEventBegin(output_event, ierr); CHKERRQ(ierr)

    if (self%output_filename /= "") then
       call vec_sequence_view_hdf5(self%fluid, &
            self%output_fluid_field_indices, "/cell_fields", time_index, &
            time, self%hdf5_viewer)
       if (self%num_sources > 0) then
          call vec_sequence_view_hdf5(self%source, &
               self%output_source_field_indices, "/source_fields", time_index, &
               time, self%hdf5_viewer)
       end if
    end if

    call PetscLogEventEnd(output_event, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_output

!------------------------------------------------------------------------

  subroutine flow_simulation_boundary_residuals(self, y, lhs, residual, err)
    !! Computes residual terms for boundary ghost cells.

    use dm_utils_module, only: global_section_offset, global_vec_section, &
         dm_get_end_interior_cell
    use cell_module, only: cell_type

    class(flow_simulation_type), intent(in out) :: self
    Vec, intent(in) :: y !! primary variables
    Vec, intent(in) :: lhs !! initial LHS vector
    Vec, intent(in out) :: residual !! residual vector
    PetscErrorCode, intent(out) :: err !! error code
    ! Locals:
    PetscInt :: c, np, nc
    PetscInt :: start_cell, end_cell, end_interior_cell
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

    call global_vec_section(self%current_fluid, fluid_section)
    call VecGetArrayReadF90(self%current_fluid, fluid_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%rock, rock_section)
    call VecGetArrayReadF90(self%rock, rock_array, ierr); CHKERRQ(ierr)

    call cell%init(nc, self%eos%num_phases)
    call DMPlexGetHeightStratum(self%mesh%dm, 0, start_cell, end_cell, ierr)
    CHKERRQ(ierr)
    end_interior_cell = dm_get_end_interior_cell(self%mesh%dm, end_cell)

    do c = end_interior_cell, end_cell - 1

       if (self%mesh%ghost_cell(c) < 0) then

          lhs_offset = global_section_offset(lhs_section, c, &
               self%solution_range_start)
          cell_lhs => lhs_array(lhs_offset : lhs_offset + np - 1)
          cell_residual => residual_array(lhs_offset : lhs_offset + np - 1)

          fluid_offset = global_section_offset(fluid_section, c, &
               self%fluid_range_start)
          rock_offset = global_section_offset(rock_section, c, &
               self%rock_range_start)

          call cell%rock%assign(rock_array, rock_offset)
          call cell%fluid%assign(fluid_array, fluid_offset)

          cell_residual = cell%balance(np) - cell_lhs

       end if

    end do

    call cell%destroy()
    call VecRestoreArrayReadF90(self%current_fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(lhs, lhs_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(residual, residual_array, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_boundary_residuals

!------------------------------------------------------------------------

  subroutine flow_simulation_get_dof(self, dof_total, dof_min, dof_max, &
       dof_imbalance)
    !! Returns total simulation degrees of freedom, and minimum and
    !! maximum over processes. Also returns the imbalance measure
    !! defined by Kumar et al. (1994).

    use dm_utils_module, only: dm_cell_counts

    class(flow_simulation_type), intent(in) :: self
    PetscInt, intent(out) :: dof_total, dof_min, dof_max
    PetscReal, intent(out) :: dof_imbalance
    ! Locals:
    PetscInt :: cells_total, cells_min, cells_max

    call dm_cell_counts(self%mesh%dm, cells_total, cells_min, cells_max)
    dof_total = cells_total * self%eos%num_primary_variables
    dof_min = cells_min * self%eos%num_primary_variables
    dof_max = cells_max * self%eos%num_primary_variables

    dof_imbalance = dble(dof_max - dof_min) / dble(dof_min)

  end subroutine flow_simulation_get_dof

!------------------------------------------------------------------------

end module flow_simulation_module
