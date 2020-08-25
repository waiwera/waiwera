!   Copyright 2020 University of Auckland.

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

module tracer_module
  !! Module for tracers.

#include <petsc/finclude/petsc.h>

  use petsc

  implicit none
  private

  PetscInt, parameter, public :: max_tracer_name_length = 32

  type, public :: tracer_type
     !! Type for passive tracer.
     private
     character(max_tracer_name_length), public :: name !! Tracer name
     PetscInt, public :: phase_index !! Tracer phase index
     PetscReal, public :: decay_constant !! Decay constant (1/s)
     PetscReal, public :: activation !! Activation energy (J/mol)
  end type tracer_type

  public :: setup_tracers, create_tracer_vector

contains

!------------------------------------------------------------------------

  subroutine setup_tracers(json, eos, tracers, logfile, err)

    !! Sets up tracers from JSON input.

    use fson
    use fson_value_m, only : TYPE_ARRAY, TYPE_OBJECT
    use fson_mpi_module
    use eos_module, only: eos_type, max_phase_name_length
    use kinds_module
    use logfile_module

    type(fson_value), pointer, intent(in) :: json
    class(eos_type), intent(in) :: eos
    type(tracer_type), allocatable, intent(out) :: tracers(:) !! Array of tracer objects
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    PetscErrorCode, intent(out) :: err
    ! Locals:
    type(fson_value), pointer :: tracer_json, traceri_json
    PetscInt :: tracer_json_type, num_tracers, i
    character(max_tracer_name_length) :: default_name
    character(max_phase_name_length) :: phase_name
    character(32) :: tracer_str
    PetscReal, parameter :: default_decay_constant = 0._dp
    PetscReal, parameter :: default_activation_energy = 0._dp

    err = 0

    if (fson_has_mpi(json, "tracer")) then

       call fson_get_mpi(json, "tracer", tracer_json)
       tracer_json_type = fson_type_mpi(json, "tracer")
       select case (tracer_json_type)
       case (TYPE_OBJECT)
          num_tracers = 1
          traceri_json => tracer_json
       case (TYPE_ARRAY)
          num_tracers = fson_value_count_mpi(json, "tracer")
          traceri_json => fson_value_children_mpi(tracer_json)
       end select
       allocate(tracers(num_tracers))

       do i = 1, num_tracers
          write(tracer_str, '(a, i0, a)') 'tracer(', i, ')'
          write(default_name, '(a, i0)') 'tracer_', i - 1
          call fson_get_mpi(traceri_json, "name", default_name, &
               tracers(i)%name, logfile, tracer_str // ".name")
          call fson_get_mpi(traceri_json, "phase", eos%default_tracer_phase, &
               phase_name, logfile, tracer_str // ".phase")
          tracers(i)%phase_index = eos%phase_index(phase_name)
          if (tracers(i)%phase_index < 0) then
             if (present(logfile)) then
                call logfile%write(LOG_LEVEL_ERR, "input", "unrecognised_phase", &
                     str_key = "name", str_value = phase_name)
             end if
             err = 1
             exit
          end if
          call fson_get_mpi(traceri_json, "decay", default_decay_constant, &
               tracers(i)%decay_constant, logfile, tracer_str // ".decay")
          call fson_get_mpi(traceri_json, "activation", default_activation_energy, &
               tracers(i)%activation, logfile, tracer_str // ".activation")
          traceri_json => fson_value_next_mpi(traceri_json)
       end do
    else
       num_tracers = 0
       allocate(tracers(num_tracers))
    end if

  end subroutine setup_tracers

!------------------------------------------------------------------------

  subroutine create_tracer_vector(dm, tracers, tracer_vector, range_start)
    !! Create tracer solution vector.

    use dm_utils_module, only: dm_set_data_layout, global_vec_range_start

    DM, intent(in) :: dm
    type(tracer_type), intent(in) :: tracers(:)
    Vec, intent(in out) :: tracer_vector
    PetscInt, intent(out) :: range_start
    ! Locals:
    DM :: dm_tracer
    PetscInt :: dim
    PetscInt, allocatable :: tracer_num_components(:), tracer_dim(:)
    character(max_tracer_name_length), allocatable :: tracer_names(:)
    ISLocalToGlobalMapping :: l2g
    PetscErrorCode :: ierr

    call DMClone(dm, dm_tracer, ierr); CHKERRQ(ierr)
    call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)

    associate(num_tracers => size(tracers))

      allocate(tracer_num_components(num_tracers))
      tracer_num_components = 1
      allocate(tracer_dim(num_tracers))
      tracer_dim = dim
      allocate(tracer_names(num_tracers))
      tracer_names = tracers%name

      call dm_set_data_layout(dm_tracer, tracer_num_components, &
           tracer_dim, tracer_names)
      call DMCreateGlobalVector(dm_tracer, tracer_vector, ierr); CHKERRQ(ierr)
      call PetscObjectSetName(tracer_vector, "tracer", ierr); CHKERRQ(ierr)
      call global_vec_range_start(tracer_vector, range_start)

      call DMGetLocalToGlobalMapping(dm_tracer, l2g, ierr); CHKERRQ(ierr)
      call VecSetLocalToGlobalMapping(tracer_vector, l2g, ierr); CHKERRQ(ierr)

      call DMDestroy(dm_tracer, ierr); CHKERRQ(ierr)
      deallocate(tracer_num_components, tracer_dim, tracer_names)

    end associate

  end subroutine create_tracer_vector

!------------------------------------------------------------------------

end module tracer_module
