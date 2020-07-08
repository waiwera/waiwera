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

  PetscInt, parameter :: max_tracer_name_length = 32
  PetscInt, parameter, public :: tracer_phase_index = 1 ! assume liquid phase tracers

  type, public :: tracer_type
     !! Type for passive tracer.
     private
     character(max_tracer_name_length), public :: name !! Tracer name
     PetscReal, public :: decay !! Decay rate (1/s)
  end type tracer_type

  public :: setup_tracer

contains

!------------------------------------------------------------------------

  subroutine setup_tracer(json, dm, auxiliary, aux_solution, range_start, &
       tracers)

    !! Sets up tracers from JSON input.

    use fson
    use fson_mpi_module
    use kinds_module
    use dm_utils_module, only: dm_set_data_layout, global_vec_range_start

    type(fson_value), pointer, intent(in) :: json
    DM, intent(in out) :: dm
    PetscBool, intent(out) :: auxiliary
    Vec, intent(in out) :: aux_solution
    PetscInt, intent(out) :: range_start
    type(tracer_type), allocatable, intent(out) :: tracers(:)
    ! Locals:
    DM :: dm_aux
    PetscInt :: dim, i, num_tracers
    PetscInt, allocatable :: tracer_num_components(:), tracer_dim(:)
    character(max_tracer_name_length), allocatable :: tracer_names(:)
    PetscReal, parameter :: default_decay_rate = 0._dp
    PetscErrorCode :: ierr

    if (fson_has_mpi(json, "tracer")) then

       auxiliary = PETSC_TRUE
       call DMClone(dm, dm_aux, ierr); CHKERRQ(ierr)
       call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)

       ! TODO get tracer parameters
       num_tracers = 1
       allocate(tracers(num_tracers))

       allocate(tracer_num_components(num_tracers))
       tracer_num_components = 1
       allocate(tracer_dim(num_tracers))
       tracer_dim = dim
       allocate(tracer_names(num_tracers))

       tracers(1) = tracer_type("tracer", default_decay_rate)
       do i = 1, num_tracers
          tracer_names(i) = tracers(i)%name
       end do

       call dm_set_data_layout(dm_aux, tracer_num_components, &
         tracer_dim, tracer_names)
       call DMCreateGlobalVector(dm_aux, aux_solution, ierr); CHKERRQ(ierr)
       call PetscObjectSetName(aux_solution, "tracer", ierr); CHKERRQ(ierr)
       call global_vec_range_start(aux_solution, range_start)
       call DMDestroy(dm_aux, ierr); CHKERRQ(ierr)

       deallocate(tracer_names)
    else
       auxiliary = PETSC_FALSE
       num_tracers = 0
       allocate(tracers(num_tracers))
    end if

  end subroutine setup_tracer

!------------------------------------------------------------------------

end module tracer_module
