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

module profiling_module
  !! Module for profiling code via PETSc log events.

  implicit none
  private

#include <petsc/finclude/petscsys.h>

  PetscClassId, public ::  log_class
  PetscLogEvent, public :: simulation_init_event
  PetscLogEvent, public :: fluid_init_event
  PetscLogEvent, public :: fluid_properties_event, fluid_transitions_event
  PetscLogEvent, public :: cell_balances_event
  PetscLogEvent, public :: cell_inflows_event, sources_event
  PetscLogEvent, public :: output_event

  public :: init_profiling

contains

!------------------------------------------------------------------------

  subroutine init_profiling()
    !! Initialize code profiling, registering PETSc class ID and all
    !! log events needed anywhere in the code.

    ! Locals:
    PetscErrorCode :: ierr

    call PetscClassIdRegister("waiwera", log_class, ierr); CHKERRQ(ierr)

    ! Register log events:
    call PetscLogEventRegister("sim_init", log_class, simulation_init_event, ierr)
    call PetscLogEventRegister("fluid_init", log_class, fluid_init_event, ierr)
    call PetscLogEventRegister("fluid_props", log_class, &
         fluid_properties_event, ierr)
    call PetscLogEventRegister("fluid_trans", log_class, &
         fluid_transitions_event, ierr)
    call PetscLogEventRegister("cell_balances", log_class, cell_balances_event, ierr)
    call PetscLogEventRegister("output", log_class, output_event, ierr)
    call PetscLogEventRegister("cell_inflows", log_class, &
         cell_inflows_event, ierr)
    call PetscLogEventRegister("sources", log_class, &
         sources_event, ierr)

  end subroutine init_profiling

!------------------------------------------------------------------------

end module profiling_module
