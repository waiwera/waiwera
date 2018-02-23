!   Copyright 2018 University of Auckland.

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

module hdf5io_module
  !! Module for HDF5 input/output.

#include <petsc/finclude/petsc.h>

  use petsc

  implicit none
  private

  PetscInt, parameter, public :: max_field_name_length = 128

  public :: vec_view_fields_hdf5, vec_load_fields_hdf5

contains

!------------------------------------------------------------------------

  subroutine dm_time_view_hdf5(dm, time_index, time, viewer)
    !! Views time sequence values to HDF5 viewer. Based on
    !! DMSequenceView_HDF5().

    DM, intent(in) :: dm
    PetscInt, intent(in) :: time_index
    PetscScalar, intent(in) :: time
    PetscViewer, intent(in out) :: viewer
    ! Locals:
    Vec :: times
    PetscInt :: localsize
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr); CHKERRQ(ierr)

    if (time_index >= 0) then

       if (rank == 0) then
          localsize = 1
       else
          localsize = 0
       end if
       call VecCreateMPI(PETSC_COMM_WORLD, localsize, 1, times, ierr)
       CHKERRQ(ierr)
       call VecSetBlockSize(times, 1, ierr); CHKERRQ(ierr)
       call PetscObjectSetName(times, "time", ierr); CHKERRQ(ierr)

       if (rank == 0) then
          call VecSetValue(times, 0, time, INSERT_VALUES, ierr)
          CHKERRQ(ierr)
       end if
       call VecAssemblyBegin(times, ierr); CHKERRQ(ierr)
       call VecAssemblyEnd(times, ierr); CHKERRQ(ierr)
       call PetscViewerHDF5PushGroup(viewer, "/", ierr); CHKERRQ(ierr)
       call PetscViewerHDF5SetTimestep(viewer, time_index, ierr)
       CHKERRQ(ierr)
       call VecView(times, viewer, ierr); CHKERRQ(ierr)
       call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
       call VecDestroy(times, ierr); CHKERRQ(ierr)

    end if
    
  end subroutine dm_time_view_hdf5

!------------------------------------------------------------------------

  subroutine vec_view_fields_hdf5(v, field_indices, field_group, viewer)
    !! Views specified fields of vector v to specified group in HDF5
    !! viewer. Based on VecView_Plex_Local_HDF5_Internal(). Currently
    !! only for scalar fields.

    use dm_utils_module, only: section_get_field_vector

    Vec, intent(in) :: v !! Vector to view
    PetscInt, intent(in) :: field_indices(:) !! Indices of fields to view
    character(*), intent(in) :: field_group !! Group to write to in HDF5 file
    PetscViewer :: viewer !! HDF5 viewer
    ! Locals:
    PetscInt, parameter :: max_name_length = 64
    DM :: dm
    PetscSection :: section, global_section
    Vec :: subv
    IS :: index_set
    PetscInt :: time_index, i, f
    PetscReal :: time
    character(max_name_length) :: name
    character(max_field_name_length) :: field_name
    PetscErrorCode :: ierr

    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMGetDefaultSection(dm, section, ierr); CHKERRQ(ierr)
    call DMGetOutputSequenceNumber(dm, time_index, time, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5SetTimestep(viewer, time_index, ierr); CHKERRQ(ierr)    
    call dm_time_view_hdf5(dm, time_index, time, viewer)
    call DMGetDefaultGlobalSection(dm, global_section, ierr); CHKERRQ(ierr)
    call PetscObjectGetName(v, name, ierr); CHKERRQ(ierr)

    do i = 1, size(field_indices)
       f = field_indices(i)
       call PetscSectionGetFieldName(section, f, field_name, ierr)
       CHKERRQ(ierr)
       if (field_name /= "") then
          call PetscViewerHDF5PushGroup(viewer, field_group, ierr)
          CHKERRQ(ierr)
          call section_get_field_vector(section, global_section, v, &
               f, index_set, subv)
          call PetscObjectSetName(subv, &
               trim(name) // "_" // trim(field_name), ierr); CHKERRQ(ierr)
          call VecView(subv, viewer, ierr); CHKERRQ(ierr)
          call VecRestoreSubVector(v, index_set, subv, ierr); CHKERRQ(ierr)
          call ISDestroy(index_set, ierr); CHKERRQ(ierr)
          call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
       end if
    end do

  end subroutine vec_view_fields_hdf5

!------------------------------------------------------------------------

  subroutine vec_load_fields_hdf5(v, field_indices, field_group, viewer)
    !! Loads specified fields of vector v from specified group in HDF5 file.

    use dm_utils_module, only: section_get_field_vector

    Vec, intent(in) :: v !! Vector to load
    PetscInt, intent(in) :: field_indices(:) !! Indices of fields to load
    character(*), intent(in) :: field_group !! Group to read from in HDF5 file
    PetscViewer :: viewer !! HDF5 viewer
    ! Locals:
    DM :: dm
    PetscSection :: section, global_section
    Vec :: subv
    IS :: index_set
    PetscInt, parameter :: max_name_length = 64
    character(max_name_length) :: name
    character(max_field_name_length) :: field_name
    PetscInt :: i, f
    PetscErrorCode :: ierr

    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMGetDefaultSection(dm, section, ierr); CHKERRQ(ierr)
    call DMGetDefaultGlobalSection(dm, global_section, ierr); CHKERRQ(ierr)
    call PetscObjectGetName(v, name, ierr); CHKERRQ(ierr)

    do i = 1, size(field_indices)
       f = field_indices(i)
       call PetscSectionGetFieldName(section, f, field_name, ierr)
       CHKERRQ(ierr)
       if (field_name /= "") then
          call PetscViewerHDF5PushGroup(viewer, field_group, ierr)
          CHKERRQ(ierr)
          call section_get_field_vector(section, global_section, v, &
               f, index_set, subv)
          call PetscObjectSetName(subv, &
               trim(name) // "_" // trim(field_name), ierr); CHKERRQ(ierr)
          call VecLoad(subv, viewer, ierr); CHKERRQ(ierr)
          call VecRestoreSubVector(v, index_set, subv, ierr); CHKERRQ(ierr)
          call ISDestroy(index_set, ierr); CHKERRQ(ierr)
          call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
       end if
    end do
    
  end subroutine vec_load_fields_hdf5

!------------------------------------------------------------------------

end module hdf5io_module
