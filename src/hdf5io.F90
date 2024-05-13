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

  public :: get_hdf5_time, vec_view_fields_hdf5, vec_load_fields_hdf5, &
       vec_sequence_view_hdf5

contains

!------------------------------------------------------------------------

  subroutine get_hdf5_time(time_index, viewer, t)
    !! Reads time from HDF5 viewer at specified time index.

    PetscInt, intent(in) :: time_index
    PetscViewer, intent(in) :: viewer
    PetscReal, intent(out) :: t
    ! Locals:
    Vec :: times
    PetscInt :: localsize
    PetscMPIInt :: rank
    PetscReal, pointer, contiguous :: times_array(:)
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
       call PetscViewerHDF5PushGroup(viewer, "/", ierr); CHKERRQ(ierr)
       call PetscViewerHDF5SetTimestep(viewer, time_index, ierr)
       CHKERRQ(ierr)
       call VecLoad(times, viewer, ierr); CHKERRQ(ierr)
       call VecGetArrayReadF90(times, times_array, ierr); CHKERRQ(ierr)
       if (rank == 0) t = times_array(1)
       call MPI_bcast(t, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)
       call VecRestoreArrayReadF90(times, times_array, ierr); CHKERRQ(ierr)
       call VecDestroy(times, ierr); CHKERRQ(ierr)
       call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
    end if

  end subroutine get_hdf5_time

!------------------------------------------------------------------------

  subroutine dm_time_view_hdf5(dm, time_index, time, viewer)
    !! Views time sequence values to HDF5 viewer. Based on
    !! DMSequenceView_HDF5().

    DM, intent(in) :: dm
    PetscInt, intent(in) :: time_index
    PetscScalar, intent(in) :: time
    PetscViewer, intent(in) :: viewer
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

  subroutine set_hdf5_field_type_attribute(viewer, field_path, dm, f)
    !! Sets HDF5 file field type attribute (if it exists) for the
    !! specified field index.

    use dm_utils_module, only: dm_global_cell_field_dof

    PetscViewer, intent(in) :: viewer !! HDF5 viewer
    character(*), intent(in) :: field_path !! Field path in HDF5 file
    DM, intent(in) :: dm !! DM
    PetscInt, intent(in) :: f !! Field index
    ! Locals:
    PetscBool :: has_field_type
    PetscInt :: dim, field_dof
    character(6) :: field_type
    PetscErrorCode :: ierr

    call DMGetDimension(dm, dim, ierr); CHKERRQ(ierr)

    call PetscViewerHDF5HasAttribute(viewer, field_path, &
         "vector_field_type", has_field_type, ierr); CHKERRQ(ierr)
    if (has_field_type) then
       field_dof = dm_global_cell_field_dof(dm, f)
       if (field_dof == dim) then
          field_type = "vector"
       else
          field_type = "scalar"
       end if
       call PetscViewerHDF5WriteAttribute(viewer, field_path, &
            "vector_field_type", PETSC_STRING, field_type, ierr)
       CHKERRQ(ierr)
    end if

  end subroutine set_hdf5_field_type_attribute

!------------------------------------------------------------------------

  subroutine vec_view_fields_hdf5(v, field_indices, field_group, viewer)
    !! Views specified fields of vector v to specified group in HDF5
    !! viewer. Based on VecView_Plex_Local_HDF5_Internal().

    use dm_utils_module, only: get_field_subvector

    Vec, intent(in) :: v !! Vector to view
    PetscInt, intent(in) :: field_indices(:) !! Indices of fields to view
    character(*), intent(in) :: field_group !! Group to write to in HDF5 file
    PetscViewer, intent(in) :: viewer !! HDF5 viewer
    ! Locals:
    PetscInt, parameter :: max_name_length = 64
    character(max_name_length) :: vector_name
    DM :: dm
    PetscSection :: section
    Vec :: sub_v
    IS :: index_set
    PetscInt :: time_index, i, f
    PetscReal :: time
    character(max_field_name_length) :: field_name
    character(:), allocatable :: subvector_name, field_path
    PetscBool :: is_timestepping
    PetscErrorCode :: ierr

    call PetscObjectGetName(v, vector_name, ierr); CHKERRQ(ierr)
    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMGetSection(dm, section, ierr); CHKERRQ(ierr)

    call PetscViewerHDF5IsTimestepping(viewer, is_timestepping, ierr); CHKERRQ(ierr)
    if (is_timestepping) then
       call DMGetOutputSequenceNumber(dm, time_index, time, ierr); CHKERRQ(ierr)
       call PetscViewerHDF5SetTimestep(viewer, time_index, ierr); CHKERRQ(ierr)
       call dm_time_view_hdf5(dm, time_index, time, viewer)
    end if

    do i = 1, size(field_indices)
       f = field_indices(i)
       call PetscSectionGetFieldName(section, f, field_name, ierr)
       CHKERRQ(ierr)
       if (field_name /= "") then
          call PetscViewerHDF5PushGroup(viewer, field_group, ierr)
          CHKERRQ(ierr)
          call get_field_subvector(v, f, index_set, sub_v)
          subvector_name = trim(vector_name) // "_" // trim(field_name)
          call PetscObjectSetName(sub_v, subvector_name, ierr); CHKERRQ(ierr)
          call VecView(sub_v, viewer, ierr); CHKERRQ(ierr)
          call VecRestoreSubVector(v, index_set, sub_v, ierr); CHKERRQ(ierr)
          call ISDestroy(index_set, ierr); CHKERRQ(ierr)
          call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
          field_path = field_group // "/" // subvector_name
          call set_hdf5_field_type_attribute(viewer, field_path, dm, f)
       end if
    end do

  end subroutine vec_view_fields_hdf5

!------------------------------------------------------------------------

  subroutine vec_load_fields_hdf5(v, field_indices, field_group, viewer, &
       cell_index, file_cell_index)
    !! Loads specified fields of vector v from specified group in HDF5
    !! file, and reorders them according to the specified cell_index
    !! ordering.

    use dm_utils_module, only: get_field_subvector, vec_reorder

    Vec, intent(in) :: v !! Vector to load
    PetscInt, intent(in) :: field_indices(:) !! Indices of fields to load
    character(*), intent(in) :: field_group !! Group to read from in HDF5 file
    PetscViewer, intent(in) :: viewer !! HDF5 viewer
    IS, intent(in) :: cell_index !! Desired cell index ordering
    IS, intent(in) :: file_cell_index !! Cell index ordering in HDF5 file
    ! Locals:
    DM :: dm
    PetscSection :: section
    Vec :: sub_v
    IS :: index_set
    PetscInt, parameter :: max_name_length = 64
    character(max_name_length) :: vector_name
    character(max_field_name_length) :: field_name
    character(:), allocatable :: subvector_name
    PetscInt :: i, f, time_index
    PetscReal :: time
    PetscBool :: has_field
    PetscErrorCode :: ierr

    call PetscObjectGetName(v, vector_name, ierr); CHKERRQ(ierr)
    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMGetSection(dm, section, ierr); CHKERRQ(ierr)

    call DMGetOutputSequenceNumber(dm, time_index, time, ierr); CHKERRQ(ierr)
    call PetscViewerHDF5SetTimestep(viewer, time_index, ierr); CHKERRQ(ierr)

    do i = 1, size(field_indices)
       f = field_indices(i)
       call PetscSectionGetFieldName(section, f, field_name, ierr)
       CHKERRQ(ierr)
       if (field_name /= "") then
          call PetscViewerHDF5PushGroup(viewer, field_group, ierr)
          CHKERRQ(ierr)
          subvector_name = trim(vector_name) // "_" // trim(field_name)
          call get_field_subvector(v, f, index_set, sub_v)
          call PetscObjectSetName(sub_v, subvector_name, ierr); CHKERRQ(ierr)
          call PetscViewerHDF5HasObject(viewer, sub_v, has_field, ierr); CHKERRQ(ierr)
          if (has_field) then
             call VecLoad(sub_v, viewer, ierr); CHKERRQ(ierr)
             call vec_reorder(sub_v, file_cell_index, cell_index)
          end if
          call VecRestoreSubVector(v, index_set, sub_v, ierr); CHKERRQ(ierr)
          call ISDestroy(index_set, ierr); CHKERRQ(ierr)
          call PetscViewerHDF5PopGroup(viewer, ierr); CHKERRQ(ierr)
       end if
    end do

  end subroutine vec_load_fields_hdf5

!------------------------------------------------------------------------

  subroutine vec_sequence_view_hdf5(v, field_indices, field_group, &
       time_index, time, viewer)
    !! Views a vector v to the specified HDF5 viewer, also setting the
    !! vector's DM output sequence number.

    Vec, intent(in) :: v
    PetscInt, intent(in) :: field_indices(:)
    character(*), intent(in) :: field_group
    PetscInt, intent(in) :: time_index
    PetscReal, intent(in) :: time
    PetscViewer, intent(in out) :: viewer
    ! Locals:
    DM :: dm
    PetscErrorCode :: ierr

    call VecGetDM(v, dm, ierr); CHKERRQ(ierr)
    call DMSetOutputSequenceNumber(dm, time_index, time, &
         ierr); CHKERRQ(ierr)
    call vec_view_fields_hdf5(v, field_indices, field_group, viewer)

  end subroutine vec_sequence_view_hdf5

!------------------------------------------------------------------------

end module hdf5io_module
