module zone_test

  ! Tests for zone module

#include <petsc/finclude/petsc.h>

  use petsc
  use kinds_module
  use fruit
  use fson
  use zone_module
  use fson_mpi_module

  implicit none
  private

  public :: test_get_zone_type, test_cell_array

contains

!------------------------------------------------------------------------

  subroutine test_get_zone_type

    ! get_zone_type

    type(fson_value), pointer :: json
    PetscInt :: zone_type
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '[1, 2, 3]')
    zone_type = get_zone_type(json)
    if (rank == 0) then
       call assert_equals(ZONE_TYPE_CELL_ARRAY, zone_type, "array")
    end if
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"type": "array", "cells": [1, 2, 3]}')
    zone_type = get_zone_type(json)
    if (rank == 0) then
       call assert_equals(ZONE_TYPE_CELL_ARRAY, zone_type, "type array")
    end if
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"cells": [1, 2, 3]}')
    zone_type = get_zone_type(json)
    if (rank == 0) then
       call assert_equals(ZONE_TYPE_CELL_ARRAY, zone_type, "cells")
    end if
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"type": "foo"}')
    zone_type = get_zone_type(json)
    if (rank == 0) then
       call assert_equals(-1, zone_type, "unknown type")
    end if
    call fson_destroy_mpi(json)

  end subroutine test_get_zone_type

!------------------------------------------------------------------------

  subroutine test_cell_array
    ! cell array

    type(fson_value), pointer :: json
    type(zone_cell_array_type) :: zone
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '[1, 2, 3]')
    call zone%init("zone1", 1, json)
    if (rank == 0) then
       call assert_equals("zone1", zone%name, 'array name')
       call assert_equals([1, 2, 3], zone%cells, 3, 'array cells')
    end if
    call zone%destroy()
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"cells": [1, 2, 3]}')
    call zone%init("zone1", 1, json)
    if (rank == 0) then
       call assert_equals([1, 2, 3], zone%cells, 3, 'cells cells')
    end if
    call zone%destroy()
    call fson_destroy_mpi(json)

  end subroutine test_cell_array

!------------------------------------------------------------------------

end module zone_test
