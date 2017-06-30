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

    json => fson_parse_mpi(str = '{"name": "zone1", "type": "array", "cells": [1, 2, 3]}')
    if (rank == 0) then
       zone_type = get_zone_type(json)
       call assert_equals(ZONE_TYPE_CELL_ARRAY, zone_type, "type array")
    end if
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"name": "zone1", "cells": [1, 2, 3]}')
    if (rank == 0) then
       zone_type = get_zone_type(json)
       call assert_equals(ZONE_TYPE_CELL_ARRAY, zone_type, "cells")
    end if
    call fson_destroy_mpi(json)

    json => fson_parse_mpi(str = '{"name": "zone1", "type": "foo"}')
    if (rank == 0) then
       zone_type = get_zone_type(json)
       call assert_equals(-1, zone_type, "unknown type")
    end if
    call fson_destroy_mpi(json)

  end subroutine test_get_zone_type

!------------------------------------------------------------------------

  subroutine test_cell_array
    ! cell array

    type(fson_value), pointer :: json
    type(zone_cell_array_type) :: zone
    PetscErrorCode :: err
    PetscMPIInt :: rank
    PetscErrorCode :: ierr

    call MPI_comm_rank(PETSC_COMM_WORLD, rank, ierr)

    json => fson_parse_mpi(str = '{"name": "zone1", "cells": [1, 2, 3]}')
    if (rank == 0) then
       call zone%init(json, err = err)
       call assert_equals(0, err, 'array err')
       call assert_equals("zone1", zone%name, 'cells name')
       call assert_equals([1, 2, 3], zone%cells, 3, 'cells cells')
       call zone%destroy()
    end if
    call fson_destroy_mpi(json)

  end subroutine test_cell_array

!------------------------------------------------------------------------

end module zone_test
