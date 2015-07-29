module boundary_module
  !! Module for implementation of simulation boundary conditions.

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  character(len = 16), public :: open_boundary_label_name = "open boundary"

  public :: setup_boundaries

contains

!------------------------------------------------------------------------

  subroutine setup_boundaries(json, num_primary, dm, bcs)
    !! Sets up labels identifying boundaries of the mesh, and returns
    !! array of boundary condition values for each boundary.

    use fson
    use fson_mpi_module

    type(fson_value), pointer, intent(in) :: json
    PetscInt, intent(in) :: num_primary
    DM, intent(in) :: dm
    PetscReal, allocatable, intent(out) :: bcs(:,:)
    ! Locals:
    PetscErrorCode :: ierr
    PetscBool :: has_label
    type(fson_value), pointer :: boundaries, bdy
    PetscInt :: num_boundaries, num_faces, ibdy, iface, f
    PetscInt, allocatable :: default_faces(:)
    PetscInt, allocatable :: faces(:)
    PetscInt :: default_region = 1
    PetscInt :: region
    PetscReal :: default_values(1) = [1.0e5]
    PetscReal, allocatable :: values(:)

    default_faces = [PetscInt::] ! empty integer array

    call DMPlexHasLabel(dm, open_boundary_label_name, has_label, &
         ierr); CHKERRQ(ierr)

    if (.not.(has_label)) then

       call DMPlexCreateLabel(dm, open_boundary_label_name, &
            ierr); CHKERRQ(ierr)

       if (fson_has_mpi(json, "boundaries")) then
          call fson_get_mpi(json, "boundaries", boundaries)
          num_boundaries = fson_value_count_mpi(boundaries, ".")
          allocate(bcs(num_primary + 1, num_boundaries))
          do ibdy = 1, num_boundaries
             bdy => fson_value_get_mpi(boundaries, ibdy)
             call fson_get_mpi(bdy, "faces", default_faces, faces)
             num_faces = size(faces)
             do iface = 1, num_faces
                f = faces(iface)
                call DMPlexSetLabelValue(dm, open_boundary_label_name, &
                     f, ibdy, ierr); CHKERRQ(ierr)
             end do
             call fson_get_mpi(bdy, "region", default_region, region)
             call fson_get_mpi(bdy, "values", default_values, values)
             bcs(1, ibdy) = dble(region)
             bcs(2 : num_primary + 1, ibdy) = values(1 : num_primary)
          end do
       end if

    end if

  end subroutine setup_boundaries

!------------------------------------------------------------------------

end module boundary_module
