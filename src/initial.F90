module initial_module
  !! Module for setting up initial conditions.

  use kinds_module
  use fson
  use fson_mpi_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  public :: setup_initial

contains

!------------------------------------------------------------------------

  subroutine setup_initial(json, mesh, num_primary, t, y, rock_vector)
    !! Initializes time t and a Vec y with initial conditions read
    !! from JSON input 'initial'.  Conditions may be specified as a
    !! constant value or as an array. The array may contain a complete
    !! set of initial conditions for all cells, or if a shorter array is
    !! given, this is repeated over initial conditions vector.

    use fson_value_m, only : TYPE_REAL, TYPE_INTEGER, TYPE_ARRAY
    use mesh_module
    use boundary_module, only: open_boundary_label_name
    use rock_module
    use dm_utils_module, only: vec_section, section_offset

    type(fson_value), pointer, intent(in) :: json
    type(mesh_type), intent(in) :: mesh
    PetscInt, intent(in) :: num_primary
    PetscReal, intent(out) :: t
    Vec, intent(in out) :: y, rock_vector
    ! Locals:
    PetscErrorCode :: ierr
    PetscReal :: const_initial_value
    PetscInt :: int_const_initial_value
    PetscInt, allocatable :: indices(:)
    PetscReal, allocatable :: initial_input(:), initial_data(:)
    PetscInt :: i, np, count
    PetscBool :: const
    DMLabel :: bdy_label
    IS :: bdy_IS
    PetscReal, pointer:: y_array(:), rock_array(:), cell_primary(:)
    PetscSection :: y_section, rock_section
    PetscInt :: ibdy, f, iface, num_faces
    PetscInt :: rock_offsets(2), rock_dof, y_offset
    PetscInt, pointer :: bdy_faces(:), cells(:)
    type(rock_type) :: rock
    PetscReal, parameter :: default_start_time = 0.0_dp
    PetscReal, parameter :: default_initial_value = 0.0_dp

    call fson_get_mpi(json, "time.start", default_start_time, t)

    call VecGetSize(y, count, ierr); CHKERRQ(ierr)
    const = .true.

    if (fson_has_mpi(json, "initial")) then

       select case (fson_type_mpi(json, "initial"))
          case (TYPE_REAL)
             call fson_get_mpi(json, "initial", val = const_initial_value)
          case (TYPE_INTEGER)
             call fson_get_mpi(json, "initial", val = int_const_initial_value)
             const_initial_value = real(int_const_initial_value)
          case (TYPE_ARRAY)
             const = .false.
             call fson_get_mpi(json, "initial", val = initial_input)
             np = size(initial_input)
             if (np >= count) then
                initial_data = initial_input(1:count)
             else ! repeat input over array:
                do i = 1, np
                   initial_data(i:count:np) = initial_input(i)
                end do
             end if
             deallocate(initial_input)
       end select
    else
       const_initial_value = default_initial_value
    end if

    if (const) then
       call VecSet(y, const_initial_value, ierr); CHKERRQ(ierr)
    else
       allocate(indices(count))
       do i = 1, count
          indices(i) = i-1
       end do
       call VecSetValues(y, count, indices, &
            initial_data, INSERT_VALUES, ierr); CHKERRQ(ierr)
       call VecAssemblyBegin(y, ierr); CHKERRQ(ierr)
       call VecAssemblyEnd(y, ierr); CHKERRQ(ierr)
       deallocate(indices, initial_data)
    end if

    ! Boundary conditions:
    call VecGetArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call vec_section(y, y_section)
    call VecGetArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)
    call vec_section(rock_vector, rock_section)
    call DMPlexGetLabel(mesh%dm, open_boundary_label_name, &
         bdy_label, ierr); CHKERRQ(ierr)
    do ibdy = 1, size(mesh%bcs, 1)
       call DMPlexGetStratumIS(mesh%dm, open_boundary_label_name, &
            ibdy, bdy_IS, ierr); CHKERRQ(ierr)
       if (bdy_IS /= 0) then
          call ISGetIndicesF90(bdy_IS, bdy_faces, ierr); CHKERRQ(ierr)
          num_faces = size(bdy_faces)
          do iface = 1, num_faces
             f = bdy_faces(iface)
             call DMPlexGetSupport(mesh%dm, f, cells, ierr); CHKERRQ(ierr)
             call section_offset(y_section, cells(2), y_offset, ierr)
             do i = 1, 2
                call section_offset(rock_section, cells(i), &
                     rock_offsets(i), ierr); CHKERRQ(ierr)
             end do
             ! Set primary variables:
             cell_primary => y_array(y_offset : y_offset + num_primary - 1)
             cell_primary = mesh%bcs(2: num_primary + 1, ibdy)
             ! Copy rock type data from interior cell to boundary ghost cell:
             rock_dof = rock%dof()
             rock_array(rock_offsets(2) : rock_offsets(2) + rock_dof) = &
                  rock_array(rock_offsets(1) : rock_offsets(1) + rock_dof)
          end do
       end if
    end do
    call VecRestoreArrayF90(y, y_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(rock_vector, rock_array, ierr); CHKERRQ(ierr)

  end subroutine setup_initial

!------------------------------------------------------------------------

end module
