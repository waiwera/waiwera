module flow_simulation_module
  !! Module for high-level representation of a flow simulation ODE.

  use ode_module
  use mesh_module
  use thermodynamics_module
  use eos_module

  implicit none

  private

#include <petsc/finclude/petsc.h90>

  PetscInt, parameter, public :: max_title_length = 120

  type, public, extends(ode_type) :: flow_simulation_type
     !! Simulation type.
     private
     character(max_title_length), public :: title
     Vec, public :: rock
     Vec, public :: fluid
     class(thermodynamics_type), allocatable, public :: thermo
     class(eos_type), allocatable, public :: eos
     PetscReal, public :: gravity
   contains
     private
     procedure, public :: init => flow_simulation_init
     procedure, public :: destroy => flow_simulation_destroy
     procedure, public :: lhs => flow_simulation_cell_balances
     procedure, public :: rhs => flow_simulation_cell_inflows
     procedure, public :: pre_eval => flow_simulation_fluid_properties
     procedure, public :: output => flow_simulation_output
  end type flow_simulation_type

contains

!------------------------------------------------------------------------

  subroutine flow_simulation_init(self, json)
    !! Initializes a flow simulation using data from the specified JSON object.

    use kinds_module
    use mpi_module
    use fson
    use fson_mpi_module
    use thermodynamics_setup_module, only: setup_thermodynamics
    use eos_setup_module, only: setup_eos
    use initial_module, only: setup_initial
    use fluid_module, only: setup_fluid_vector
    use rock_module, only: setup_rock_vector, setup_rocktype_labels

    class(flow_simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    character(len = max_title_length), parameter :: default_title = ""
    PetscReal, parameter :: default_gravity = 9.8_dp
    PetscErrorCode :: ierr

    call fson_get_mpi(json, "title", default_title, self%title)
    call setup_thermodynamics(json, self%thermo)
    call setup_eos(json, self%thermo, self%eos)
    call self%mesh%init(json)
    call setup_rocktype_labels(json, self%mesh%dm)
    call self%mesh%setup_boundaries(self%eos%num_primary_variables, json)
    call self%mesh%configure(self%eos%primary_variable_names)
    call DMCreateGlobalVector(self%mesh%dm, self%solution, ierr)
    CHKERRQ(ierr)
    call PetscObjectSetName(self%solution, "primary", ierr); CHKERRQ(ierr)
    call setup_rock_vector(json, self%mesh%dm, self%rock)
    call setup_initial(json, self%mesh, self%eos%num_primary_variables, &
         self%time, self%solution, self%rock)
    call setup_fluid_vector(self%mesh%dm, self%eos%num_phases, &
         self%eos%num_components, self%fluid)
    call fson_get_mpi(json, "gravity", default_gravity, self%gravity)

  end subroutine flow_simulation_init

!------------------------------------------------------------------------

  subroutine flow_simulation_destroy(self)
    !! Destroys the simulation.

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call VecDestroy(self%solution, ierr); CHKERRQ(ierr)
    call VecDestroy(self%fluid, ierr); CHKERRQ(ierr)
    call VecDestroy(self%rock, ierr); CHKERRQ(ierr)
    call self%mesh%destroy()
    call self%thermo%destroy()
    deallocate(self%thermo)
    deallocate(self%eos)

  end subroutine flow_simulation_destroy

!------------------------------------------------------------------------

  subroutine flow_simulation_cell_balances(self, t, y, lhs)
    !! Computes mass and energy balance for each cell, for the given
    !! primary thermodynamic variables and time.

    use dm_utils_module, only: global_section_offset, global_vec_section
    use cell_module, only: cell_type

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time (s)
    Vec, intent(in) :: y !! global primary variables vector
    Vec, intent(out) :: lhs
    ! Locals:
    PetscInt :: c, ghost, np, nc
    PetscSection :: fluid_section, rock_section, lhs_section
    PetscInt :: fluid_offset, rock_offset, lhs_offset
    PetscReal, pointer :: fluid_array(:), rock_array(:), lhs_array(:)
    PetscReal, pointer :: balance(:)
    type(cell_type) :: cell
    DMLabel :: ghost_label
    PetscErrorCode :: ierr

    np = self%eos%num_primary_variables
    nc = self%eos%num_components

    call global_vec_section(lhs, lhs_section)
    call VecGetArrayF90(lhs, lhs_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%fluid, fluid_section)
    call VecGetArrayReadF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%rock, rock_section)
    call VecGetArrayReadF90(self%rock, rock_array, ierr); CHKERRQ(ierr)

    call cell%init(nc, self%eos%num_phases)

    call DMPlexGetLabel(self%mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    do c = self%mesh%start_cell, self%mesh%end_cell - 1

       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then

          call global_section_offset(lhs_section, c, lhs_offset, ierr)
          CHKERRQ(ierr)
          balance => lhs_array(lhs_offset : lhs_offset + np - 1)

          call global_section_offset(fluid_section, c, fluid_offset, ierr)
          CHKERRQ(ierr)
          call global_section_offset(rock_section, c, rock_offset, ierr)
          CHKERRQ(ierr)

          call cell%assign( &
               rock_data = rock_array, rock_offset = rock_offset, &
               fluid_data = fluid_array, fluid_offset = fluid_offset)

          balance = cell%balance(self%eos%isothermal)

       end if

    end do

    call cell%destroy()
    nullify(balance)
    call VecRestoreArrayReadF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(lhs, lhs_array, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_cell_balances

!------------------------------------------------------------------------

  subroutine flow_simulation_cell_inflows(self, t, y, rhs)
    !! Computes net inflow (per unit volume) into each cell, from
    !! flows through faces and source terms, for the given primary
    !! thermodynamic variables and time.

    use kinds_module
    use dm_utils_module
    use face_module, only: face_type

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time
    Vec, intent(in) :: y !! global primary variables vector
    Vec, intent(out) :: rhs
    ! Locals:
    PetscInt :: f, ghost_cell, ghost_face, i, np
    Vec :: local_fluid, local_rock
    PetscReal, pointer :: rhs_array(:)
    PetscReal, pointer :: cell_geom_array(:), face_geom_array(:)
    PetscReal, pointer :: fluid_array(:), rock_array(:)
    PetscSection :: rhs_section, rock_section, fluid_section
    PetscSection :: cell_geom_section, face_geom_section
    type(face_type) :: face
    PetscInt :: face_geom_offset, cell_geom_offsets(2)
    PetscInt :: rock_offsets(2), fluid_offsets(2), rhs_offsets(2)
    DMLabel :: ghost_label
    PetscInt, pointer :: cells(:)
    PetscReal, pointer :: inflow(:)
    PetscReal, allocatable :: face_flow(:)
    PetscReal, parameter :: flux_sign(2) = [-1._dp, 1._dp]
    PetscReal, allocatable :: primary(:)
    PetscErrorCode :: ierr

    np = self%eos%num_primary_variables
    allocate(face_flow(np), primary(np))

    call global_vec_section(rhs, rhs_section)
    call VecGetArrayF90(rhs, rhs_array, ierr); CHKERRQ(ierr)
    rhs_array = 0._dp

    call local_vec_section(self%mesh%cell_geom, cell_geom_section)
    call VecGetArrayReadF90(self%mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)
    call local_vec_section(self%mesh%face_geom, face_geom_section)
    call VecGetArrayReadF90(self%mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)

    call global_to_local_vec_section(self%fluid, local_fluid, fluid_section)
    call VecGetArrayReadF90(local_fluid, fluid_array, ierr); CHKERRQ(ierr)

    call global_to_local_vec_section(self%rock, local_rock, rock_section)
    call VecGetArrayReadF90(local_rock, rock_array, ierr); CHKERRQ(ierr)

    call face%init(self%eos%num_components, self%eos%num_phases)

    call DMPlexGetLabel(self%mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    do f = self%mesh%start_face, self%mesh%end_face - 1

       call DMLabelGetValue(ghost_label, f, ghost_face, ierr); CHKERRQ(ierr)
       if (ghost_face < 0) then

          call section_offset(face_geom_section, f, face_geom_offset, &
               ierr); CHKERRQ(ierr)

          call DMPlexGetSupport(self%mesh%dm, f, cells, ierr); CHKERRQ(ierr)
          do i = 1, 2
             call section_offset(cell_geom_section, cells(i), &
                  cell_geom_offsets(i), ierr); CHKERRQ(ierr)
             call section_offset(fluid_section, cells(i), &
                  fluid_offsets(i), ierr); CHKERRQ(ierr)
             call section_offset(rock_section, cells(i), &
                  rock_offsets(i), ierr); CHKERRQ(ierr)
             call global_section_offset(rhs_section, cells(i), &
                  rhs_offsets(i), ierr); CHKERRQ(ierr)
          end do

          call face%assign(face_geom_array, face_geom_offset, &
               cell_geom_array, cell_geom_offsets, &
               rock_array, rock_offsets, fluid_array, fluid_offsets)

          face_flow = face%flux(self%eos%isothermal, self%gravity) * &
               face%area

          do i = 1, 2
             call DMLabelGetValue(ghost_label, cells(i), ghost_cell, &
                  ierr); CHKERRQ(ierr)
             if ((ghost_cell < 0) .and. &
                  (cells(i) <= self%mesh%end_interior_cell - 1)) then
                inflow => rhs_array(rhs_offsets(i) : rhs_offsets(i) + np - 1)
                inflow = inflow + flux_sign(i) * face_flow / &
                     face%cell(i)%volume
             end if
          end do

       end if
    end do

    call face%destroy()
    nullify(inflow)
    call VecRestoreArrayReadF90(self%mesh%face_geom, face_geom_array, ierr)
    CHKERRQ(ierr)

    ! TODO: source term loop will go here

    call VecRestoreArrayReadF90(local_rock, rock_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(local_fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(self%mesh%cell_geom, cell_geom_array, ierr)
    CHKERRQ(ierr)
    call VecRestoreArrayF90(rhs, rhs_array, ierr); CHKERRQ(ierr)
    call restore_dm_local_vec(local_fluid)
    call restore_dm_local_vec(local_rock)
    deallocate(face_flow, primary)

  end subroutine flow_simulation_cell_inflows

!------------------------------------------------------------------------

  subroutine flow_simulation_fluid_properties(self, t, y)
    !! Computes fluid properties in all cells, based on the current time
    !! and primary thermodynamic variables.

    use dm_utils_module, only: global_section_offset, global_vec_section
    use fluid_module, only: fluid_type
    use mpi_module

    class(flow_simulation_type), intent(in out) :: self
    PetscReal, intent(in) :: t !! time
    Vec, intent(in) :: y !! global primary variables vector
    ! Locals:
    PetscInt :: c, region, np, nc, ghost
    PetscSection :: y_section, fluid_section
    PetscInt :: y_offset, fluid_offset
    PetscReal, pointer :: y_array(:), cell_primary(:)
    PetscReal, pointer :: fluid_array(:)
    type(fluid_type) :: fluid
    DMLabel :: ghost_label
    PetscErrorCode :: ierr

    np = self%eos%num_primary_variables
    nc = self%eos%num_components

    call global_vec_section(y, y_section)
    call VecGetArrayReadF90(y, y_array, ierr); CHKERRQ(ierr)

    call global_vec_section(self%fluid, fluid_section)
    call VecGetArrayF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)

    call fluid%init(nc, self%eos%num_phases)

    call DMPlexGetLabel(self%mesh%dm, "ghost", ghost_label, ierr)
    CHKERRQ(ierr)

    do c = self%mesh%start_cell, self%mesh%end_cell - 1

       call DMLabelGetValue(ghost_label, c, ghost, ierr); CHKERRQ(ierr)
       if (ghost < 0) then

          call global_section_offset(y_section, c, y_offset, ierr)
          CHKERRQ(ierr)
          cell_primary => y_array(y_offset : y_offset + np - 1)

          call global_section_offset(fluid_section, c, fluid_offset, ierr)
          CHKERRQ(ierr)

          call fluid%assign(fluid_array, fluid_offset)
          region = nint(fluid%region)

          call self%eos%fluid_properties(region, cell_primary, fluid)

       end if

    end do

    call VecRestoreArrayF90(self%fluid, fluid_array, ierr); CHKERRQ(ierr)
    call VecRestoreArrayReadF90(y, y_array, ierr); CHKERRQ(ierr)
    call fluid%destroy()

  end subroutine flow_simulation_fluid_properties

!------------------------------------------------------------------------

  subroutine flow_simulation_output(self)
    !! Output from flow simulation.

    use mpi_module

    class(flow_simulation_type), intent(in out) :: self
    ! Locals:
    PetscViewer :: viewer
    PetscErrorCode :: ierr

    call VecView(self%solution, PETSC_VIEWER_STDOUT_WORLD, ierr)
    CHKERRQ(ierr)

    call PetscViewerCreate(mpi%comm, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerSetType(viewer, PETSCVIEWERVTK, ierr); CHKERRQ(ierr)
    call PetscViewerFileSetName(viewer, "solution.vtu", ierr); CHKERRQ(ierr)
    call VecView(self%solution, viewer, ierr); CHKERRQ(ierr)
    call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)

  end subroutine flow_simulation_output

!------------------------------------------------------------------------

end module flow_simulation_module
