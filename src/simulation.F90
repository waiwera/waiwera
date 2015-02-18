module simulation_module
  !! Module for high-level representation of a simulation.

  use mpi_module
  use kinds_module
  use mesh_module
  use timestepping_module
  use thermodynamics_module
  use IAPWS_module
  use IFC67_module
  use eos_w_module
  use eos_module
  use fson

  implicit none

  private

#include <petsc-finclude/petscsys.h>
#include <petsc-finclude/petscdef.h>
#include <petsc-finclude/petscdm.h>

  integer, parameter, public :: max_filename_length = 200
  integer, parameter, public :: max_title_length = 120

  type, public :: simulation_type
     !! Simulation type.
     private
     type(mesh_type) :: mesh
     type(timestepper_type) :: timestepper
     class(thermodynamics_type), pointer :: thermo
     class(eos_type), pointer :: eos
     character(max_filename_length) :: input_filename
     character(max_title_length), public :: title
   contains
     private
     procedure :: read_title => simulation_read_title
     procedure :: read_thermodynamics => simulation_read_thermodynamics
     procedure :: read_eos => simulation_read_eos
     procedure, public :: init => simulation_init
     procedure, public :: run => simulation_run
     procedure, public :: destroy => simulation_destroy
  end type simulation_type

contains

!------------------------------------------------------------------------

  subroutine simulation_read_title(self, json)
    !! Reads simulation title from JSON input file.
    !! If not present, a default value is assigned.

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    type(fson_value), pointer :: title
    character(len = max_title_length), parameter :: default_title = ""
    PetscErrorCode :: ierr

    if (rank == input_rank) then
       call fson_get(json, "title", title)
       if (associated(title)) then
          call fson_get(title, ".", self%title)
       else
          self%title = default_title
       end if
    end if

    ! Broadcast to all processors:
    call MPI_bcast(self%title, max_title_length, MPI_CHARACTER, &
         input_rank, comm, ierr)

  end subroutine simulation_read_title

!------------------------------------------------------------------------

  subroutine simulation_read_thermodynamics(self, json)
    !! Reads simulation thermodynamic formulation from JSON input file.
    !! If not present, a default value is assigned.

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscErrorCode :: ierr
    type(fson_value), pointer :: thermo
    integer, parameter :: max_thermo_ID_length = 8
    character(max_thermo_ID_length) :: thermo_ID
    character(max_thermo_ID_length), parameter :: &
         default_thermo_ID = "IAPWS"

    if (rank == input_rank) then
       call fson_get(json, "thermodynamics", thermo)
       if (associated(thermo)) then
          call fson_get(thermo, ".", thermo_ID)
       else
          thermo_ID = default_thermo_ID
       end if
    end if

    ! Broadcast to all processors:
    call MPI_bcast(thermo_ID, max_thermo_ID_length, &
         MPI_CHARACTER, input_rank, comm, ierr)

    select case (thermo_ID)
    case ("IFC67")
       self%thermo => IFC67
    case default
       self%thermo => IAPWS
    end select

    call self%thermo%init()

  end subroutine simulation_read_thermodynamics

!------------------------------------------------------------------------

  subroutine simulation_read_eos(self, json)
    !! Reads simulation equation of state from JSON input file.
    !! If not present, a default value is assigned.

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscErrorCode :: ierr
    type(fson_value), pointer :: eos
    integer, parameter :: max_eos_ID_length = 12
    character(max_eos_ID_length) :: eos_ID
    character(max_eos_ID_length), parameter :: &
         default_eos_ID = "W"

    if (rank == input_rank) then
       call fson_get(json, "eos", eos)
       if (associated(eos)) then
          call fson_get(eos, ".", eos_ID)
       else
          eos_ID = default_eos_ID
       end if
    end if

    ! Broadcast to all processors:
    call MPI_bcast(eos_ID, max_eos_ID_length, &
         MPI_CHARACTER, input_rank, comm, ierr)

    select case (eos_ID)
    case ("EW")
       self%eos => eos_w  ! change to eos_ew when it's ready
    case default
       self%eos => eos_w
    end select

    call self%eos%init(self%thermo)

  end subroutine simulation_read_eos

!------------------------------------------------------------------------

  subroutine simulation_init(self, filename)
    !! Initializes a simulation using data from the input file with 
    !! specified name.

    class(simulation_type), intent(in out) :: self
    character(*), intent(in) :: filename !! Input file name
    ! Locals:
    type(fson_value), pointer :: json
    PetscErrorCode :: ierr
    PetscInt :: dof

    self%input_filename = filename

    if (rank == input_rank) then
       json => fson_parse(filename)
    end if

    call self%read_title(json)

    call self%read_thermodynamics(json)
    call self%read_eos(json)

    dof = self%eos%num_primary
    call self%mesh%init(json, dof)

    if (rank == input_rank) then
       call fson_destroy(json)
    end if

  end subroutine simulation_init

!------------------------------------------------------------------------

  subroutine simulation_run(self)
    !! Runs the simulation.

    class(simulation_type), intent(in out) :: self

    ! call self%timestepper%run()

    ! maybe some final output?

  end subroutine simulation_run

!------------------------------------------------------------------------

  subroutine simulation_destroy(self)
    !! Destroys the simulation.

    class(simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call self%mesh%destroy()
    call self%thermo%destroy()
    ! call self%timestepper%destroy()

  end subroutine simulation_destroy

!------------------------------------------------------------------------

end module simulation_module
