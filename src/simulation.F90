module simulation_module
  !! Module for high-level representation of a simulation.

  use kinds_module
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

  integer, parameter :: max_title_length = 120

  type, public :: simulation_type
     !! Simulation type.
     private
     type(timestepper_type) :: timestepper
     class(thermodynamics_type), pointer :: thermo
     class(eos_type), pointer :: eos
     MPI_Comm :: comm
     PetscMPIInt :: rank, io_rank = 0
     DM :: dm
     character(max_title_length), public :: title
   contains
     private
     procedure :: read_title => simulation_read_title
     procedure, public :: init => simulation_init
     procedure, public :: run => simulation_run
     procedure, public :: destroy => simulation_destroy
  end type simulation_type

contains

!------------------------------------------------------------------------

  subroutine simulation_read_title(self, json)
    !! Reads simulation title from JSON input file.

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscErrorCode :: ierr

    if (self%rank == self%io_rank) then
       call fson_get(json, "title", self%title)
    end if

    ! Broadcast to all processors:
    call MPI_bcast(self%title, max_title_length, MPI_CHARACTER, &
         self%io_rank, self%comm, ierr)

  end subroutine simulation_read_title
    
!------------------------------------------------------------------------

  subroutine simulation_init(self, filename, comm)
    !! Initializes a simulation using data from the input file with 
    !! specified name.

    class(simulation_type), intent(in out) :: self
    character(*), intent(in) :: filename !! Input file name
    MPI_Comm, intent(in) :: comm !! MPI communicator
    ! Locals:
    type(fson_value), pointer :: json
    PetscErrorCode :: ierr

    self%comm = comm
    call MPI_comm_rank(self%comm, self%rank, ierr); CHKERRQ(ierr)

    if (self%rank == self%io_rank) then
       json => fson_parse(filename)
    end if

    call self%read_title(json)

    self%thermo => IAPWS

    call self%thermo%init()

    self%eos => eos_w
    call self%eos%init(self%thermo)

    if (self%rank == self%io_rank) then
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

    ! call self%timestepper%destroy()
    call self%thermo%destroy()
    ! call DMDestroy(self%dm, ierr); CHKERRQ(ierr)
    ! etc.

  end subroutine simulation_destroy

!------------------------------------------------------------------------

end module simulation_module
