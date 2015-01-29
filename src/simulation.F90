module simulation_module

  ! Type for representing a simulation.

  use kinds_module
  use timestepping_module
  use thermodynamics_module

  implicit none

  private

#include <petsc-finclude/petscsys.h>
#include <petsc-finclude/petscdef.h>
#include <petsc-finclude/petscdm.h>

  type, public :: simulation_type
     private
     type(timestepper_type) :: timestepper
     ! type(thermodynamics_type) :: thermo
     ! type(eos_type) :: eos  ... when EOS type defined
     ! output ?
     DM :: dm
   contains
     private
     procedure, public :: init => simulation_init
     procedure, public :: run => simulation_run
     procedure, public :: destroy => simulation_destroy
  end type simulation_type

contains

!------------------------------------------------------------------------

  subroutine simulation_init(self, filename)

    ! Initializes a simulation using data from the input file with 
    ! specified name.

    class(simulation_type), intent(in out) :: self
    character(*), intent(in) :: filename

    ! Open input file, read data and initialize simulation
    ! including self%timestepper, self%thermo, self%dm etc.

  end subroutine simulation_init

!------------------------------------------------------------------------

  subroutine simulation_run(self)

    ! Runs simulation.

    class(simulation_type), intent(in out) :: self

    ! call self%timestepper%run()

    ! maybe some final output?

  end subroutine simulation_run

!------------------------------------------------------------------------

  subroutine simulation_destroy(self)

    class(simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    ! call self%timestepper%destroy()
    ! call self%thermo%destroy()
    ! call DMDestroy(self%dm, ierr); CHKERRQ(ierr)
    ! etc.

  end subroutine simulation_destroy

!------------------------------------------------------------------------

end module simulation_module
