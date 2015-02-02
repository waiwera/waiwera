module simulation_module

  ! Type for representing a simulation.

  use kinds_module
  use timestepping_module
  use thermodynamics_module
  use IAPWS_module
  use IFC67_module
  use eos_module
  use eos_w_module

  implicit none

  private

#include <petsc-finclude/petscsys.h>
#include <petsc-finclude/petscdef.h>
#include <petsc-finclude/petscdm.h>

  type, public :: simulation_type
     private
     type(timestepper_type) :: timestepper
     class(thermodynamics_type), pointer :: thermo
     class(eos_type), pointer :: eos
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
    ! including self%timestepper, self%dm etc.
    
    self%thermo => IAPWS

    call self%thermo%init()

    self%eos => eos_w
    call self%eos%init(self%thermo)

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
    call self%thermo%destroy()
    ! call DMDestroy(self%dm, ierr); CHKERRQ(ierr)
    ! etc.

  end subroutine simulation_destroy

!------------------------------------------------------------------------

end module simulation_module
