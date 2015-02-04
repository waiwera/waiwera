module eos_w_module

  ! Isothermal pure water equation of state

  use kinds_module
  use eos_module
  use thermodynamics_module

  implicit none
  private

  type, public, extends(eos_type) :: eos_w_type
     private
     real(dp), public :: temperature = 20._dp
   contains
     private
     procedure :: transition => eos_w_transition
     procedure, public :: init => eos_w_init
     procedure, public :: check_primary => eos_w_check_primary
     procedure, public :: secondary => eos_w_secondary
  end type eos_w_type

  type(eos_w_type), public, target :: eos_w

contains

!------------------------------------------------------------------------

  subroutine eos_w_init(self, thermo)
    ! Initialise isothermal pure water EOS

    class(eos_w_type), intent(in out) :: self
    class(thermodynamics_type), intent(in), target :: thermo

    self%num_primary = 1
    self%num_secondary = 3

    self%thermo => thermo

  end subroutine eos_w_init

!------------------------------------------------------------------------
  
  subroutine eos_w_transition(self, region1, region2, primary)
    ! Perform transitions between thermodynamic regions for isothermal
    ! pure water 

    class(eos_w_type), intent(in) :: self
    integer, intent(in) :: region1, region2
    real(dp), intent(in out), target :: primary(self%num_primary)

    ! no transitions needed for isothermal pure water
    continue

  end subroutine eos_w_transition

!------------------------------------------------------------------------

  subroutine eos_w_check_primary(self, region, primary)
    ! Check primary variables for current region and make
    ! transition if needed for isothermal pure water

    class(eos_w_type), intent(in) :: self
    integer, intent(in out) :: region
    real(dp), intent(in out), target :: primary(self%num_primary)

    ! no checks needed for isothermal pure water
    continue

  end subroutine eos_w_check_primary

!------------------------------------------------------------------------

  subroutine eos_w_secondary(self, region, primary, secondary)

    ! Calculate secondary variables from region and primary variables

    class(eos_w_type), intent(in out) :: self
    integer, intent(in) :: region
    real(dp), intent(in), target :: primary(self%num_primary)
    real(dp), intent(out), target :: secondary(self%num_secondary)
    ! Locals:
    integer :: err
    real(dp), pointer :: pressure, density, viscosity

    pressure => primary(1)
    density => secondary(1)
    viscosity => secondary(3)

    call self%thermo%region(1)%ptr%properties([pressure, self%temperature],\
    secondary(1:2), err)
    call self%thermo%region(1)%ptr%viscosity(self%temperature, pressure, \
    density, viscosity)

  end subroutine eos_w_secondary

!------------------------------------------------------------------------

end module eos_w_module
