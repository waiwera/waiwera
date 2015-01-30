module eos_module

  ! Abstract type for equation of state (EOS) objects

  use kinds_module

  implicit none
  private

  type, public, abstract :: eos_type
     private
     integer, public :: num_primary, num_secondary
   contains
     private
     procedure(eos_transition_procedure), deferred :: transition
     procedure(eos_check_primary_procedure), public, deferred :: check_primary
     procedure(eos_secondary_procedure), public, deferred :: secondary
  end type eos_type

  abstract interface

     subroutine eos_transition_procedure(self, &
          region1, region2, primary)
       ! Perform transitions between thermodynamic regions
       import :: eos_type, dp
       implicit none
       class(eos_type), intent(in out) :: self
       integer, intent(in) :: region1, region2
       real(dp), intent(in out) :: primary(self%num_primary)
     end subroutine eos_transition_procedure

     subroutine eos_check_primary_procedure(self, region, primary)
       ! Check primary variables for current region and make
       ! transition if needed
       import :: eos_type, dp
       implicit none
       class(eos_type), intent(in out) :: self
       integer, intent(in out) :: region
       real(dp), intent(in out) :: primary(self%num_primary)
     end subroutine eos_check_primary_procedure

     subroutine eos_secondary_procedure(self, region, &
          primary, secondary)
       ! Calculate secondary variables from region and primary variables
       import :: eos_type, dp
       implicit none
       class(eos_type), intent(in out) :: self
       integer, intent(in) :: region
       real(dp), intent(in) :: primary(self%num_primary)
       real(dp), intent(out) :: secondary(self%num_secondary)
     end subroutine eos_secondary_procedure

  end interface

end module eos_module
