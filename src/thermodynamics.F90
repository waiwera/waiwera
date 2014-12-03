module thermodynamics_module

  ! Thermodynamics constants and abstract types

  use kinds_module

  implicit none
  private

!------------------------------------------------------------------------
  ! Physical constants
!------------------------------------------------------------------------

  real(dp), parameter, public :: rconst     = 0.461526e3_dp     ! Gas constant
  real(dp), parameter, public :: tc_k       = 273.15_dp         ! Conversion from Celsius to Kelvin
  real(dp), parameter, public :: tcriticalk = 647.096_dp        ! Critical temperature (Kelvin)
  real(dp), parameter, public :: tcritical  = tcriticalk - tc_k ! Critical temperature (Celcius)
  real(dp), parameter, public :: dcritical  = 322.0_dp          ! Critical density (kg/m3)
  real(dp), parameter, public :: pcritical  = 22.064e6_dp       ! Critical pressure (Pa)

!------------------------------------------------------------------------
  ! Thermodynamic region type
!------------------------------------------------------------------------

  type, public, abstract :: region_type
     contains
       private
       procedure(region_init), public, deferred :: init
       procedure(region_destroy), public, deferred :: destroy
       procedure(region_properties), public, deferred :: properties
  end type region_type

  ! Pointer to region:
  type, public :: pregion_type
     class(region_type), pointer, public :: ptr
   contains
     procedure, public :: set => pregion_set
  end type pregion_type

  ! Thermodynamics type:
  type, public, abstract :: thermodynamics_type
     private
     integer, public :: num_regions
     type(pregion_type), allocatable, public :: region(:)
   contains
     private
     procedure(thermodynamics_init_procedure), public, deferred :: init
     procedure(thermodynamics_destroy_procedure), public, deferred :: destroy
  end type thermodynamics_type

  abstract interface

     subroutine region_init(self)
       import :: region_type
       implicit none
       class(region_type), intent(in out) :: self
     end subroutine region_init

     subroutine region_destroy(self)
       import :: region_type
       implicit none
       class(region_type), intent(in out) :: self
     end subroutine region_destroy

     subroutine region_properties(self, param, props, err)
       import :: region_type, dp
       implicit none
       class(region_type), intent(in out) :: self
       real(dp), intent(in), target :: param(:)
       real(dp), intent(out) :: props(:)
       integer, intent(out) :: err
     end subroutine region_properties

     subroutine thermodynamics_init_procedure(self)
       import :: thermodynamics_type
       implicit none
       class(thermodynamics_type), intent(in out) :: self
     end subroutine thermodynamics_init_procedure

     subroutine thermodynamics_destroy_procedure(self)
       import :: thermodynamics_type
       implicit none
       class(thermodynamics_type), intent(in out) :: self
     end subroutine thermodynamics_destroy_procedure

  end interface

!------------------------------------------------------------------------

contains

!------------------------------------------------------------------------
  ! Region pointers
!------------------------------------------------------------------------

  subroutine pregion_set(self, tgt)
    
    ! Sets a region pointer. This is just a workaround to give
    ! tgt the 'target' attribute, which can't always be done as part of
    ! its declaration, e.g. if it's a component of a derived type.

    implicit none
    class(pregion_type), intent(in out) :: self
    class(region_type), target, intent(in) :: tgt

    self%ptr => tgt

  end subroutine pregion_set

!------------------------------------------------------------------------

end module thermodynamics_module
