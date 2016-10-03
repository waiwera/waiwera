module source_control_module
  !! Module for source controls- for controlling source parameters (e.g. flow rate, enthalpy) over time.

  use kinds_module
  use eos_module, only: max_phase_name_length
  use list_module
  use source_module
  use interpolation_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  character(max_phase_name_length), parameter, public :: default_limiter_type = "total"
  PetscReal, parameter, public :: default_limiter_limit = 1._dp
  PetscReal, parameter, public :: default_deliverability_productivity_index = 1.e-11_dp
  PetscReal, parameter, public :: default_deliverability_bottomhole_pressure = 1.e5_dp

  type, public, abstract :: source_control_type
     !! Abstract type for source control, controlling source
     !! parameters over time, for one or more sources.
     private
     type(list_type), public :: sources
   contains
     procedure(source_control_init_procedure), public, deferred :: init
     procedure(source_control_destroy_procedure), public, deferred :: destroy
     procedure(source_control_update_procedure), public, deferred :: update
  end type source_control_type

  type, abstract, extends(source_control_type) :: source_control_table_type
     !! Controls a source parameter (e.g. rate or enthalpy) via a
     !! table of values vs. time.
     private
     type(interpolation_table_type), public :: table !! Table of values vs. time
   contains
     private
     procedure, public :: init => source_control_table_init
     procedure, public :: destroy => source_control_table_destroy
  end type source_control_table_type

  type, public, extends(source_control_table_type) :: source_control_rate_table_type
     !! Controls source rate via a table of values vs. time.
   contains
     private
     procedure, public :: update => source_control_rate_table_update
  end type source_control_rate_table_type

  type, public, extends(source_control_table_type) :: source_control_enthalpy_table_type
     !! Controls source injection enthalpy via a table of values vs. time.
   contains
     private
     procedure, public :: update => source_control_enthalpy_table_update
  end type source_control_enthalpy_table_type

  type, public, extends(source_control_type) :: source_control_deliverability_type
     !! Controls a source on deliverability.
     private
     PetscReal, public :: productivity_index !! Productivity index
     PetscReal, public :: bottomhole_pressure
   contains
     procedure, public :: init => source_control_deliverability_init
     procedure, public :: destroy => source_control_deliverability_destroy
     procedure, public :: update => source_control_deliverability_update
  end type source_control_deliverability_type

  type, public, extends(source_control_type) :: source_control_limiter_type
     !! Limits flow in a particular phase (or total flow) through a source.
     private
     PetscInt, public :: phase !! Which phase is limited (0 for total)
     PetscReal, public :: limit !! Flow limit
   contains
     private
     procedure :: rate_scale => source_control_limiter_rate_scale
     procedure, public :: init => source_control_limiter_init
     procedure, public :: destroy => source_control_limiter_destroy
     procedure, public :: update => source_control_limiter_update
  end type source_control_limiter_type

  abstract interface

     subroutine source_control_init_procedure(self)
       !! Initialises a source control object.
       import :: source_control_type
       class(source_control_type), intent(in out) :: self
     end subroutine source_control_init_procedure

     subroutine source_control_destroy_procedure(self)
       !! Destroys a source control object.
       import :: source_control_type
       class(source_control_type), intent(in out) :: self
     end subroutine source_control_destroy_procedure

     subroutine source_control_update_procedure(self, t, interval, &
          fluid_data, fluid_section)
       !! Updates sources at the specified time.
       import :: source_control_type
       class(source_control_type), intent(in out) :: self
       PetscReal, intent(in) :: t, interval(2)
       PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
       PetscSection, intent(in) :: fluid_section
     end subroutine source_control_update_procedure

  end interface

contains

!------------------------------------------------------------------------
! Table source control:
!------------------------------------------------------------------------

  subroutine source_control_table_init(self)
    !! Initialises source_control_table object.

    class(source_control_table_type), intent(in out) :: self

    call self%sources%init()

  end subroutine source_control_table_init

!------------------------------------------------------------------------

  subroutine source_control_table_destroy(self)
    !! Destroys source_control_table_type object.

    class(source_control_table_type), intent(in out) :: self

    call self%table%destroy()
    call self%sources%destroy()

  end subroutine source_control_table_destroy

!------------------------------------------------------------------------
! Rate table source control:
!------------------------------------------------------------------------

  subroutine source_control_rate_table_update(self, t, interval, &
       fluid_data, fluid_section)
    !! Update flow rate for source_control_rate_table_type.

    class(source_control_rate_table_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    ! Locals:
    PetscReal :: rate

    rate = self%table%average(interval)
    call self%sources%traverse(source_control_rate_table_update_iterator)

  contains

    subroutine source_control_rate_table_update_iterator(node, stopped)
      !! Sets source rate at a list node.
      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      select type (source => node%data)
      type is (source_type)
         source%rate = rate
      end select
      stopped = PETSC_FALSE
    end subroutine source_control_rate_table_update_iterator

  end subroutine source_control_rate_table_update

!------------------------------------------------------------------------
! Enthalpy table source control:
!------------------------------------------------------------------------

  subroutine source_control_enthalpy_table_update(self, t, interval, &
       fluid_data, fluid_section)
    !! Update injection enthalpy for source_control_enthalpy_table_type.

    class(source_control_enthalpy_table_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    ! Locals:
    PetscReal :: enthalpy

    enthalpy = self%table%average(interval)
    call self%sources%traverse(source_control_enthalpy_table_update_iterator)

  contains

    subroutine source_control_enthalpy_table_update_iterator(node, stopped)
      !! Sets source injection enthalpy at a list node.
      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      select type (source => node%data)
      type is (source_type)
         source%injection_enthalpy = enthalpy
      end select
      stopped = PETSC_FALSE
    end subroutine source_control_enthalpy_table_update_iterator

  end subroutine source_control_enthalpy_table_update

!------------------------------------------------------------------------
! Deliverability source control:
!------------------------------------------------------------------------

  subroutine source_control_deliverability_init(self)
    !! Initialises source_control_deliverability object.

    class(source_control_deliverability_type), intent(in out) :: self

    call self%sources%init()

  end subroutine source_control_deliverability_init

!------------------------------------------------------------------------

  subroutine source_control_deliverability_destroy(self)
    !! Destroys source_control_deliverability object.

    class(source_control_deliverability_type), intent(in out) :: self

    call self%sources%destroy()

  end subroutine source_control_deliverability_destroy

!------------------------------------------------------------------------

  subroutine source_control_deliverability_update(self, t, interval, &
       fluid_data, fluid_section)
    !! Update flow rate for source_control_deliverability_type.

    class(source_control_deliverability_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section

    call self%sources%traverse(source_control_deliverability_update_iterator)

  contains

    subroutine source_control_deliverability_update_iterator(node, stopped)
      !! Applies deliverability control at a list node.
      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped

      select type (source => node%data)
      type is (source_type)

         call source%update_fluid(fluid_data, fluid_section)

         source%rate = self%productivity_index * (source%fluid%pressure - &
              self%bottomhole_pressure)

      end select

      stopped = PETSC_FALSE

    end subroutine source_control_deliverability_update_iterator

  end subroutine source_control_deliverability_update

!------------------------------------------------------------------------
! Limiter source control:
!------------------------------------------------------------------------

  subroutine source_control_limiter_init(self)
    !! Initialises source_control_limiter object.

    class(source_control_limiter_type), intent(in out) :: self

    call self%sources%init()

  end subroutine source_control_limiter_init

!------------------------------------------------------------------------

  subroutine source_control_limiter_destroy(self)
    !! Destroys source_control_limiter_type object.

    class(source_control_limiter_type), intent(in out) :: self

    call self%sources%destroy()

  end subroutine source_control_limiter_destroy

!------------------------------------------------------------------------

  PetscReal function source_control_limiter_rate_scale(self, rate) &
       result(scale)
    !! Returns factor by which flow should be scaled to avoid
    !! exceededing limit.

    class(source_control_limiter_type), intent(in) :: self
    PetscReal, intent(in) :: rate
    ! Locals:
    PetscReal, parameter :: small = 1.e-6_dp

    if ((rate > self%limit) .and. (rate > small)) then
       scale = self%limit / rate
    else
       scale = 1._dp
    end if

  end function source_control_limiter_rate_scale

!------------------------------------------------------------------------

  subroutine source_control_limiter_update(self, t, interval, &
       fluid_data, fluid_section)
    !! Update flow rate for source_control_limiter_type.

    class(source_control_limiter_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section

    call self%sources%traverse(source_control_limiter_update_iterator)

  contains

    subroutine source_control_limiter_update_iterator(node, stopped)
      !! Limits source rate at a list node.

      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscReal :: rate, scale

      select type (source => node%data)
      type is (source_type)

         if (self%phase == 0) then ! Limit total flow:
            rate = abs(source%rate)
         else
            call source%update_fluid(fluid_data, fluid_section)
            rate = source%phase_flow_fractions(self%phase) * &
                 abs(source%rate)
         end if

         scale = self%rate_scale(rate)
         source%rate = source%rate * scale

      end select

      stopped = PETSC_FALSE

    end subroutine source_control_limiter_update_iterator

  end subroutine source_control_limiter_update

!------------------------------------------------------------------------

end module source_control_module
