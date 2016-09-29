module source_control_module
  !! Module for source controls- for controlling source parameters (e.g. flow rate, enthalpy) over time.

  use kinds_module
  use list_module
  use source_module
  use interpolation_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

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

     subroutine source_control_update_procedure(self, t, interval)
       !! Updates sources at the specified time.
       import :: source_control_type
       class(source_control_type), intent(in out) :: self
       PetscReal, intent(in) :: t, interval(2)
     end subroutine source_control_update_procedure

  end interface

contains

!------------------------------------------------------------------------
! Source control table:
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
! Source control rate table:
!------------------------------------------------------------------------

  subroutine source_control_rate_table_update(self, t, interval)
    !! Update flow rate for source_control_rate_table_type.

    class(source_control_rate_table_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
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
! Source control enthalpy table:
!------------------------------------------------------------------------

  subroutine source_control_enthalpy_table_update(self, t, interval)
    !! Update injection enthalpy for source_control_enthalpy_table_type.

    class(source_control_enthalpy_table_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
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

end module source_control_module
