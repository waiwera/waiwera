module source_control_module
  !! Module for source controls- for controlling source parameters (e.g. flow rate, enthalpy) over time.

  use kinds_module
  use eos_module, only: max_phase_name_length
  use thermodynamics_module
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
   contains
     procedure(source_control_destroy_procedure), public, deferred :: destroy
     procedure(source_control_update_procedure), public, deferred :: update
  end type source_control_type

  type, abstract, extends(source_control_type) :: source_control_table_type
     !! Controls a source parameter (e.g. rate or enthalpy) via a
     !! table of values vs. time.
     private
     type(list_type), public :: sources
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
     type(list_type), public :: sources
     PetscReal, public :: productivity_index !! Productivity index
     PetscReal, public :: bottomhole_pressure
   contains
     procedure, public :: init => source_control_deliverability_init
     procedure, public :: destroy => source_control_deliverability_destroy
     procedure, public :: update => source_control_deliverability_update
  end type source_control_deliverability_type

  type, public, extends(source_control_type) :: source_control_separator_type
     !! Takes flow from a single source and outputs water and steam flow.
     private
     type(source_type), pointer :: source
     class(thermodynamics_type), pointer :: thermo
     PetscReal, public :: separator_pressure !! Separator pressure
     PetscReal :: water_enthalpy !! Enthalpy of water at separator pressure
     PetscReal :: steam_enthalpy !! Enthalpy of steam at separator pressure
     PetscReal, public :: water_flow_rate !! Output separated water flow rate
     PetscReal, public :: steam_flow_rate !! Output separated steam flow rate
   contains
     procedure, public :: init => source_control_separator_init
     procedure, public :: destroy => source_control_separator_destroy
     procedure, public :: update => source_control_separator_update
     procedure :: steam_fraction => source_control_separator_steam_fraction
  end type source_control_separator_type

  type, public, extends(source_control_type) :: source_control_limiter_type
     !! Limits flow in a particular phase (or total flow) through a source.
     private
     type(list_type), public :: sources
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

  subroutine source_control_table_init(self, data, interpolation_type, &
       averaging_type, sources)
    !! Initialises source_control_table object.

    class(source_control_table_type), intent(in out) :: self
    PetscReal, intent(in) :: data(:,:)
    PetscInt, intent(in) :: interpolation_type
    PetscInt, intent(in) :: averaging_type
    type(list_type), intent(in out) :: sources

    call self%sources%init()
    call self%table%init(data, interpolation_type, averaging_type)
    call self%sources%add(sources)

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

  subroutine source_control_deliverability_init(self, productivity_index, &
       bottomhole_pressure, sources)
    !! Initialises source_control_deliverability object.

    class(source_control_deliverability_type), intent(in out) :: self
    PetscReal, intent(in) :: productivity_index
    PetscReal, intent(in) :: bottomhole_pressure
    type(list_type), intent(in out) :: sources

    call self%sources%init()
    call self%sources%add(sources)
    self%productivity_index = productivity_index
    self%bottomhole_pressure = bottomhole_pressure

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
! Separator source control:
!------------------------------------------------------------------------

  subroutine source_control_separator_init(self)
    !! Initialises source_control_separator object.

    class(source_control_separator_type), intent(in out) :: self
    ! Locals:
    PetscReal :: saturation_temperature
    PetscReal :: params(2), water_props(2), steam_props(2)
    PetscErrorCode :: err

    call self%thermo%saturation%temperature(self%separator_pressure, &
         saturation_temperature, err)
    params = [self%separator_pressure, saturation_temperature]
    call self%thermo%water%properties(params, water_props, err)
    call self%thermo%steam%properties(params, steam_props, err)

    associate(water_density => water_props(1), &
         water_internal_energy => water_props(2), &
         steam_density => steam_props(1), &
         steam_internal_energy => steam_props(2))
      self%water_enthalpy = water_internal_energy + &
           self%separator_pressure / water_density
      self%steam_enthalpy = steam_internal_energy + &
           self%separator_pressure / steam_density
    end associate

  end subroutine source_control_separator_init

!------------------------------------------------------------------------

  subroutine source_control_separator_destroy(self)
    !! Destroys source_control_separator object.

    class(source_control_separator_type), intent(in out) :: self

    self%source => null()
    self%thermo => null()

  end subroutine source_control_separator_destroy

!------------------------------------------------------------------------

  PetscReal function source_control_separator_steam_fraction(self, &
       enthalpy) result(sf)
    !! Returns steam fraction for given enthalpy, based on the
    !! separator pressure.

    class(source_control_separator_type), intent(in) :: self
    PetscReal, intent(in) :: enthalpy

    if (enthalpy <= self%water_enthalpy) then
       sf = 0._dp
    else if (enthalpy <= self%steam_enthalpy) then
       sf = (enthalpy - self%water_enthalpy) / &
            (self%steam_enthalpy - self%water_enthalpy)
    else
       sf = 1._dp
    end if

  end function source_control_separator_steam_fraction

!------------------------------------------------------------------------

  subroutine source_control_separator_update(self, t, interval, &
       fluid_data, fluid_section)
    !! Update separated water and steam flow rates for
    !! source_control_separator_type.

    class(source_control_separator_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: fluid_data(:)
    PetscSection, intent(in) :: fluid_section
    ! Locals:
    PetscReal :: phase_flow_fractions(self%source%fluid%num_phases)
    PetscReal :: enthalpy, steam_fraction

    call self%source%update_fluid(fluid_data, fluid_section)
    phase_flow_fractions = self%source%fluid%phase_flow_fractions()
    enthalpy = self%source%fluid%specific_enthalpy(phase_flow_fractions)

    steam_fraction = self%steam_fraction(enthalpy)
    self%water_flow_rate = (1._dp - steam_fraction) * self%source%rate
    self%steam_flow_rate = steam_fraction * self%source%rate

  end subroutine source_control_separator_update

!------------------------------------------------------------------------
! Limiter source control:
!------------------------------------------------------------------------

  subroutine source_control_limiter_init(self, phase, limit, sources)
    !! Initialises source_control_limiter object.

    class(source_control_limiter_type), intent(in out) :: self
    PetscInt, intent(in) :: phase
    PetscReal, intent(in) :: limit
    type(list_type), intent(in out) :: sources

    call self%sources%init()
    call self%sources%add(sources)
    self%phase = phase
    self%limit = limit

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
      PetscReal, allocatable :: phase_flow_fractions(:)

      select type (source => node%data)
      type is (source_type)

         if (self%phase == 0) then ! Limit total flow:
            rate = abs(source%rate)
         else
            call source%update_fluid(fluid_data, fluid_section)
            allocate(phase_flow_fractions(source%fluid%num_phases))
            phase_flow_fractions = source%fluid%phase_flow_fractions()
            rate = phase_flow_fractions(self%phase) * &
                 abs(source%rate)
            deallocate(phase_flow_fractions)
         end if

         scale = self%rate_scale(rate)
         source%rate = source%rate * scale

      end select

      stopped = PETSC_FALSE

    end subroutine source_control_limiter_update_iterator

  end subroutine source_control_limiter_update

!------------------------------------------------------------------------

end module source_control_module
