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

  PetscReal, parameter, public :: default_source_control_separator_pressure = 0.55e6_dp
  PetscInt, parameter, public :: max_limiter_type_length = 5
  character(max_limiter_type_length), parameter, public :: &
       default_source_control_limiter_type_str = "total"
  PetscInt, parameter, public :: SRC_CONTROL_LIMITER_TYPE_TOTAL = 1, &
       SRC_CONTROL_LIMITER_TYPE_WATER = 2, SRC_CONTROL_LIMITER_TYPE_STEAM = 3
  PetscReal, parameter, public :: default_source_control_limiter_limit = 1._dp
  PetscReal, parameter, public :: default_deliverability_productivity_index = 1.e-11_dp
  PetscReal, parameter, public :: default_deliverability_reference_pressure = 1.e5_dp
  PetscInt, parameter, public :: SRC_DIRECTION_PRODUCTION = 1, &
       SRC_DIRECTION_INJECTION = 2, SRC_DIRECTION_BOTH = 3
  PetscInt, parameter, public :: default_source_direction = SRC_DIRECTION_BOTH

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
     type(interpolation_table_type), public :: productivity_index !! Productivity index vs. time
     PetscReal, public :: reference_pressure
     PetscInt, public :: direction
   contains
     procedure, public :: init => source_control_deliverability_init
     procedure, public :: destroy => source_control_deliverability_destroy
     procedure, public :: update => source_control_deliverability_update
     procedure, public :: set_reference_pressure_initial => &
          source_control_deliverability_set_reference_pressure_initial
     procedure, public :: calculate_PI_from_rate => &
          source_control_deliverability_calculate_PI_from_rate
     procedure, public :: calculate_PI_from_recharge => &
          source_control_deliverability_calculate_PI_from_recharge
  end type source_control_deliverability_type

  type, public, extends(source_control_type) :: source_control_separator_type
     !! Takes flow from a single source and outputs water and steam flow.
     private
     type(source_type), pointer :: source
     class(thermodynamics_type), pointer :: thermo
     PetscReal, public :: separator_pressure !! Separator pressure
     PetscReal :: water_enthalpy !! Enthalpy of water at separator pressure
     PetscReal :: steam_enthalpy !! Enthalpy of steam at separator pressure
     PetscReal, public :: steam_fraction  !! Output steam fraction
     PetscReal, public :: water_flow_rate !! Output separated water flow rate
     PetscReal, public :: steam_flow_rate !! Output separated steam flow rate
   contains
     procedure, public :: init => source_control_separator_init
     procedure, public :: destroy => source_control_separator_destroy
     procedure, public :: update => source_control_separator_update
     procedure :: update_steam_fraction => source_control_separator_update_steam_fraction
  end type source_control_separator_type

  type, public, extends(source_control_type) :: source_control_limiter_type
     !! Limits total flow, or separated steam or water flow (output
     !! from a separator) through a source.
     private
     type(list_type), public :: sources
     PetscReal, pointer, public :: rate !! Pointer to flow rate variable
     PetscInt, public :: type !! Type of limiter (water, steam or total)
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
          local_fluid_data, local_fluid_section)
       !! Updates sources at the specified time.
       import :: source_control_type
       class(source_control_type), intent(in out) :: self
       PetscReal, intent(in) :: t, interval(2)
       PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
       PetscSection, intent(in) :: local_fluid_section
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
       local_fluid_data, local_fluid_section)
    !! Update flow rate for source_control_rate_table_type.

    class(source_control_rate_table_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
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
       local_fluid_data, local_fluid_section)
    !! Update injection enthalpy for source_control_enthalpy_table_type.

    class(source_control_enthalpy_table_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
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

  subroutine source_control_deliverability_init(self, productivity_data, &
       interpolation_type, averaging_type, reference_pressure, direction, &
       sources)
    !! Initialises source_control_deliverability object.

    class(source_control_deliverability_type), intent(in out) :: self
    PetscReal, intent(in) :: productivity_data(:,:)
    PetscInt, intent(in) :: interpolation_type, averaging_type
    PetscReal, intent(in) :: reference_pressure
    PetscInt, intent(in) :: direction
    type(list_type), intent(in out) :: sources

    call self%sources%init()
    call self%sources%add(sources)
    call self%productivity_index%init(productivity_data, &
         interpolation_type, averaging_type)
    self%reference_pressure = reference_pressure
    self%direction = direction

  end subroutine source_control_deliverability_init

!------------------------------------------------------------------------

  subroutine source_control_deliverability_destroy(self)
    !! Destroys source_control_deliverability object.

    class(source_control_deliverability_type), intent(in out) :: self

    call self%sources%destroy()

  end subroutine source_control_deliverability_destroy

!------------------------------------------------------------------------

  subroutine source_control_deliverability_update(self, t, interval, &
       local_fluid_data, local_fluid_section)
    !! Update flow rate for source_control_deliverability_type.

    class(source_control_deliverability_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section

    call self%sources%traverse(source_control_deliverability_update_iterator)

  contains

    subroutine source_control_deliverability_update_iterator(node, stopped)
      !! Applies deliverability control at a list node.
      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped
      ! Locals:
      PetscInt :: p, phases
      PetscReal :: pressure_difference, productivity_index
      PetscBool :: flowing

      select type (source => node%data)
      type is (source_type)

         call source%update_fluid(local_fluid_data, local_fluid_section)
         pressure_difference = source%fluid%pressure - self%reference_pressure

         source%rate = 0._dp

         select case (self%direction)
         case (SRC_DIRECTION_PRODUCTION)
            flowing = (pressure_difference > 0._dp)
         case (SRC_DIRECTION_INJECTION)
            flowing = (pressure_difference < 0._dp)
         case default
            flowing = PETSC_TRUE
         end select

         if (flowing) then
            productivity_index = self%productivity_index%average(interval)
            phases = nint(source%fluid%phase_composition)
            do p = 1, source%fluid%num_phases
               if (btest(phases, p - 1)) then
                  source%rate = source%rate - productivity_index * &
                       source%fluid%phase(p)%mobility() * pressure_difference
               end if
            end do
         end if

      end select

      stopped = PETSC_FALSE

    end subroutine source_control_deliverability_update_iterator

  end subroutine source_control_deliverability_update

!------------------------------------------------------------------------

  subroutine source_control_deliverability_set_reference_pressure_initial(self, &
       global_fluid_data, global_fluid_section, fluid_range_start)
    !! Sets reference pressure for deliverability control to be the
    !! initial fluid pressure.

    use dm_utils_module, only: global_section_offset

    class(source_control_deliverability_type), intent(in out) :: self
    PetscReal, pointer, contiguous, intent(in) :: global_fluid_data(:)
    PetscSection, intent(in) :: global_fluid_section
    PetscInt, intent(in) :: fluid_range_start
    ! Locals:
    type(list_node_type), pointer :: node
    PetscInt :: c, fluid_offset
    PetscErrorCode :: ierr

    if (self%sources%count == 1) then
       node => self%sources%head
       select type (source => node%data)
       type is (source_type)

          c = source%cell_index
          call global_section_offset(global_fluid_section, c, &
               fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)

          call source%fluid%assign(global_fluid_data, fluid_offset)

          self%reference_pressure = source%fluid%pressure

       end select
    end if

  end subroutine source_control_deliverability_set_reference_pressure_initial

!------------------------------------------------------------------------

  subroutine source_control_deliverability_calculate_PI_from_rate(&
       self, initial_rate, global_fluid_data, global_fluid_section, &
       fluid_range_start)
    !! Calculates productivity index for deliverability control, from
    !! specified initial flow rate. This only works if the control has
    !! exactly one source, otherwise the correct productivity index
    !! would not be well-defined. If the productivity index can't be
    !! calculated, it is left at its initial value.

    use dm_utils_module, only: global_section_offset

    class(source_control_deliverability_type), intent(in out) :: self
    PetscReal, intent(in) :: initial_rate
    PetscReal, pointer, contiguous, intent(in) :: global_fluid_data(:)
    PetscSection, intent(in) :: global_fluid_section
    PetscInt, intent(in) :: fluid_range_start
    ! Locals:
    type(list_node_type), pointer :: node
    PetscInt :: c, fluid_offset
    PetscReal, allocatable :: phase_mobilities(:)
    PetscReal :: factor
    PetscErrorCode :: ierr
    PetscReal, parameter :: tol = 1.e-9_dp

    if (self%sources%count == 1) then
       node => self%sources%head
       select type (source => node%data)
       type is (source_type)

          c = source%cell_index
          call global_section_offset(global_fluid_section, c, &
               fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)

          call source%fluid%assign(global_fluid_data, fluid_offset)
          allocate(phase_mobilities(source%fluid%num_phases))
          phase_mobilities = source%fluid%phase_mobilities()

          factor = sum(phase_mobilities) * &
               (source%fluid%pressure - self%reference_pressure)

          if (abs(factor) > tol) then
             self%productivity_index%val(1) = abs(initial_rate) / factor
          end if

          deallocate(phase_mobilities)

       end select
    end if

  end subroutine source_control_deliverability_calculate_PI_from_rate

!------------------------------------------------------------------------

  subroutine source_control_deliverability_calculate_PI_from_recharge(&
       self, recharge_coefficient, global_fluid_data, global_fluid_section, &
       fluid_range_start)
    !! Calculates productivity index for deliverability control, from
    !! specified recharge coefficient. This only works if the control has
    !! exactly one source, otherwise the correct productivity index
    !! would not be well-defined. If the productivity index can't be
    !! calculated, it is left at its initial value.

    use dm_utils_module, only: global_section_offset

    class(source_control_deliverability_type), intent(in out) :: self
    PetscReal, intent(in) :: recharge_coefficient
    PetscReal, pointer, contiguous, intent(in) :: global_fluid_data(:)
    PetscSection, intent(in) :: global_fluid_section
    PetscInt, intent(in) :: fluid_range_start
    ! Locals:
    type(list_node_type), pointer :: node
    PetscInt :: c, fluid_offset
    PetscReal, allocatable :: liquid_mobility
    PetscErrorCode :: ierr
    PetscReal, parameter :: tol = 1.e-9_dp

    if (self%sources%count == 1) then
       node => self%sources%head
       select type (source => node%data)
       type is (source_type)

          c = source%cell_index
          call global_section_offset(global_fluid_section, c, &
               fluid_range_start, fluid_offset, ierr); CHKERRQ(ierr)

          call source%fluid%assign(global_fluid_data, fluid_offset)

          liquid_mobility = source%fluid%phase(1)%mobility()

          if (liquid_mobility > tol) then
             self%productivity_index%val(1) = -recharge_coefficient / liquid_mobility
          end if

       end select
    end if

  end subroutine source_control_deliverability_calculate_PI_from_recharge

!------------------------------------------------------------------------
! Separator source control:
!------------------------------------------------------------------------

  subroutine source_control_separator_init(self, source, &
       thermo, separator_pressure)
    !! Initialises source_control_separator object.

    class(source_control_separator_type), intent(in out) :: self
    type(source_type), target, intent(in) :: source
    class(thermodynamics_type), target, intent(in) :: thermo
    PetscReal, intent(in) :: separator_pressure
    ! Locals:
    PetscReal :: saturation_temperature
    PetscReal :: params(2), water_props(2), steam_props(2)
    PetscErrorCode :: err

    self%source => source
    self%thermo => thermo
    self%separator_pressure = separator_pressure

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

  subroutine source_control_separator_update_steam_fraction(self, &
       enthalpy)
    !! Updates steam fraction for given enthalpy, based on the
    !! separator pressure.

    class(source_control_separator_type), intent(in out) :: self
    PetscReal, intent(in) :: enthalpy

    if (enthalpy <= self%water_enthalpy) then
       self%steam_fraction = 0._dp
    else if (enthalpy <= self%steam_enthalpy) then
       self%steam_fraction = (enthalpy - self%water_enthalpy) / &
            (self%steam_enthalpy - self%water_enthalpy)
    else
       self%steam_fraction = 1._dp
    end if

  end subroutine source_control_separator_update_steam_fraction

!------------------------------------------------------------------------

  subroutine source_control_separator_update(self, t, interval, &
       local_fluid_data, local_fluid_section)
    !! Update separated water and steam flow rates for
    !! source_control_separator_type.

    class(source_control_separator_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    ! Locals:
    PetscReal :: phase_flow_fractions(self%source%fluid%num_phases)
    PetscReal :: enthalpy

    associate(np => size(self%source%flow))

      if ((self%source%rate < 0._dp) .and. &
           (self%source%production_component < np)) then

         call self%source%update_fluid(local_fluid_data, local_fluid_section)
         phase_flow_fractions = self%source%fluid%phase_flow_fractions()
         enthalpy = self%source%fluid%specific_enthalpy(phase_flow_fractions)

         call self%update_steam_fraction(enthalpy)
         self%water_flow_rate = (1._dp - self%steam_fraction) * self%source%rate
         self%steam_flow_rate = self%steam_fraction * self%source%rate

      else

         self%water_flow_rate = 0._dp
         self%steam_flow_rate = 0._dp

      end if

    end associate

  end subroutine source_control_separator_update

!------------------------------------------------------------------------
! Limiter source control:
!------------------------------------------------------------------------

  subroutine source_control_limiter_init(self, limiter_type, rate, &
       limit, sources)
    !! Initialises source_control_limiter object.

    class(source_control_limiter_type), intent(in out) :: self
    PetscInt, intent(in) :: limiter_type
    PetscReal, target, intent(in) :: rate
    PetscReal, intent(in) :: limit
    type(list_type), intent(in out) :: sources

    call self%sources%init()
    call self%sources%add(sources)
    self%type = limiter_type
    self%rate => rate
    self%limit = limit

  end subroutine source_control_limiter_init

!------------------------------------------------------------------------

  subroutine source_control_limiter_destroy(self)
    !! Destroys source_control_limiter_type object.

    class(source_control_limiter_type), intent(in out) :: self

    self%rate => null()
    call self%sources%destroy()

  end subroutine source_control_limiter_destroy

!------------------------------------------------------------------------

  PetscReal function source_control_limiter_rate_scale(self) &
       result(scale)
    !! Returns factor by which flow should be scaled to avoid
    !! exceededing limit.

    class(source_control_limiter_type), intent(in) :: self
    ! Locals:
    PetscReal :: abs_rate
    PetscReal, parameter :: small = 1.e-6_dp

    abs_rate = abs(self%rate)

    if ((abs_rate > self%limit) .and. (abs_rate > small)) then
       scale = self%limit / abs_rate
    else
       scale = 1._dp
    end if

  end function source_control_limiter_rate_scale

!------------------------------------------------------------------------

  subroutine source_control_limiter_update(self, t, interval, &
       local_fluid_data, local_fluid_section)
    !! Update flow rate for source_control_limiter_type.

    class(source_control_limiter_type), intent(in out) :: self
    PetscReal, intent(in) :: t, interval(2)
    PetscReal, pointer, contiguous, intent(in) :: local_fluid_data(:)
    PetscSection, intent(in) :: local_fluid_section
    PetscReal :: scale

    scale = self%rate_scale()
    call self%sources%traverse(source_control_limiter_update_iterator)

  contains

    subroutine source_control_limiter_update_iterator(node, stopped)
      !! Limits source rate at a list node.

      type(list_node_type), pointer, intent(in out)  :: node
      PetscBool, intent(out) :: stopped

      select type (source => node%data)
      type is (source_type)
         source%rate = source%rate * scale
      end select

      stopped = PETSC_FALSE

    end subroutine source_control_limiter_update_iterator

  end subroutine source_control_limiter_update

!------------------------------------------------------------------------

end module source_control_module
