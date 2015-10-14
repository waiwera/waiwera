module cell_module
  !! Defines types for accessing local quantities defined on a cell- geometry, rock
  !! and fluid properties.
  !! The components of these types all point to values in arrays obtained from
  !! global parallel vectors.

  use kinds_module
  use rock_module
  use fluid_module

  implicit none
  private

#include <petsc/finclude/petscdef.h>

  type cell_type
     !! Type for accessing local cell geometry, rock and
     !! fluid properties.
     private
     PetscReal, pointer, public :: volume       !! cell volume
     PetscReal, pointer, public :: centroid(:)  !! cell centroid
     type(rock_type), public :: rock   !! rock properties
     type(fluid_type), public :: fluid !! fluid properties
   contains
     private
     procedure, public :: init => cell_init
     procedure, public :: assign => cell_assign
     procedure, public :: destroy => cell_destroy
     procedure, public :: dof => cell_dof
     procedure, public :: balance => cell_balance
  end type cell_type

  PetscInt, parameter :: num_cell_variables = 2
  PetscInt, parameter, public :: &
       cell_variable_num_components(num_cell_variables) = &
       [1, 3]

  public :: cell_type

contains

!------------------------------------------------------------------------

  subroutine cell_init(self, num_components, num_phases)
    !! Initialises a cell.

    class(cell_type), intent(in out) :: self
    PetscInt, intent(in) :: num_components !! Number of fluid components
    PetscInt, intent(in) :: num_phases     !! Number of fluid phases

    call self%fluid%init(num_components, num_phases)

  end subroutine cell_init

!------------------------------------------------------------------------

  subroutine cell_assign(self, geom_data, geom_offset, &
       rock_data, rock_offset, fluid_data, fluid_offset)
    !! Assigns pointers in a cell to values from the specified data arrays,
    !! starting from the given offset.

    class(cell_type), intent(in out) :: self
    PetscReal, target, intent(in), optional :: geom_data(:)  !! array with geometry data
    PetscInt, intent(in), optional  :: geom_offset  !! geometry array offset for this cell
    PetscReal, target, intent(in), optional :: rock_data(:)  !! array with rock data
    PetscInt, intent(in), optional  :: rock_offset  !! rock array offset for this cell
    PetscReal, target, intent(in), optional :: fluid_data(:)  !! array with fluid data
    PetscInt, intent(in), optional  :: fluid_offset  !! fluid array offset for this cell

    if ((present(geom_data)) .and. ((present(geom_offset)))) then
       self%centroid => geom_data(geom_offset: geom_offset + 2)
       self%volume => geom_data(geom_offset + 3)
    end if

    if ((present(rock_data)) .and. ((present(rock_offset)))) then
       call self%rock%assign(rock_data, rock_offset)
    end if
    if ((present(fluid_data)) .and. ((present(fluid_offset)))) then
       call self%fluid%assign(fluid_data, fluid_offset)
    end if

  end subroutine cell_assign

!------------------------------------------------------------------------

  subroutine cell_destroy(self)
    !! Destroys a cell.

    class(cell_type), intent(in out) :: self

    nullify(self%volume)
    nullify(self%centroid)
    call self%fluid%destroy()
    call self%rock%destroy()
    
  end subroutine cell_destroy

!------------------------------------------------------------------------

  PetscInt function cell_dof(self)
    !! Returns number of degrees of freedom in a cell object.

    class(cell_type), intent(in) :: self

    cell_dof = sum(cell_variable_num_components)

  end function cell_dof

!------------------------------------------------------------------------

  function cell_balance(self, eos) result(balance)
    !! Returns array containing mass balance (per unit volume) for each
    !! mass component in the cell, and energy balance for non-isothermal
    !! simulations.

    use eos_module, only: eos_type

    class(cell_type), intent(in) :: self
    class(eos_type), intent(in) :: eos
    PetscReal :: balance(eos%num_primary_variables)
    ! Locals:
    PetscInt :: nc
    PetscReal :: er, ef

    nc = self%fluid%num_components

    ! Mass balances:
    balance(1: nc) = self%rock%porosity * self%fluid%component_density(eos)

    if (.not. eos%isothermal) then
       ! Energy balance:
       er = self%rock%energy(self%fluid%temperature)
       ef = self%fluid%energy(eos)
       balance(eos%num_primary_variables) = self%rock%porosity * ef + &
            (1._dp - self%rock%porosity) * er
    end if

  end function cell_balance

!------------------------------------------------------------------------
  
end module cell_module
