module simulation_module
  !! Module for high-level representation of a simulation.

  use mpi_module
  use kinds_module
  use mesh_module
  use timestepping_module
  use thermodynamics_module
  use IAPWS_module
  use IFC67_module
  use eos_w_module
  use eos_module
  use fson
  use fson_mpi_module

  implicit none

  private

#include <petsc-finclude/petscsys.h>
#include <petsc-finclude/petscdef.h>
#include <petsc-finclude/petscvec.h>

  integer, parameter, public :: max_filename_length = 200
  integer, parameter, public :: max_title_length = 120

  type, public :: simulation_type
     !! Simulation type.
     private
     type(mesh_type) :: mesh
     Vec :: initial
     type(timestepper_type) :: timestepper
     class(thermodynamics_type), pointer :: thermo
     class(eos_type), pointer :: eos
     character(max_filename_length) :: input_filename
     character(max_title_length), public :: title
   contains
     private
     procedure :: read_title => simulation_read_title
     procedure :: read_thermodynamics => simulation_read_thermodynamics
     procedure :: read_eos => simulation_read_eos
     procedure :: read_initial => simulation_read_initial
     procedure :: read_timestepping => simulation_read_timestepping
     procedure, public :: init => simulation_init
     procedure, public :: run => simulation_run
     procedure, public :: destroy => simulation_destroy
  end type simulation_type

contains

!------------------------------------------------------------------------

  subroutine simulation_read_title(self, json)
    !! Reads simulation title from JSON input file.
    !! If not present, a default value is assigned.

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    character(len = max_title_length), parameter :: default_title = ""

    call fson_get_mpi(json, "title", default_title, self%title)

  end subroutine simulation_read_title

!------------------------------------------------------------------------

  subroutine simulation_read_thermodynamics(self, json)
    !! Reads simulation thermodynamic formulation from JSON input file.
    !! If not present, a default value is assigned.

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscErrorCode :: ierr
    integer, parameter :: max_thermo_ID_length = 8
    character(max_thermo_ID_length) :: thermo_ID
    character(max_thermo_ID_length), parameter :: &
         default_thermo_ID = "IAPWS"

    call fson_get_mpi(json, "thermodynamics", default_thermo_ID, thermo_ID)

    select case (thermo_ID)
    case ("IFC67")
       self%thermo => IFC67
    case default
       self%thermo => IAPWS
    end select

    call self%thermo%init()

  end subroutine simulation_read_thermodynamics

!------------------------------------------------------------------------

  subroutine simulation_read_eos(self, json)
    !! Reads simulation equation of state from JSON input file.
    !! If not present, a default value is assigned.

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    PetscErrorCode :: ierr
    integer, parameter :: max_eos_ID_length = 12
    character(max_eos_ID_length) :: eos_ID
    character(max_eos_ID_length), parameter :: &
         default_eos_ID = "W"

    call fson_get_mpi(json, "eos", default_eos_ID, eos_ID)

    select case (eos_ID)
    case ("EW")
       self%eos => eos_w  ! change to eos_ew when it's ready
    case default
       self%eos => eos_w
    end select

    call self%eos%init(self%thermo)

  end subroutine simulation_read_eos

!------------------------------------------------------------------------

  subroutine simulation_read_initial(self, json)
    !! Reads initial conditions from JSON input. These may be specified
    !! as a constant value or as an array. The array may be of length
    !! equal to the number of primary variables (to apply the same initial
    !! conditions in each cell) or as a complete specification of 
    !! initial conditions for all cells.

    use fson_value_m

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    type(fson_value), pointer :: initial
    PetscErrorCode :: ierr
    real(dp) :: const_initial_value
    integer :: int_const_initial_value
    integer, allocatable :: indices(:)
    real(dp), allocatable :: initial_input(:), initial_data(:)
    integer :: i, np, count
    logical :: const
    real(dp), parameter :: default_initial_value = 0.0_dp

    call DMCreateGlobalVector(self%mesh%dm, self%initial, ierr); CHKERRQ(ierr)
    call VecGetSize(self%initial, count, ierr); CHKERRQ(ierr)

    if (mpi%rank == mpi%input_rank) then
       
       const = .true.
       call fson_get(json, "initial", initial)
       if (associated(initial)) then

           ! if (initial%value_type == TYPE_STRING) then
              ! interpret as filename and read from separate file,
              ! then continue as below
           ! end if

          if (initial%value_type == TYPE_REAL) then

             call fson_get(initial, ".", const_initial_value)

          else if (initial%value_type == TYPE_INTEGER) then

             call fson_get(initial, ".", int_const_initial_value)
             const_initial_value = real(int_const_initial_value)

          else if (initial%value_type == TYPE_ARRAY) then

             const = .false.
             call fson_get(initial, ".", initial_input)
             np = size(initial_input)
             allocate(initial_data(count))
             if (np == count) then
                initial_data = initial_input
             else
                ! repeat input over array:
                do i = 1, np
                   initial_data(i:count:np) = initial_input(i)
                end do
             end if
             deallocate(initial_input)

          end if
       else
          const_initial_value = default_initial_value
       end if
    end if

    call MPI_bcast(const, 1, MPI_LOGICAL, mpi%input_rank, mpi%comm, ierr)
    
    if (const) then

       call MPI_bcast(const_initial_value, 1, MPI_DOUBLE_PRECISION, &
            mpi%input_rank, mpi%comm, ierr)
       call VecSet(self%initial, const_initial_value, ierr); CHKERRQ(ierr)

    else

       if (mpi%rank /= mpi%input_rank) then
          allocate(initial_data(count))
       end if
       call MPI_bcast(initial_data, count, MPI_DOUBLE_PRECISION, &
            mpi%input_rank, mpi%comm, ierr)
       allocate(indices(count))
       do i = 1, count
          indices(i) = i-1
       end do
       call VecSetValues(self%initial, count, indices, &
            initial_data, INSERT_VALUES, ierr); CHKERRQ(ierr)
       call VecAssemblyBegin(self%initial, ierr); CHKERRQ(ierr)
       call VecAssemblyEnd(self%initial, ierr); CHKERRQ(ierr)
       deallocate(indices, initial_data)

    end if

    call VecView(self%initial, PETSC_VIEWER_STDOUT_WORLD, ierr); CHKERRQ(ierr)

  end subroutine simulation_read_initial

!------------------------------------------------------------------------

  subroutine simulation_read_timestepping(self, json)
    !! Reads time stepping data from JSON input.

    class(simulation_type), intent(in out) :: self
    type(fson_value), pointer, intent(in) :: json
    ! Locals:
    type(fson_value), pointer :: time
    PetscReal, parameter :: default_initial_time = 0.0_dp
    PetscReal :: initial_time

    if (mpi%rank == mpi%input_rank) then
        call fson_get(json, "time", time)
        if (associated(time)) then
           call fson_get_default(time, "start", default_initial_time, &
                initial_time)
        end if
       
    end if

    ! call self%timestepper%init(method, self%mesh%dm, simulation_balance, &
    !      simulation_inflow, initial_time, self%initial, 

  end subroutine simulation_read_timestepping

!------------------------------------------------------------------------

  subroutine simulation_init(self, filename)
    !! Initializes a simulation using data from the input file with 
    !! specified name.

    class(simulation_type), intent(in out) :: self
    character(*), intent(in) :: filename !! Input file name
    ! Locals:
    type(fson_value), pointer :: json
    PetscErrorCode :: ierr
    PetscInt :: dof

    self%input_filename = filename

    if (mpi%rank == mpi%input_rank) then
       json => fson_parse(filename)
    end if

    call self%read_title(json)

    call self%read_thermodynamics(json)
    call self%read_eos(json)

    dof = self%eos%num_primary
    call self%mesh%init(json, dof)

    call self%read_initial(json)

    call self%read_timestepping(json)

    if (mpi%rank == mpi%input_rank) then
       call fson_destroy(json)
    end if

  end subroutine simulation_init

!------------------------------------------------------------------------

  subroutine simulation_run(self)
    !! Runs the simulation.

    class(simulation_type), intent(in out) :: self

    ! call self%timestepper%run()

    ! maybe some final output?

  end subroutine simulation_run

!------------------------------------------------------------------------

  subroutine simulation_destroy(self)
    !! Destroys the simulation.

    class(simulation_type), intent(in out) :: self
    ! Locals:
    PetscErrorCode :: ierr

    call VecDestroy(self%initial, ierr); CHKERRQ(ierr)
    call self%mesh%destroy()
    call self%thermo%destroy()
    ! call self%timestepper%destroy()

  end subroutine simulation_destroy

!------------------------------------------------------------------------

  subroutine simulation_balance(t, primary, balance)
    !! Computes mass and energy balance for each cell, for the given
    !! primary thermodynamic variables and time.
    PetscReal, intent(in) :: t !! time (s)
    Vec, intent(in) :: primary
    Vec, intent(out) :: balance
    ! Locals:
    PetscErrorCode :: ierr

    ! Dummy identity function:
    call VecCopy(primary, balance, ierr); CHKERRQ(ierr)

  end subroutine simulation_balance

!------------------------------------------------------------------------

  subroutine simulation_inflow(t, primary, inflow)
    !! Computes net inflow into each cell, from flows through faces and
    !! source terms, for the given primary thermodynamic variables and
    !! time.

    PetscReal, intent(in) :: t
    Vec, intent(in) :: primary
    Vec, intent(out) :: inflow
    ! Locals:
    PetscErrorCode :: ierr

    ! Dummy zero inflows:
    call VecSet(inflow, 0.0_dp, ierr); CHKERRQ(ierr)

  end subroutine simulation_inflow

!------------------------------------------------------------------------

end module simulation_module
