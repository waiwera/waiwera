module source_module
 !! Module for handling sinks and sources.

  use kinds_module
  use list_module

  implicit none
  private

#include <petsc/finclude/petsc.h90>

  type, public :: source_type
     !! Type for mass / energy source, applying specified values of
     !! generation to each equation in a particular cell at the
     !! current time.
     private
     PetscReal, allocatable, public :: flow(:) !! Generation flow rate for each equation
     PetscInt, public :: cell_natural_index !! Natural index of cell the source is in
     PetscInt, public :: cell_index !! Local index of cell the source is in
   contains
     private
     procedure, public :: init => source_init
     procedure, public :: destroy => source_destroy
     procedure, public :: total_flow => source_total_flow
  end type source_type

  type, public, abstract :: source_control_type
     !! Abstract type for mass / energy source control, controlling
     !! generation values over time in one or more sources.
     private
     type(list_type), public :: sources
   contains
     procedure(source_control_init_procedure), public, deferred :: init
     procedure(source_control_destroy_procedure), public, deferred :: destroy
     procedure(source_control_update_procedure), public, deferred :: update
  end type source_control_type

  abstract interface

     subroutine source_control_init_procedure(self)
       !! Initialises source control object
       import :: source_control_type
       class(source_control_type), intent(in out) :: self
     end subroutine source_control_init_procedure

     subroutine source_control_destroy_procedure(self)
       !! Destroys source control object
       import :: source_control_type
       class(source_control_type), intent(in out) :: self
     end subroutine source_control_destroy_procedure

     subroutine source_control_update_procedure(self, time)
       !! Updates sources at the specified time.
       import :: source_control_type
       class(source_control_type), intent(in out) :: self
       PetscReal, intent(in) :: time
     end subroutine source_control_update_procedure

  end interface

  public :: setup_sources

contains

!------------------------------------------------------------------------

  subroutine source_init(self, num_primary_variables, &
       cell_natural_index, cell_index)
    !! Initialises a source object.

    class(source_type), intent(in out) :: self
    PetscInt, intent(in) :: num_primary_variables
    PetscInt, intent(in) :: cell_natural_index
    PetscInt, intent(in) :: cell_index

    allocate(self%flow(num_primary_variables))
    self%flow = 0._dp
    self%cell_natural_index = cell_natural_index
    self%cell_index = cell_index

  end subroutine source_init

!------------------------------------------------------------------------

  subroutine source_total_flow(self, fluid, isothermal, q)
    !! Returns array with total flow through the source, including any
    !! energy production from mass production.

    use fluid_module, only: fluid_type

    class(source_type), intent(in) :: self
    type(fluid_type), intent(in) :: fluid
    PetscBool, intent(in) :: isothermal
    PetscReal, intent(out) :: q(:)
    ! Locals:
    PetscInt :: p, np, phases, c
    PetscReal :: flow_fractions(fluid%num_phases), hc

    q = self%flow
    if (.not. isothermal) then
       np = size(self%flow)
       associate (qenergy => q(np))

         phases = nint(fluid%phase_composition)
         flow_fractions = fluid%flow_fractions()

         do c = 1, fluid%num_components
            associate(f => self%flow(c))
              if (f < 0._dp) then
                 hc = 0._dp
                 do p = 1, fluid%num_phases
                    if (btest(phases, p - 1)) then
                       associate(phase => fluid%phase(p))
                         hc = hc + flow_fractions(p) * &
                              phase%specific_enthalpy * &
                              phase%mass_fraction(c)
                       end associate
                    end if
                 end do
                 qenergy = qenergy + f * hc
              end if
            end associate
         end do

       end associate
    end if

  end subroutine source_total_flow

!------------------------------------------------------------------------

  subroutine source_destroy(self)
    !! Destroys a source object.
    class(source_type), intent(in out) :: self

    deallocate(self%flow)

  end subroutine source_destroy

!------------------------------------------------------------------------

  subroutine setup_sources(json, dm, np, isothermal, sources_list, logfile)
    !! Sets up list of sinks and sources.

    use kinds_module
    use fson
    use fson_mpi_module
    use logfile_module
    use mesh_module, only: cell_order_label_name

    type(fson_value), pointer, intent(in) :: json !! JSON file object
    DM, intent(in) :: dm !! Mesh DM 
    PetscInt, intent(in) :: np !! Number of primary variables
    PetscBool, intent(in) :: isothermal !! Whether EOS is isothermal
    type(list_type), intent(in out) :: sources_list !! List of sources
    type(logfile_type), intent(in out), optional :: logfile !! Logfile for log output
    ! Locals:
    PetscErrorCode :: ierr
    PetscInt :: c, isrc, i, component
    type(fson_value), pointer :: sources, src
    PetscInt :: num_sources, cell, num_cells, ghost
    type(source_type), pointer :: s
    PetscInt, pointer :: cells(:)
    PetscReal :: q, enthalpy
    PetscBool :: mass_inject
    character(:), allocatable :: name
    IS :: cell_IS
    DMLabel :: ghost_label
    character(len=64) :: srcstr
    character(len=12) :: istr
    PetscInt, parameter ::  default_component = 1
    PetscReal, parameter :: default_enthalpy = 83.9e3

    call sources_list%init(delete_deallocates = PETSC_TRUE)

    if (fson_has_mpi(json, "source")) then

       call DMGetLabel(dm, "ghost", ghost_label, ierr); CHKERRQ(ierr)

       call fson_get_mpi(json, "source", sources)
       num_sources = fson_value_count_mpi(sources, ".")

       do isrc = 1, num_sources

          write(istr, '(i0)') isrc - 1
          srcstr = 'source[' // trim(istr) // '].'
          src => fson_value_get_mpi(sources, isrc)
          if (allocated(name)) deallocate(name)
          if (fson_has_mpi(src, "name")) then
             call fson_get_mpi(src, "name", val = name)
          end if
          call fson_get_mpi(src, "cell", val = cell)
          call fson_get_mpi(src, "component", default_component, &
               component, logfile, trim(srcstr) // "component")
          call fson_get_mpi(src, "value", val = q)
          mass_inject = ((.not.(isothermal)) .and. (component < np) &
               .and. (q > 0._dp))
          if (mass_inject) then
             call fson_get_mpi(src, "enthalpy", default_enthalpy, &
                  enthalpy, logfile, trim(srcstr) // "enthalpy")
          else
             enthalpy = 0._dp
          end if

          call DMGetStratumIS(dm, cell_order_label_name, &
               cell, cell_IS, ierr); CHKERRQ(ierr)

          if (cell_IS /= 0) then

             call ISGetIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)
             num_cells = size(cells)
             do i = 1, num_cells
                c = cells(i)
                call DMLabelGetValue(ghost_label, c, ghost, ierr)
                CHKERRQ(ierr)
                if (ghost < 0) then
                   allocate(s)
                   call s%init(np, cell, c)
                   s%flow(component) = q
                   if (mass_inject) then
                      s%flow(np) = s%flow(np) + enthalpy * q
                   end if
                   if (allocated(name)) then
                      call sources_list%append(s, name)
                   else
                      call sources_list%append(s)
                   end if
                end if
             end do
             call ISRestoreIndicesF90(cell_IS, cells, ierr); CHKERRQ(ierr)

          end if
          call ISDestroy(cell_IS, ierr); CHKERRQ(ierr)

       end do

    else
       call logfile%write(LOG_LEVEL_INFO, "input", "no_sources")
    end if

  end subroutine setup_sources

!------------------------------------------------------------------------

end module source_module
