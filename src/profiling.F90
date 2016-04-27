module profiling_module
  !! Module for profiling code via PETSc log events.

  implicit none
  private

#include <petsc/finclude/petscsys.h>

  PetscClassId, public ::  log_class
  PetscLogEvent, public :: simulation_init_event
  PetscLogEvent, public :: fluid_init_event
  PetscLogEvent, public :: fluid_properties_event, fluid_transitions_event
  PetscLogEvent, public :: lhs_fn_event, rhs_fn_event
  PetscLogEvent, public :: output_event
  PetscLogEvent, public :: thermo_region1_properties_event, &
       thermo_region2_properties_event, thermo_region3_properties_event
  PetscLogEvent, public :: thermo_viscosity_event
  PetscLogEvent, public :: thermo_sat_pressure_event, &
       thermo_sat_temperature_event
  PetscLogEvent, public :: eos_bulk_properties_event, &
       eos_phase_properties_event
  PetscLogEvent, public :: powertable_compute_event
  PetscLogEvent, public :: assign_pointers_event
  PetscLogEvent, public :: face_flux_event
  PetscLogEvent, public :: cell_inflows_event, sources_event

  public :: init_profiling

contains

!------------------------------------------------------------------------

  subroutine init_profiling()
    !! Initialize code profiling.

    ! Locals:
    PetscErrorCode :: ierr

    call PetscClassIdRegister("supermodel", log_class, ierr); CHKERRQ(ierr)

    ! Register log events:
    call PetscLogEventRegister("sim_init", log_class, simulation_init_event, ierr)
    call PetscLogEventRegister("fluid_init", log_class, fluid_init_event, ierr)
    call PetscLogEventRegister("fluid_props", log_class, &
         fluid_properties_event, ierr)
    call PetscLogEventRegister("fluid_trans", log_class, &
         fluid_transitions_event, ierr)
    call PetscLogEventRegister("lhs_function", log_class, lhs_fn_event, ierr)
    call PetscLogEventRegister("rhs_function", log_class, rhs_fn_event, ierr)
    call PetscLogEventRegister("output", log_class, output_event, ierr)
    call PetscLogEventRegister("thermo_props1", log_class, &
         thermo_region1_properties_event, ierr)
    call PetscLogEventRegister("thermo_props2", log_class, &
         thermo_region2_properties_event, ierr)
    call PetscLogEventRegister("thermo_visc", log_class, &
         thermo_viscosity_event, ierr)
    call PetscLogEventRegister("thermo_sat_p", log_class, &
         thermo_sat_pressure_event, ierr)
    call PetscLogEventRegister("thermo_sat_t", log_class, &
         thermo_sat_temperature_event, ierr)
    call PetscLogEventRegister("eos_bulk_props", log_class, &
         eos_bulk_properties_event, ierr)
    call PetscLogEventRegister("eos_phase_props", log_class, &
         eos_phase_properties_event, ierr)
    call PetscLogEventRegister("powertable_comp", log_class, &
         powertable_compute_event, ierr)
    call PetscLogEventRegister("assign_ptrs", log_class, &
         assign_pointers_event, ierr)
    call PetscLogEventRegister("face_flux", log_class, &
         face_flux_event, ierr)
    call PetscLogEventRegister("cell_inflows", log_class, &
         cell_inflows_event, ierr)
    call PetscLogEventRegister("sources", log_class, &
         sources_event, ierr)

  end subroutine init_profiling

!------------------------------------------------------------------------

end module profiling_module
