.. index:: output; HDF5, HDF5
.. _setup_output:

*****************
Simulation output
*****************

Waiwera outputs the simulation results (not log messages, which are written to the :ref:`setup_logfile`) to an output file in the `HDF5 <https://portal.hdfgroup.org/display/HDF5/HDF5>`_ file format. HDF5 is a binary file format and data model designed for efficiently storing and handling large, complex datasets. A variety of software `tools <https://portal.hdfgroup.org/display/HDF5/Libraries+and+Tools+Reference>`_ are available for managing, viewing and analysing HDF5 files.

Waiwera simulation results consist mainly of:

* selected fluid properties (e.g. pressures, temperatures) in each cell
* selected flux properties (e.g. mass fluxes) through each face
* selected source properties (e.g. flow rates, enthalpies) at each source (see :ref:`source_terms`)
* selected mesh geometry properties (e.g. cell centroids and volumes, and face areas)

which are written to the HDF5 file at specified times.

Simulation output can be controlled via the **"output"** value in the Waiwera JSON input file. This can be specified as a boolean value, so that it is possible (though unusual) to disable simulation output by setting it to ``false``. In most cases, however, it is specified as an object.

This object has a **"filename"** string value for specifying the filename of the simulation output. If this is not specified, then a default filename will be used, formed from the Waiwera JSON input filename but with the extension changed to ".h5". The other values in the "output" object control the times at which results are output, and which fluid and source properties are written.

.. note::
   **JSON object**: simulation output

   **JSON path**: output

   +-------------+--------------+----------------------+-------------------------------+
   |**name**     |**type**      |**default**           |**value**                      |
   +-------------+--------------+----------------------+-------------------------------+
   |"filename"   |string        |input filename with   |simulation output filename     |
   |             |              |extension changed from|                               |
   |             |              |".json" to ".h5"      |                               |
   +-------------+--------------+----------------------+-------------------------------+
   |"initial"    |boolean       |true                  |whether initial conditions     |
   |             |              |                      |are included in output         |
   |             |              |                      |                               |
   |             |              |                      |                               |
   +-------------+--------------+----------------------+-------------------------------+
   |"final"      |boolean       |true                  |whether final results are      |
   |             |              |                      |always included in output      |
   +-------------+--------------+----------------------+-------------------------------+
   |"frequency"  |integer       |1                     |number of time steps between   |
   |             |              |                      |:ref:`regular_output`          |
   |             |              |                      |                               |
   +-------------+--------------+----------------------+-------------------------------+
   |"checkpoint" |object        |{}                    |parameters for                 |
   |             |              |                      |:ref:`checkpoint_output`       |
   |             |              |                      |                               |
   +-------------+--------------+----------------------+-------------------------------+
   |"fields"     |object        |depends on EOS        |fields to output (see          |
   |             |              |                      |:ref:`output_fields`)          |
   |             |              |                      |                               |
   +-------------+--------------+----------------------+-------------------------------+
   |"jacobian"   |boolean |     |false                 |whether to output Jacobian     |
   |             |object        |                      |matrix                         |
   |             |              |                      |                               |
   +-------------+--------------+----------------------+-------------------------------+

.. index:: output; regular
.. _regular_output:

Regular output
==============

By default, results are output at every time step. However, if this amount of output is not needed, results may be output less frequently using the **"frequency"** integer value. This sets the number of time steps between regular output, so for example:

.. code-block:: json

   {"output": {"frequency": 5}}

gives output only every fifth time step.

Setting the "frequency" value to zero disables regular output. This is usually desirable only in conjunction with "checkpoint" output (see :ref:`checkpoint_output`), or if only the final simulation results are required (see :ref:`initial_and_final_output`), e.g. for :ref:`steady_state` using adaptive time-stepping.

.. index:: output; initial, output; final
.. _initial_and_final_output:

Initial and final output
========================

The **"initial"** and **"final"** boolean values control whether results are output at the start and end of the simulation respectively. Both are set to ``true`` by default, so that the initial conditions are written to the output file, as well as the results after the final time step (regardless of whether this would have been written anyway).

For example:

.. code-block:: json

   {"output": {"frequency": 0, "initial": false, "final": true}}

disables regular and initial output, but retains final output (suitable for a steady-state simulation).

.. index:: output; checkpoint, checkpoints
.. _checkpoint_output:

Output at specified times
=========================

Results can also be output at specified "checkpoint" times, as well as (or instead of) :ref:`regular_output`. Checkpoint output is written at the specified times, regardless of the time step sizes being used. At the start of each time step, a check is carried out to see if a checkpoint time would be passed using the current time step size. If so, then the time step size is reduced to hit the checkpoint time exactly.

Checkpoint output is specified using the **"checkpoint"** value, which is an object.

.. note::
   **JSON object**: checkpoint output

   **JSON path**: output.checkpoint

   +------------+------------+------------+--------------------------+
   |**name**    |**type**    |**default** |**value**                 |
   +------------+------------+------------+--------------------------+
   |"time"      |array       |[]          |checkpoint times          |
   |            |            |            |                          |
   +------------+------------+------------+--------------------------+
   |"step"      |array       |[]          |intervals between         |
   |            |            |            |checkpoint times          |
   |            |            |            |                          |
   +------------+------------+------------+--------------------------+
   |"repeat"    |integer |   |1           |how many times to repeat  |
   |            |boolean     |            |checkpoint sequence       |
   +------------+------------+------------+--------------------------+
   |"tolerance" |number      |0.1         |non-dimensional tolerance |
   |            |            |            |for detecting checkpoint  |
   |            |            |            |times                     |
   +------------+------------+------------+--------------------------+

Checkpoint times can be specified directly using the **"time"** array value. Alternatively, the intervals between checkpoint times can be specified via the **"step"** array value. In this case, the first checkpoint time is equal to the simulation start time as specified in the "time.start" value (see :ref:`time_stepping`), plus the first interval specified in the "step" array.

The specified sequence of checkpoint times (or intervals) can be repeated using the **"repeat"** value. This may be either an integer, in which case the checkpoint sequence will be repeated the specified number of times, or a boolean value. Setting it to ``true`` means the checkpoint sequence will be repeated indefinitely, until the simulation stops. Setting it to ``false`` has the same effect as setting it to 1 (i.e. the sequence is done once only, and not repeated again).

Note that if the checkpoint times are specified via the "time" array, and are repeated, then the *pattern* of times (i.e. the intervals between them) is repeated rather than the absolute times themselves (which would make no sense).

For example:

.. code-block:: json

   {"output": {"checkpoint": {"time": [1000, 2000, 3000]}}}

specifies a simple sequence of three checkpoint times. This could also be specified using steps:

.. code-block:: json

   {"output": {"checkpoint": {"step": [1000, 1000, 1000]}}}

or more simply (as the steps are all equal) using repeated steps:

.. code-block:: json

   {"output": {"checkpoint": {"step": [1000], "repeat": 3}}}

It could also be done using repeated times:

.. code-block:: json

   {"output": {"checkpoint": {"time": [1000], "repeat": 3}}}

Checkpoints every 1000 s for the entire simulation could be specified by:

.. code-block:: json

   {"output": {"checkpoint": {"time": [1000], "repeat": true}}}

The **"tolerance"** value specifies a tolerance :math:`\epsilon` for detecting when the time-stepping algorithm has hit a checkpoint. This is a non-dimensional (i.e. relative) tolerance, with the absolute tolerance given by this value multiplied by the current time step size :math:`\Delta t^n` (see :ref:`time_stepping_methods`). Specifically, the next checkpoint time :math:`t_c` will be hit (and the time step size altered to :math:`t_c - t^n`)  in the current time step if:

.. math::

   t^n + (1 + \epsilon) \Delta t^n \ge t_c

This tolerance :math:`\epsilon` is necessary for two reasons. Firstly, with no tolerance, detecting checkpoints would in some situations (e.g. when a checkpoint coincides nearly exactly with a simulated time :math:`t^n`) be subject to rounding errors, and therefore unreliable.

Secondly, the tolerance can give better time-stepping behaviour if a time step happens to fall just short of a checkpoint time. Without the tolerance, the time step would be completed, and the size of the following time step would have to be reduced to a very small value to hit the checkpoint. With the tolerance, the time step size can instead be increased slightly so that it hits the checkpoint, with no need for a subsequent reduction. This is the reason the default tolerance is relatively large (10%), larger than what would otherwise be needed simply to avoid rounding error issues.

.. index:: output; fields
.. _output_fields:

Output fields
=============

The main simulation results consist of fluid, flux and source properties, or "fields", output for each cell, face and source. It is possible to control which fields are output using the **"output.fields"** value. This is an object, with five values, **"fluid"**, **flux**, **"source"**, **"cell_geometry"** and **"face_geometry"**, specifying the fluid, flux, source, cell geometry and face geometry output fields respectively.

.. note::
   **JSON object**: output fields

   **JSON path**: output.fields

   +---------------+---------------+--------------+--------------+
   |**name**       |**type**       |**default**   |**value**     |
   +---------------+---------------+--------------+--------------+
   |"fluid"        |array | string |depends on    |fluid output  |
   |               |               |EOS           |fields        |
   +---------------+---------------+--------------+--------------+
   |"flux"         |array | string |[]            |flux output   |
   |               |               |              |fields        |
   |               |               |              |              |
   +---------------+---------------+--------------+--------------+
   |"source"       |array | string |["component", |source output |
   |               |               |"rate",       |fields        |
   |               |               |"enthalpy"]   |              |
   +---------------+---------------+--------------+--------------+
   |"cell_geometry"|array | string |["centroid",  |cell geometry |
   |               |               |"volume"]     |fields        |
   |               |               |              |              |
   +---------------+---------------+--------------+--------------+
   |"face_geometry"|array | string |["area"]      |face geometry |
   |               |               |              |fields        |
   |               |               |              |              |
   +---------------+---------------+--------------+--------------+

Each of these values can be specified as an array of strings, containing the field names. Alternatively, they can be set to the single string value **"all"**, in which case all available fields will be output.

.. index:: output; fluid
.. _output_fluid_fields:

Fluid fields
------------

The fluid fields available for output are of two types: "bulk" fields and "phase" fields. The latter are properties of particular fluid phases (e.g. liquid or vapour) whereas the former pertain to bulk properties of the fluid mixture as a whole.

The available bulk fluid fields are:

+---------------------------+-----------------------------+
|**field name**             |**value**                    |
+---------------------------+-----------------------------+
|"pressure"                 |fluid pressure (Pa)          |
|                           |                             |
+---------------------------+-----------------------------+
|"temperature"              |fluid temperature            |
|                           |(:math:`^{\circ}`\ C)        |
+---------------------------+-----------------------------+
|"region"                   |thermodynamic region         |
|                           |                             |
+---------------------------+-----------------------------+
|"phases"                   |fluid phase composition      |
|                           |                             |
+---------------------------+-----------------------------+
|`component_name` +         |partial pressures of mass    |
|"_partial_pressure"        |components (Pa)              |
+---------------------------+-----------------------------+

There is a partial pressure field for each mass component in the :ref:`eos` module being used. For example, for the :ref:`water_air_energy_eos` EOS, the mass component names are "water" and "air", so the corresponding partial pressure fluid field names are "water_partial_pressure" and "air_partial_pressure".

The available fluid phase fields are:

+--------------------------+-------------------------+
|**field name**            |**value**                |
+--------------------------+-------------------------+
|"density"                 |phase density (kg/m\     |
|                          |:sup:`3`)                |
+--------------------------+-------------------------+
|"viscosity"               |phase dynamic viscosity  |
|                          |(Pa s)                   |
+--------------------------+-------------------------+
|"saturation"              |phase saturation         |
+--------------------------+-------------------------+
|"relative_permeability"   |phase relative           |
|                          |permeability             |
+--------------------------+-------------------------+
|"capillary_pressure"      |phase capillary pressure |
|                          |(Pa)                     |
+--------------------------+-------------------------+
|"specific_enthalpy"       |phase enthalpy (J/kg)    |
+--------------------------+-------------------------+
|"internal_energy"         |phase internal energy    |
|                          |(J/kg)                   |
+--------------------------+-------------------------+
|`component_name` +        |phase component mass     |
|"_mass_fraction"          |fraction                 |
+--------------------------+-------------------------+

The name of each fluid phase field is also prepended by the phase name (and an underscore). Hence, for example, for the "liquid" phase, the field name for the saturation is "liquid_saturation".

In each phase, there is a mass fraction field for each mass component in the EOS module being used. For example, for the :ref:`water_air_energy_eos` EOS, the field name for the mass fraction of air in the "vapour" phase is "vapour_air_mass_fraction".

Each EOS module has a default set of output fluid fields, listed in the documentation for each :ref:`eos`.

Fluid fields for restarting
---------------------------

The Waiwera HDF5 output files can be used to provide initial conditions for restarting a subsequent simulation (see :ref:`restarting`). To make sure this is always possible, the fluid output fields must contain the fields corresponding to the thermodynamic :ref:`primary_variables` for the :ref:`eos` being used. Note that primary variable fields for all possible :ref:`thermodynamic_regions` must be included.

For example, for the :ref:`water_energy_eos` EOS, the fluid output fields must include "pressure", "temperature" and "vapour_saturation".

If the necessary primary variable fields are not specified in the "output.fields.fluid" array, Waiwera will automatically add them. 

.. index:: output; tracers

Tracer fields
-------------

If tracers are being simulated (see :ref:`setup_tracers`), then an output field for cell tracer mass fraction is automatically included for each tracer (there is usually little point in simulating a tracer unless it is going to be output). The field name is the same as the tracer name.

.. index:: output; flux
.. _output_flux_fields:

Flux fields
-----------

Flux fields may be output for any of the fluid mass or energy components, or the fluid phases:

+--------------------------+-------------------------+
|**field name**            |**value**                |
+--------------------------+-------------------------+
|`component_name`          |mass or energy component |
|                          |flux (kg/m\ :sup:`2`/s or|
|                          |W/m\ :sup:`2`)           |
+--------------------------+-------------------------+
|`phase_name`              |phase flux (kg/m\        |
|                          |:sup:`2`/s)              |
+--------------------------+-------------------------+

There is a flux field for each mass or energy component in the :ref:`eos` module being used. For example, for the :ref:`water_air_energy_eos` EOS, the mass component names are "water" and "air", so the corresponding mass component flux field names are also simply "water" and "air". As this EOS is non-isothermal, it also has an energy component, so there is an additional "energy" flux component field name.

There is also a flux field for each fluid phase in the :ref:`eos` module. For the :ref:`water_air_energy_eos` EOS, for example, the phases are "liquid" and "vapour", so the corresponding flux field names are also "liquid" and "vapour".

By default (regardless of the equation of state), no flux fields are output. Note that for a typical mesh, particularly in 3-D, there are significantly more faces than cells, so specifying flux output will considerably increase output file sizes.

.. index:: output; sources
.. _output_source_fields:

Source fields
-------------

The available source output fields are:

+-----------------------+-------------------------------+
|**name**               |**value**                      |
+-----------------------+-------------------------------+
|"component"            |mass or energy component       |
+-----------------------+-------------------------------+
|"rate"                 |flow rate (kg/s or J/s)        |
+-----------------------+-------------------------------+
|"enthalpy"             |enthalpy (J/kg)                |
+-----------------------+-------------------------------+
|`component_name` +     |mass or energy component flow  |
|"_flow"                |(kg/s or J/s)                  |
+-----------------------+-------------------------------+
|"source_index"         |index of source in input       |
+-----------------------+-------------------------------+
|"local_source_index"   |index of source on local       |
|                       |processor                      |
+-----------------------+-------------------------------+
|"natural_cell_index"   |cell index of source           |
+-----------------------+-------------------------------+
|"local_cell_index"     |cell index of source on local  |
|                       |processor                      |
+-----------------------+-------------------------------+
|"injection_enthalpy"   |enthalpy applied for injection |
|                       |(J/kg)                         |
+-----------------------+-------------------------------+
|"injection_component"  |component for injection        |
+-----------------------+-------------------------------+
|"production_component" |component for production       |
+-----------------------+-------------------------------+
|`tracer_name` + "_flow"|tracer flow rate (kg/s)        |
+-----------------------+-------------------------------+

There is a mass component flow field for each mass component in the :ref:`eos` module being used. For example, for the :ref:`water_air_energy_eos` EOS, there will be two mass component flow fields, "water_flow" and "air_flow".

For non-isothermal EOS modules there is also a "heat_flow" field, for flow in the energy component.

If tracers are being simulated (see :ref:`setup_tracers`), then there is an additional flow field for each tracer, with "_flow" appended to the tracer name. (Note that tracer flow rates at sources are not output by default.)

Regardless of the :ref:`eos`, the default source output fields are ["component", "rate", "enthalpy"].

.. index:: output; mesh geometry
.. _output_geometry_fields:

Mesh geometry fields
--------------------

Selected geometric properties of the :ref:`simulation_mesh` can optionally be output. These can be useful for post-processing purposes.

The available cell geometry fields are:

+--------------------------+-------------------------+
|**field name**            |**value**                |
+--------------------------+-------------------------+
|"centroid"                |array containing         |
|                          |coordinates (m) of cell  |
|                          |centroid                 |
+--------------------------+-------------------------+
|"volume"                  |cell volume (m\ :sup:`3`)|
|                          |                         |
+--------------------------+-------------------------+

The available face geometry fields are:

+--------------------------+-------------------------+
|**field name**            |**value**                |
+--------------------------+-------------------------+
|"area"                    |face area (m\ :sup:`2`)  |
|                          |                         |
|                          |                         |
+--------------------------+-------------------------+
|"distance"                |array containing distance|
|                          |(m) from face to cell    |
|                          |centroids on either side |
+--------------------------+-------------------------+
|"distance12"              |distance (m) between cell|
|                          |centroids on either side |
|                          |of the face              |
+--------------------------+-------------------------+
|"normal"                  |unit normal vector to    |
|                          |face (m)                 |
|                          |                         |
+--------------------------+-------------------------+
|"gravity_normal"          |dot product of unit      |
|                          |normal with gravity      |
|                          |vector (m/s\ :sup:`2`)   |
+--------------------------+-------------------------+
|"centroid"                |array containing         |
|                          |coordinates (m) of face  |
|                          |centroid                 |
+--------------------------+-------------------------+
|"permeability_direction"  |:ref:`rock_permeability` |
|                          |direction assigned to    |
|                          |face                     |
|                          |                         |
+--------------------------+-------------------------+

Note that if no :ref:`output_flux_fields` are specified (the default), then no face geometry fields will be output either.

Examples
--------

In the following example, the water / energy EOS is specified, with the default fluid output fields plus the densities of both the liquid and vapour phases:

.. code-block:: json

   {"eos": {"name": "we"},
    "output": {"fields": {
                  "fluid": ["pressure", "temperature", "vapour_saturation",
                            "liquid_density", "vapour_density"]}}}

Because Waiwera will automatically add all primary variable fluid fields (in this case, "pressure", "temperature" and "vapour_saturation") if they are not specified, the following JSON input would have the same effect:

.. code-block:: json

   {"eos": {"name": "we"},
    "output": {"fields": {
                  "fluid": ["liquid_density", "vapour_density"]}}}

The next example specifies the water / air / energy EOS, with source output fields of enthalpy plus the separate flows in the two mass components (water and air):

.. code-block:: json

   {"eos": {"name": "wae"},
    "output": {"fields": {
                  "source": ["enthalpy", "water_flow", "air_flow"]}}}

In this example all available fluid fields will be output:

.. code-block:: json

   {"eos": {"name": "we"},
    "output": {"fields": {"fluid": "all"}}}

The next example specifies the water / air / energy EOS, with default fluid output fields, but also requires water and air mass fluxes to be output:

.. code-block:: json

   {"eos": {"name": "wae"},
    "output": {"fields": {"flux": ["water", "air"]}}}

Here the liquid and vapour mass fluxes are output instead of component fluxes:

.. code-block:: json

   {"eos": {"name": "wae"},
    "output": {"fields": {"flux": ["liquid", "vapour"]}}}

The following example defines two tracers named "T1" and "T2" and specifies that their flow rates should be included in the source output (along with the fluid flow rate and enthalpy):

.. code-block:: json

   {"tracer": [{"name": "T1"}, {"name": "T2"}],
    "output": {"fields": {
                  "source": ["rate", "enthalpy", "T1_flow", "T2_flow"]}}}

In this example, the output mesh geometry fields are specified, so that only cell volumes (and no face geometry fields) are output:

.. code-block:: json

   {"output": {"fields": {
                  "cell_geometry": ["volume"],
                  "face_geometry": []}}}

Jacobian matrix output
======================

For some applications (e.g. inverse modelling, uncertainty quantification) it may be useful to output not only the simulation results but also the Jacobian matrix (see :ref:`jacobian_matrix`).

The Jacobian matrix is output not to the main simulation output HDF5 file (PETSc does not yet support this) but to an additional binary file. This file is in the native PETSc binary matrix file format and may be read into other codes via the PETSc library. As well as the main binary file, and accompanying text file with the ``*.info`` extension is also created. Note that the ordering of rows and columns in the output Jacobian matrix corresponds to that of the cells in the main simulation output (see :ref:`index_datasets`).

Jacobian output can be specified in the JSON input file via the **"jacobian"** value. Setting this value to the boolean ``false`` (the default) disables Jacobian output. Setting it to ``true`` enables Jacobian output, and writes the Jacobian to a binary matrix file with the same file name as the main HDF5 simulation output, but with file extension changed to ``*.jac``.

The **"jacobian"** value can also be specifed as an object containing a **"filename"** string value, which can be used to override the default file name for the binary Jacobian matrix file. Setting the file name in this way implicitly enables Jacobian output.

If Jacobian output is enabled, the Jacobian matrix is written whenever the main simulation results are written (see :ref:`initial_and_final_output` and :ref:`checkpoint_output`). For example, for a steady-state simulation in which results are written only at the end of the run, the Jacobian matrix will similarly be written only once to the binary matrix file.

Examples
--------

This example enables Jacobian output and uses the default Jacobian file name, based on the main simulation HDF5 output file name (which is itself not specified here, and is hence based on the input file name):

.. code-block:: json

    {"output": {"jacobian": true}}

This example also enables Jacobian output and specifies the Jacobian file name explicitly:

.. code-block:: json

    {"output": {"jacobian": {"filename": "model.jac"}}}

In this example, the simulation output file name is specified, and Jacobian output enabled, which will be written to a file called ``production.jac``:

.. code-block:: json

    {"output": {"filename": "production.h5", "jacobian": true}}
