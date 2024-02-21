.. index:: simulation; sources, sources
.. _source_terms:

************
Source terms
************

A **source** injects or produces (i.e. extracts) mass or energy at a particular **flow rate** (kg/s or J/s) in a given cell. By convention, the flow rate is positive for injection, and negative for production. (Here we use the term "source" to include "sinks", which are just considered to be sources with negative flow rate.)

Sources may be used to represent physical sources or sinks such as wells or springs. They may also be used to implement :ref:`neumann_boundary_conditions`, i.e. to model boundaries where flow rates of mass or energy are prescribed.

Sources are set up in the Waiwera JSON input file via the **"source"** value. This is an array of objects. Each object in the array contains a source specification, which can set up a single source, or multiple sources with similar parameters in different cells. The flow rate for each source can be specified using its **"rate"** value.

.. admonition:: JSON input
                
   **JSON object**: source specification

   **JSON path**: source[`index`]

   +-----------------------+----------------+------------+-------------------------+
   |**name**               |**type**        |**default** |**value**                |
   +-----------------------+----------------+------------+-------------------------+
   |"name"                 |string          |""          |optional source name     |
   |                       |                |            |                         |
   +-----------------------+----------------+------------+-------------------------+
   |"component"            |integer | string|"water"     |mass or energy component |
   |                       |                |(injection) |                         |
   |                       |                || 0         |                         |
   |                       |                |(production)|                         |
   +-----------------------+----------------+------------+-------------------------+
   |"production_component" |integer | string|"energy" if |mass or energy component |
   |                       |                |"component" |for :ref:`mixed_flow`    |
   |                       |                |= "energy", |production               |
   |                       |                |otherwise 0 |                         |
   +-----------------------+----------------+------------+-------------------------+
   |"rate"                 |number | array ||0           |flow rate (kg/s or J/s)  |
   |                       |object          |            |                         |
   +-----------------------+----------------+------------+-------------------------+
   |"enthalpy"             |number | array ||83.9×10\    |injection enthalpy (J/kg)|
   |                       |object          |:sup:`3`    |                         |
   |                       |                |J/kg        |                         |
   +-----------------------+----------------+------------+-------------------------+
   |"tracer"               |number | array ||0           |tracer injection rate    |
   |                       |object          |            |(kg/s)                   |
   +-----------------------+----------------+------------+-------------------------+
   |"cell" | "cells"       |integer | array |[]          |indices of cells with    |
   |                       || ``null``      |            |this source specification|
   +-----------------------+----------------+------------+-------------------------+
   |"zones"                |array           |[]          |:ref:`mesh_zones`        |
   |                       |                |            |containing cells with    |
   |                       |                |            |this source specification|
   +-----------------------+----------------+------------+-------------------------+
   |"interpolation"        |string          |"linear"    |interpolation method for |
   |                       |                |            |data tables              |
   +-----------------------+----------------+------------+-------------------------+
   |"averaging"            |string          |"integrate" |averaging method for data|
   |                       |                |            |tables                   |
   +-----------------------+----------------+------------+-------------------------+
   |"separator"            |boolean | object|``false``   |separator properties (see|
   |                       |                |            |:ref:`source_separators`)|
   +-----------------------+----------------+------------+-------------------------+
   |"deliverability"       |object          |{}          |:ref:`deliverability`    |
   |                       |                |            |source control           |
   +-----------------------+----------------+------------+-------------------------+
   |"recharge"             |object          |{}          |:ref:`recharge` source   |
   |                       |                |            |control                  |
   +-----------------------+----------------+------------+-------------------------+
   |"injectivity"          |object          |{}          |:ref:`injectivity`       |
   |                       |                |            |source control           |
   +-----------------------+----------------+------------+-------------------------+
   |"limiter"              |object          |{}          |:ref:`limiter` source    |
   |                       |                |            |control                  |
   +-----------------------+----------------+------------+-------------------------+
   |"direction"            |string          |"both"      |:ref:`direction` source  |
   |                       |                |            |control                  |
   +-----------------------+----------------+------------+-------------------------+
   |"factor"               |number | array ||{}          |:ref:`factor` source     |
   |                       |object          |            |control                  |
   +-----------------------+----------------+------------+-------------------------+

.. index:: sources; names
.. _source_names:

Source names
============

A source can be given a name using its **"name"** value, which can be an arbitrary string.

It is not required to give a source a name, unless a source network (see :ref:`source_networks`) is defined and the source is referenced by other objects (e.g. groups or reinjectors) in the source network. This referencing is done on the basis of names.

Source cells and zones
======================

Each source specification object has a **"cell"** value which can be used to specify a single cell index. There is also a **"cells"** value which can be either an integer (in which case it works exactly the same way as the "cell" value) or an array of cell indices, if multiple sources are to be set up in different cells but with the same parameters. (Note that the "cells" value must not contain duplicate cells, otherwise the "source_index" array in the output (see :ref:`index_datasets`) cannot be determined correctly.)

There is also a **"zones"** array value which can be used to specify :ref:`mesh_zones`, instead of (or as well as) cells containing the sources. When this is used, a source will be added to each cell in any of the specified zones.

There are situations in which it is useful to be able to specify a source that is not associated with any cells in the model (see :ref:`reinjector_outside`). This can be done by specifying either the "cell" or "cells" value as ``null``, or by omitting all of the "cell", "cells" and "zones" values.

.. index:: sources; injection
.. _injection:

Injection
=========

Specifying a positive **"rate"** value gives an injection source. It is also possible to specify time-dependent injection rates using a "table" source control (see :ref:`table_source_controls`), i.e. a rank-2 array of times and rates instead of a single fixed scalar flow rate.

When injecting mass, each source can inject only one mass component at a time. Depending on the :ref:`eos` (EOS) module being used, there may be multiple mass components being simulated (see :ref:`conservation_equations`). Hence, it is necessary to specify which component is being injected, via the **"component"** value. Components can be referred to either by name (string) or by index (1-based), according to the list of component names for the EOS being used. If no component is specified, the default "water" will be used. 

When mass is injected, a corresponding amount of energy is also automatically injected, according to the enthalpy (J/kg) of the injected fluid. If the mass flow rate is :math:`q`, then the corresponding energy flow rate is :math:`hq`, where :math:`h` is the enthalpy. Hence, for mass injection it is necessary to specify the injection enthalpy as well as the flow rate, via the **"enthalpy"** value. If no enthalpy is specified, a default value of 83.9 kJ/kg will be used (corresponding approximately to injecting water of temperature 20\ :math:`^{\circ}`\ C).

For example:

.. code-block:: json

   {"source": [
     {"cell": 331, "component": "water", "enthalpy": 334.9e3, "rate": 2.5}
   ]}

specifies a source injecting water with enthalpy 334.9 kJ/kg into a single cell, at a fixed rate of 2.5 kg/s.

It is also possible for a source to inject heat only (rather than mass), simply by setting the "component" value to "energy" (or to the index of the energy component, which is :math:`C+1` if the number of mass components in the EOS is :math:`C`). In this case, as no mass is being injected, there is no need to specify an enthalpy.

For example:

.. code-block:: json

   {"source": [
     {"cells": [99, 103, 231], "component": "energy", "rate": 1e3}
   ]}

specifies 1 kW heat sources in three different cells.

The following example shows how time-dependent injection may be specified using a rank-2 array of times and flow rates for the **"rate"** value. For more detail, see :ref:`table_source_controls`.

.. code-block:: json

   {"source": [
     {"cell": 313, "enthalpy": 350e3,
      "rate": [[0, 3.5], [3600, 2.4], [7200, 3.3]]}
   ]}

.. index:: sources; production

Production
==========

Specifying a negative **"rate"** value gives a production source. It is also possible to specify time-dependent production rates using a "table" source control (see :ref:`table_source_controls`), i.e. a rank-2 array of times and rates instead of a single fixed scalar flow rate.

For production, it is possible to specify a mass component to be produced (again via the "component" value), in which case only that component will be extracted from the cell. However, it is more usual to produce all mass components present in the cell. This can be done by either not specifying the "component" value, or setting it to zero.

Whenever mass is produced (either one component or all of them) the associated energy is also produced, according to the enthalpy of the fluid being extracted. However this enthalpy is a function of the thermodynamic conditions in the cell and does not need to be specified.

The JSON input for a production source is the same as for injection, apart from the flow rate being negative, and the absence of the "enthalpy" value. For example:

.. code-block:: json

   {"source": [
     {"cell": 313, "rate": -2.5}
   ]}

specifies a source producing all mass components at a fixed rate of 2.5 kg/s.

As for injection, it is also possible to produce heat only, rather than mass (e.g. to simulate a heat exchanger). For example:

.. code-block:: json

   {"source": [
     {"cells": [99, 103, 231], "component": "energy", "rate": -1e3}
   ]}

specifies three sources each extracting 1 kW of heat.

The following example shows how a time-dependent production rate can be specified using a rank-2 array of times and flow rates for the **"rate"** value. For more detail, see :ref:`table_source_controls`.

.. code-block:: json

   {"source": [
     {"cell": 313, "rate": [[0, -2.5], [3600, -2.8], [7200, -3.2]]}
   ]}

.. index:: sources; mixed flow
.. _mixed_flow:

Mixed flow
==========

The flow rate in a source may vary with time (see :ref:`source_controls`), and while it is uncommon, by default there is nothing to prevent a source from switching between production and injection during a simulation. (It is possible to limit the flow direction using a :ref:`direction` source control.)

For mixed-flow sources, it is possible to specify the production component independently of the injection component (determined by the "component" value) if desired, so that a source may inject one component and produce a different one. This can be done by specifying the **"production_component"** value. If not specified, by default it is given the value "energy" if the "component" value is also "energy". If the "component" value specifies a mass component, then "production_component" takes the default value of zero (i.e. produce all mass components).

Note that it is not necessary to set the "production_component" value except in this special case of mixed-flow sources with different components for production and injection. In all other cases, setting the "component" value by itself is sufficient.

.. index:: tracers; sources, sources; tracer

Tracer injection
================

If tracers are being simulated (see :ref:`setup_tracers`), then for :ref:`injection` sources it is possible to specify the tracer injection rates via the **"tracer"** value. This can be either:

- a scalar, representing a constant value to be applied to all tracers defined in the simulation
- a rank-1 array of numbers, with one constant value for each tracer
- a rank-2 array of numbers, representing a table of tracer injection rates vs. time (to be applied to all tracers)
- an object, with tracer names as keys and corresponding number or rank-2 array values, representing constant or time-dependent tracer injection rates

For example:

.. code-block:: json

   {"source": [
     {"cell": 331, "component": "water",
      "enthalpy": 350e3, "rate": 2.5, "tracer": 1e-6}
   ]}

specifies a source injecting water with enthalpy 350 kJ/kg at a constant rate of 2.5 kg/s, with tracer injected at :math:`10^{-6}` kg/s. In the following example:

.. code-block:: json

   {"source": [
     {"cell": 331, "component": "water",
      "enthalpy": 350e3, "rate": 2.5, "tracer": [1e-6, 1e-5, 0]}
   ]}

constant injection rates are specified for three tracers, the third one being zero. In the following example:

.. code-block:: json

   {"source": [
     {"cell": 331, "component": "water",
      "enthalpy": 350e3, "rate": 2.5,
       "tracer": [[0, 1e-5], [3600, 1e-6], [9600, 5e-7]]}
   ]}

a time-dependent tracer injection rate is specified, with values provided for three times (see :ref:`table_source_controls`).

Here is an example of specifying tracer injection using an object to refer to individual tracers by name:

.. code-block:: json

   {"source": [
     {"cell": 331, "component": "water",
      "enthalpy": 350e3, "rate": 2.5,
      "interpolation": "step",
       "tracer": {
         "T1": [[0, 1e-6], [3600, 0]],
         "T2": [[0, 0], [3600, 1e-5], [7200, 0]]}}
   ]}

In this case, it is assumed that tracers with names "T1" and "T2" have been defined separately in the input JSON file (see :ref:`setup_tracers`). For this source, tracer "T1" is injected at :math:`10^{-6}` kg/s for the first hour, after which tracer "T2" is injected at :math:`10^{-5}` kg/s for the second hour. Any tracers not included in this type of source specification will be given the default injection rate of zero.

.. index:: sources; separators, separators
.. _source_separators:

Separators
==========

A source may optionally have a **separator** which calculates separated water and steam flows from the total flow. These separated flows may then be included in the simulation output. They can also be used by certain kinds of source controls (e.g. a steam :ref:`limiter`) which require the separated flows as input. A separator also calculates the separated water and steam enthalpies.

The separated steam (:math:`q_s`) and water (:math:`q_w`) flows are calculated from the source flow rate :math:`q` as follows:

.. math::

   q_s & = f q \\
   q_w & = (1 - f) q

where :math:`f` is the steam fraction, calculated from:

.. math::

   f = \begin{cases}
   0 & h \le h_w \\
   \frac{h - h_w}{h_s - h_w} & h_w < h \le h_s \\
   1 & h > h_s
   \end{cases}

where the reference steam and water enthalpies :math:`h_s`, :math:`h_w` are calculated from their respective internal energies (:math:`U_s`, :math:`U_w`) and densities (:math:`\rho_s`, :math:`\rho_w`), and the **separator pressure** :math:`P_0` (a specified parameter) as follows:

.. math::

   h_s & = U_s + P_0 / \rho_s \\
   h_w & = U_w + P_0 / \rho_w \\

Multi-stage separators (e.g. two-stage flash) are also used, which essentially consist of several separators chained together, the separated water output of one stage being fed into the input of the next. Typically the first stage has a high separator pressure and produces high-pressure steam, the next has lower separator pressure and produces lower-pressure steam, and so on. For multi-stage separators, the calculated separated steam flow represents the total from all stages, and the steam fraction is the ratio of the total separated steam flow to the total mass flow. Similarly, the separated steam enthalpy is the (mass-flow weighted) combined enthalpy from all stages. The separated water flow rate and enthalpy are taken from the final stage.

Separator properties for a source can be specified via its **"separator"** value. This can be either a Boolean value or an object. Setting it to ``true`` specifies a separator with default separator pressure, while setting it to ``false`` (the default) means no separator is used (and any output values such as separated steam flow will be zero).

Setting it as an object containing a **"pressure"** value allows the separator pressure to be specified. This can be either a single number (for a single-stage separator) or a rank-1 array of stage separator pressures (for a multi-stage separator). Any number of stages may be used.

.. admonition:: JSON input

   **JSON object**: source separator

   **JSON path**: source[`index`]["separator"]

   +---------------+-----------------+--------------+---------------------+
   |**name**       |**type**         |**default**   |**value**            |
   +---------------+-----------------+--------------+---------------------+
   |"pressure"     |number | array   |0.55 MPa      |separator pressure   |
   |               |                 |              |:math:`P_0` (Pa), or |
   |               |                 |              |stage separator      |
   |               |                 |              |pressures            |
   +---------------+-----------------+--------------+---------------------+

For example, here a fixed-rate production source is defined with a default separator:

.. code-block:: json

   {"source": [
     {"cell": 53, "rate": -6.2, "separator": true}
   ]}

Here the separator is given a separator pressure of 50 bar:

.. code-block:: json

   {"source": [
     {"cell": 53, "rate": -6.2, "separator": {"pressure": 50e5}}
   ]}

.. index:: sources; controls, source controls
.. _source_controls:

Source controls
===============

In many cases, it is necessary to simulate sources with flow rates (and possibly other quantities such as enthalpy or tracer flow rates, for injection) that vary with time. To do this, a variety of different "source controls" may be added to a source, depending on what type of time variation is needed.

These may be straight-forward controls in which the time variation is simply prescribed, or dynamic controls which vary flow rates in response to fluid conditions in the cell or other factors. Most types of controls may be combined together to simulate more complex source behaviour (see :ref:`combining_source_controls`).

.. index:: source controls; table
.. _table_source_controls:

Tables
------

The simplest type of time variation results from flow rates or other quantities (e.g. enthalpy, tracer injection rate) being prescribed in the form of tables of values vs. time.

In the JSON input for a source specification, this can be achieved simply by specifying these values as rank-2 arrays (rather than numbers). These arrays are treated as :ref:`interpolation_tables` to enable Waiwera to compute the quantity at any time, and compute average values over the time step. The associated **"interpolation"** and **"averaging"** JSON values control the details of how these processes are carried out. (Note that the same interpolation and averaging parameters apply to different tables in the same source.)

For example:

.. code-block:: json

   {"source": [
     {"cell": 313, "rate": [[0, -2.5], [3600, -2.8], [7200, -3.2]],
      "interpolation": "step"}
   ]}

specifies a source with time-varying flow rate, defined by tabulated points at three times (0, 1 hour and 2 hours). Step (i.e. piecewise constant) interpolation is used. Since an explicit "averaging" value is not specified, the default (integration) is used.

The following example has an injection source with both flow rate and enthalpy varying piecewise-linearly with time:

.. code-block:: json

   {"source": [
     {"cell": 300,
      "rate": [[0, 1.7], [3600, 1.9], [7200, 1.6]],
      "enthalpy": [[0, 83.9e3], [1800, 98.1e3], [3600, 101.2e3], [4800, 88.7e3]],
      "interpolation": "linear"}
   ]}

Note that the tabulated flow rate and enthalpy data need not be specified at the same times.

The flow rate and / or enthalpy can equivalently be specified not as arrays but as objects containing a **"time"** array value, for example:

.. code-block:: json

   {"source": [
     {"cell": 313,
      "rate": {"time": [[0, -2.5], [3600, -2.8], [7200, -3.2]]},
      "interpolation": "step"}
   ]}

This alternative syntax is generally not needed, but is provided for consistency with other data that may be specified as tables in which the independent variable can either be time or another quantity.

.. index:: source controls; deliverability
.. _deliverability:

Deliverability
--------------

The "deliverability" source control dynamically changes the flow rate in a production source, according to the difference between the pressure in the cell and a reference pressure. This control is typically used for wells, in which case the reference pressure represents a wellbore pressure.

The total mass flow rate :math:`q` (kg/s) is given by:

.. math::
   :label: deliverability

   q = - \alpha \sum_p { \frac{k^r_p \rho_p}{\mu_p} (P - P_0)}

where the sum is taken over all phases present. The quantity :math:`\alpha` is a specified "productivity index", :math:`P` is the pressure and :math:`P_0` is the reference pressure. The quantities :math:`k^r_p`, :math:`\rho_p` and :math:`\mu_p` are respectively the phase relative permeability, density and viscosity of the fluid in the cell.

In the Waiwera JSON input file, a deliverability control is added to a source specification via its **"deliverability"** value.

.. admonition:: JSON input

   **JSON object**: deliverability source control

   **JSON path**: source[`index`]["deliverability"]

   +---------------+-----------------+--------------+---------------------+
   |**name**       |**type**         |**default**   |**value**            |
   +---------------+-----------------+--------------+---------------------+
   |"pressure"     |number | array | |10\ :sup:`5`  |reference pressure   |
   |               |object | string  |Pa            |:math:`P_0` (Pa)     |
   |               |                 |              |                     |
   |               |                 |              |                     |
   +---------------+-----------------+--------------+---------------------+
   |"productivity" |number | array | |calculated    |productivity index   |
   |               |object           |from initial  |:math:`\alpha` (m\   |
   |               |                 |rate (if      |:sup:`3`)            |
   |               |                 |specified),   |                     |
   |               |                 |otherwise 10\ |                     |
   |               |                 |:sup:`-11` m\ |                     |
   |               |                 |:sup:`3`      |                     |
   +---------------+-----------------+--------------+---------------------+
   |"threshold"    |number           |undefined     |threshold pressure   |
   |               |                 |              |(Pa)                 |
   +---------------+-----------------+--------------+---------------------+

Within a deliverability object, the reference pressure :math:`P_0` is specified via the **"pressure"** value, which may be given as:

* a constant number
* a rank-2 array representing an interpolation table (see :ref:`interpolation_tables`) of reference pressure vs. time
* an object, containing a **"time"** array value (equivalent to specifying the reference pressure itself as an array)
* an object containing an **"enthalpy"** array value, representing an interpolation table of values vs. flowing enthalpy, rather than time
* a string with value "initial", in which case the reference pressure is set equal to the pressure in the source cell at the start of the simulation

Similarly, the productivity index :math:`\alpha` is specified via the **"productivity"** value, which may be given as:

* a constant number
* a rank-2 array representing an interpolation table of productivity index vs. time
* an object, containing a **"time"** array value (equivalent to specifying the productivity index itself as an array)

If the productivity index is not specified, but an initial flow rate is specified instead via the source specification's **"rate"** value, then the productivity index will be calculated (using equation :eq:`deliverability`) to match the given flow rate. If the flow rate is not specified either, then a default value will be used.
   
The deliverability **"threshold"** value gives the option of switching on the deliverability control only when the pressure drops below the specified threshold pressure, and deactivating it again if the pressure rises back over the threshold. This option can be used, for example, for history matching simulations in which measured flow rates are specified for a well, but the model permeability is insufficient to maintain the specified flow rates without the pressure dropping towards zero, stalling the simulation.In such cases, using the "threshold" option causes the measured flow rates to be treated effectively as a target, with the well switching to deliverability if the target cannot be met. When the threshold is used, the productivity index is calculated automatically from the flow rate as the pressure drops below the threshold pressure, so that the flow rate remains consistent as the deliverability control switches on. The deliverability control will also switch off if the flow rate it computes is lower (i.e. more negative) than the specified flow rate (which can occur, for example, if the specified flow rate is time-dependent and reduces suddenly to zero).

When a deliverability control is used to model a production well, normally the flow rate should be limited to production only (i.e. if the pressure drops below the reference pressure, the well will not flow), by using a direction control (see :ref:`direction`).

For example, the source below has the simplest possible type of deliverability control, in which both the reference pressure (2 bar) and productivity index (10\ :sup:`-12` m\ :sup:`3`) are constant:
:

.. code-block:: json

   {"source": [{"cell": 10,
                "deliverability": {"pressure": 2e5, "productivity": 1e-12}}
              ]}

This source has a time-varying reference pressure as well as time-varying productivity index:

.. code-block:: json

   {"source": [{"cell": 10,
                "deliverability": {"pressure": [[0, 2.5e5],
                                                [1.5e4, 2.4e5],
                                                [4.1e4, 2.2e5]],
                                   "productivity": [[0, 1e-11],
                                                    [1.5e4, 3e-12],
                                                    [4.1e4, 1.2e-12]]}}
              ]}

This source has a constant productivity index, but an enthalpy-dependent reference pressure, decreasing from 25 bar at low enthalpies to 15 bar at 2000 kJ/kg:

.. code-block:: json

   {"source": [{"cell": 10,
                "deliverability": {
                  "productivity": 2.2e-11,
                  "pressure": {"enthalpy": [[0, 25e5],
                                            [1000e3, 25e5],
                                            [2000e3, 15e5]]}
                }}]}

This source also has an enthalpy-dependent reference pressure, and has its productivity index calculated from a specified initial flow rate of -3.2 kg/s:

.. code-block:: json

   {"source": [{"cell": 10,
                "rate": -3.2,
                "deliverability": {
                  "pressure": {"enthalpy": [[0, 25e5],
                                            [1000e3, 25e5],
                                            [2000e3, 15e5]]}
                }}]}

This source has a table of specified flow rates vs. time, but switches to deliverability if the pressure drops below the threshold value of 2 bar:

.. code-block:: json

   {"source": [
     {"cell": 313, "rate": [[0, -2.5], [3600, -2.8], [7200, -3.2]],
      "deliverability": {"pressure": 1e5, "productivity": 1e-12, "threshold": 2e5}}
   ]}

.. index:: source controls; recharge
.. _recharge:

Recharge
--------

Like the deliverability source control, the "recharge" control also dynamically controls the source flow rate based on the difference between the pressure and a reference pressure. However, the relationship between flow rate :math:`q` and pressure difference is via a simple proportionality constant, called the "recharge coefficient":

.. math::

   q = -\beta (P - P_0)

where :math:`P` is the pressure, :math:`P_0` is the reference pressure and :math:`\beta` is the recharge coefficient.

Recharge controls are most commonly used to implement boundary conditions, for example at the side boundaries of a transient reservoir model, where it may be necessary to allow inflow or outflow as the pressures in the interior change.

In the Waiwera JSON input file, a recharge control is added to a source specification via its **"recharge"** value.

.. admonition:: JSON input

   **JSON object**: recharge source control

   **JSON path**: source[`index`]["recharge"]

   +--------------+------------+------------+-------------------+
   |**name**      |**type**    |**default** |**value**          |
   +--------------+------------+------------+-------------------+
   |"pressure"    |number |    |10\ :sup:`5`|reference pressure |
   |              |array |     |Pa          |:math:`P_0` (Pa)   |
   |              |object |    |            |                   |
   |              |string      |            |                   |
   +--------------+------------+------------+-------------------+
   |"coefficient" |number |    |10\         |recharge           |
   |              |array |     |:sup:`-2`   |coefficient        |
   |              |object      |m.s         |:math:`\beta` (m.s)|
   |              |            |            |                   |
   +--------------+------------+------------+-------------------+

Within a recharge object, the reference pressure :math:`P_0` is specified via the **"pressure"** value, which may be given as:

* a constant number
* a rank-2 array representing an interpolation table (see :ref:`interpolation_tables`) of reference pressure vs. time
* an object, containing a **"time"** array value (equivalent to specifying the reference pressure itself as an array)
* an object containing an **"enthalpy"** array value, representing an interpolation table of values vs. flowing enthalpy, rather than time
* a string with value "initial", in which case the reference pressure is set equal to the pressure in the source cell at the start of the simulation

Similarly, the recharge coefficient :math:`\beta` is specified via the **"coefficient"** value, which may be given as:

* a constant number
* a rank-2 array representing an interpolation table of productivity index vs. time
* an object, containing a **"time"** array value (equivalent to specifying the productivity index itself as an array)

For example, the source below has a recharge control with reference pressure set to the pressure at the start of the simulation, and a recharge coefficient of 10\ :sup:`-3` m.s:

.. code-block:: json

   {"source": [
     {"cell": 200, "recharge": {"pressure": "initial", "coefficient": 1e-3}}
   ]}

.. index:: source controls; injectivity
.. _injectivity:

Injectivity
-----------

The injectivity source control is typically used for injection wells which inject fluid against a specified pressure, i.e. the injection rate decreases as the fluid pressure increases.

The injectivity source control is in fact exactly the same as the :ref:`recharge` control, and the input and options are identical, except that in the Waiwera JSON input file, an injectivity control is added to a source specification via its **"injectivity"** value.

.. admonition:: JSON input

   **JSON object**: injectivity source control

   **JSON path**: source[`index`]["injectivity"]

   +--------------+------------+------------+-------------------+
   |**name**      |**type**    |**default** |**value**          |
   +--------------+------------+------------+-------------------+
   |"pressure"    |number |    |10\ :sup:`5`|reference pressure |
   |              |array |     |Pa          |:math:`P_0` (Pa)   |
   |              |object |    |            |                   |
   |              |string      |            |                   |
   +--------------+------------+------------+-------------------+
   |"coefficient" |number |    |10\         |injectivity        |
   |              |array |     |:sup:`-2`   |coefficient        |
   |              |object      |m.s         |:math:`\beta` (m.s)|
   |              |            |            |                   |
   +--------------+------------+------------+-------------------+

Usually it is not desirable for the injection rate to go negative if the fluid pressure increases above the reference pressure. Hence, like the :ref:`deliverability` source control, the injectivity control is usually used in conjunction with a direction control (see :ref:`direction`). For example:

.. code-block:: json

   {"source": [
     {"cell": 200,
      "direction": "injection",
      "injectivity": {"pressure": 65e5, "coefficient": 1e-4}
      }
   ]}

.. index:: source controls; limiter
.. _limiter:

Limiter
-------

In some situations it is necessary to limit the flow rate of a source, so that it cannot exceed a prescribed maximum value -- for example, when a well has a prescribed maximum flow rate to comply with regulations. In the simplest case the limit applies to the total flow, but in other situations the source output may be passed through a separator, and the limit is set on either separated steam or water. It is also possible for a limiter to limit multiple different flow types at once, with different limits set for total flow and/or separated water and steam flows.

A limiter may be added to a source in the Waiwera JSON input file by specifying the **"limiter"** value in that source. This value is an object, which contains values corresponding to the flow types being limited (e.g. "total" or "steam"). These values are positive and apply to the absolute value of the flow rate. They can be either single constant limits or a rank-2 arrays representing an interpolation table of limit values vs. time.

When "water" or "steam" limit values are specified, a separator (see :ref:`source_separators`) is used to compute the flow rates of separated steam and water. The separator pressure can be specified  via the **"separator.pressure"** value for the source.

If time-dependent limits are specified, it is possible to specify the interpolation and averaging type for the interpolation tables via the **"interpolation"** and **"averaging"** values respectively (see :ref:`interpolation_tables`). Note that the limit value for each simulation time step will be taken from the average value of the specified limit table over the time step.

.. admonition:: JSON input

   **JSON object**: limiter source control

   **JSON path**: source[`index`]["limiter"]

   +---------------------+------------+------------+-------------------------+
   |**name**             |**type**    |**default** |**value**                |
   +---------------------+------------+------------+-------------------------+
   |"total"              |number |    |no limit    |total flow rate limit    |
   |                     |array       |            |(kg/s)                   |
   +---------------------+------------+------------+-------------------------+
   |"water"              |number |    |no limit    |separated water flow rate|
   |                     |array       |            |limit (kg/s)             |
   +---------------------+------------+------------+-------------------------+
   |"steam"              |number |    |no limit    |separated steam flow rate|
   |                     |array       |            |limit (kg/s)             |
   +---------------------+------------+------------+-------------------------+
   |"interpolation"      |string      |"linear"    |interpolation method for |
   |                     |            |            |limit tables             |
   +---------------------+------------+------------+-------------------------+
   |"averaging"          |string      |"integrate" |averaging method for     |
   |                     |            |            |limit tables             |
   +---------------------+------------+------------+-------------------------+

The example below specifies a source on deliverability, with a simple limit of 5.1 kg/s on the total flow rate. (Because it is the total flow being limited, a separator is not needed.) 

.. code-block:: json

   {"source": [
     {"cell": 100,
      "deliverability": {"pressure": 2e5, "productivity": 1e-12},
      "limiter": {"total": 5.1}}
   ]}

Here is the same source but with a limit of 3.5 kg/s on the steam flow, and a separator is defined with separator pressure set at 50 bar:

.. code-block:: json

   {"source": [
     {"cell": 100,
      "deliverability": {"pressure": 2e5, "productivity": 1e-12},
      "separator": {"pressure": 50e5},
      "limiter": {"steam": 3.5}}
   ]}

Here the above source is modified with a time-dependent steam limit, reducing stepwise from 3.5 kg/s at time zero to 2.5 kg/s after 1e6 seconds:

.. code-block:: json

   {"source": [
     {"cell": 100,
      "deliverability": {"pressure": 2e5, "productivity": 1e-12},
      "separator": {"pressure": 50e5},
      "limiter": {"steam": [[0, 3.5], [1e6, 2.5]],
                  "interpolation": "step"}}
   ]}

Here is the same source but with constant limits applied to both total flow and steam flow:

.. code-block:: json

   {"source": [
     {"cell": 100,
      "deliverability": {"pressure": 2e5, "productivity": 1e-12},
      "separator": {"pressure": 50e5},
      "limiter": {"total": 5.1, "steam": 3.5}}
   ]}

Limiters can also be specified using an alternative (older) syntax, in which a **"type"** string value is used to specify the limit type ("total", "water" or "steam"). The flow rate limit is then set via a **"limit"** value. This value can again be specified as a constant number limit or an array of limit values vs. time.

Note that this older syntax is less flexible and cannot be used to specify multiple limits. It is retained for backwards compatibility.

With this alternative syntax, when the "type" value is "water" or "steam" the separator pressure may be specified within the limiter specification (again for backwards compatibility), although this is now deprecated: it should usually be specified  via the **"separator.pressure"** value for the source instead.

.. admonition:: JSON input

   **JSON object**: limiter source control (alternative syntax)

   **JSON path**: source[`index`]["limiter"]

   +---------------------+------------+------------+-------------------------+
   |**name**             |**type**    |**default** |**value**                |
   +---------------------+------------+------------+-------------------------+
   |"type"               |string      |"total"     |limiter type             |
   |                     |            |            |("total" | "water"       |
   |                     |            |            || "steam")               |
   |                     |            |            |                         |
   +---------------------+------------+------------+-------------------------+
   |"limit"              |number |    |1 kg/s      |flow rate limit          |
   |                     |array       |            |(kg/s)                   |
   +---------------------+------------+------------+-------------------------+
   |"interpolation"      |string      |"linear"    |interpolation method for |
   |                     |            |            |limit table              |
   +---------------------+------------+------------+-------------------------+
   |"averaging"          |string      |"integrate" |averaging method for     |
   |                     |            |            |limit table              |
   +---------------------+------------+------------+-------------------------+
   |"separator_pressure" |number      |0.55 MPa    |separator pressure       |
   |                     |            |            |(Pa)                     |
   +---------------------+------------+------------+-------------------------+

.. index:: source controls; direction
.. _direction:

Direction
---------

As mentioned above (see :ref:`mixed_flow`), it is possible for a source's flow rate to change sign during a simulation. The flow rate in a specified rate table may contain both positive and negative flow rates, although this is not common (it could potentially be used e.g. for a production well which is shut in, and later used as an reinjection well). Deliverability and recharge source controls may give flow rates that change sign, if the pressure drops below (or rises above) the reference pressure.

The flow rate may be limited to a particular direction by using a "direction" source control, via the **"direction"** value of the source. This is a simple string value which may be set to "production" or "out" if the flow rate should always remain negative, or to "injection" or "in" if the flow rate should always remain positive.

With this control applied, flow rates are set to zero if they would otherwise flow in the direction opposite to that specified. Setting the limiter value to "both" is equivalent to not specifying a limiter -- both directions are allowed.

For example:

.. code-block:: json

   {"source": [
     {"cell": 200, "recharge": {"pressure": "initial", "coefficient": 1e-3},
      "direction": "in"
     }
   ]}

specifies a recharge source that can only flow into the model, not out. A direction control can be added to a well on deliverability as follows, to ensure it stops flowing if the pressure drops below the reference pressure:

.. code-block:: json

   {"source": [{"cell": 10,
                "deliverability": {"pressure": 2e5, "productivity": 1e-12},
                "direction": "production"}
              ]}

.. index:: source controls; factor
.. _factor:

Factor
------

In some situations it can be useful to apply a scale factor to the flow rate, particularly if the flow rate is not prescribed but is computed using a dynamic control such as :ref:`deliverability`. Multiplying the flow rate by a factor might be used to simulate changes in well performance over time, e.g. from scaling or makeovers, or to shut in a well on deliverability at a particular time.

A factor control can be added to a source via its **"factor"** value. This can take several forms:

* a simple number, to apply a constant scale factor to the flow rate
* a rank-2 array representing an interpolation table (see :ref:`interpolation_tables`) of scale factor vs. time, to apply a time-dependent scale factor
* an object, containing a **"time"** array value, as well as optional **"interpolation"** and **"averaging"** values (see :ref:`interpolation_tables`)

Specifying the "factor" value as an object allows it to have its own parameters for interpolation and averaging, separate from those used to interpolate or average the source flow rate and enthalpy. This can be useful if, for example, a well uses linear interpolation for flow rate, but a step interpolation is more appropriate for the factor control, to simulate shutting the well in at a particular time.

For example:

.. code-block:: json

   {"source": [{"cell": 10,
                "deliverability": {"pressure": 2e5, "productivity": 1e-12},
                "direction": "production",
                "factor": [[0, 1],
                           [3.15576e7, 0.95],
                           [6.31152e7, 0.73],
                           [9.46728e7, 0.89]]}
              ]}

specifies a production well on deliverability, with a declining scale factor applied over the first three years of production. Here no parameters are specified for interpolation or averaging, so the defaults (linear interpolation, integration averaging) are used for both flow rates and the scale factor.

The following example uses step interpolation to simulate shutting in a deliverability well at time 10\ :sup:`8` seconds:

.. code-block:: json

   {"source": [{"cell": 10,
                "deliverability": {"pressure": 2e5, "productivity": 1e-12},
                "direction": "production",
                "factor": {"time": [[0, 1], [1e8, 0]], "interpolation": "step"}}
              ]}

.. index:: source controls; combining
.. _combining_source_controls:

Combining source controls
-------------------------

As we have seen in some of the examples above, it is possible to use different source controls together on one source, to simulate more complex behaviour. In fact, in principle it is possible to use any combination of source controls together on the same source.

However, some of these combinations are more useful than others. There is no point in having multiple controls that independently assign different flow rates to the same source, for example, a deliverability control and a recharge control.

Waiwera applies controls to a source in a pre-defined order -- in fact, the same order they have been described here. (The order in which they are specified in the JSON input file is not important.) So, for example, if a source did have both a deliverability control and a recharge control, the flow rate computed by the deliverability control would be overridden by the flow rate computed by the recharge control. Controls which do not compute a flow rate (e.g. limiters, direction and factor controls), but only modify flow rates computed by other controls, are applied last.

Interactions between sources
============================

For simulations in which sources interact (e.g. for groupings of wells, or reinjection), these interactions are defined not as part of the individual source specifications, but rather by means of a "source network" -- see :ref:`source_networks`.
