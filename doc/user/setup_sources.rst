.. _source_terms:

************
Source terms
************

A **source** injects or produces (i.e. extracts) mass or energy at a particular **flow rate** (kg/s or J/s) in a given cell. By convention, the flow rate is positive for injection, and negative for production. (Here we use the term "source" to include "sinks", which are just considered to be sources with negative flow rate.)

Sources may be used to represent physical sources or sinks such as wells or springs. They may also be used to implement :ref:`neumann_boundary_conditions`, i.e. model boundaries where flow rates of mass or energy are prescribed.

Sources are set up in the Waiwera JSON input file via the **"source"** value. This is an array of objects. Each object in the array contains a source specification, which can set up a single source, or multiple sources with similar parameters in different cells.

.. note::
   **JSON object**: source specification

   +-----------------------+----------------+------------+-------------------------+
   |**name**               |**type**        |**default** |**value**                |
   +-----------------------+----------------+------------+-------------------------+
   |"name"                 |string          |""          |optional source name     |
   |                       |                |            |                         |
   +-----------------------+----------------+------------+-------------------------+
   |"component"            |number | string |"water"     |mass or energy component |
   |                       |                |(injection) |                         |
   |                       |                || 0         |                         |
   |                       |                |(production)|                         |
   +-----------------------+----------------+------------+-------------------------+
   |"production_component" |number | string |0           |mass or energy component |
   |                       |                |            |for production           |
   +-----------------------+----------------+------------+-------------------------+
   |"rate"                 |number | array ||0           |flow rate (kg/s or J/s)  |
   |                       |object          |            |                         |
   +-----------------------+----------------+------------+-------------------------+
   |"enthalpy"             |number | array ||83.9×10\    |injection enthalpy (J/kg)|
   |                       |object          |:sup:`3`    |                         |
   |                       |                |J/kg        |                         |
   +-----------------------+----------------+------------+-------------------------+
   |"cell" | "cells"       |number | array  |[]          |indices of cells with    |
   |                       |                |            |this source specification|
   +-----------------------+----------------+------------+-------------------------+
   |"interpolation"        |string          |"linear"    |interpolation method for |
   |                       |                |            |data tables              |
   +-----------------------+----------------+------------+-------------------------+
   |"averaging"            |string          |"integrate" |averaging method for data|
   |                       |                |            |tables                   |
   +-----------------------+----------------+------------+-------------------------+
   |"deliverability"       |object          |{}          |:ref:`deliverability`    |
   |                       |                |            |source control           |
   +-----------------------+----------------+------------+-------------------------+
   |"recharge"             |object          |{}          |:ref:`recharge` source   |
   |                       |                |            |control                  |
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

Source cells
============

Each source specification object has a **"cell"** value which can be used to specify a single cell index. There is also a **"cells"** value which can be either a number (in which case it works exactly the same way as the "cell" value) or an array of cell indices, if multiple sources are to be set up in different cells but with the same parameters.

Injection
=========

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

Production
==========

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

Source controls
===============

In many cases, it is necessary to simulate sources with flow rates (and possibly enthalpies, for injection) that vary with time. To do this, a variety of different "source controls" may be added to a source, depending on what type of time variation is needed.

These may be straight-forward controls in which the time variation is simply prescribed, or dynamic controls which vary flow rates in response to fluid conditions in the cell, or other factors. Most types of controls may be combined together to simulate more complex source behaviour (see :ref:`combining_source_controls`).

Tables
------

The simplest type of time variation results from flow rates and / or injection enthalpies being prescribed in the form of tables of values vs. time.

In the JSON input for a source specification, this can be achieved simply by specifying the "rate" and / or "enthalpy" values as rank-2 arrays (rather than numbers). These are treated as :ref:`interpolation_tables` to enable Waiwera to compute the flow rate and / or enthalpy at any time, and compute average values over the time step. The associated **"interpolation"** and **"averaging"** JSON values control the details of how these processes are carried out. (Note that the same interpolation and averaging parameters apply to both flow rate and enthalpy.)

For example:

.. code-block:: json

   {"source": [
     {"cell": 313, "rate": [[0, -2.5], [3600, -2.8], [7200, -3.2]], "interpolation": "step"}
   ]}

specifies a source with time-varying flow rate, defined by tabulated points at three times. Step (i.e. piecewise constant) interpolation is used. Since an explicit "averaging" value is not specified, the default (integration) is used.

The following example has a production source with both flow rate and enthalpy varying piecewise-linearly with time:

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

.. _deliverability:

Deliverability
--------------

The "deliverability" source control dynamically changes the flow rate in a production source, according to the difference between the pressure in the cell and a reference pressure. This control is typically used for wells, and the reference pressure represents a wellbore pressure.

The total mass flow rate :math:`q` (kg/s) is given by:

.. math::
   :label: deliverability

   q = - \alpha \sum_p { \frac{k_r^p \rho_p}{\mu_p} (P - P_0)}

where the sum is taken over all phases present. The quantity :math:`\alpha` is a specified "productivity index", :math:`P` is the pressure and :math:`P_0` is the reference pressure. The quantities :math:`k_r^p`, :math:`\rho_p` and :math:`\mu_p` are respectively the phase relative permeability, density and viscosity of the fluid in the cell.

In the Waiwera JSON input file, a deliverability control is added to a source specification via its **"deliverability"** value.

.. note::
   **JSON object**: deliverability source control

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
   
The deliverability **"threshold"** value gives the option of switching on the deliverability control only when the pressure drops below the specified threshold pressure, and deactivating it again if the pressure rises back over the threshold. This option can be used, for example, for history matching simulations in which measured flow rates are specified for a well, but the model permeability is insufficient to maintain the specified flow rates without the pressure dropping towards zero, stalling the simulation. In such cases, using the "threshold" option causes the measured flow rates to be treated effectively as a target, with the well switching to deliverability if the target cannot be met.

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
                "deliverability": {"pressure": [[0, 2.5e5], [1.5e4, 2.4e5], [4.1e4, 2.2e5]],
                                   "productivity": [[0, 1e-11], [1.5e4, 3e-12], [4.1e4, 1.2e-12]]}}
              ]}

This source has a constant productivity index, but an enthalpy-dependent reference pressure, decreasing from 25 bar at low enthalpies to 15 bar at 2000 kJ/kg:

.. code-block:: json

   {"source": [{"cell": 10,
                "deliverability": {
                  "productivity": 2.2e-11,
                  "pressure": {"enthalpy": [[0, 25e5], [1000e3, 25e5], [2000e3, 15e5]]}
                }}]}

This source also has an enthalpy-dependent reference pressure, and has its productivity index calculated from a specified initial flow rate of -3.2 kg/s:

.. code-block:: json

   {"source": [{"cell": 10,
                "rate": -3.2,
                "deliverability": {
                  "pressure": {"enthalpy": [[0, 25e5], [1000e3, 25e5], [2000e3, 15e5]]}
                }}]}

This source has a table of specified flow rates vs. time, but switches to deliverability if the pressure drops below the threshold value of 2 bar:

.. code-block:: json

   {"source": [
     {"cell": 313, "rate": [[0, -2.5], [3600, -2.8], [7200, -3.2]],
      "deliverability": {"pressure": 1e5, "productivity": 1e-12, "threshold": 2e5}}
   ]}

.. _recharge:

Recharge
--------

Like the deliverability source control, the "recharge" control also dynamically controls the source flow rate based on the difference between the pressure and a reference pressure. However, the relationship between flow rate :math:`q` and pressure difference is via a simple proportionality constant, called the "recharge coefficient":

.. math::

   q = -\beta (P - P_0)

where :math:`P` is the pressure, :math:`P_0` is the reference pressure and :math:`\beta` is the recharge coefficient.

Recharge controls are most commonly used to implement boundary conditions, for example at the side boundaries of a transient reservoir model, where it may be necessary to allow inflow or outflow as the pressures in the interior change.

In the Waiwera JSON input file, a recharge control is added to a source specification via its **"recharge"** value.

.. note::
   **JSON object**: recharge source control

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

.. _limiter:

Limiter
-------

In some situations it is necessary to limit the flow rate of a source, so that it cannot exceed a prescribed maximum value -- for example, when a well has a prescribed maximum flow rate to comply with regulations. In the simplest case the limit applies to the total flow, but in other situations the source output may be passed through a separator, and the limit is set on either separated steam or water.

A limiter may be added to a source in the Waiwera JSON input file by specifying the **"limiter"** value in that source. This value is an object, which has a **"type"** string value specifying whether the limit is set on total flow, separated water flow or steam flow. The flow rate limit is set via the **"limit"** value. Note that this value is positive and applies to the absolute value of the flow rate.

.. note::
   **JSON object**: limiter source control

   +---------------------+------------+------------+------------------+
   |**name**             |**type**    |**default** |**value**         |
   +---------------------+------------+------------+------------------+
   |"type"               |string      |"total"     |limiter type      |
   |                     |            |            |("total" | "water"|
   |                     |            |            || "steam")        |
   |                     |            |            |                  |
   +---------------------+------------+------------+------------------+
   |"limit"              |number      |1 kg/s      |flow rate limit   |
   |                     |            |            |(kg/s)            |
   +---------------------+------------+------------+------------------+
   |"separator_pressure" |number      |55×10\      |separator pressure|
   |                     |            |:sup:`5` Pa |:math:`P_0` (Pa)  |
   +---------------------+------------+------------+------------------+

When the "type" value is "water" or "steam", a simple separator is simulated to compute the flow rates of separated steam (:math:`q_s`) and water (:math:`q_w`) from the source flow rate :math:`q` and fluid composition:

.. math::

   q_s & = f q \\
   q_w & = (1 - f) q

where :math:`f` is the steam fraction, calculated from:

.. math::

   f = \begin{cases}
   0 & h \le h_w \\
   \frac{h - h_w}{h_s - h_w} & hw < h \le h_s \\
   1 & h > h_s
   \end{cases}

where the steam and water enthalpies :math:`h_s`, :math:`h_w` are calculated from their respective internal energies (:math:`U_s`, :math:`U_w`) and densities (:math:`\rho_s`, :math:`\rho_w`), and the separator pressure :math:`P_0` (specified via the **"separator_pressure"** value), as follows:

.. math::

   h_s & = U_s + P_0 / \rho_s \\
   h_w & = U_w + P_0 / \rho_w \\

The example below specifies a source on deliverability, with a simple limit of 5.1 kg/s on the total flow rate. (Because it is the total flow being limited, there is no need to specify a separator pressure.) 

.. code-block:: json

   {"source": [
     {"cell": 100,
      "deliverability": {"pressure": 2e5, "productivity": 1e-12},
      "limiter": {"limit": 5.1}}
   ]}

Here is the same source but with a limit of 3.5 kg/s on the steam flow, and the separator pressure set at 50 bar:

.. code-block:: json

   {"source": [
     {"cell": 100,
      "deliverability": {"pressure": 2e5, "productivity": 1e-12},
      "limiter": {"limit": 3.5, "type": "steam", "separator_pressure": 50e5}}
   ]}

.. _direction:

Direction
---------

By default, there is nothing to prevent a source from switching between production and injection during a simulation. The flow rate in a specified rate table may contain both positive and negative flow rates, although this is not common. Deliverability and recharge source controls may give flow rates that change sign, if the pressure drops below (or rises above) the reference pressure. In the case of recharge this may happen naturally if, for example, pressures in a reservoir drop during production and rise again after production ceases.

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
                "factor": [[0, 1], [3.15576e7, 0.95], [6.31152e7, 0.73], [9.46728e7, 0.89]]}
              ]}

specifies a production well on deliverability, with a declining scale factor applied over the first three years of production. Here no parameters are specified for interpolation or averaging, so the defaults (linear interpolation, integration averaging) are used for both flow rates and the scale factor.

The following example uses step interpolation to simulate shutting in a deliverability well at time 10\ :sup:`8` seconds:

.. code-block:: json

   {"source": [{"cell": 10,
                "deliverability": {"pressure": 2e5, "productivity": 1e-12},
                "direction": "production",
                "factor": {"time": [[0, 1], [1e8, 0]], "interpolation": "step"}}
              ]}

.. _combining_source_controls:

Combining source controls
-------------------------

As we have seen in some of the examples above, it is possible to use different source controls together on one source, to simulate more complex behaviour. In fact, in principle it is possible to use any combination of source controls together on the same source.

However, some of these combinations are more useful than others. There is no point having multiple controls that independently assign different flow rates to the same source, for example, a deliverability control and a recharge control.

Waiwera applies controls to a source in a pre-defined order -- in fact, the same order they have been described here. So, for example, if a source did have both a deliverability control and a recharge control, the flow rate computed by the deliverability control would be overridden by the flow rate computed by the recharge control. Controls which do not compute a flow rate (e.g. limiters, direction and factor controls), but only modify flow rates computed by other controls, are applied last.

.. source may change between injection and production (production_component value?)

