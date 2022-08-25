.. index:: simulation; source networks, source networks
.. _source_networks:

***************
Source networks
***************

In many simulations, some sources (see :ref:`source_terms`) may not act entirely independently, but instead may interact in some way. Examples of this include:

- multi-feed wells, with sources for individual feedzones grouped together into a well, which behaves (from the above-ground perspective) in many ways as a single unit
- borefields, with production wells grouped together and limits or targets imposed on the group
- reinjection, with some of the output from a group of production wells distributed amongst a set of injection wells

To model the interactions between sources in a Waiwera simulation, a **source network** may be defined. A source network has two parts:

- **groups**, defining groupings of production sources
- **reinjectors**, defining the distribution of produced fluid (from sources or groups) amongst injection sources

In the Waiwera JSON input file a source network may be set up using the **network** value. This is an object containing two array values, **"group"** and **"reinject"**, defining source network :ref:`source_network_groups` and :ref:`source_network_reinjectors` respectively.

.. note::
   **JSON object**: source network

   **JSON path**: network

   +-----------------------+----------------+------------+-------------------------+
   |**name**               |**type**        |**default** |**value**                |
   +-----------------------+----------------+------------+-------------------------+
   |"group"                |array           |[]          |source groups            |
   |                       |                |            |                         |
   +-----------------------+----------------+------------+-------------------------+
   |"reinject"             |array           |[]          |source reinjectors       |
   |                       |                |            |                         |
   +-----------------------+----------------+------------+-------------------------+

.. index:: source networks; groups, source groups, groups
.. _source_network_groups:

Groups
======

A **group** takes multiple inputs and combines them into a single output. The inputs may be production sources (i.e. with negative rate), or other groups, enabling a nested hierarchy of sub-groups and groups to be defined if needed.

The mass flow rate at the output is simply the sum of the mass flow rates from the inputs. The enthalpy at the output is a weighted sum of the input enthalpies (weighted by mass flow rate). Hence, the total input and output mass and energy flow rates are equal.

A group also calculates output mass flow rates and enthalpies for separated water and steam (i.e. liquid and vapour phases) from two-phase inputs. How this is done depends on whether the group has its own separator or not  (see :ref:`group_separator`).

Source groups are set up in the Waiwera JSON input file via the **"network.group"** value. This is an array of objects. Each object in the array specifies a single group.

.. note::
   **JSON object**: source group

   **JSON path**: network.group[`index`]

   +-----------------------+----------------+------------+---------------------------+
   |**name**               |**type**        |**default** |**value**                  |
   +-----------------------+----------------+------------+---------------------------+
   |"name"                 |string          |""          |:ref:`group_name`          |
   |                       |                |            |                           |
   +-----------------------+----------------+------------+---------------------------+
   |"in"                   |array           |[]          |:ref:`group_inputs`        |
   |                       |                |            |                           |
   +-----------------------+----------------+------------+---------------------------+
   |"limiter"              |object          |{}          |:ref:`group_limiter`       |
   |                       |                |            |                           |
   +-----------------------+----------------+------------+---------------------------+
   |"scaling"              |string          |"uniform"   |:ref:`group_scaling` type  |
   |                       |                |            |                           |
   |                       |                |            |                           |
   +-----------------------+----------------+------------+---------------------------+
   |"separator"            |boolean | object|``false``   |:ref:`group_separator`     |
   |                       |                |            |                           |
   +-----------------------+----------------+------------+---------------------------+

.. index:: source groups; name
.. _group_name:

Group name
----------

A group can be given a name using its **"name"** value, which can be an arbitrary string.

Although the name is optional, references between objects (groups or reinjectors) in the source network are done by name. So if a group is to be referenced by another object, for example to include a group in another group, or feed its output to a reinjector, it must have a name.

The name must also be unique, that is, no other sources, groups or reinjectors may have the same name.

.. index:: source groups; inputs
.. _group_inputs:

Group inputs
------------

The **"in"** value defines the inputs for the group. This is an array of strings, listing the names of all the inputs. The inputs may be production sources, other groups, or a mix of both.

Note that groups can be defined in arbitrary order, so a group can include other groups in its inputs even if they are defined later in the group array.

The order of the inputs in a group is also usually arbitrary, but can be significant if for example progressive :ref:`group_scaling` is used.

A source or group cannot belong to (i.e. be specified as an input for) more than one group.

Example:

.. code-block:: json

   {"source": [
     {"name": "prod 1", "cell": 331, "rate": -2.5},
     {"name": "prod 2", "cell": 306,
      "deliverability": {"pressure": 5e5, "productivity": 1e-11}}
   ],
   "network": {
     "group": [{"name": "production", "in": ["prod 1", "prod 2"]}]
   }}

defines two production sources (one with fixed rate and the other on deliverability) and a source network group named "production" containing the two sources.

In the following example, four production sources are defined. A group "g1" is defined which contains two sub-groups, "sub1" and "sub2". These sub-groups each contain two of the sources:

.. code-block:: json

   {"source": [
     {"name": "p1", "cell": 331, "rate": -2.5},
     {"name": "p2", "cell": 307, "rate": -1.7},
     {"name": "p3", "cell": 351, "rate": -3.2},
     {"name": "p4", "cell": 301, "rate": -4.3}
   ],
   "network": {
     "group": [
        {"name": "g1", "in": ["sub1", "sub2"]},
        {"name": "sub1", "in": ["p1", "p2"]},
        {"name": "sub2", "in": ["p3", "p4"]}
     ]
   }}

.. index:: source groups; limiters
.. _group_limiter:

Group limiter
-------------

A group limiter is very similar to a :ref:`limiter` source control, which controls the flow rate in a source so that a specified maximum flow rate is not exceeded. A group limiter does the same thing for a group. The difference is that whereas a limiter source control controls the source's own flow rate, a group limiter controls the flow rates of its inputs.

As for a limiter source control, a group limiter may specify limits on total flow, separated water flow or separated steam flow, or any combination of these, and the specified limits can be constant or time-dependent. The JSON input for a group limiter is specified using the group's **"limiter"** value and is exactly the same as for a limiter source control.

.. note::
   **JSON object**: group limiter

   **JSON path**: network.group[`index`]["limiter"]

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
   |                     |            |            |limit arrays             |
   +---------------------+------------+------------+-------------------------+
   |"averaging"          |string      |"integrate" |averaging method for     |
   |                     |            |            |limit arrays             |
   +---------------------+------------+------------+-------------------------+

If any of the output flow rates (total, water or steam) exceeds its specified corresponding limit, the flow rates in some or all of the group inputs are reduced so that the output flow rate equals the limit. Exactly how this is done depends on the :ref:`group_scaling` policy.

If a group with a limiter has other groups in its inputs, then the limiter is applied recursively - that is, the limiter imposes limits on the input groups (according to its scaling policy), which in turn impose limits on their own inputs. This continues down until the level of individual sources is reached, at which point the flow rates in these sources are adjusted to meet the higher-level group limits.

Note that in more complex situations where there are higher-level groups with limiters, and they have lower-level groups with their own limiters as inputs, conflicts between the limiters can occur, and it may not be possible for all the groups to meet their specified limits exactly. Some of the groups may have output flow rates that are below their specified limits. In such cases, higher-level groups are given priority - that is, the highest-level groups will meet their limits exactly while the lower-level groups may not.

Example:

.. code-block:: json

   {"network": {
     "group": [{"name": "group 1",
                "in": ["rx 321", "rx 410", "rx 50"],
                "limiter": {"steam": 45, "total": 82.5}}]
   }}

Here a group is defined with three inputs and constant limits imposed on both separated steam flow and total mass flow.

.. index:: source groups; scaling
.. _group_scaling:

Group scaling
-------------

When a :ref:`group_limiter` is used to specify maximum flow rates in the output of a group, the limiter scales down the flow rates in its inputs so that the limit is not exceeded. Since a group may have many inputs, there are many possible ways to scale the input flow rates to meet the output limit.

Exactly how the inputs are scaled by a limiter is determined by the group's scaling policy, specified in the Waiwera JSON input via the group's **"scaling"** value. This is a string which can have one of two values, "uniform" (the default) or "progressive".

When **uniform scaling** is used, all inputs are scaled uniformly by the same factor. Hence, for example, if a group has several inputs and a total flow rate limit which is exceeded by a factor of 2, then all inputs would have their flow rates scaled by a factor of 0.5. Uniform scaling is appropriate, for example, for multi-feed wells, in which case if a flow rate limit is imposed on the well as a whole, then flow rates in all feedzones are scaled by the same amount.

When **progressive scaling** is used, inputs are scaled progressively, starting from the end of the input list. That is, the last input is scaled first, and if this is sufficient to meet the output limit, no other inputs are scaled. If it is not sufficient, then the last input's flow rate is set to zero, and the second-to-last input is scaled. This process continues with as many inputs as necessary until the limit is satisfied. Clearly, when progressive scaling is used, the order of the inputs is significant.

Progressive scaling can be used, for example, to simulate groups of make-up wells for a production scenario in which a (possibly time-dependent) limit on total mass flow is imposed. At any given time in the simulation, the group limiter will ensure that the right number of make-up wells are flowing to get the correct total flow. This may change during the simulation if, for example, the production wells are on deliverability and the reservoir pressures decline, so that more wells are needed to make up the total.

Example:

.. code-block:: json

   {"network": {
     "group": [{"name": "group 1",
                "in": ["rx 321", "rx 410", "rx 50"],
                "limiter": {"steam": 45, "total": 82.5},
                "scaling": "progressive"}]
   }}

Here a group is defined with three inputs, limits on steam and total flow and progressive scaling.

Note that water or steam limiters with progressive scaling should not be used for groups which have their own :ref:`group_separator`. This is because the inputs are scaled by different factors, and if they have different enthalpies, then the output enthalpy will be changed by the progressive scaling process. This would also alter the proportion of water to steam flow coming out of a group separator, resulting in a complex non-linear feedback loop. Hence, groups with water or steam limiters and their own separators (e.g. representing multi-feed wells) should use uniform scaling. Conversely, groups with water or steam limiters and progressive scaling (e.g. make-up well groups) should not have their own separator.
   
.. index:: source groups; separators, separators
.. _group_separator:

Group separator
---------------

If a source network group does not have its own separator (which is the default), then the separated water and steam flows at the group output are calculated simply by summing the separated water and steam flows from its inputs. This assumes that the inputs have their own separators, so that separated output flows have been computed for them. This is appropriate, for example, for groups of wells which each have their own individual separators.

For some types of groups, it is more appropriate for the group to have its own separator. For example, a multi-feed well would normally have a single separator taking fluid from all its feedzones (rather than individual separators for each feedzone).

In this case, a group separator can be defined in the Waiwera JSON input using the group's **"separator"** value. It is defined in exactly the same way as a separator for a source (see :ref:`source_separators`): it can be a Boolean value or an object containing a "pressure" value. Multi-stage separators can be used by specifying the separator pressure as an array.

.. note::
   **JSON object**: group separator

   **JSON path**: network.group[`index`]["separator"]

   +---------------+-----------------+--------------+---------------------+
   |**name**       |**type**         |**default**   |**value**            |
   +---------------+-----------------+--------------+---------------------+
   |"pressure"     |number | array   |0.55 MPa      |separator pressure   |
   |               |                 |              |:math:`P_0` (Pa), or |
   |               |                 |              |stage separator      |
   |               |                 |              |pressures            |
   +---------------+-----------------+--------------+---------------------+

Example:

.. code-block:: json

   {"network": {
     "group": [{"name": "w 210",
                "in": ["f 210a", "f 210b"],
                "limiter": {"steam": 12.5},
                "scaling": "uniform",
                "separator": {"pressure": [1.45e6, 0.55e6]}}]
   }}

Here a group is defined with two inputs, a steam limiter, uniform scaling and a two-stage separator.

.. index:: source groups; file output
.. _group_file_output:

File output from groups
-----------------------

Like sources, source network groups also generate output datasets (e.g. flow rates and enthalpies) in the Waiwera HDF5 output file. These datasets are contained in the same HDF5 group as the source output datasets. The output fields for network groups can be specified in the Waiwera JSON input file (see :ref:`output_fields`).

Note that, as for production sources, source groups always have negative (or zero) flow rates, as they represent fluid being removed from the model.

.. index:: source networks; reinjectors, reinjectors
.. _source_network_reinjectors:

Reinjectors
===========

A **reinjector** is in some ways the opposite of a group: instead of taking multiple input flows and combining them into a single output, it takes a single input flow and distributes it amongst multiple outputs. The input (see :ref:`reinjector_input`) can be a production source, a group or another reinjector. The outputs can be injection sources or other reinjectors.

Each reinjector output can deliver fluid from one of two categories: separated water or separated steam (not two-phase fluid). The outputs are specified in two arrays, **water** and **steam** (see :ref:`reinjector_outputs`). For geothemal reservoir models, reinjected fluid in the "steam" category is usually not in fact steam but steam condensate, after the steam from the production wells has passed through a power plant and cooling facilities.

A reinjector also has an **overflow** which handles any fluid left over after the outputs have been processed (see :ref:`reinjector_overflow`).

Reinjectors are set up in the Waiwera JSON input file via the **"network.reinject"** value. This is an array of objects. Each object in the array specifies a single reinjector.

.. note::
   **JSON object**: source reinjector

   **JSON path**: network.reinject[`index`]

   +-----------------------+----------------+------------+----------------------------------+
   |**name**               |**type**        |**default** |**value**                         |
   +-----------------------+----------------+------------+----------------------------------+
   |"name"                 |string          |""          |:ref:`reinjector_name`            |
   |                       |                |            |                                  |
   +-----------------------+----------------+------------+----------------------------------+
   |"in"                   |string          |""          |:ref:`reinjector_input`           |
   |                       |                |            |                                  |
   +-----------------------+----------------+------------+----------------------------------+
   |"water"                |array           |[]          |separated water outputs (see      |
   |                       |                |            |:ref:`reinjector_outputs`)        |
   +-----------------------+----------------+------------+----------------------------------+
   |"steam"                |array           |[]          |separated steam outputs (see      |
   |                       |                |            |:ref:`reinjector_outputs`)        |
   |                       |                |            |                                  |
   +-----------------------+----------------+------------+----------------------------------+
   |"overflow"             |string | object |{}          |:ref:`reinjector_overflow`        |
   |                       |                |            |                                  |
   +-----------------------+----------------+------------+----------------------------------+

.. index:: reinjectors; name
.. _reinjector_name:

Reinjector name
---------------

Like a group, a reinjector can be given a name using its **"name"** value, which can be an arbitrary string. This is optional, unless the reinjector is to be referenced by another reinjector (for example, if it serves as an output or overflow from another reinjector), in which case a name will be needed.

The name must also be unique, that is, no other sources, groups or reinjectors may have the same name.

.. index:: reinjectors; input
.. _reinjector_input:

Reinjector input
----------------

If a reinjector is taking its input from a production source or group, this is specified using the reinjector's  **"in"** value. This is string value containing the name of the appropriate source or group.

If the input is a source, it must have its own separator, so that input separated water and steam flow rates are available to the reinjector. Similarly, if the input is a group, it should either have its own separator, or all its inputs should have separated water and steam flows available.

If a reinjector's input comes not from a source or group but from one of the outputs or the overflow from another reinjector, this is not specified using the "in" value. This is because reinjectors can have multiple outputs, so it would be necessary to specify not only the reinjector name but also which output (or the overflow) was being used. Instead, the appropriate reinjector output or overflow specifies the receiving reinjector using its own "out" value (see :ref:`reinjector_outputs`). The receiving reinjector's "in" value does not need to be specified.

The order of reinjectors in the Waiwera JSON input file is not significant, and a reinjector output can be assigned to another reinjector that appears later in the reinjector list.

.. index:: reinjectors; outputs
.. _reinjector_outputs:

Reinjector outputs
------------------

Reinjector outputs are divided into two categories: separated water and separated steam (or steam condensate). For each category, a list of outputs is specified (see below).

Within each category, the corresponding input fluid is distributed **progressively** amongst the list of outputs - that is, the first output is assigned a flow rate, then the second output, and so on until either all the input flow has been distributed, or all the outputs have been assigned. If all the input flow is distributed to the outputs, then any remaining outputs in the list will have zero flow (and the overflow for that category will be zero). Otherwise, there will be a non-zero overflow (see :ref:`reinjector_overflow`).

Hence, the order of reinjector outputs is significant (for the same reason that the order of group inputs is significant when progressive scaling is used).

The reinjector outputs are specified in the Waiwera JSON input using two values, **"water"** for separated water outputs and **"steam"** for separated steam (or steam condensate) outputs. These are both arrays of objects. Each object represents a single separated water or steam output and its values are listed below.

.. note::
   **JSON object**: reinjector output

   **JSON path**: network.reinject[`index`]["water"][`index`] **or** network.reinject[`index`]["steam"][`index`]

   +---------------------+------------+------------+-------------------------+
   |**name**             |**type**    |**default** |**value**                |
   +---------------------+------------+------------+-------------------------+
   |"out"                |string      |""          |name of output source or |
   |                     |            |            |reinjector               |
   +---------------------+------------+------------+-------------------------+
   |"rate"               |number |    |see below   |output flow rate (kg/s)  |
   |                     |array       |            |                         |
   +---------------------+------------+------------+-------------------------+
   |"proportion"         |number |    |see below   |proportion of input flow |
   |                     |array       |            |                         |
   +---------------------+------------+------------+-------------------------+
   |"enthalpy"           |number |    |see below   |output enthalpy (J/kg)   |
   |                     |array       |            |                         |
   +---------------------+------------+------------+-------------------------+
   |"interpolation"      |string      |"linear"    |interpolation method for |
   |                     |            |            |arrays                   |
   +---------------------+------------+------------+-------------------------+
   |"averaging"          |string      |"integrate" |averaging method for     |
   |                     |            |            |arrays                   |
   +---------------------+------------+------------+-------------------------+

The **"out"** string value specifies the name of the output source or reinjector. Note that a source or reinjector cannot be specified as an output for more than one reinjector. If multiple reinjector outputs or reinjectors need to inject fluid to the same cell, individual injection sources (with the same cell) must be defined for each one.

A reinjector output can specify its flow rate using either the **"rate"** or the **"proportion"** value. The "rate" value specifies an absolute flow rate, whereas the "proportion" value specifies the output flow rate as a proportion (between zero and one) of the input flow rate for the corresponding phase (separated water or steam).

If neither a rate nor a proportion is specified, the flow rate for the output will be set to the capacity of the receiving source or reinjector, if this has been defined. If the output is an injection source, this means that a flow rate has been specified for that source, or it has a source control computing its flow rate (e.g. an :ref:`injectivity`). If the output is another reinjector, its capacity is determined by summing the capacities of its outputs. If any of its outputs are also reinjectors, then its capacity is calculated recursively.

It is also possible to specify a flow rate for both the reinjector output and its receiving source. In this case, the injection rate used will be the minimum of the two values.

If there is no flow rate specified in either the reinjector output (via the "rate" or "proportion" values) or the receiving source, then the source is treated as if it has effectively infinite capacity, and the injection rate is set to the balance of the input flow not already reinjected by previous outputs in the list. Clearly, this means any subsequent reinjector outputs in the list will have zero flow.

The **"enthalpy"** of the output can also be specified. If it is not specified, it is set equal to the enthalpy of the input, for the corresponding category (separated water or steam). In most cases the enthalpy should be specified, particularly for steam condensate outputs which will have a much lower (liquid) enthalpy than the produced steam.

The "rate", "proportion" and "enthalpy" values can all be specified either as fixed constants or rank-2 arrays of time-dependent values. If array values are used, the interpolation and averaging types can be set via the **"interpolation"** and **"averaging"** values (see :ref:`interpolation_tables`).

Examples:

.. code-block:: json

   {"source": [
     {"name": "p1", "cell": 234, "rate": -2.4, "separator": {"pressure": 5e5}},
     {"name": "p2", "cell": 275, "rate": -3.1, "separator": {"pressure": 6e5}},
     {"name": "i1", "cell": 658},
     {"name": "i2", "cell": 697}
    ],
    "network": {
       "group": [{"name": "g1", "in": ["p1", "p2"]}],
       "reinject": [{"name": "r1", "in": "g1",
                     "water": [
                       {"out": "i1", "rate": 1.5, "enthalpy": 85e3},
                       {"out": "i2", "proportion": 0.3, "enthalpy": 85e3}
                     ]}]
   }}

Here there are two production wells ("p1" and "p2") in a group called "g1", and a reinjector "r1" taking its input from the group "g1". The reinjector has two separated water outputs, going to the two injection wells, "i1" and "i2". The first output has a fixed rate of 1.5 kg/s while the second output's flow is 30\% of the input separated water flow. Both outputs have a specified enthalpy of 85 kJ/kg.

.. code-block:: json

   {"source": [
     {"name": "p1", "cell": 234, "rate": -2.4, "separator": {"pressure": 5e5}},
     {"name": "p2", "cell": 275, "rate": -3.1, "separator": {"pressure": 6e5}},
     {"name": "i1", "cell": 658},
     {"name": "i2", "cell": 697},
     {"name": "i3", "cell": 602, "injectivity": {"pressure": 6e5, "coefficient": 1e-6},
                                 "limiter": {"total": 5.5}},
     {"name": "i4", "cell": 633}
    ],
    "network": {
       "group": [{"name": "g1", "in": ["p1", "p2"]}],
       "reinject": [{"name": "r1", "in": "g1",
                     "water": [
                       {"out": "i1", "proportion": 0.2, "enthalpy": 85e3},
                       {"out": "r2", "proportion": 0.5, "enthalpy": 85e3},
                       {"out": "i2", "proportion": 0.3, "enthalpy": 85e3}
                     ]},
                    {"name": "r2",
                     "water": [
                       {"out": "i3"},
                       {"out": "i4", "rate": 2.4}
                     ]}]
   }}

Here there are again two production wells and a group "g1", which feeds into the reinjector "r1". This has three outputs, the second of which is another reinjector, "r2". Note that the reinjector "r2" does not need to specify its input, as this is defined in the corresponding output of "r1". The reinjector "r2" outputs to two injection wells, the first of which has its flow rate determined by an injectivity control. There is no flow rate defined in this reinjector output.

.. index:: reinjectors; overflow
.. _reinjector_overflow:

Reinjector overflow
-------------------

If there is more fluid entering a reinjector than can be delivered through its outputs (i.e. the input separated water or steam (condensate) flow is greater than the sum of its corresponding outputs), then there will be some left over. This is handled by the reinjector's "overflow".

Unlike :ref:`reinjector_outputs`, the reinjector overflow handles both separated water and steam overflows (i.e. there are not separate overflows for water and steam). The overflow flow rates can be monitored or post-processed via their values in the Waiwera output (see :ref:`reinjector_file_output`).

Overflows can also be directed to an injection source or another reinjector using the reinjector's **"overflow"** value. This is a string or an object with a string **"out"** value, in either case containing the name of the source or reinjector.

Example:

.. code-block:: json

   {"source": [
     {"name": "p1", "cell": 234,
      "deliverability": {"pressure": 4e5, "productivity": 1e-11},
      "separator": {"pressure": 5e5}},
     {"name": "p2", "cell": 275,
      "deliverability": {"pressure": 3e5, "productivity": 5e-12},
      "separator": {"pressure": 6e5}},
     {"name": "i1", "cell": 658},
     {"name": "i2", "cell": 697},
     {"name": "i3", "cell": 605},
     {"name": "i4", "cell": 618}
    ],
    "network": {
       "group": [{"name": "g1", "in": ["p1", "p2"]}],
       "reinject": [{"name": "r1", "in": "g1",
                     "water": [
                       {"out": "i1", "rate": 1.5, "enthalpy": 85e3},
                       {"out": "i2", "rate": 0.9, "enthalpy": 85e3}
                     ],
                     "overflow": "re2"},
                    {"name": "re2",
                     "water": [
                       {"out": "i3", "proportion": 0.5, "enthalpy": 85e3},
                       {"out": "i4", "proportion": 0.5, "enthalpy": 85e3}
                     ]}]
   }}

In this example, two production wells on deliverability feed into a group "g1" and this in turn feeds into a reinjector "r1", which has two fixed-rate outputs into the injection wells "i1" and "i2". Any water flow left over after injection into these two wells is directed to a second reinjector, "r2", which divides the overflow equally between the injection wells "i3" and "i4".

.. index:: reinjectors; file output
.. _reinjector_file_output:

File output from reinjectors
----------------------------

Like source network groups, reinjectors generate output datasets in the Waiwera HDF5 output file, in the same HDF5 group as the source output datasets. The output fields for reinjectors can be specified in the Waiwera JSON input file (see :ref:`output_fields`).

Note that, as for injection sources, reinjectors always have positive (or zero) flow rates, as they represent fluid being added to the model.
