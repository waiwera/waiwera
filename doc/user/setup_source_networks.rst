.. index:: simulation; source networks, source networks
.. _source_networks:

***************
Source networks
***************

In many simulations, some sources (see :ref:`source_terms`) may not act entirely independently, but instead may interact in some way. Examples of this include:

- multi-feed wells, with sources for individual feedzones grouped together into a well, which behaves in many ways as a single unit
- borefields, with production wells grouped together and limits or targets imposed on the group
- reinjection, with some of the output from a group of production wells distributed amongst a set of injection wells

To model the interactions between sources in a Waiwera simulation, a **source network** may be defined. A source network has two parts:

- **groups**, defining groupings of production sources
- **reinjectors**, defining the distribution of produced fluid amongst injection sources

In the Waiwera JSON input file a source network may be set up using the **network** value. This is an object containing two array values, **"group"** and **"reinject"**, defining source groups and reinjectors respectively.

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

.. index:: source networks; groups, source groups
.. _source_network_groups:

Groups
======

A **group** takes multiple inputs and combines them into a single output. The inputs may be production sources (i.e. with negative rate), or other groups, enabling a nested hierarchy of sub-groups and groups to be defined if needed.

The mass flow rate at the output is simply the sum of the mass flow rates from the inputs. The enthalpy at the output is a weighted sum of the input enthalpies (weighted by mass flow rate). Hence, the total input and output mass and energy flow rates are equal.

A group also calculates output mass flow rates and enthalpies for separated water and steam (i.e. liquid and vapour phases) from two-phase inputs. How this is done depends on whether the group has its own separator or not  (see :ref:`group_separators`).

Source groups are set up in the Waiwera JSON input file via the **"network.group"** value. This is an array of objects. Each object in the array specifies a single group.

.. note::
   **JSON object**: source group

   **JSON path**: network.group[`index`]

   +-----------------------+----------------+------------+-------------------------+
   |**name**               |**type**        |**default** |**value**                |
   +-----------------------+----------------+------------+-------------------------+
   |"name"                 |string          |""          |optional group name      |
   |                       |                |            |                         |
   +-----------------------+----------------+------------+-------------------------+
   |"in"                   |array           |[]          |names of inputs          |
   |                       |                |            |                         |
   +-----------------------+----------------+------------+-------------------------+
   |"limiter"              |object          |{}          |                         |
   |                       |                |            |                         |
   +-----------------------+----------------+------------+-------------------------+
   |"scaling"              |string          |"uniform"   |scaling type             |
   |                       |                |            |                         |
   |                       |                |            |                         |
   +-----------------------+----------------+------------+-------------------------+
   |"separator"            |boolean | object|``false``   |separator                |
   |                       |                |            |                         |
   +-----------------------+----------------+------------+-------------------------+


.. index:: source groups; name
.. _group_name:

Name
----

.. index:: source groups; inputs
.. _group_inputs:

Inputs
------


.. index:: source groups; limiter
.. _group_limiter:

Limiters
--------

.. index:: source groups; scaling
.. _group_scaling:

Scaling
-------

.. index:: source groups; separator, separators
.. _group_separators:

Separators
----------

.. index:: source groups; output
.. _group_output:

Output from groups
------------------

.. index:: source networks; reinjectors
.. _source_network_reinjectors:

Reinjectors
===========

