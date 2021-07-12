.. index:: tracers
.. _setup_tracers:

*******
Tracers
*******

Parameters related to the configuration of tracers (see :ref:`simulating_tracers`) are specified in the Waiwera JSON input file via the **"tracer"** value. This value is an object, and may contain the values shown below. Multiple tracers can be simulated simultaneously by specifying the **tracer** value as an array of these objects.

.. note::

   **JSON object**: tracer parameters

   **JSON path**: tracer (or tracer[index], if there are multiple tracers)

   +-----------------------+--------------+----------------------+-----------------------+
   |**name**               |**type**      |**default**           |**value**              |
   +-----------------------+--------------+----------------------+-----------------------+
   |"name"                 |string        |"tracer", or          |tracer name            |
   |                       |              |"tracer_i" where i is |                       |
   |                       |              |the tracer index      |                       |
   +-----------------------+--------------+----------------------+-----------------------+
   |"phase"                |string        |depends on EOS, but   |tracer phase           |
   |                       |              |usually "liquid"      |                       |
   +-----------------------+--------------+----------------------+-----------------------+
   |"diffusion"            |number        |0                     |diffusion coefficient  |
   |                       |              |                      |(m\ :sup:`2`/s)        |
   +-----------------------+--------------+----------------------+-----------------------+
   |"decay"                |number        |0                     |constant decay rate    |
   |                       |              |                      |(1/s)                  |
   +-----------------------+--------------+----------------------+-----------------------+
   |"activation"           |number        |0                     |activation energy      |
   |                       |              |                      |(J/mol)                |
   +-----------------------+--------------+----------------------+-----------------------+

The **"name"** value is a string identifying each tracer. If it is not specified, the default name "tracer" is used if there is a single tracer specified. If an array of tracers is specified, the default name for each is  based on its (zero-based) index, e.g. "tracer_0", "tracer_1", etc.

Each tracer is assumed to be specific to a particular fluid phase (see :ref:`tracer_terms`). This phase can be specified via the **"phase"** string value. If not specified, a default value is used, which is specific to the :ref:`eos`. However, for most equations of state, the default tracer phase is "liquid".

As well as being advected by the movement of fluid, tracers may also be transported by diffusion, which moves tracer from areas of high concentration to areas of low concentration. The strength of the diffusion effect is determined by a coefficient defined via the **"diffusion"** value. For the details of how this diffusion coefficient is defined and used, see :ref:`tracer_terms`.

Each tracer can be given a decay rate (see :ref:`tracer_eqns`), which can be either constant or temperature-dependent (see :ref:`tracer_temp_decay`). A constant decay rate can be set directly via the **"decay"** value.

Temperature-dependent decay is introduced by specifying the activation energy for the tracer (:math:`E_0` in equation :eq:`arrhenius`), via the **"activation"** value. In this case the "decay" value specifies the constant decay rate :math:`\alpha^0` in the same equation.

Examples
========

.. code-block:: json

  {"eos": "we",
   "tracer": {"name": "T1"}
  }

Here a single non-decaying liquid-phase tracer named "T1" is defined. The phase is not specified, so it is determined by the default tracer phase for the "we" equation of state, which is "liquid".

The following example defines a single non-decaying vapour phase tracer, with the phase specified explicitly:

.. code-block:: json

  {"tracer": {"name": "vtracer", "phase": "vapour"}
  }

Here a single liquid-phase tracer is defined, with a constant decay rate of :math:`10^{-6} s^{-1}`:

.. code-block:: json

  {"tracer": {"name": "T2", "phase": "liquid", "decay": 1e-6}
  }

The example below defines a single liquid-phase tracer, with a temperature-dependent decay rate determined by the activation energy :math:`E_0` = 2 kJ/mol:

.. code-block:: json

  {"tracer": {"name": "T3", "phase": "liquid", "decay": 1e-6, "activation": 2e3}
  }

This example defines three tracers with various properties:

.. code-block:: json

  {"tracer": [
              {"name": "T1", "phase": "liquid", "diffusion": 1e-6},
              {"name": "T2", "phase": "vapour", "decay": 2e-7, "diffusion": 1.5e-6},
              {"name": "T3", "phase": "liquid", "decay": 1e-6, "activation": 1850}
             ]
  }
