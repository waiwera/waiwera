.. index:: simulation; gravity, gravity

*******
Gravity
*******

The face fluxes in the mass and energy balance equations (see :ref:`function_evaluations`) include a term involving the gravity vector (:math:`\mathbf{g}`). The magnitude (and direction) of this acceleration due to gravity for the simulation can be specified in the Waiwera JSON input file via the **"gravity"** value.

.. note::

   +-----------------+-----------------------------------------------------------------------+
   | **JSON value**: | "gravity"                                                             |
   +-----------------+-----------------------------------------------------------------------+
   | **type**:       | array | number | ``null``                                             |
   +-----------------+-----------------------------------------------------------------------+
   | **default**:    | [0, 0] (2-D meshes) | [0, 0, -9.8] m/s\ :superscript:`2` (3-D meshes) |
   +-----------------+-----------------------------------------------------------------------+
   | **specifies**:  | acceleration due to gravity (m/s\ :superscript:`2`)                   |
   +-----------------+-----------------------------------------------------------------------+

The gravity vector is a 2-D vector for 2-D simulations, and a 3-D vector for 3-D simulations. The entire vector can be specified in full, by using an array value for gravity in the JSON input.

For example:

.. code-block:: json

  {"gravity": [0, 0, -9.81]}

In many cases, however, the gravity vector is aligned with the last coordinate axis of the mesh, but in the opposite direction (e.g. in the :math:`-y` direction for 2-D vertical slice meshes, and in the :math:`-z` direction for 3-D meshes). In this case only the `magnitude` of the gravity vector need be specified, by using a number value for gravity in the JSON input.

For example:

.. code-block:: json

  {"gravity": 9.81}

sets the gravity vector to [0, -9.81] m/s\ :superscript:`2` for a 2-D simulation, or [0, 0, -9.81] m/s\ :superscript:`2` for a 3-D simulation.

Not specifying any gravity value, or specifying it as ``null`` causes a default value to be used. For 2-D simulations the default is [0, 0], i.e. no gravity, which is appropriate for horizontal layer models. For 3-D simulations the default is [0, 0, -9.8] m/s\ :superscript:`2`.
