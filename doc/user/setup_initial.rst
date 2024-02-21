.. index:: simulation; initial conditions, initial conditions

.. _initial_conditions:

******************
Initial conditions
******************

Before the simulation can run, the :ref:`primary_variables` in each cell at the start time must be specified. In addition, the :ref:`thermodynamic_regions` in the cells must also be specified, so that the initial primary variables can be interpreted correctly.

There are two ways of doing this. However, both ways specify the initial conditions via the **"initial"** value in the Waiwera JSON input file, and in both cases this value is an object (although some of its values are different in each case).

.. admonition:: JSON input

   **JSON object**: initial conditions

   **JSON path**: initial

   +------------+----------------+----------------+-------------------------------+
   |**name**    |**type**        |**default**     |**value**                      |
   +------------+----------------+----------------+-------------------------------+
   |"primary"   |number | array  |(`no default`)  |initial                        |
   |            |                |                |:ref:`primary_variables`       |
   |            |                |                |                               |
   +------------+----------------+----------------+-------------------------------+
   |"region"    |integer | array |depends on EOS  |initial                        |
   |            |                |                |:ref:`thermodynamic_regions`   |
   +------------+----------------+----------------+-------------------------------+
   |"filename"  |string          | (`no default`) |name of Waiwera output HDF5    |
   |            |                |                |file to restart from           |
   +------------+----------------+----------------+-------------------------------+
   |"index"     |integer         |0               |results index (zero-based)     |
   |            |                |                |in file to restart from        |
   +------------+----------------+----------------+-------------------------------+
   |"minc"      |boolean         |false           |whether initial conditions     |
   |            |                |                |include values for MINC cells  |
   |            |                |                |(see                           |
   |            |                |                |:ref:`minc_initial_conditions`)|
   +------------+----------------+----------------+-------------------------------+
   |"tracer"    |number | array  |0               |:ref:`initial_tracer`          |
   +------------+----------------+----------------+-------------------------------+

.. index:: initial conditions; in JSON input file

Specifying initial conditions in the JSON input file
----------------------------------------------------

It is possible to specify initial conditions directly in the Waiwera JSON input file. Normally this is practical only for relatively small models.

In this case the primary variables for each cell are specified via the **"primary"** value in the "initial" object. This value can be either:

* a single number, which is assigned to initial conditions in all cells (this only really makes sense for the "w" EOS (:ref:`eos`), as all others have more than one primary variable per cell)
* a rank-1 array of numbers, of length equal to the number of primary variables per cell (which depends on the EOS), representing a set of uniform primary variables to be assigned to every cell
* a rank-2 array of numbers, one row for each cell (in cell index order) and one column for each primary variable, representing the full specification of all primary variables in all cells

The thermodynamic region for each cell is specified via the **"region"** value. This value can be either:

* a single integer, assigned to all cells
* an array of integers, one value per cell (in cell index order)

If the "region" value is not specified, the default thermodynamic region for the EOS is used.

For example:

.. code-block:: json

   {"eos": {"name": "w", "temperature": 18.0},
    "initial": {"primary": 1.5e5}
   }

specifies a uniform initial pressure of 1.5 bar in an isothermal (18 :math:`^{\circ}`\ C) water simulation. (Here there is no need to specify the thermodynamic region, as the "w" EOS is limited to liquid water anyway.) In the example below, the same initial conditions are specified in a non-isothermal water / energy simulation:

.. code-block:: json

   {"eos": {"name": "we"},
    "initial": {"primary": [1.5e5, 18.0], "region": 1}
   }

In the next example, different primary variables are specified in each cell (in a very small four-cell model), but the thermodynamic region is uniform (liquid water):

.. code-block:: json

   {"eos": {"name": "we"},
    "initial": {"primary": [[1.50e5, 18.0],
                            [1.46e5, 18.2],
                            [1.43e5, 18.3],
                            [1.41e5, 18.6]],
    "region": 1}
   }

The next example demonstrates different primary variables and regions being assigned to the cells, the first two cells being two-phase and the other two liquid:

.. code-block:: json

   {"eos": {"name": "we"},
    "initial": {"primary": [[1.01e5, 0.5],
                            [1.01e5, 0.1],
                            [1.1e5,  100],
                            [1.2e5,   98]],
    "region": [4, 4, 1, 1]}
   }

.. index:: initial conditions; restarting
.. _restarting:

Restarting from a previous output file
--------------------------------------

Waiwera saves its results to an HDF5 output file (see :ref:`setup_output`), and a new simulation can be restarted directly from the output of a previous run, using it as initial conditions.

In this case, the "initial" object in the JSON input file for the restarted simulation takes a different form.

.. index::
   pair: MINC; initial conditions

The filename of the output from the previous simulation is specified using the **"filename"** value. In general, an output file may contain results for more than one time. The new simulation can be restarted from any of the results in the previous output file. The index of the desired set of results can be specified using the **"index"** value, which defaults to zero.

Restarting from a previous output file will read both the primary variables and the thermodynamic regions from the file. Clearly, the output file should contain results for the same number of cells as the restarted simulation (except in the case of restarting a MINC simulation from single-porosity initial conditions -- see :ref:`minc_initial_conditions`).

Initial conditions read from an HDF5 output file via the **filename** value will override any initial conditions specified directly in the JSON input file (via the **primary** or **region** values).

Generally the previous output file should have been generated using the same :ref:`eos` used by the restarted simulation. However, this is not strictly necessary, as long as the output file contains results for all the primary variables of the EOS used in the restarted run.

For example:

.. code-block:: json

   {"initial": {"filename": "previous_run.h5", "index": 99}}

restarts a simulation from a Waiwera HDF5 output file named "previous_run.h5", starting from the set of results in the file with zero-based index 99.

Note that if *both* the **"primary"** and **"filename"** values are specified, the simulation will be restarted from the results in the specified previous output file, and the **"primary"** values will be ignored.

.. index:: initial conditions; default

Default initial conditions
--------------------------

If no **"initial"** value is present in the Waiwera JSON input file, default initial conditions will be assigned. A warning message to that effect will be written to the logfile (see :ref:`setup_logfile`).

In this case, the default primary variables (and thermodynamic region) for the :ref:`eos` being used will be assigned to all cells.

.. index:: initial conditions; tracer, tracers; initial conditions
.. _initial_tracer:

Tracer initial conditions
-------------------------

If tracers are being simulated (see :ref:`setup_tracers`) then it is possible to set initial conditions for the tracer mass fractions, via the **"tracer"** value. This can be either a scalar, to be applied to all tracers defined in the simulation, or an array, with one value for each tracer. In either case, the specified tracer initial conditions will be applied over all cells in the mesh (note there is no provision for specifying different tracer initial conditions in different cells via the **"tracer"** value, as this is usually not needed).

If no initial conditions for tracer are specified, then all tracer mass fractions will be initialised to zero. For most problems, this is what is desired.

Tracer initial conditions can be specified via the **"tracer"** value regardless of whether the fluid initial conditions are specified using the **"primary"** and **"region"** values, or by :ref:`restarting`. If restarting from an output file which does not contain any tracer results (e.g. a steady-state solution), initial conditions from the **"tracer"** value will be applied. However, if restarting from an output file which does contain tracer results, these will override the **"tracer"** value.
