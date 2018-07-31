.. _setup_logfile:

**********
Log output
**********

Waiwera writes a log file with diagnostic information about the simulation as it runs. This contains log messages about a variety of different conditions and events, including:

* software and compiler versions
* number of processors being used
* wall clock start and end times, and duration of the simulation
* default input values applied
* summary information of important simulation options
* diagnostic information about each time step
* diagnostic information about the non-linear solution process at each time step
* phase transitions
* time step size reductions

Note that the log file does not contain simulation results. These are written out to the main :ref:`setup_output` results file, in a different format that is more appropriate for large quantities of numerical data.

YAML format
===========

The log file is written out in `YAML <http://yaml.org/>`_ format. Somewhat like `JSON <https://www.json.org/>`_ (which is used for :ref:`waiwera_input`) this is a lightweight, text-based data interchange format, designed to be easily human-readable and also appropriate for parsing by scripts and other programs. Like JSON, YAML supports hierarchical organisation of data and representation of arrays and objects as well as simpler data types like strings and numbers.

The reason YAML is used for the log output instead of JSON lies in the way it structures and delimits the data in the file. YAML does not rely so heavily on start and end markers (e.g. brackets) to delimit data, but makes more use of structural indenting (i.e. indents are important, which is not the case in JSON).

This means that if the simulation is stopped prematurely it will not leave partly-written data structures (with missing end delimiters) in the log file which cannot be parsed.

For more details on YAML, see the `YAML website <http://yaml.org/>`_.

Log message structure
=====================

All messages written to the Waiwera log file have the same structure, with four parts:

1) **level**: the log message level, which can be either:

   * **"info"**, for messages with information on normal functioning of the simulation
   * **"warn"**, for warnings (conditions that are unusual but will not cause the simulation to stop)
   * **"err"**, for errors, which will cause the simulation to stop
2) **source**: which part of the code has generated the message
3) **event**: the event that the message is about
4) **data**: supplementary data providing additional detail on the event

Log message format
==================

Each log message is formatted in the YAML file as a single line in the following format:

.. code-block:: yaml

   - [level, source, event, data]

In YAML lines starting with a dash (-) denote array items, so the lines in the log file collectively represent an array of messages. Arrays can also be written in YAML in abbreviated one-line form inside pairs of square brackets (as in JSON), so each message also represents an array, with four items (level, source, event and data).

The "level", "source" and "event" items are simple strings. The "data" item is an object (i.e. a collection of named items) which can be written in YAML (as in JSON) inside braces, with paired names and their corresponding values (see :ref:`json_data_structures`).

In general, log messages are written to the file on consecutive lines, but blank lines are inserted after some groups of messages (e.g. messages for a particular time step) to make them more easily human-readable. These blank lines are ignored when the file is parsed using a script or other program.

Example log messages
====================

The following log message:

.. code-block:: yaml

   - [info, run, start, {"num_processors": 16}]

is an information message at the start of the run, showing that the simulation will be run on 16 processors. The "data" object in this message has only one value, "num_processors".

Once the simulation has begun, typical log output for a single time step may look something like this:

.. code-block:: yaml

   - [info, timestep, start, {"count": 14, "size": 0.819200E+09}]
   - [info, nonlinear_solver, iteration, {"count": 1, "cell": 29, "equation": 3, "residual": 0.549950E+00}]
   - [info, nonlinear_solver, iteration, {"count": 2, "cell": 29, "equation": 2, "residual": 0.847490E-01}]
   - [info, nonlinear_solver, iteration, {"count": 3, "cell": 28, "equation": 2, "residual": 0.225553E-03}]
   - [info, nonlinear_solver, end, {"iterations": 3, "converged": T, "reason": "function_relative"}]
   - [info, timestep, end, {"tries": 1, "size": 0.819200E+09, "time": 0.163830E+10, "status": "increase"}]

These are all information messages (no warnings or errors). First there is a message at the start of the timestep, showing the current timestep count and size. This is followed by three non-linear solver iteration messages, showing the iteration count, together with the size of the maximum (non-dimensionalised) residual (see :ref:`nonlinear_solution`) and the cell index and equation (i.e. component) number with the largest residual.

A cell going through a phase transition might generate a log message like this:

.. code-block:: yaml

   - [info, fluid, transition, {"cell": 109, "old_region": 1, "new_region": 4, "new_primary": [0.319733E+06, 0.100000E-05, 0.318214E+06]}]

Here the cell with index 109 is undergoing a phase transition from thermodynamic region 1 to 4 (see :ref:`thermodynamic_regions`), i.e. liquid to two-phase; in other words, the fluid in the cell is boiling. The data object in the message also gives the array of new primary variables after the phase transition. In this example the :ref:`eos` is water / CO\ :sub:`2` / energy, so the primary variables in region 4 are pressure, vapour saturation and CO\ :sub:`2` partial pressure.

In the following log messages:

.. code-block:: yaml

   - [warn, nonlinear_solver, end, {"iterations": 8, "converged": F, "reason": "max_iterations"}]
   - [warn, timestep, reduction, {"new_size": 0.838861E+10}]

the first message is a warning showing that the non-linear solver (see :ref:`nonlinear_solution`) has reached the maximum allowable number of iterations without converging, and is therefore stopping. The second message is also a warning, showing the the time step size is being reduced (see :ref:`time_step_reductions`), and the time step will be re-tried with the new size shown.

Controlling log output
======================

Log output is enabled by default, with the filename of the log file formed from the filename of the JSON input file, but with the extension changed from ".json" to ".yaml". However, log output can be controlled by setting the **"logfile"** value in the JSON input file.

The "logfile" value can take a boolean value and be used simply to turn log output on or off, for example:

.. code-block:: json

   {"logfile": false}

Alternatively, the "logfile" value can be specified as an object, with a **"filename"** string value for specifying the filename. It also has an **"echo"** boolean value  which controls whether log output is echoed to the console display as the simulation runs.

.. note::
   **JSON object**: log output

   **JSON path**: logfile

   +------------+------------+-----------------------+------------------+
   |**name**    |**type**    |**default**            |**value**         |
   +------------+------------+-----------------------+------------------+
   |"filename"  |string      |input filename with    |log filename      |
   |            |            |extension changed from |                  |
   |            |            |".json" to ".yaml"     |                  |
   +------------+------------+-----------------------+------------------+
   |"format"    |object      |{"max_num_length": 12, |number formatting |
   |            |            |"num_real_digits": 6}  |parameters        |
   +------------+------------+-----------------------+------------------+
   |"echo"      |boolean     |true                   |whether log output|
   |            |            |                       |is echoed to      |
   |            |            |                       |console           |
   +------------+------------+-----------------------+------------------+

The **"format"** object value controls the formatting of numerical data in the log output. Its **"max_num_length"** integer value specifies the maximum length (in characters) of a number, and its **"num_real_digits"** integer value specifies the number of digits after the decimal point in floating point numbers.

For example:

.. code-block:: json

   {"logfile": {"filename": "foo.yaml", "echo": false,
                "format": {"max_num_length": 14}}}

specifies log output to file "foo.yaml", without echoing log messages to console output, and with numerical values allowed to take up to 14 characters.
