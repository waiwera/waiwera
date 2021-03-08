.. index:: scripting; setup
.. _setup_script:

************************************
Setting up simulations using scripts
************************************

Because many scripting languages are able to parse JSON files (see :ref:`json_script`), it is possible to write scripts to set up Waiwera JSON input files. This can simplify some aspects of simulation setup (for example, setting up boundary conditions, or large numbers of sources).

A script for setting up a simulation must usually use some form of library for interacting with the simulation mesh, for example to identify the correct cells for particular boundary conditions or sources. The appropriate mesh library depends on the type of mesh being used (see :ref:`mesh_formats`).

One option for creating meshes for use with Waiwera (see :ref:`creating_meshes`) is `Layermesh <https://github.com/acroucher/layermesh>`_, a `Python <https://www.python.org/>`_ library for creating and manipulating meshes with a layer/column structure. By importing the Layermesh module together with Python's built-in JSON module, both the mesh and all other aspects of a Waiwera simulation can be set up from a Python script. (This is similar to using the `PyTOUGH <https://github.com/acroucher/PyTOUGH>`_ library to set up a TOUGH2 simulation.)

Example
=======

The example Python script below sets up a simple steady-state non-isothermal 3-D geothermal simulation and saves it to an input JSON file. The simulation uses the default :ref:`water_energy_eos` equation of state.

A simulation mesh is created with 10 × 1 km cells in the `x`- and `y`-directions, and 25 layers of varying thicknesses. Scattered surface elevation data are fitted to the mesh to represent a spatially-varying water table. The mesh is written to a Layermesh HDF5 file (which can be used for post-processing) and an ExodusII mesh file for use with Waiwera.

Uniform initial conditions are prescribed for the start of the steady-state run, and atmospheric boundary conditions are applied at the top surface of the model. Two rock types are created, one for the reservoir and another for cap-rock near the top of the model. An upflow of hot water (enthalpy 1300 kJ/kg) is introduced around the centre of the model.

Time-stepping parameters are specified for the simulation, which stops when the time reaches 10\ :sup:`16` seconds, approximating a steady state. Simulation output for the initial state and intermediate time steps is disabled, so that the solution is only output at the end of the steady-state run.

Finally, the simulation is output to a JSON file (with a specified indent size and sorted keys).

After the simulation is run, the Layermesh library can also be used to plot the results (see :ref:`output_script`).

.. literalinclude:: demo_setup.py
   :language: Python
