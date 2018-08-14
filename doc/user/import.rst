.. _importing:

***********************
Importing Waiwera input
***********************

.. index:: JSON; importing, importing; JSON

Importing from a TOUGH2 input data file
=======================================

It is possible to import a Waiwera model from the input for an existing TOUGH2 [Prue04]_ model. A small Python module is included in the Waiwera ``utils/`` directory to automate this process. The module, ``t2data_json.py``, makes use of the `PyTOUGH <https://github.com/acroucher/PyTOUGH>`_ [WeCr12]_ library for handling TOUGH2 simulations via Python scripts, so this must be installed before the import module can be used. (Eventually it is planned to incorporate the ``t2data_json.py`` module into PyTOUGH itself.)

Also, the ``t2data_json.py`` module needs to be on the Python path of your system, so Python can find it. The ``utils/`` directory can be added to the Python path, or a copy of or link to the module can be added to a directory already on your Python path.

The module defines a class ``t2data_export_json``, which is a subclass of PyTOUGH's ``t2data`` class for representing the main TOUGH2 input data file. As well as all the usual functionality of the ``t2data`` class, the extended class ``t2data_export_json`` has an additional ``json()`` method which returns a JSON representation of the TOUGH2 model, in the form expected by Waiwera. This JSON representation (in the form of a Python dictionary) can be written out to a Waiwera JSON input file using the ``json`` module built in to Python.

The simulation mesh for the TOUGH2 simulation must be provided in the form of a MULgrid geometry file (see [Crou18]_), before the ``t2data_json`` module can be used. If a MULgrid geometry file is not available for the TOUGH2 model, PyTOUGH includes functionality for reverse-engineering such a geometry file from the TOUGH2 input data file, if the mesh is rectangular.

The ``json()`` method takes the following parameters (most of which are optional and have their defaults shown):

.. code-block:: python

   json(geo, mesh_filename, atmos_volume = 1.e25, incons = None,
             eos = None, bdy_incons = None, mesh_coords = 'xyz')

* **geo**: ``mulgrid``

  The ``mulgrid`` geometry object (defined in PyTOUGH's ``mulgrids`` module) representing the simulation mesh.
* **mesh_filename**: string

  The filename of the mesh file (e.g. ExodusII or GMSH mesh --  see :ref:`mesh_formats`) for the Waiwera simulation.

* **atmos_volume**: float

  Volume above which blocks in the TOUGH2 model are identified as boundary condition (e.g. atmosphere) blocks, rather than being part of the simulation mesh.

* **incons**: string | ``t2incon`` | ``None``

  Initial conditions for the TOUGH2 model. If specified as a string, this should be the filename of the initial conditions file; if a ``t2incon`` object, this should be the contents of such a file, represented by PyTOUGH's ``t2incon`` class. If ``None`` is specified, then default initial conditions will be applied from the TOUGH2 model's "PARAM" section.

* **eos**: string | integer | ``None``

  The equation of state used for the TOUGH2 simulation. For AUTOUGH2 simulations, this can generally be set to ``None``, and the EOS will be read from the "SIMUL" or "MULTI" sections in the TOUGH2 input file. Otherwise, it can be specified as an integer corresponding to the EOS number (1 being pure water, 2 being water / CO\ :sub:`2` etc.) or as a string corresponding to the AUTOUGH2 EOS names (EOS1 being 'EW', EOS2 being 'EWC' etc.). Note that only EOS modules 1 -- 4 (i.e. 'W', 'EW', 'EWC', 'EWA' and 'EWAV' in terms of AUTOUGH2 EOS names) are supported.

* **bdy_incons**: ``t2incon`` | ``None``

  Initial conditions from which boundary conditions are to be derived. Usually this will be the same as the ``t2incon`` object specifying the initial conditions (via the ``incons`` parameter). If set to ``None``, then default boundary conditions will be applied from the TOUGH2 model's default initial conditions in the "PARAM" section. Faces on which to apply boundary conditions are identified by the presence of connections to blocks with either zero or large volume (above the volume specified by the ``atmos_volume`` parameter).

* **mesh_coords**: string

  String representing the coordinate system to be used in the Waiwera model. 3-D Cartesian meshes are identified as 'xyz'. 2-D Cartesian meshes may be identified as either 'xy', 'xz', or 'yz' (depending on orientation), while 2-D radial meshes are identified as 'rz'.

Limitations of the import process
=================================

* Only certain equations of state are supported: EOS1 ('W' or 'EW'), EOS2 ('EWC'), EOS3 ('EWA') and EOS4 ('EWAV').
* Some generator types (mostly AUTOUGH2 types) are not supported: 'CO2', 'DMAK', 'FEED', 'FINJ', 'HLOS', 'IMAK', 'MAKE', 'PINJ', 'POWR', 'RINJ', 'TMAK', 'TOST', 'VOL.', 'WBRE', 'WFLO', 'XINJ', 'XIN2'.
* Only some relative permeability curves are supported: linear, Pickens, Corey, Grant, fully mobile, van Genuchten, and the AUTOUGH2 tri-linear relative permeability curve (type 19).
* Only some capillary pressure functions are supported: linear, van Genuchten and zero.

.. index:: mesh; importing, importing; mesh

Importing meshes
================

Simulation meshes in the MULgrid geometry format may be imported into ExodusII or GMSH format, suitable for Waiwera, using the PyTOUGH ``mulgrid`` object's ``write_mesh()`` method. For more details, refer to the PyTOUGH user guide [Crou18]_.

The ``write_mesh()`` method relies on the `meshio <https://pypi.org/project/meshio/>`_ Python library, so this must be installed first. The ``meshio`` library can also be used to import other types of mesh files into formats suitable for Waiwera.

Example
=======

The following Python script reads in a MULgrid geometry file, together with a TOUGH2 input data file and initial conditions file. It exports the mesh to a GMSH file and uses the ``json()`` method of the ``t2data_export_json`` class to generate a JSON representation of the simulation. An additional value is added to the JSON representation before it is written out to a JSON file.

.. code-block:: python

   # import 2-D vertical slice model from TOUGH2 to Waiwera

   from t2data_json import *
   from t2incons import *
   import json
   json.encoder.FLOAT_REPR = lambda o: format(o, '0.12g')

   geo = mulgrid('geometry.dat')
   mesh_filename = 'geometry.msh'
   geo.write_mesh(mesh_filename, dimension = 2, slice = 'x')

   dat = t2data_export_json('model.dat')
   inc = t2incon('model.incon')

   jsondata = dat.json(geo, mesh_filename, mesh_coords = 'xz', incons = inc,
                        bdy_incons = inc)
   jsondata['mesh']['thickness'] = 10

   json.dump(jsondata, file('model.json', 'w'), indent = 2, sort_keys = True)

After the necessary Python modules have been imported at the start of the script, a line is added to change the way floating point values are written out by the ``json`` module, essentially suppressing unnecessary digits after the decimal point. This is optional but can make the resulting JSON files easier to read.

In this example the mesh is a 2-D Cartesian vertical slice mesh, so the ``dimension`` parameter of the ``write_mesh()`` method is set to 2, and the ``slice`` parameter is set to 'x'. The ``mesh_coords`` parameter is also set to 'xz' in the ``json()`` method call, to define the orientation of the mesh.

After the ``jsondata`` dictionary variable has been created, containing a JSON representation of the simulation, this variable is edited to add a value for the 2-D mesh thickness. Other aspects of the simulation could be altered as desired here by editing the ``jsondata`` variable. 

The last line of the script writes the Waiwera input JSON file, "model.json". Two optional parameters are added here, one to specify the size of the indenting, and the other to sort the keys in the JSON file. Since the ``jsondata`` variable is a Python dictionary, and dictionary variables have no implied order, the keys could otherwise be written out in arbitrary order. Setting the ``sort_keys`` parameter ensures the keys are written out in alphabetical order, which can make particular keys easier to find in the file.

.. [Crou18] Croucher, A. (2018). "PyTOUGH user's guide", version 1.5.1, University of Auckland, July 2018. Available from the `PyTOUGH website <https://github.com/acroucher/PyTOUGH>`_.
.. [Prue04] Pruess, K. (2004). "The TOUGH codes -- a family of simulation tools for multiphase flow and transport processes in permeable media". Vadose Zone Journal, 3(3), 738 -- 746.
.. [WeCr12] Wellmann, J.F., Croucher, A.E. and Regenauer-Lieb, K. (2012). "Python scripting libraries for subsurface fluid and heat flow simulations with TOUGH2 and SHEMAT". Computers & Geosciences, 43 (197-206).
