.. _importing:

***********************
Importing Waiwera input
***********************

.. index:: JSON; importing, importing; JSON

Importing from a TOUGH2 input data file
=======================================

It is possible to import a Waiwera model from the input for an existing TOUGH2 [Prue04]_ model, using `PyTOUGH <https://github.com/acroucher/PyTOUGH>`_ [WeCr12]_, a library for handling TOUGH2 simulations via Python scripts.

PyTOUGH (since version 1.5.2) defines a ``t2data`` Python class for representing the contents of the main TOUGH2 input data file. This class has a ``json()`` method which returns a JSON representation of the TOUGH2 model, in the form expected by Waiwera. This JSON representation (in the form of a Python dictionary) can be written out to a Waiwera JSON input file using the ``json`` module built in to Python.

The simulation mesh for the TOUGH2 simulation must be provided in the form of a MULgrid geometry file (see [Crou18]_), before the ``json()`` method can be used. If a MULgrid geometry file is not available for the TOUGH2 model, PyTOUGH includes functionality for reverse-engineering such a geometry file from the TOUGH2 input data file, if the mesh is rectangular.

Full details of the parameters to be passed into the ``json()`` method can be found in the PyTOUGH user guide.

Limitations of the import process
=================================

* Only certain equations of state are supported: EOS1 ('W' or 'EW'), EOS2 ('EWC'), EOS3 ('EWA') and EOS4 ('EWAV').
* Some generator types (mostly AUTOUGH2 types) are not supported: 'CO2', 'DMAK', 'FEED', 'FINJ', 'HLOS', 'IMAK', 'MAKE', 'PINJ', 'POWR', 'RINJ', 'TMAK', 'TOST', 'VOL.', 'WBRE', 'WFLO', 'XINJ', 'XIN2'.
* Only some relative permeability curves are supported: linear, Pickens, Corey, Grant, fully mobile, van Genuchten, and the AUTOUGH2 tri-linear relative permeability curve (type 19).
* Only some capillary pressure functions are supported: linear, van Genuchten and zero.
* If the model has side boundary conditions applied (via horizontal connections to large volume blocks), the boundary blocks must have their centres defined, so that the appropriate normal vectors can be determined for specifying the boundary faces.

.. index:: mesh; importing, importing; mesh

Importing meshes
================

Simulation meshes in the MULgrid geometry format may be imported into ExodusII or GMSH format, suitable for Waiwera, using the PyTOUGH ``mulgrid`` object's ``write_mesh()`` method. For more details, refer to the PyTOUGH user guide [Crou18]_.

The ``write_mesh()`` method relies on the `meshio <https://pypi.org/project/meshio/>`_ Python library, so this must be installed first. The ``meshio`` library can also be used to import other types of mesh files into formats suitable for Waiwera.

Example
=======

The following Python script reads in a MULgrid geometry file, together with a TOUGH2 input data file and initial conditions file. It exports the mesh to a GMSH file and uses the ``json()`` method of the ``t2data`` class to generate a JSON representation of the simulation. An additional value is added to the JSON representation before it is written out to a JSON file.

.. code-block:: python

   # import 2-D vertical slice model from TOUGH2 to Waiwera

   from t2data import *
   from t2incons import *
   import json
   json.encoder.FLOAT_REPR = lambda o: format(o, '0.12g')

   geo = mulgrid('geometry.dat')
   mesh_filename = 'geometry.msh'
   geo.write_mesh(mesh_filename, dimension = 2, slice = 'x')

   dat = t2data('model.dat')
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
