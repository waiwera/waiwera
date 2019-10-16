.. index:: mesh, simulation; mesh

***************
Simulation mesh
***************

Waiwera uses a :ref:`finite_volume_discretisation` of the simulation domain to solve the mass and energy conservation equations within it. Such a discretisation strictly requires only the specification of the grid cell volumes, together with the areas of the mesh faces and their distances to the centres of the cells on either side, as these are the only mesh quantities used in the solution of the equations.

However, Waiwera requires a full **geometric mesh** to be given for the simulation, with specification of the positions of all the cell vertices, and definition of the vertex list for each cell. This is exactly the same kind of mesh definition used for e.g. finite element models. The finite volume discretisation is generated internally by Waiwera from the geometric mesh. In practice it is often more convenient to generate a geometric mesh than a pure finite volume discretisation. A geometric mesh is also necessary for post-processing of simulation results.

Mesh files
==========

The simulation mesh is not stored in the main Waiwera JSON input file. Instead, it is stored in an auxiliary file. The JSON input file only specifies the location of this auxiliary mesh file.

.. _specifying_mesh:

Specifying the mesh
===================

All mesh-related data in the Waiwera input JSON file is contained in the **"mesh"** value. This can be either:

* a string value specifying the mesh filename
* an object, with a **"filename"** string value (and other optional values as well)

Examples:

.. code-block:: json

  {"mesh": "my_mesh.msh"}

or:

.. code-block:: json

  {"mesh": {"filename": "my_mesh.msh"}}

Note that the mesh filename (specified in either of the above two ways) is a required input parameter in the simulation (there is no default).

.. note::

   **JSON object**: mesh

   **JSON path**: mesh

   +--------------------+-----------+-------------------+-----------------------------------------+
   |**name**            |**type**   |**default**        |**value**                                |
   +--------------------+-----------+-------------------+-----------------------------------------+
   |"filename"          |string     |(`no default`)     |mesh filename                            |
   |                    |           |                   |                                         |
   +--------------------+-----------+-------------------+-----------------------------------------+
   |"radial"            |boolean    |``false``          |whether 2-D mesh is radial (see          |
   |                    |           |                   |:ref:`radial_mesh`)                      |
   |                    |           |                   |                                         |
   +--------------------+-----------+-------------------+-----------------------------------------+
   |"thickness"         |number     |1.0 m              |thickness (m) for :ref:`2d_mesh`         |
   |                    |           |                   |                                         |
   +--------------------+-----------+-------------------+-----------------------------------------+
   |"permeability_angle"|number     |0\ :math:`^{\circ}`|horizontal rotation angle (\             |
   |                    |           |                   |:math:`^{\circ}`) of permeability axes   |
   |                    |           |                   |(see :ref:`rock_types`)                  |
   |                    |           |                   |                                         |
   +--------------------+-----------+-------------------+-----------------------------------------+
   |"faces"             |array      |[]                 |face properties (see                     |
   |                    |           |                   |:ref:`rock_permeability`)                |
   +--------------------+-----------+-------------------+-----------------------------------------+
   |"minc"              |object |   |[]                 |MINC treatment of fractured media (see   |
   |                    |array      |                   |:ref:`minc`)                             |
   |                    |           |                   |                                         |
   |                    |           |                   |                                         |
   +--------------------+-----------+-------------------+-----------------------------------------+
   |"rebalance"         |boolean    |``true`` for MINC, |whether to rebalance MINC mesh           |
   |                    |           |otherwise ``false``|                                         |
   +--------------------+-----------+-------------------+-----------------------------------------+
   |"zones"             |object     |{}                 |definitions of :ref:`mesh_zones`         |
   |                    |           |                   |                                         |
   +--------------------+-----------+-------------------+-----------------------------------------+

.. index:: mesh; formats
.. _mesh_formats:

Mesh formats
============

Mesh handling in Waiwera makes use of the `PETSc <https://www.mcs.anl.gov/petsc/>`_ library -- more specifically, it uses PETSc's **DMPlex** unstructured mesh functionality. As a result, Waiwera meshes may be stored in any of the file formats supported by DMPlex. Currently these formats include:

* ExodusII [ScYa95]_
* `GMSH <http://gmsh.info/>`_
* Salome `MED <http://www.salome-platform.org/user-section/about/med>`_

Meshes in these formats may be generated by a variety of stand-alone mesh generation software packages, e.g. `CUBIT <https://cubit.sandia.gov/>`_ or `GMSH <http://gmsh.info/>`_. It is also possible to convert between various mesh formats (including ExodusII, GMSH, Salome MED and VTK) using the `meshio <https://pypi.org/project/meshio/>`_ Python library. Alternatively, meshes may be imported from TOUGH2 models (see :ref:`importing`).

.. [ScYa95] Schoof, L.A. and Yarberry, V.R. (1995). "ExodusII: a finite element data model". Technical Report SAND92-2137, Sandia National Laboratories, Albuquerque, New Mexico, USA.

.. index:: mesh; coordinate systems

Mesh coordinate systems
=======================

.. index:: mesh; 3-D Cartesian

3-D Cartesian meshes
--------------------

This is the default mesh type. PETSc's DMPlex implementation supports most common 3-D element types such as 8-node hexahedral ("brick") and 4-node tetrahedral elements, but does not currently support 6-node prism (or "wedge") elements.

.. index:: mesh; 2-D Cartesian
.. _2d_mesh:

2-D Cartesian meshes
--------------------

For 2-D problems (e.g. horizontal or vertical slice models), a mesh file containing a 2-D mesh can be used. In this case, the mesh thickness can be specified by the **"mesh.thickness"** value.

For example:

.. code-block:: json

  {"mesh": {"filename": "2D.msh", "thickness": 35.0}}

If the thickness is not specified, a default value of 1.0 m is assumed.

.. index:: mesh; 2-D radial
.. _radial_mesh:

2-D radial meshes
-----------------

For radial problems, a 2-D mesh file can also be used, and the **"mesh.radial"** boolean value should be set to true.

For example:

.. code-block:: json

  {"mesh": {"filename": "cylindrical.msh", "radial": true}}

In this case, the mesh is interpreted as being in :math:`r-z` (cylindrical) coordinates.

.. index:: mesh; orthogonality
.. _mesh_orthogonality:

Mesh orthogonality
==================

Because Waiwera uses a :ref:`finite_volume_discretisation` to solve the mass and energy conservation equations, with a two-point flux approximation for evaluating the pressure and temperature gradients at the mesh faces (see :ref:`function_evaluations`), the mesh must satisfy the "orthogonality criterion", i.e. the line joining any two cell centres must be orthogonal to to face between them. In other words, the angle :math:`\theta` in :numref:`orthog_fig` should always be 90\ :math:`^{\circ}`.

.. _orthog_fig:
.. figure:: orthogonality.*
           :scale: 67 %
           :align: center

           Mesh orthogonality

Care must therefore be taken to ensure the mesh satisfies this criterion. Some simpler kinds of mesh satisfy it trivially, e.g. all regular or irregular rectangular meshes.

.. index:: mesh; partitioning
.. _mesh_partitioning:

Mesh partitioning
=================

When running Waiwera in parallel, the mesh is "partitioned" so that each parallel process contains only part of the mesh. The mesh partitioning algorithm attempts to balance the computational load between the different processes, while also making the interfaces between the partitions as small as possible, so that the minimum amount of data need be communicated between partitions during the solution process.

Waiwera uses the mesh partitioning algorithms provided by PETSc. By default, the Chaco partitioner is used.

.. index:: mesh; rebalancing
.. _mesh_rebalancing:

Mesh rebalancing
================

When MINC is used to simulate flows in fractured media (see :ref:`minc`), MINC matrix rock cells are added to the mesh. This process is carried out in parallel. If MINC is applied only to some parts of the mesh, this may result in some parallel processes having significantly more cells than others, which degrades load balancing and reduces parallel scaling performance.

If the **"mesh.rebalance"** value is set to true (the default for MINC meshes), then Waiwera will rebalance the mesh after adding MINC cells, to restore optimal load balancing.

.. index:: mesh; zones, zones
.. _mesh_zones:

Mesh zones
==========

It is possible to define named "zones" on the Waiwera mesh, to facilitate assigning different properties (e.g. rock properties) to different parts of it, without having to specify individual cells. Once a zone has been defined, rock or other properties can be defined on all cells in the zone simply by referring to the appropriate zone name.

Zones can be re-used for different purposes, and some types of zone specification are purely geometrical (not relying on cell indexing) and can therefore be re-used for different meshes.

Mesh zones are defined in the **"mesh.zones"** value in the Waiwera JSON input file. This value in an object, containing pairs of zone names and their corresponding zone definitions.

The available types of zones are as follows.

.. index:: zones; cell array

Cell array type
---------------

In this type of zone, cells in the zone are explicitly identified by their cell indices. In this case the zone definition can be either:

* an array of integers (the cell indices)
* an object with a **"cells"** value, containing the array of cell indices

All cell indices are zero-based (i.e. start from zero) and refer to the cell indices in the serial mesh (i.e. before partitioning, if the simulation is run in parallel).

If the zone definition is an object, it can optionally also contain a **"type"** string value, set to "array", to make the zone type more explicit.

.. note::

   **JSON object**: cell array zone

   **JSON path**: mesh.zones[`name`]

   +----------+----------+-----------+--------------+
   |**name**  |**type**  |**default**|**value**     |
   +----------+----------+-----------+--------------+
   |"cells"   |array     |[]         |cell indices  |
   +----------+----------+-----------+--------------+
   |"type"    |string    |"array"    |zone type     |
   +----------+----------+-----------+--------------+

Example:

.. code-block:: json

  {"mesh": {"filename": "model.msh",
            "zones": {"bottom": [0, 4, 5, 6, 7],
                      "cap": {"cells": [10, 12, 14]},
                      "feed": {"type": "array", "cells": [100, 253, 342]}
                      }
   }}

Here three zones are defined. The "bottom" zone is defined as an array of cell indices. The "cap" zone is defined as an object containing a "cells" array value. The "feed" zone is also an object, with its "type" value explicitly set to "array".

.. index:: zones; box

Box type
--------

In this type of zone, a "box" is defined by coordinate ranges, and only the cells with centres inside that box belong to the zone.

The zone definition is an object, with between one and three coordinate ranges named **"x"**, **"y"** and **"z"** (or **"r"** for radial meshes). Each coordinate range is a two-element array of numbers containing the minimum and maximum of the range.

If any of these coordinate ranges is absent, there is assumed to be no limitation on that coordinate.

The zone definition can optionally also contain a **"type"** string value, set to "box", to make the zone type more explicit.

.. note::

   **JSON object**: box zone

   **JSON path**: mesh.zones[`name`]

   +----------+----------+-----------+---------------------+
   |**name**  |**type**  |**default**|**value**            |
   +----------+----------+-----------+---------------------+
   |"x"       |array     |[]         |:math:`x`-coordinate |
   |          |          |           |range                |
   |          |          |           |                     |
   +----------+----------+-----------+---------------------+
   |"y"       |array     |[]         |:math:`y`-coordinate |
   |          |          |           |range                |
   +----------+----------+-----------+---------------------+
   |"z"       |array     |[]         |:math:`z`-coordinate |
   |          |          |           |range                |
   +----------+----------+-----------+---------------------+
   |"r"       |array     |[]         |:math:`r`-coordinate |
   |          |          |           |range                |
   +----------+----------+-----------+---------------------+
   |"type"    |string    |"box"      |zone type            |
   +----------+----------+-----------+---------------------+

Examples:

.. code-block:: json

  {"mesh": {"filename": "model.msh",
            "zones": {"production": {"x": [0, 1000],
                                     "y": [500, 900],
                                     "z": [-1000, -200]},
                                     "NE": {"type": "box",
                                            "x": [1000, 4000],
                                            "y": [1200, 3300]},
                      "basement": {"z": [-3000, -2000]}}
   }}

Here the "production" zone is defined by coordinate ranges in all three coordinates. The "NE" zone (explicitly given the "box" type) restricts only the :math:`x` and :math:`y` coordinates, with no limit on :math:`z`. The "basement" zone contains all cells in the model with elevations between -3000 and -2000 m.

.. code-block:: json

  {"mesh": {"filename": "welltest.msh", "radial": true,
            "zones": {"skin": {"r": [0.1, 0.6]}}
   }}

Here a "skin" zone is defined in a radial model, containing all cells with centre radii between 0.1 and 0.6 m of the origin.

A zone covering the entire mesh can be constructed by specifying a box zone with no coordinate ranges, as in the zone named "all" below:

.. code-block:: json

  {"mesh": {"filename": "model.msh",
            "zones": {"all": {"type": "box"}}
   }}

.. index:: zones; combining

Combination type
----------------

This type forms a new zone by combining other zones together. Zones may be combined by:

* "adding" zones together (union)
* "subtracting" zones from one another (complement)
* "multiplying" zones together (intersection)

The zone definition is an object, with between one and three values named **"+"**, **"-"** and **"*"**, corresponding to the zone combination operations listed above. Each value can be either a single zone name (string), an array of zone names, or ``null``.

The zone definition can optionally also contain a **"type"** string value, set to "combine", to make the zone type more explicit.

.. note::

   **JSON object**: combination zone

   **JSON path**: mesh.zones[`name`]

   +----------+----------------+-----------+-------------+
   |**name**  |**type**        |**default**|**value**    |
   |          |                |           |             |
   +----------+----------------+-----------+-------------+
   |"+"       |string | array  |[]         |zones to add |
   |          |                |           |             |
   +----------+----------------+-----------+-------------+
   |"-"       |string | array  |[]         |zones to     |
   |          |                |           |subtract     |
   +----------+----------------+-----------+-------------+
   |"*"       |string | array  |[]         |zones to     |
   |          |                |           |multiply     |
   +----------+----------------+-----------+-------------+
   |"type"    |string          |"combine"  |zone type    |
   +----------+----------------+-----------+-------------+


Combination zones do not need to be defined in any particular order with respect to the other zones. They may refer to zones defined further down in the Waiwera JSON input file.

Examples:

.. code-block:: json

  {"mesh": {"filename": "model.msh",
            "zones": {"NE": {"x": [1000, 4000], "y": [1200, 3300]},
                      "NW": {"x": [0, 1000], "y": [1200, 3300]},
                      "N": {"+": ["NE", "NW"]},
                      "N basement": {"*": ["N", "basement"]},
                      "basement": {"z": [-3000, -2000]}}
            }}

Here the "NE" and "NW" zones are added to produce an "N" zone, the union of the two. The "N basement" zone consists of all cells in the "N" zone with elevations between -3000 and -2000 m.

In the next example, the "production matrix" zone consists of all cells in the "production" zone but not in either of the "fracture1" or "fracture2" zones:

.. code-block:: json

  {"mesh": {"filename": "model.msh",
            "zones": {"production": {"x": [0, 1000], "y": [0, 1000]},
                      "fracture1": {"z": [-2000, -1900]},
                      "fracture2": {"z": [-2400, -2350]},
                      "production matrix": {"+": "production",
                                            "-": ["fracture1", "fracture2"]}
                      }
            }}
            
A zone covering the entire mesh can be constructed as the complement of a ``null`` zone, as in the zone named "all" below:

.. code-block:: json

  {"mesh": {"filename": "model.msh",
            "zones": {"all": {"-": null}}
            }}
