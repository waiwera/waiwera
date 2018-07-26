.. index:: simulation; boundary conditions, boundary conditions

*******************
Boundary conditions
*******************

.. index:: boundary conditions; Dirichlet

Dirichlet boundary conditions
=============================

"Dirichlet" boundary conditions prescribe the solution (in this case, the thermodynamic :ref:`primary_variables`) on a boundary face. These can be used, for example, to apply atmospheric boundary conditions at the ground surface of a model.

In the Waiwera JSON input file, Dirichlet boundary conditions are set up via the **"boundaries"** value, which is an array of objects. Each object contains a boundary condition specification, which may apply to one or more mesh faces.

In each boundary condition, the primary variables are specified via the **"primary"** array value. The length and default value of this array both depend on which :ref:`eos` (EOS) is being used. For example, for the :ref:`water_energy_eos` EOS, there are two primary variables in the array, representing pressure and temperature for two-phase conditions, or pressure and vapour saturation for two-phase conditions.

As well as the primary variables, the thermodynamic region (see :ref:`thermodynamic_regions`) must be specified as part of the boundary condition, via the **"region"** integer value. This is needed so that the primary variables can be interpreted correctly for the desired phase conditions on the boundary. If it is not specified, a default region is used, which generally depends on the EOS (but usually corresponds to liquid conditions).

.. note::

   **JSON object**: boundary condition

   **JSON path**: boundaries[`index`]

   +------------+---------------+------------+-------------------------+
   |**name**    |**type**       |**default** |**value**                |
   +------------+---------------+------------+-------------------------+
   |"primary"   |array          |depends on  |prescribed primary       |
   |            |               |EOS         |variables                |
   |            |               |            |                         |
   +------------+---------------+------------+-------------------------+
   |"region"    |integer        |depends on  |prescribed thermodynamic |
   |            |               |EOS         |region                   |
   |            |               |            |                         |
   +------------+---------------+------------+-------------------------+
   |"faces"     |object | array |{}          |faces on which to apply  |
   |            |               |            |boundary condition       |
   |            |               |            |                         |
   +------------+---------------+------------+-------------------------+

The faces on which to apply the boundary condition are specified via the **"faces"** value. This can be either a single object or an array of objects, depending on whether the boundary condition is applied to one face specification or several different faces specifications.

In either case, the faces are specified by a combination of cell indices and an outward normal vector for the boundary. The cell indices, specified by the **"cells"** array value, are the indices of the cells just inside the boundary. Because it is an array value, multiple faces may be contained in each face specification. In general, specifying the boundary cells is not sufficient to specify the faces uniquely, because on edges of the mesh these cells may have multiple boundary faces. Hence, it is also necessary to specify the outward normal vector of the boundary, via the **"normal"** array value.

.. figure:: boundary_conditions.*
           :scale: 50 %
           :align: center

           Specifying boundary faces by cells and normal vector

.. note::

   **JSON object**: boundary condition face specification

   **JSON path**: boundaries[`index`]["faces"]

   +------------+------------+------------+-----------------------+
   |**name**    |**type**    |**default** |**value**              |
   +------------+------------+------------+-----------------------+
   |"cells"     |array       |[]          |boundary cell indices  |
   +------------+------------+------------+-----------------------+
   |"normal"    |array       |[0, 0, 1]   |outward normal vector  |
   |            |            |            |of boundary            |
   +------------+------------+------------+-----------------------+

For example:

.. code-block:: json

  {"mesh": {"filename": "my_mesh.exo"},
   "eos": {"name": "we"},
   "boundaries": [
     {"primary": [1e5, 18], "region": 1,
      "faces": {"cells": [101, 102, 103, 104], "normal": [0, 0, 1]}}
   ]}

specifies a liquid water boundary condition on four faces at the top surface of the mesh (normal vector pointing up). Because the water / energy EOS is used, which has liquid water primary variables of pressure and temperature, the boundary condition sets a pressure 1 bar and temperature 18\ :math:`^{\circ}`\ C.

In this example:

.. code-block:: json

  {"mesh": {"filename": "my_mesh.exo"},
   "eos": {"name": "wae"},
   "boundaries": [
     {"primary": [1e5, 18, 0.99e5], "region": 1,
      "faces": {"cells": [101, 102, 103, 104], "normal": [0, 0, 1]}},
     {"primary": [2e5, 0.5, 0], "region": 4,
     "faces": {"cells": [200, 201, 202, 203], "normal": [1, 0, 0]}}
   ]}

two boundaries are set up, the first again on four faces at the top surface of the mesh and with liquid conditions. Because the water, air and energy EOS is used, the region 1 boundary conditions specify pressure, temperature and air partial pressure. The second boundary sets two-phase conditions (pressure, vapour saturation and air partial pressure) with zero partial pressure of air on a horizontal boundary in the positive :math:`x`\ -direction.

.. index:: boundary conditions; Neumann
.. _neumann_boundary_conditions:

Neumann boundary conditions
===========================

"Neumann" boundary conditions prescribe the mass or energy flux through a boundary face. For example, in a geothermal reservoir model, Neumann boundary conditions may be used to specify basal mass and energy fluxes at the bottom boundary of the model.

In the finite volume framework (see :ref:`finite_volume_discretisation`), a specified flux through a boundary face (which would otherwise be zero) is formally identical to adding a source term to the cell just inside the boundary. In either case, a term is simply added to the right-hand side of the discretised conservation equations for that cell.

Hence, it is not necessary to provide a separate mechanism for implementing Neumann boundary conditions, as they can always be implemented using equivalent source terms instead (see :ref:`source_terms`).
