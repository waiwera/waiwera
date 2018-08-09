.. index:: MINC, fractured media, mesh; MINC
.. _minc:

************************************
Modelling fractured media using MINC
************************************

When highly fractured rocks are simulated, the assumption that the rock acts purely as a porous medium within each cell may no longer hold. There may be significant flow along fractures, combined with storage of fluid within the matrix rock between fractures.

One way of simulating this kind of fracture / matrix flow is by using the method known as Multiple INteracting Continua, or MINC [PrNa85]_. In a MINC model, the mesh is altered so that each cell in a fractured medium now represents the fracture within the cell, and one or more additional cells are added "inside" this cell to represent the rock matrix.

If the matrix is represented by just a single cell inside each fracture cell, this corresponds to the "dual porosity" approximation. In general, the matrix is represented by multiple cells nested inside each other. Flow is permitted between:

* adjacent fracture cells
* each fracture cell and the first matrix cell inside it
* each matrix cell and the next matrix cell inside itself

but not between matrix cells that are inside different fracture cells. Flow in the rock matrix is locally one-dimensional. Generally permeabilities are low in the matrix cells, but high in the fracture cells.

.. figure:: minc.*
           :scale: 67 %
           :align: center

           Schematic of MINC cell connectivity (fracture cells and three-level nested matrix cells)

Once the additional MINC matrix cells have been added to the mesh, they are treated exactly the same as any other cells, in terms of how the mass and energy conservation equations are solved in them.

Note that it is not necessary to add the additional MINC cells to the Waiwera simulation mesh explicitly. The original single-porosity mesh is still used, and Waiwera will add the MINC cells internally to the PETSc DMPlex mesh (see :ref:`mesh_formats`). The only change to the input for a MINC model lies in the specification of the MINC parameters.

MINC parameters may be specified in the Waiwera JSON input file via the **"mesh.minc"** value. This is an object with two values, **"geometry"** and **"rock"**. Note that is also possible to set up multiple MINC zones in the mesh, each with different MINC parameters, simply by specifying the "mesh.minc" value as an array of these objects.

.. note::
   **JSON object**: MINC parameters

   **JSON path**: mesh.minc (or mesh.minc[`index`], if there are multiple MINC zones)

   +-----------+---------------+----------------+---------------------------+
   |**name**   |**type**       |**default**     |**value**                  |
   +-----------+---------------+----------------+---------------------------+
   |"geometry" |object         |see below       |:ref:`minc_geometry`       |
   |           |               |                |parameters                 |
   +-----------+---------------+----------------+---------------------------+
   |"rock"     |object | array |see below       |:ref:`minc_rock_properties`|
   |           |               |                |                           |
   +-----------+---------------+----------------+---------------------------+

.. [PrNa85] Pruess, K. and Narasimhan, T.N. (1985). "A practical method for modeling fluid and heat flow in fractured porous media". Soc. Pet. Eng. J. 25, 14 -- 26.

.. index:: MINC; geometry
.. _minc_geometry:

MINC geometry
=============

As the MINC matrix cells are treated the same as any of the other cells in the :ref:`finite_volume_discretisation`, it is necessary to assign appropriate volumes to them. Similarly, areas must be assigned to the additional mesh faces introduced by MINC, as well as connection distances from each face to the cells on either side.

For conventional non-MINC cells these geometric parameters are computed simply from the mesh geometry. For MINC cells they must be computed based on some assumptions about the geometry of the fractures in the rock medium.

The MINC method makes the simplifying assumption that thermodynamic conditions in the rock matrix cells depend primarily on the distance from the nearest fracture. The MINC cells are therefore constructed so that the interfaces between them are parallel to the nearest fracture. There are various ways of doing this. The simplest and most common way proceeds by idealising the matrix cells as sets of identical cubes, and this is the approach adopted in Waiwera's implementation of the MINC method. More details can be found in [PrNa85]_.

The main parameters controlling the MINC geometry are:

* the **fracture volume fraction** :math:`f`: the fraction of the original cell volume occupied by fractures (typically small)
* the **matrix volume fractions** :math:`m_i`: the fraction of the original cell volume occupied by each level of rock matrix
* the number of sets of perpendicular **fracture planes** (1, 2 or 3) enclosing the matrix rock
* the **spacing** between fracture planes

The **"minc.geometry"** value in the "mesh" object contains two values, **"fracture"** and **"matrix"**, specifying the MINC fracture and matrix geometry parameters respectively.

.. note::
   **JSON object**: MINC geometry parameters

   **JSON path**: mesh.minc.geometry (or mesh.minc[`index`].geometry)

   +----------------+----------------+----------------+-----------------------+
   |**name**        |**type**        |**default**     |**value**              |
   +----------------+----------------+----------------+-----------------------+
   |"fracture"      |object          |see below       |MINC fracture geometry |
   |                |                |                |parameters             |
   |                |                |                |                       |
   +----------------+----------------+----------------+-----------------------+
   |"matrix"        |object          |see below       |MINC matrix geometry   |
   |                |                |                |parameters             |
   +----------------+----------------+----------------+-----------------------+

Volume fractions
----------------

The fracture and matrix volume fractions may be specified using the **"fracture.volume"** and **"matrix.volume"** values in the "mesh.minc" object. Because the volume fractions must sum to one (i.e. :math:`f + \sum_i{m_i} = 1`), it is not usually necessary to specify both of these values. If either one is omitted, the default value it takes is set so that the above equation is satisfied.

If there are multiple MINC matrix levels, the "matrix.volume" value in the "mesh.minc" object should be specified as an array of matrix volume fractions :math:`m_i`. If it is specified as a single number, there is just one matrix level (i.e. a "dual-porosity" model). The size of this value implicitly determines the number of matrix levels.

If both fracture and matrix volume fractions are specified, but they do not sum to one, then they will be scaled so that they do sum to one. (This means, for example, that volume fractions can be specified as percentage values instead of decimal fractions.) Note that this is not possible unless both fracture and matrix volume fractions are specified.

Fracture spacings and planes
----------------------------

The fracture spacing is specified using the the **"fracture.spacing"** value in the "mesh.minc" object. This can be either a single number or an array, depending on how many sets of fracture planes have been specified via the **"fracture.planes"** value (1, 2 or 3). If it is specified as an array, then the different sets of fracture planes may have different spacings. If there are multiple sets of fracture planes, but the spacing is specified as a number, then that value will be applied uniformly to all sets of fracture planes.

Fracture connection distance
----------------------------

The **"fracture.connection"** value in the "mesh.minc" object specifies the distance from each MINC fracture cell to the face connecting it to its first-level matrix rock cell. This is zero by default, but can optionally be set to a small positive value if this improves numerical performance.

.. note::
   **JSON object**: MINC fracture geometry parameters

   **JSON path**: mesh.minc.geometry.fracture (or mesh.minc[`index`].geometry.fracture)

   +----------------+----------------+----------------+-----------------------------+
   |**name**        |**type**        |**default**     |**value**                    |
   +----------------+----------------+----------------+-----------------------------+
   |"volume"        |number          |:math:`1 -      |fracture volume fraction     |
   |                |                |\sum_i{m_i}`    |:math:`f`                    |
   |                |                |                |                             |
   +----------------+----------------+----------------+-----------------------------+
   |"planes"        |integer         |1               |number of fracture planes    |
   +----------------+----------------+----------------+-----------------------------+
   |"spacing"       |number | array  |50 m            |fracture spacings (m)        |
   +----------------+----------------+----------------+-----------------------------+
   |"connection"    |number          |0               |fracture connection distance |
   |                |                |                |(m)                          |
   +----------------+----------------+----------------+-----------------------------+

.. note::
   **JSON object**: MINC matrix geometry parameters

   **JSON path**: mesh.minc.geometry.matrix (or mesh.minc[`index`].geometry.matrix)

   +---------+---------------+---------------------+----------------+
   |**name** |**type**       |**default**          |**value**       |
   +---------+---------------+---------------------+----------------+
   |"volume" |number | array |:math:`1 - f` (if    |matrix volume   |
   |         |               |:math:`f` specified),|fractions       |
   |         |               |otherwise 0.9        |:math:`m_i`     |
   |         |               |                     |                |
   |         |               |                     |                |
   +---------+---------------+---------------------+----------------+

.. index:: MINC; rock properties
.. _minc_rock_properties:

MINC rock properties
====================

The **"rock"** value in the "mesh.minc" object specifies rock properties for the fracture and matrix rocks in the MINC zone, as well as which parts of the mesh these rock properties are assigned to. In the simplest case, the "rock" value is a single object, although for more flexibility it can also be specified as an array of objects, with different rock properties assigned to different parts of the mesh zone.

.. note::
   **JSON object**: MINC rock properties

   **JSON path**: mesh.minc.rock, mesh.minc[`index`].rock, mesh.minc.rock[`index`] or mesh.minc[`index1`].rock[`index2`]

   +----------------+----------------+----------------+-----------------+
   |**name**        |**type**        |**default**     |**value**        |
   +----------------+----------------+----------------+-----------------+
   |"fracture"      |object          |(no default)    |fracture rock    |
   |                |                |                |type             |
   +----------------+----------------+----------------+-----------------+
   |"matrix"        |object          |(no default)    |matrix rock type |
   +----------------+----------------+----------------+-----------------+
   |"cells"         |array           |[]              |cell indices for |
   |                |                |                |this MINC zone   |
   +----------------+----------------+----------------+-----------------+
   |"zones"         |string | array  |[]              |mesh zones for   |
   |                |                |                |this MINC zone   |
   +----------------+----------------+----------------+-----------------+
   |"types"         |string | array  |[]              |rock types for   |
   |                |                |                |this MINC zone   |
   +----------------+----------------+----------------+-----------------+

.. index:: MINC; zones
.. _minc_zone_extent:

Defining the extent of MINC zones
---------------------------------

A MINC zone can cover all or only part of the simulation mesh. Because the MINC process adds cells to the computational mesh, it increases the computational cost of the simulation, so in many cases it is applied only to the parts of the mesh where it is necessary.

The parts of the mesh covered by a MINC zone can be specified in the "mesh.minc.rock" object in three ways, by specifying:

* **"cells"**: an array of cell indices for individual cells that belong to the MINC zone
* **"zones"**: one or more :ref:`mesh_zones` defining the extent of the MINC zone
* **"types"**: one or more :ref:`rock_types` which cover the MINC zone

The "zones" and "types" values can be either single strings or arrays of strings, containing the names of the mesh zones or rock types defining the extent of the MINC zone. It is possible (though not usual) to define the extent of a MINC zone by a combination of cells, zones and rock types.

MINC rock types
---------------

The **"fracture"** and **"matrix"** values in the "mesh.minc.rock" object define the rock properties in the MINC zone's fracture and matrix cells, respectively. Currently these MINC rock properties must be defined via :ref:`rock_types`, so both the "fracture" and "matrix" values are objects containing a single string value, **"type"**, specifying the appropriate rock type name. A rock type with the specified name must be defined in the main "rock.types" array (`not` the "mesh.minc.rock.types" array), where all the other rock types for single-porosity cells are also defined.

Normally a rock type definition includes the specification of which cells are assigned that rock type (see :ref:`rock_type_cells_and_zones`). For MINC fracture and matrix rock types this is not necessary as the MINC rock assignments are defined in the "minc.rock" value (see :ref:`minc_zone_extent`). Hence, the "cells" and "zones" values in MINC rock types do not need to be specified.

Default MINC rock properties
----------------------------

The fracture and matrix rock types are "based" on the original single-porosity rock type for the cell, in the sense that if any fracture or matrix rock properties are not specified, they are given default values from the original single-porosity rock type (if this exists -- otherwise the standard default values for rock properties are applied, as for any other rock type). The only exception to this is the MINC matrix rock porosity, which is given a special default value as described below.

MINC matrix porosities
----------------------

When the MINC method is used, an original single-porosity cell is replaced by a MINC fracture cell and one or more MINC matrix cells. In some cases it is necessary to ensure that the MINC cells have the same total void fraction as the original single-porosity cells.

For example, if a natural-state model is run, followed by a transient (e.g. production) model run that uses the natural-state results as initial conditions, then it is usual to run the natural-state model in single-porosity, even if the subsequent transient model is run using MINC. (This is because at steady state each MINC fracture cell should be in thermodynamic equilibrium with its matrix cells, giving the same results as for a single-porosity run.) In this case it is necessary to make sure the total void fractions in the transient MINC model are consistent with those in the original natural-state single-porosity cells, so that the two models have the same total pore volume.

The total void fraction in a MINC cell is consistent with the original porosity :math:`\phi` if:

.. math::

   f \phi_f + (1 - f) \phi_m = \phi

where :math:`\phi_f` and :math:`\phi_m` are the fracture and matrix porosities (and, as above, :math:`f` is the fracture volume fraction). This can be ensured by setting the matrix porosity to:

.. math::
   :label: minc_matrix_porosity

   \phi_m = \frac{\phi - f \phi_f}{1 - f}

While it is possible to calculate and set the porosities of the matrix rock types manually, Waiwera will do this automatically for any MINC matrix rock types for which the porosity :math:`\phi_m` is simply not specified in the JSON input file. That is, the porosity given by equation :eq:`minc_matrix_porosity` is the default for MINC matrix rock types.

Example
=======

In the example below, a model is set up with MINC applied to a central "production" zone (see :ref:`mesh_zones`), with single porosity outside of that zone. The default MINC geometry parameters are used except that the fracture spacing is set to 45 m, and three matrix rock levels are used inside each cell, occupying 15%, 30% and 50% of the cell volumes respectively. The fracture volume fraction is not specified, so by default the fractures make up the remaining 5% of the volume. The fracture rock type has a high permeability (10\ :sup:`-12` m\ :sup:`2`) and porosity (0.5), while the matrix rock type has low permeability (10\ :sup:`-16` m\ :sup:`2`) and porosity (0.05).

.. code-block:: json

   {"mesh": {"filename": "my_mesh.msh",
             "zones": {"production": {"x": [-500, 500], "y": [-500, 500], "z": [-1000, -200]},
                       "outer": {"-": "production"}},
             "minc": {"geometry": {"fracture": {"spacing": 45},
                                   "matrix": {"volume": [0.15, 0.3, 0.5]}},
                      "rock": {"fracture": {"type": "fracture"},
                               "matrix": {"type": "matrix"},
                               "zones": "production"
                     }}
            },
   "rock": {"types": [{"name": "formation",
                       "permeability": [1e-14, 1e-14, 1e-15],
                       "zones": ["outer"]},
                       {"name": "fracture",
                        "permeability": [1e-12, 1e-12, 1e-12],
                        "porosity": 0.5},
                       {"name": "matrix",
                        "permeability": [1e-16, 1e-16, 1e-16],
                        "porosity": 0.05}]}
   }

.. index:: MINC; initial conditions, initial conditions; MINC

.. _minc_initial_conditions:

MINC initial conditions
=======================

When specifying :ref:`initial_conditions` for a MINC simulation, it is possible either to specify initial conditions for all cells, including the MINC fracture and matrix cells, or to specify initial conditions only in the fracture cells (and cells outside of any MINC zones).

In the latter case, the initial conditions for the fracture cells are also applied to the matrix cells inside them. Hence, the matrix cells are initially in thermodynamic equilibrium with the fracture cells. This option can be used, for example, when restarting a transient simulation using the results from a steady-state single-porosity run as initial conditions (see :ref:`restarting`).

In the Waiwera JSON input file, the **"initial.minc"** boolean value specifies whether the initial conditions contain data for the MINC matrix cells. If it is set to ``false`` (the default), then the initial conditions are assumed not to contain data for the MINC matrix cells, and the initial conditions for the fracture cells will be extended into the matrix cells. If it is set to ``true``, the initial conditions are assumed to contain data for all MINC cells (e.g. when restarting one MINC simulation from a previous MINC run).

.. index:: MINC; output, output; MINC

MINC output
===========

The output from MINC simulations is much the same as for a single-porosity simulation (see :ref:`setup_output`), except that results are also output for the MINC matrix cells. These are output after the results for all the single-porosity and fracture cells. They are ordered first by MINC matrix level and second by the order of the corresponding fracture cell. So after the single-porosity and fracture cell results, first come all the first-level MINC matrix cell results (ordered by fracture cell index), then the second-level MINC matrix cell results, etc. This is true also if there are multiple MINC zones.

:ref:`setup_logfile` for MINC simulations is also very similar to that for single-porosity simulations. The main difference is that when cell indices are given (e.g. for phase transitions, or the location of maximum residuals) these always refer to single-porosity or fracture cells, and the MINC matrix level is also given. Level zero indicates single-porosity or fracture cells, and subsequent levels indicate the matrix level.
