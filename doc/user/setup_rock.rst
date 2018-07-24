.. index:: rock

***************
Rock parameters
***************

Parameters related to the details of the rock media for the simulation are specified in the Waiwera JSON input file via the **"rock"** value.

.. note::

   +-----------------+-----------------+
   | **JSON value**: | "rock"          |
   +-----------------+-----------------+
   | **type**:       | object          |
   +-----------------+-----------------+
   | **specifies**:  | rock parameters |
   +-----------------+-----------------+

The rock value is an object, and may contain the values shown below.

.. note::

   **JSON object**: rock

   +-----------------------+--------------+----------------------+-----------------------+
   |**name**               |**type**      |**default**           |**value**              |
   +-----------------------+--------------+----------------------+-----------------------+
   |"types"                |array         |[]                    |:ref:`rock_types`      |
   +-----------------------+--------------+----------------------+-----------------------+
   |"relative_permeability"|object        |``{"type": "linear"}``|:ref:`relative_perm`   |
   |                       |              |                      |                       |
   +-----------------------+--------------+----------------------+-----------------------+
   |"capillary_pressure"   |object | null |``{"type": "zero"}``  |:ref:`capillarity`     |
   |                       |              |                      |                       |
   +-----------------------+--------------+----------------------+-----------------------+

.. index:: simulation; rock properties, rock; properties

Rock properties
===============

Each cell in the Waiwera simulation mesh contains rock with particular physical properties (e.g. permeability, porosity etc.). The rock properties may potentially be different in each cell.

.. index:: rock; types
.. _rock_types:

Rock types
----------

In many cases, rock properties may be assumed approximately constant over specified parts of the simulation mesh. These may, for example, correspond to lithologic units. To facilitate setting up this kind of rock property distribution, "rock types" may be defined which specify a set of rock properties to be assigned to a given list of cells and/or :ref:`mesh_zones`.

Rock types may be specified in the Waiwera JSON input file via the **"rock.types"** value. This is an array of objects, each object corresponding to a rock type definition.

.. note::

   **JSON object**: rock type

    +------------------+--------------+----------------------+----------------------------------+
    |**name**          |**type**      |**default**           |**value**                         |
    +------------------+--------------+----------------------+----------------------------------+
    |"name"            |string        |""                    |optional rock type name           |
    |                  |              |                      |                                  |
    +------------------+--------------+----------------------+----------------------------------+
    |"permeability"    |array         |[10\ :sup:`-13`, 10\  |permeability (\                   |
    |                  |              |:sup:`-13`, 10\       |m\ :sup:`2`\ )                    |
    |                  |              |:sup:`-13`] m\        |                                  |
    |                  |              |:sup:`2`              |                                  |
    +------------------+--------------+----------------------+----------------------------------+
    |"wet_conductivity"|number        |2.5 W / (m            |heat conductivity of fully        |
    |                  |              |:math:`^{\circ}`\ C)  |saturated rock (W / (m            |
    |                  |              |                      |:math:`^{\circ}`\ C))             |
    +------------------+--------------+----------------------+----------------------------------+
    |"dry_conductivity"|number        |= wet conductivity    |heat conductivity of unsaturated  |
    |                  |              |                      |rock (W / (m :math:`^{\circ}`\ C))|
    +------------------+--------------+----------------------+----------------------------------+
    |"porosity"        |number        |0.1                   |porosity                          |
    +------------------+--------------+----------------------+----------------------------------+
    |"density"         |number        |2200 kg / m\ :sup:`3` |grain density (kg / m\ :sup:`3`)  |
    +------------------+--------------+----------------------+----------------------------------+
    |"specific_heat"   |number        |1000 J / (kg          |specific heat (J / (kg            |
    |                  |              |:math:`^{\circ}`\ C)  |:math:`^{\circ}`\ C))             |
    +------------------+--------------+----------------------+----------------------------------+
    |"cells"           |array         |[]                    |cells with this rock type         |
    +------------------+--------------+----------------------+----------------------------------+
    |"zones"           |string | array|[]                    |:ref:`mesh_zones` with this rock  |
    |                  |              |                      |type                              |
    +------------------+--------------+----------------------+----------------------------------+

.. index:: rock; permeability
.. _rock_permeability:

Rock permeability
-----------------
Permeability is in many simulations the most influential rock property. In the Waiwera JSON input, permeability is specified as a 2-D array for 2-D simulations, or as a 3-D array for 3-D simulations, to allow for anisotropy.

In the mass and energy balance equations, permeability appears only in the face flux terms (see :ref:`function_evaluations`), where the value at each mesh face is determined by harmonic weighting of the cell values on either side of the face. The scalar effective permeability normal to the face is chosen from the permeability array according to the **permeability direction** assigned to that face. By default, these directions are chosen according to the **permeability axes** of the mesh. These axes are, in turn, aligned by default with the mesh coordinate axes, so that the elements of the permeability array are associated with the :math:`x`, :math:`y` and :math:`z` axes (in a Cartesian mesh). For faces which are not perfectly aligned with any permeability axis (e.g. in non-rectangular, unstructured meshes) the axis most closely aligned with the face normal vector is used to determine the default permeability direction.

However, the mesh permeability axes can be rotated in the horizontal plane (for example, to align them with fault planes, or with the principal axes of a mesh that has itself been rotated) by specifying the **"mesh.permeability_angle"** value (see :ref:`specifying_mesh`). In this case, the permeability axes are rotated by the specified angle anti-clockwise from the :math:`x` axis.

For example:

.. code-block:: json

  {"mesh": {"filename": "my_mesh.msh",
            "permeability_angle": 30}}

rotates the permeability axes anti-clockwise in the horizontal plane by 30\ :math:`^{\circ}`.

.. index:: mesh; faces

In addition, individual face permeability directions can be manually overridden, via the **"mesh.faces"** value (see :ref:`specifying_mesh`). This value is an array of objects, each of which has a **"cells"** array value containing the indices of the two cells on either side of the face. There is also a **"permeability_direction"** value which specifies the integer permeability direction for the face, overriding the default value.

.. note::

   **JSON object**: mesh face
   
   +------------------------+----------+-----------+----------------------------+
   |**name**                |**type**  |**default**|**value**                   |
   +------------------------+----------+-----------+----------------------------+
   |"cells"                 |array     |[]         |cell indices                |
   +------------------------+----------+-----------+----------------------------+
   |"permeability_direction"|integer   |1          |face permeability direction |
   +------------------------+----------+-----------+----------------------------+

For example:
 
.. code-block:: json

  {"mesh": {"filename": "my_mesh.msh",
            "faces": [
              {"cells": [99, 100],  "permeability_direction": 2},
              {"cells": [152,  9],  "permeability_direction": 1},
              {"cells": [205, 288], "permeability_direction": 1},
            ]}
  }

overrides the permeability directions for three faces in the mesh, leaving all others at their default values.

Rock type cells and zones
-------------------------

The cells assigned to each rock type can be specified explicitly using the **"cells"** value, an array of integer cell indices.

If :ref:`mesh_zones` have been defined, then zones can also be assigned to the rock type, via the rock type **"zones"** value. This can be either a string specifying a single zone name, or an array of strings, specifying multiple zone names.

It is possible to specify both cells and zones for a rock type, in which case all cells identified either explicitly or via zones are assigned to the rock type.

If there are any cells not assigned to any rock type, they will be given default rock properties (as shown above under :ref:`rock_types`).

.. index:: rock; relative permeability, relative permeability
.. _relative_perm:

Relative permeability curves
============================

Relative permeability curves are a way of adapting Darcy's law to multi-phase flow. When multiple phases are present in a cell, the flow in a given phase may be inhibited by the presence of the other phases. For example, in two-phase flow, the presence of the vapour phase (steam) may reduce the effective permeability for the liquid water phase.

In the equation for mass flux through the cell faces (see :ref:`function_evaluations`) the relative permeability :math:`k_r^p` is a factor applied the rock permeability :math:`k` for phase :math:`p`. The relative permeability curves describe how :math:`k_r^p` for each phase varies as a function of saturation.

A variety of different types of relative permeability curves have been proposed. Waiwera offers several of these, and allows the desired curves to be specified in the JSON input file via the **"rock.relative_permeability"** value. This value is an object, containing a **"type"** string value which selects the type of curves, along with other parameters which depend on the curve type.

The different types of relative permeability curves available in Waiwera are described below.

.. index:: relative permeability; fully mobile

Fully mobile
------------

This type of relative permeability curve maintains full mobility for all phases, regardless of saturation (i.e. :math:`k_r^p = 1` for all phases). It can be specified by setting the **"type"** value to "fully mobile" (or "fully_mobile"). There are no other parameters.

.. note::

   **JSON object**: fully mobile relative permeability

   +----------+----------+--------------+----------------------+
   |**name**  |**type**  |**default**   |**value**             |
   +----------+----------+--------------+----------------------+
   |"type"    |string    |"fully mobile"|relative permeability |
   |          |          |              |curve type            |
   +----------+----------+--------------+----------------------+

For example:

.. code-block:: json

  {"rock": {"relative_permeability": {"type": "fully mobile"}}}

.. index:: relative permeability; linear

Linear
------

Setting the **"type"** value to "linear" selects linear relative permeability functions. Here the relative permeabilities are linear functions of saturation.

For both liquid and vapour phases, the curves vary linearly from zero to one between the specified saturation limits for that phase. Below the lower limit, the relative permeability is identically zero, and above the upper limit it is identically one. The limits are specified in the Waiwera JSON input file via the **"liquid"** and **"vapour"** array values.

.. note::

   **JSON object**: linear relative permeability

   +------------+------------+------------+----------------------------+
   |**name**    |**type**    |**default** |**value**                   |
   +------------+------------+------------+----------------------------+
   |"type"      |string      |"linear"    |relative permeability curve |
   |            |            |            |type                        |
   |            |            |            |                            |
   +------------+------------+------------+----------------------------+
   |"liquid"    |array       |[0, 1]      |liquid saturation limits    |
   +------------+------------+------------+----------------------------+
   |"vapour"    |array       |[0, 1]      |vapour saturation limits    |
   +------------+------------+------------+----------------------------+

For example:

.. code-block:: json

  {"rock": {"relative_permeability": {"type": "linear",
                                      "liquid": [0.1, 0.9],
                                      "vapour": [0.1, 0.9]}}}

specifies linear relative permeability curves for both liquid and vapour phases, with limits 0.1 and 0.9, as in the figure below:

.. figure:: relative_permeability_linear.*
           :scale: 67 %
           :align: center

           Example linear relative permeability curves

Hence, in this example, for liquid saturations below 0.1 the liquid phase is immobile, while the vapour phase is fully mobile (as the vapour saturation is above 0.9). Conversely, for liquid saturations above 0.9 the liquid phase is fully mobile but the vapour phase is immobile.

.. index:: relative permeability; Pickens

Pickens
-------

For the Pickens-type relative permeability curves, chosen by setting the **"type"** value to "pickens", the liquid relative permeability varies with liquid saturation :math:`S_1` according to a power law: :math:`k_r^1 = S_1^{\alpha}`. The exponent :math:`\alpha` can be specified in the Waiwera JSON input file via the **"power"** value. The vapour relative permeability is identically one (i.e. vapour is fully mobile; :math:`k_r^2 = 1`).

.. note::

   **JSON object**: Pickens relative permeability

   +------------+------------+------------+----------------------+
   |**name**    |**type**    |**default** |**value**             |
   +------------+------------+------------+----------------------+
   |"type"      |string      |"pickens"   |relative permeability |
   |            |            |            |curve type            |
   |            |            |            |                      |
   +------------+------------+------------+----------------------+
   | "power"    |number      |1           |exponent              |
   |            |            |            |:math:`\alpha` for    |
   |            |            |            |liquid power law      |
   +------------+------------+------------+----------------------+

For example:

.. code-block:: json

  {"rock": {"relative_permeability": {"type": "pickens", "power": 1.5}}}

specifies Pickens curves with the power-law exponent :math:`\alpha = 1.5`.

.. index:: relative permeability; Corey

Corey
-----

Corey relative permeability curves are selected by setting the **"type"** value to "corey". Here the relative permeabilities are defined as functions of an intermediate quantity :math:`S_*`:

.. math::

   S_* = \frac{S_1 - S_{lr}}{1 - S_{lr} - S_{sr}}

where :math:`S_1` is the liquid saturation, and :math:`S_{lr}` and :math:`S_{sr}` are specified constant parameters. Then, if :math:`S_2 = 1 - S_1` is the vapour saturation:

.. math::

   k_r^1 =
   \begin{cases}
   1 & S_2 < S_{sr} \\
   S_*^4 & S_{sr} \leq S_2 \leq 1 - S_{lr} \\
   0 & S_2 > 1 - S_{lr}
   \end{cases}

.. math::

   k_r^2 =
   \begin{cases}
   0 & S_2 < S_{sr} \\
   (1 - S_*)^2 (1 - S_*^2) & S_{sr} \leq S_2 \leq 1 - S_{lr} \\
   1 & S_2 > 1 - S_{lr}
   \end{cases}

The two parameters :math:`S_{lr}` and :math:`S_{sr}` are specified in the Waiwera JSON input file via the **"slr"** and **"ssr"** values in the relative permeability object.

.. note::

   **JSON object**: Corey relative permeability

   +------------+------------+------------+-------------------------+
   |**name**    |**type**    |**default** |**value**                |
   +------------+------------+------------+-------------------------+
   |"type"      |string      |"corey"     |relative permeability    |
   |            |            |            |curve type               |
   |            |            |            |                         |
   +------------+------------+------------+-------------------------+
   |"slr"       |number      |0.3         |:math:`S_{lr}` parameter |
   |            |            |            |                         |
   +------------+------------+------------+-------------------------+
   |"ssr"       |number      |0.05        |:math:`S_{sr}` parameter |
   +------------+------------+------------+-------------------------+

For example:

.. code-block:: json

  {"rock": {"relative_permeability": {"type": "corey", "slr": 0.4, "ssr": 0.1}}}

specifies Corey relative permeability curves with :math:`S_{lr} = 0.4` and :math:`S_{sr} = 0.1`.

.. index:: relative permeability; Grant

Grant
-----

For the Grant relative permeability curves, selected by setting the **"type"** value to "grant", the liquid relative permeability is the same as for Corey curves. However, the vapour relative permeability is defined as :math:`k_r^2 = 1 - k_r^1`, so the liquid and vapour relative permeabilities always sum to one.

In the Waiwera JSON input file, the **"type"** value of the relative permeability object is set to "grant". All other values are the same as for the Corey curves (though the :math:`S_{sr}` parameter has a different default value).

.. note::

   **JSON object**: Grant relative permeability

   +------------+------------+------------+-------------------------+
   |**name**    |**type**    |**default** |**value**                |
   +------------+------------+------------+-------------------------+
   |"type"      |string      |"grant"     |relative permeability    |
   |            |            |            |curve type               |
   |            |            |            |                         |
   +------------+------------+------------+-------------------------+
   |"slr"       |number      |0.3         |:math:`S_{lr}` parameter |
   |            |            |            |                         |
   +------------+------------+------------+-------------------------+
   |"ssr"       |number      |0.6         |:math:`S_{sr}` parameter |
   +------------+------------+------------+-------------------------+

.. index:: relative permeability; Van Genuchten

Van Genuchten
-------------

Setting the relative permeability **"type"** value to "van genuchten" selects the Van Genuchten curves. The liquid relative permeability curve is defined in terms of an intermediate variable :math:`S_*`:

.. math::

   S_* = \frac{S_1 - S_{lr}}{S_{ls} - S_{lr}}

where :math:`S_1` is the liquid saturation, and :math:`S_{lr}` and :math:`S_{ls}` are specified constant parameters. Then the liquid relative permeability is given by:

.. math::

   k_r^1 =
   \begin{cases}
   0 & S_* < 0 \\
   \sqrt{S_*} (1 - (1 - S_*^{1 / \lambda})^{\lambda})^2 & 0 \le S_* < 1 \\
   1 & S_* \ge 1
   \end{cases}

where :math:`\lambda` is also a specified constant parameter.

For the vapour relative permeability, there are two variations.

In the first variation, the liquid and vapour relative permeabilities are forced to sum to one, by setting :math:`k_r^2 = 1 - k_r^1`. This variation can be selected in the Waiwera JSON input file by setting the **"sum_unity"** value in the relative permeability object to ``true`` (the default).

In the second variation, the vapour relative permeability curve is defined in terms of another intermediate variable :math:`\hat{s}`:

.. math::

   \hat{s} = \frac{S_1 - S_{lr}}{1 - S_{lr} - S_{sr}}

where :math:`S_{sr}` is another specified constant parameter. Then the vapour relative permeability is given by:

.. math::

   k_r^2 = \min{((1 - \hat{s})^2 (1 - \hat{s}^2), 1)}

.. note::

   **JSON object**: Van Genuchten relative permeability

   +------------+------------+----------------+--------------------------+
   |**name**    |**type**    |**default**     |**value**                 |
   +------------+------------+----------------+--------------------------+
   |"type"      |string      |"van genuchten" |relative permeability     |
   |            |            |                |curve type                |
   |            |            |                |                          |
   +------------+------------+----------------+--------------------------+
   |"lambda"    |number      |0.45            |:math:`\lambda` parameter |
   |            |            |                |                          |
   +------------+------------+----------------+--------------------------+
   |"slr"       |number      |10\ :sup:`-3`   |:math:`S_{lr}` parameter  |
   +------------+------------+----------------+--------------------------+
   |"sls"       |number      |1               |:math:`S_{ls}` parameter  |
   +------------+------------+----------------+--------------------------+
   |"ssr"       |number      |0.6             |:math:`S_{sr}` parameter  |
   +------------+------------+----------------+--------------------------+
   |"sum_unity" |boolean     |``true``        |enforce :math:`k_r^1 +    |
   |            |            |                |k_r^2 = 1`                |
   +------------+------------+----------------+--------------------------+

The :math:`S_{sr}` parameter is used only for the second variation of the vapour relative permeability curves, and has no effect if the "sum_unity" value is ``true``.

For example:

.. code-block:: json

  {"rock": {"relative_permeability": {"type": "van genuchten", "lambda": 0.4}}}

specifies Van Genuchten relative permeability curves with :math:`\lambda = 0.4` and all other parameters left at their default values.

.. index:: relative permeability; table

Table
-----
Setting the relative permeability **"type"** value to "table" allows specification of relative permeability curves defined as general piecewise-linear tables. For each phase :math:`p`, the relative permeability curve is specified as a table of :math:`(S_p, k_r^p)` values. In the Waiwera JSON input file these tables take the form of rank-2 arrays (i.e. arrays of arrays), specified via the **"liquid"** and **"vapour"** values.

.. note::

   **JSON object**: table relative permeability

   +------------+------------+---------------+-----------------------------------+
   |**name**    |**type**    |**default**    |**value**                          |
   +------------+------------+---------------+-----------------------------------+
   |"type"      |string      |"table"        |relative permeability curve type   |
   +------------+------------+---------------+-----------------------------------+
   |"liquid"    |array       |[[0,0], [1,1]] |table of liquid relative           |
   |            |            |               |permeability :math:`k_r^1`         |
   |            |            |               |vs. liquid saturation :math:`S_1`  |
   +------------+------------+---------------+-----------------------------------+
   |"vapour"    |array       |[[0,0], [1,1]] |table of vapour relative           |
   |            |            |               |permeability :math:`k_r^2`         |
   |            |            |               |vs. vapour saturation :math:`S_2`  |
   +------------+------------+---------------+-----------------------------------+

For example:

.. code-block:: json

  {"rock": {"relative_permeability": {
     "type": "table",
     "liquid": [[0,0], [0.1, 0.01], [0.9, 0.99], [1,1]],
     "vapour": [[0,0], [0.1, 0.01], [0.9, 0.99], [1,1]]
     }}}

specifies both liquid and vapour relative permeability curves as in the figure below, with a small slope at the extremes of saturation.

.. figure:: relative_permeability_table.*
           :scale: 67 %
           :align: center

           Example table relative permeability curves

.. index:: rock; capillary pressure, capillary pressure
.. _capillarity:

Capillary pressure functions
============================

Waiwera can optionally include capillary pressure effects when calculating pressure gradients across mesh faces. For the liquid phase, the effective pressure in each cell is calculated from the sum of the fluid pressure and capillary pressure, which in turn is calculated from a specified function of saturation. These effective pressures are then used to calculate the effective pressure gradient across the mesh face. (If the saturations are the same in both cells on either side of the face, then the capillary pressures are also equal and have no effect on the calculated pressure gradient.)

As for relative permeability curves, a variety of different capillary pressure functions have been proposed, and Waiwera offers several of them. The desired capillary pressure function is specified in the Waiwera JSON input file via the **"rock.capillary_pressure"** value. This value is an object (or ``null``), containing a **"type"** string value which selects the type of function, along with other parameters which depend on the function type.

The different types of capillary pressure functions available in Waiwera are described below.

.. index:: capillary pressure; zero

Zero
----

Capillary pressure effects can be disabled by setting the **"type"** value of the capillary pressure object to "zero" (or setting the capillary pressure value to ``null``). This is the default. In this case, the capillary pressure is identically zero regardless of saturation.

.. note::

   **JSON object**: zero capillary pressure function

   +----------+----------+--------------+----------------------+
   |**name**  |**type**  |**default**   |**value**             |
   +----------+----------+--------------+----------------------+
   |"type"    |string    |"zero"        |capillary pressure    |
   |          |          |              |function type         |
   +----------+----------+--------------+----------------------+

For example:

.. code-block:: json

  {"rock": {"capillary_pressure": {"type": "zero"}}}

or

.. code-block:: json

  {"rock": {"capillary_pressure": null}}

both disable capillary pressure effects.

.. index:: capillary pressure; linear

Linear
------

Setting the capillary pressure **"type"** value to "linear" selects the linear capillary pressure function, in which capillary pressure is a linear function of liquid saturation. Lower and upper saturation limits are specified via the **"saturation_limits"** array value.

When liquid saturation is below the lower limit, the capillary pressure is fixed at :math:`-P`, where :math:`P` is a specified (positive) constant. Between the limits, the capillary pressure is linearly interpolated between :math:`-P` and zero. Above the upper limit, the capillary pressure is identically zero.

.. note::

   **JSON object**: linear capillary pressure function

   +--------------------+------------+------------+-------------------------+
   |**name**            |**type**    |**default** |**value**                |
   +--------------------+------------+------------+-------------------------+
   |"type"              |string      |"linear"    |capillary pressure       |
   |                    |            |            |function type            |
   +--------------------+------------+------------+-------------------------+
   |"saturation_limits" |array       |[0, 1]      |liquid saturation limits |
   +--------------------+------------+------------+-------------------------+
   |"pressure"          |number      |0.125×10\   |magnitude :math:`P` of   |
   |                    |            |:sup:`5` Pa |maximum capillary        |
   |                    |            |            |pressure (Pa)            |
   +--------------------+------------+------------+-------------------------+

For example:

.. code-block:: json

  {"rock": {"capillary_pressure": {"type": "linear",
                                   "saturation_limits": [0.1, 0.9],
                                   "pressure": 10.0e3}}}

gives the linear capillary pressure curve shown in the figure below.

.. figure:: capillary_linear.*
           :scale: 67 %
           :align: center

           Example linear capillary pressure function

.. index:: capillary pressure; Van Genuchten

Van Genuchten
-------------

Setting the capillary pressure **"type"** value to "van genuchten" selects the Van Genuchten capillary pressure function. The capillary pressure is defined in terms of an intermediate quantity :math:`S_*`:

.. math::

   S_* = \frac{S_1 - S_{lr}}{S_{ls} - S_{lr}}

where :math:`S_1` is the liquid saturation. and :math:`S_{lr}` and :math:`S_{ls}` are specified constant parameters. Then the capillary pressure :math:`P_c` is given by

.. math::

   P_c =
   \begin{cases}
   -P_0 & S_* < 0\\
   \min{(-P_0 (S_*^{-1 / \lambda} -1) ^ {1 - \lambda}, 0)} & 0 \le S_* < 1\\
   0 & S_* \ge 1
   \end{cases}

where :math:`P_0` and :math:`\lambda` are specified constant parameters (:math:`P_0 > 0`). An optional limit :math:`P_{max}` can be set on the magnitude of the capillary pressure determined by the above equation. If this limit is not specified, no limit is applied.

.. note::

   **JSON object**: Van Genuchten capillary pressure function

   +------------+------------+----------------+--------------------+
   |**name**    |**type**    |**default**     |**value**           |
   +------------+------------+----------------+--------------------+
   |"type"      |string      |"van genuchten" |capillary pressure  |
   |            |            |                |function type       |
   |            |            |                |                    |
   +------------+------------+----------------+--------------------+
   |"lambda"    |number      |0.45            |:math:`\lambda`     |
   |            |            |                |parameter           |
   +------------+------------+----------------+--------------------+
   |"slr"       |number      |10\ :sup:`-3`   |:math:`S_{lr}`      |
   |            |            |                |parameter           |
   +------------+------------+----------------+--------------------+
   |"sls"       |number      |1               |:math:`S_{ls}`      |
   |            |            |                |parameter           |
   +------------+------------+----------------+--------------------+
   |"P0"        |number      |0.125×10\       |:math:`P_0`         |
   |            |            |:sup:`5` Pa     |parameter (Pa)      |
   +------------+------------+----------------+--------------------+
   |"Pmax"      |number      |undefined       |:math:`P_{max}`     |
   |            |            |                |parameter           |
   +------------+------------+----------------+--------------------+

For example:

.. code-block:: json

  {"rock": {"capillary_pressure": {"type": "van genuchten", "lambda": 0.5}}}

gives the Van Genuchten capillary pressure function with :math:`\lambda = 0.5`, no :math:`P_{max}` parameter applied, and all other parameters left at their default values.

.. index:: capillary pressure; table

Table
-----

Setting the capillary pressure **"type"** value to "table" allows specification of a capillary pressure function defined by a general piecewise-linear table. The capillary pressure function is specified as a table of :math:`(S_1, P_c)` values (i.e. capillary pressure vs. liquid saturation). In the Waiwera JSON input file this table takes the form of a rank-2 array (i.e. array of arrays), specified via the **"pressure"** value.

.. note::

   **JSON object**: table capillary pressure function

   +------------+------------+---------------+--------------------+
   |**name**    |**type**    |**default**    |**value**           |
   +------------+------------+---------------+--------------------+
   |"type"      |string      |"table"        |capillary pressure  |
   |            |            |               |function type       |
   |            |            |               |                    |
   +------------+------------+---------------+--------------------+
   |"pressure"  |array       |[[0,0], [1,0]] |table of capillary  |
   |            |            |               |pressure vs. liquid |
   |            |            |               |saturation          |
   +------------+------------+---------------+--------------------+

If the table does not cover the entire liquid saturation range :math:`0 \le S_1 \le 1`, the values at the limits of the table are used outside the table range.

For example:

.. code-block:: json

  {"rock": {"capillary_pressure": {
     "type": "table",
     "pressure": [[0.1, -0.1e5], [1, 0]]
     }}}

specifies a capillary pressure function with constant value -0.1 bar for liquid saturations between zero and 0.1, decreasing linearly to zero at fully-saturated conditions (:math:`S_1 = 1`).
