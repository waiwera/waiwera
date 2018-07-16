.. _water_thermodynamics:

********************
Water thermodynamics
********************

Waiwera includes implementations of two different "thermodynamic formulations" for water, i.e. sets of equations for calculating the thermodynamic properties of pure water as functions of the :ref:`primary_variables`. These are:

* the IFC-67 formulation [IFC-67]_
* the IAPWS-97 formulation [IAPWS-97]_

The Waiwera JSON input file has a **"thermodynamics"** string value, which can be set to either "IFC67" or "IAPWS" (either upper or lower case) to select the thermodynamic formulation. If not specified, the default is "iapws".

.. note::

   +-----------------+---------------------------------------+
   | **JSON value**: | "thermodynamics"                      |
   +-----------------+---------------------------------------+
   | **type**:       | string                                |
   +-----------------+---------------------------------------+
   | **values**:     | "IFC67" | "ifc67" | "IAPWS" | "iapws" |
   +-----------------+---------------------------------------+
   | **default**:    | "iapws"                               |
   +-----------------+---------------------------------------+
   | **specifies**:  | water thermodynamics formulation      |
   +-----------------+---------------------------------------+

Example:

.. code-block:: json

  {"thermodynamics": "ifc67"}

.. [IFC-67] International Formulation Committee (1967). "A formulation of the thermodynamic properties of ordinary water substance", Düsseldorf, Germany, 1967.
.. [IAPWS-97] Wagner, W., Cooper, J.R., Dittman, A., Kijima, J., Kretzschmar, H.-J., Kruse, A., Mares, R., Oguchi, K., Sato, H., Stöcker, I., Sifner, O., Takaishi, Y., Tanishita, I., Trübenbach, J., Willkommen, Th. (2000). "The IAPWS Industrial Formulation 1997 for the thermodynamic properties of water and steam". ASME J. Eng. Gas Turbines Power 122, 150 -- 182.

Thermodynamic regions
=====================

Both IFC-67 and IAPWS-97 thermodynamic formulations divide the primary variable space into distinct **regions**, most of which represent different phase conditions.

The four thermodynamic regions used by the IAPWS-97 formulation are:

1) Liquid water
2) Vapour
3) Supercritical
4) Two-phase

This thermodynamic region numbering is used by Waiwera when reporting phase conditions (regardless of which thermodynamic formulation is used), for example when phase changes occur.

The IAPWS-97 regions are shown on a pressure-temperature diagram below. The diagram extends over the range of validitiy of the IAPWS-97 formulation (pressure :math:`\leq` 100 MPa, 0 :math:`^{\circ}`\ C :math:`\leq` temperature :math:`\leq` 800 :math:`^{\circ}`\ C).

.. figure:: iapws_regions.*
           :scale: 67 %
           :align: center

           IAPWS-97 thermodynamic regions

