.. index:: thermodynamics, simulation; thermodynamics
.. _water_thermodynamics:

********************
Water thermodynamics
********************

Waiwera includes implementations of two different "thermodynamic formulations" for water, i.e. sets of equations for calculating the thermodynamic properties of pure water as functions of the :ref:`primary_variables`. These are:

.. index:: thermodynamics; IFC-67, thermodynamics; IAPWS-97

* the IFC-67 formulation [IFC-67]_
* the IAPWS-97 formulation [IAPWS-97]_

The Waiwera JSON input file has a **"thermodynamics"** value for specifying the water thermodynamic formulation. This can be either a string containing the formulation name ("IFC67" or "IAPWS", either upper or lower case) or an object containing a **name** string value. If not specified, the default is "iapws".

.. note::
   **JSON object**: water thermodynamic formulation

   **JSON path**: thermodynamics

   +-------------+----------+-------------------+-----------------------+
   |**name**     |**type**  |**default**        |**value**              |
   +-------------+----------+-------------------+-----------------------+
   |"name"       |string    |"iapws"            |thermodynamic          |
   |             |          |                   |formulation name       |
   +-------------+----------+-------------------+-----------------------+
   |"extrapolate"|boolean   |``false``          |liquid thermodynamics  |
   |             |          |                   |extrapolation          |
   +-------------+----------+-------------------+-----------------------+

Examples:

.. code-block:: json

  {"thermodynamics": "ifc67"}

.. code-block:: json

  {"thermodynamics": {"name": "ifc67"}}

.. [IFC-67] International Formulation Committee (1967). "A formulation of the thermodynamic properties of ordinary water substance", Düsseldorf, Germany, 1967.
.. [IAPWS-97] Wagner, W., Cooper, J.R., Dittman, A., Kijima, J., Kretzschmar, H.-J., Kruse, A., Mares, R., Oguchi, K., Sato, H., Stöcker, I., Sifner, O., Takaishi, Y., Tanishita, I., Trübenbach, J., Willkommen, Th. (2000). "The IAPWS Industrial Formulation 1997 for the thermodynamic properties of water and steam". ASME J. Eng. Gas Turbines Power 122, 150 -- 182.

.. index:: thermodynamics; regions
.. _thermodynamic_regions:

Thermodynamic regions
=====================

Both IFC-67 and IAPWS-97 thermodynamic formulations divide the primary variable space into distinct **regions**, most of which represent different phase conditions.

The four thermodynamic regions used by the IAPWS-97 formulation are:

1) Liquid water
2) Vapour
3) Supercritical
4) Two-phase

This thermodynamic region numbering is used by Waiwera when reporting phase conditions (regardless of which thermodynamic formulation is used), for example when phase changes occur.

The IAPWS-97 regions are shown on a pressure-temperature diagram in :numref:`iapws_regions_plot`. The diagram extends over the range of validity of the IAPWS-97 formulation (pressure :math:`\leq` 100 MPa, 0 :math:`^{\circ}`\ C :math:`\leq` temperature :math:`\leq` 800 :math:`^{\circ}`\ C).

.. _iapws_regions_plot:
.. figure:: iapws_regions.*
           :scale: 67 %
           :align: center

           IAPWS-97 thermodynamic regions

Extrapolating liquid water thermodynamics
=========================================

The thermodynamics for liquid water (region 1) are valid up to a maximum temperature of 350 :math:`^{\circ}`\ C (for both IAPWS-97 and IFC-67 formulations). Currently Waiwera does not offer equations of state for supercritical water, so liquid temperatures over 350 :math:`^{\circ}`\ C cannot be simulated.

However, for some models temperatures may need to exceed this limit temporarily, for example, while running to steady state. In such cases it can be valid to relax this hard limit on liquid temperatures slightly in order to obtain a solution. The **"thermodynamics.extrapolate"** JSON input value can be used to activate this option. This is a Boolean value which defaults to ``false``.

Setting it to ``true`` allows the liquid water thermodynamics to be extrapolated up to a revised maximum of 360 :math:`^{\circ}`\ C. The liquid water thermodynamics are still approximately correct up to this temperature. However, it is not recommended to rely on this option for models that genuinely require output temperatures over 350 :math:`^{\circ}`\ C.

Example:

.. code-block:: json

  {"thermodynamics": {"name": "iapws", "extrapolate": true}}
