.. index:: thermodynamics, simulation; thermodynamics
.. _water_thermodynamics:

********************
Water thermodynamics
********************

Waiwera includes implementations of two different "thermodynamic formulations" for water, i.e. sets of equations for calculating the thermodynamic properties of pure water as functions of the :ref:`primary_variables`. These are:

.. index:: thermodynamics; IFC-67, thermodynamics; IAPWS-97

* the [IFC-67]_ formulation
* the [IAPWS-97]_ formulation

The Waiwera JSON input file has a **"thermodynamics"** value for specifying the water thermodynamic formulation. This can be either a string containing the formulation name ("IFC67" or "IAPWS", either upper or lower case) or an object containing a **name** string value. If not specified, the default is "iapws".

.. admonition:: JSON input

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

The thermodynamics for liquid water (region 1) are valid up to a maximum temperature of 350 :math:`^{\circ}`\ C (for both IAPWS-97 and IFC-67 formulations). Most of Waiwera's equations of state (apart from the :ref:`supercritical_eoses`) do not have support for region 3, so liquid temperatures over 350 :math:`^{\circ}`\ C cannot be simulated.

However, for some models temperatures may need to exceed this limit temporarily, for example, while running to steady state. In such cases it can be valid to relax this hard limit on liquid temperatures slightly in order to obtain a solution. The **"thermodynamics.extrapolate"** JSON input value can be used to activate this option. This is a Boolean value which defaults to ``false``.

Setting it to ``true`` allows the region 1 liquid water thermodynamics to be extrapolated up to a revised maximum of 360 :math:`^{\circ}`\ C. The liquid water thermodynamics are still approximately correct up to this temperature. However, it is not recommended to rely on this option for models that genuinely require output temperatures over 350 :math:`^{\circ}`\ C. In such cases one of the :ref:`supercritical_eoses` should be used if possible.

Example:

.. code-block:: json

  {"thermodynamics": {"name": "iapws", "extrapolate": true}}

Supercritical water thermodynamics
==================================

As the critical point of water (at approximately P = 22 MPa, T = 374 :math:`^{\circ}`\ C for IAPWS-97) is approached along the two-phase saturation line, the properties of liquid water and vapour converge, until at the critical point the two phases are indistinguishable. If the pressure and temperature are both above their critical values, there are no longer separate phases, only single-phase supercritical fluid.

.. More general discussion of liquid- and vapour-like behaviour, then Widom line, then Widom delta, pi_liq (below)

However, it is currently thought that supercritical fluid, despite being single-phase, contains at the microscopic level a mixture of liquid-like and vapour-like particles, and depending on the relative proportions of these types of particles, the supercritical fluid has more liquid-like or vapour-like behaviour. This can be described quantitatively by the "liquid-like fraction" :math:`\pi_{liq}` of the supercritical fluid, which is the proportion of liquid-like particles, taking values between 0 and 1 [Ha_et_al_2018]_.

.. Figure showing SC zone, Widom line and delta

.. Point out that SC zone is shared by regions 2 and 3

.. Widom line and delta (parameters alpha, beta)

.. Calculation of pi_liq in Waiwera

.. [Ha_et_al_2018] Ha, M.Y., Yoon, T.J., Tlusty, T., Jho, Y. and Lee, W.B. (2018). "Widom Delta of Supercritical Gas-Liquid Coexistence", J. Phys. Chem. Lett. 2018 (9), 1734 - 1738.
