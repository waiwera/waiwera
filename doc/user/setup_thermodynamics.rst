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
   |"widom"      |object    |{}                 |supercritical Widom    |
   |             |          |                   |delta parameters (IAPWS|
   |             |          |                   |only)                  |
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
           :scale: 75 %
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

.. _supercritical_thermodynamics:

Supercritical water thermodynamics
==================================

As the critical point of water (at approximately P = 22 MPa, T = 374 :math:`^{\circ}`\ C for IAPWS-97) is approached along the two-phase saturation line, the properties of liquid water and vapour converge, until at the critical point the two phases are indistinguishable. If the pressure and temperature are both above their critical values, there are no longer separate phases, only single-phase supercritical fluid.

.. _widom_line:

The Widom line
--------------

The properties of supercritical fluid can be liquid-like or vapour-like, or something in between, depending on its pressure and temperature. Several lines on the pressure-temperature diagram, based on various physical properties, and extending from the critical point, have been proposed to distinguish liquid-like and vapour-like supercritical fluid.

One of these lines is known as the **Widom line**. An expression for the location of the Widom line, based on experimental data, has been given by [Banuti_et_al_2017]_:

.. math::
   :label: widom_eqn

   P = P_c e^{A_s (T^k/T^k_c - 1)}

where :math:`P_c` is the critical pressure (Pa), :math:`T^k_c` is the critical temperature (K), :math:`T^k` is the temperature (K) and :math:`A_s` is a dimensionless constant, which for water has the value 6.479.

.. _widom_delta:

The Widom delta
---------------

The behaviour of supercritical fluid does not change suddenly from liquid-like to vapour-like as the Widom line is crossed. A transition zone around the line has been identified, known as the **Widom delta**, within which the fluid's behaviour is in between liquid-like and vapour-like [Ha_et_al_2018]_.

As yet, there do not appear to be any published expressions for the locations of the Widom delta boundaries.  In Waiwera the temperature difference between the vapour-like and liquid-like Widom delta boundaries is assumed to vary linearly with pressure, so that the boundaries :math:`T^v_{\delta}` and :math:`T^l_{\delta}` are given by:

.. math::
   :label: widom_delta_eqn

           T^v_{\delta} = T_w(P) + \frac{\alpha}{2} (P - P_0) / P_c

           T^l_{\delta} = T_w(P) - \frac{\alpha}{2} (P - P_0) / P_c

where :math:`T_w` is the Widom line temperature (:math:`^{\circ}`\ C) for the pressure :math:`P`, found by inverting equation :eq:`widom_eqn`, and :math:`\alpha` is a growth factor with default value 25 :math:`^{\circ}`\ C. :math:`P_0` is a reference pressure corresponding to a reference temperature :math:`T_0` on the saturation line, chosen slightly lower than the critical temperature :math:`T_c` (by default 0.1 :math:`^{\circ}`\ C lower), so that the Widom delta has a very small but finite width at the critical point. This avoids numerical issues with infinitely sharp transitions between liquid-like and vapour-like supercritical fluid at the critical point.

:numref:`supercritical_plot` shows the supercritical zone shaded grey on a pressure-temperature plot, with the dashed lines representing the boundaries between IAPWS-97 region 3 and regions 1 and 2. It also shows the Widom line (red) given by equation :eq:`widom_eqn` and the Widom delta boundaries (blue) given by equation :eq:`widom_delta_eqn`, with the default value of the growth parameter :math:`\alpha`. It can be seen that, for this value of :math:`\alpha`, the Widom delta lies within region 3. Most, but not all, of region 3 is supercritical. A significant proportion of region 2 is vapour-like supercritical fluid.

.. _supercritical_plot:
.. figure:: supercritical.*
           :scale: 50 %
           :align: center

           Supercritical zone (shaded grey), Widom line (red) and Widom delta boundaries (blue)

The Waiwera JSON input file has a **"thermodynamics.widom"** value for specifying Widom parameters. This contains a value **"thermodynamics.widom.delta"**, an object in which the Widom delta growth and offset parameters may be specified.

.. admonition:: JSON input

   **JSON object**: supercritical Widom parameters

   **JSON path**: thermodynamics.widom

   +-------------+----------+-------------------+-----------------------+
   |**name**     |**type**  |**default**        |**value**              |
   +-------------+----------+-------------------+-----------------------+
   |"delta"      |object    |{}                 |Widom delta parameters |
   |             |          |                   |                       |
   +-------------+----------+-------------------+-----------------------+

   **JSON object**: supercritical Widom delta parameters

   **JSON path**: thermodynamics.widom.delta

   +-------------+----------+-------------------+-----------------------+
   |**name**     |**type**  |**default**        |**value**              |
   +-------------+----------+-------------------+-----------------------+
   |"growth"     |number    |25                 |Widom delta growth     |
   |             |          |                   |factor :math:`\alpha`  |
   |             |          |                   |(:math:`^{\circ}`\ C)  |
   +-------------+----------+-------------------+-----------------------+
   |"offset"     |number    |0                  |Widom delta temperature|
   |             |          |                   |offset :math:`T_c -    |
   |             |          |                   |T_0` (:math:`^{\circ}`\|
   |             |          |                   |C)                     |
   +-------------+----------+-------------------+-----------------------+

.. _liquidlike_fraction:

Liquid-like fraction
--------------------

It is currently thought that supercritical fluid, despite being single-phase, contains at the microscopic level a mixture of molecules with liquid-like and vapour-like behaviour. The relative proportions of these two types of particles determines whether the supercritical fluid has more liquid-like or vapour-like behaviour. This can be described quantitatively by the "liquid-like fraction" :math:`\pi_{liq}` of the supercritical fluid, which is the number fraction of liquid-like particles, taking values between 0 and 1 [Ha_et_al_2018]_.

In Waiwera the liquid-like fraction is estimated (in the absence of experimental data) for a given supercritical pressure and temperature based on the Widom delta boundaries. First the non-dimensional location of the temperature within the Widom delta boundaries is calculated from:

.. math::
   :label: widom_xi

   \xi = \frac{T - T^l_{\delta}}{T^v_{\delta} - T^l_{\delta}}

Then the liquid-like fraction is calculated as:

.. math::
   :label: widom_pi_liq

   \pi_{liq} = \begin{cases}
   1 & \xi < 0 \\
   h_{00}(\xi) & 0 \le \xi \le 1 \\
   0 & \xi > 1
   \end{cases}

where :math:`h_{00} = (1 + 2 \xi) (1 - \xi)^2` is the cubic Hermite spline basis function between 1 and 0 on the unit interval, with zero gradient at each end. Using a cubic spline interpolant within the Widom delta results in a smooth and differentiable :math:`\pi_{liq}` surface over the entire supercritical zone.

.. [Banuti_et_al_2017] Banuti, D.T., Raju, M. and Ihme, M. (2017). "Similarity law for Widom lines and coexistence lines", Phys. Rev. E 95, 052120.
.. [Ha_et_al_2018] Ha, M.Y., Yoon, T.J., Tlusty, T., Jho, Y. and Lee, W.B. (2018). "Widom Delta of Supercritical Gas-Liquid Coexistence", J. Phys. Chem. Lett. 2018 (9), 1734 - 1738.
