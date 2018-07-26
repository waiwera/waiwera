.. index:: equation of state (EOS), simulation; EOS
.. _eos:

*****************
Equation of state
*****************

In different simulations there may be different combinations of mass and energy components present. In most subsurface simulations, water is present, but there may be others such as carbon dioxide or air. In non-isothermal simulations it is also necessary to model the energy component, to account for the effects of temperature on fluid properties. The set of equations describing the behaviour of the particular mixture of mass and energy components present under a range of thermodynamic conditions is referred to as the "equation of state" (EOS) module.

Given the primary thermodynamic variables, the EOS module calculates bulk fluid properties such as pressure and temperature, as well as the properties of the individual phases present, such as density, internal energy, viscosity etc. It also checks the primary variables to make sure they have not gone outside acceptable physical bounds, and handles primary variable switching when phase transitions occur.

In the Waiwera JSON input file, the **"eos"** value specifies the equation of state module to be used for the simulation.

.. note::
   **JSON object**: equation of state

   **JSON path**: eos

   +-------------+----------+-------------------+-----------------------+
   |**name**     |**type**  |**default**        |**value**              |
   +-------------+----------+-------------------+-----------------------+
   |"name"       |string    |"we"               |Abbreviated EOS module |
   |             |          |                   |name                   |
   +-------------+----------+-------------------+-----------------------+
   |"temperature"|number    |20\                |constant temperature ( |
   |             |          |:math:`^{\circ}`\ C|:math:`^{\circ}`\ C)   |
   |             |          |                   |for :ref:`water_eos`   |
   |             |          |                   |EOS                    |
   |             |          |                   |                       |
   |             |          |                   |                       |
   +-------------+----------+-------------------+-----------------------+

For example:

 .. code-block:: json

  {"eos": {"name": "wae"}}

selects the water / air / energy EOS.

The equation of state modules included in Waiwera are described below.

Water EOS modules
=================

.. index:: equation of state (EOS); water ("w")
.. _water_eos:

Water ("w")
-----------

.. note::

   * **name**: "w"
   * **component names**: ["water"]
   * **phase conditions**: liquid
   * **primary variable names**: ["pressure"]
   * **default primary variables**: [10\ :sup:`5` Pa]
   * **default region**: 1 (liquid)
   * **default output fluid fields**: ["pressure", "region"]

This is the simplest equation of state module, simulating the behaviour of isothermal, single-phase liquid water. There is only one primary thermodynamic variable: pressure.

The temperature of the simulation can be specified in the Waiwera JSON input file via the **eos.temperature** value. This is a number value, specifying the temperature in degrees Celsius (:math:`^{\circ}`\ C). Note that this value is not needed (and will be ignored) if specified for other, non-isothermal EOS modules.

For example:

 .. code-block:: json

  {"eos": {"name": "w", "temperature": 18.5}}

Fluid properties are calculated directly from the thermodynamic formulation for water (see :ref:`water_thermodynamics`), at the specified temperature.

.. index:: equation of state (EOS); water / energy ("we")
.. _water_energy_eos:

Water and energy ("we")
-----------------------

.. note::

   * **name**: "we"
   * **component names**: ["water", "energy"]
   * **phase conditions**: liquid, vapour, two-phase
   * **primary variable names**:

     * **liquid / vapour**: ["pressure", "temperature"]
     * **two-phase**: ["pressure", "vapour_saturation"]
   * **default primary variables**: [10\ :sup:`5` Pa, 20 :math:`^{\circ}`\ C]
   * **default region**: 1 (liquid)
   * **default output fluid fields**: ["pressure", "temperature", "region", "vapour_saturation"]

This is the simplest non-isothermal equation of state module, with only one mass component (water) but also including the energy component. Water may be in liquid, vapour or two-phase conditions, and may change between these states. Primary variables are pressure and temperature for single-phase conditions but switch to pressure and vapour saturation under two-phase conditions.

Fluid properties are calculated directly from the thermodynamic formulation for water (see :ref:`water_thermodynamics`).

Water / NCG EOS modules
=======================

These EOS modules simulate mixtures of water and non-condensible gases (NCGs), together with energy. They work in much the same way as the water / energy EOS ("we") apart from modifications to the fluid properties resulting from the presence of the non-condensible gas.

The primary variables for these EOS modules are as for the water / energy EOS, but with an added third variable, the partial pressure of the non-condensible gas.

.. add detail on how NCG mixture EOS modules work? - using Henry's derivative to compute energy of solution etc.

.. index:: equation of state (EOS); water / air / energy ("wae")

Water, air and energy ("wae")
-----------------------------

.. note::

   * **name**: "wae"
   * **component names**: ["water", "air", "energy"]
   * **phase conditions**: liquid, vapour, two-phase
   * **primary variable names**:

     * **liquid / vapour**: ["pressure", "temperature", "air_partial_pressure"]
     * **two-phase**: ["pressure", "vapour_saturation", "air_partial_pressure"]
   * **default primary variables**: [10\ :sup:`5` Pa, 20 :math:`^{\circ}`\ C, 0 Pa]
   * **default region**: 1 (liquid)
   * **default output fluid fields**: ["pressure", "temperature", "region", "air_partial_pressure", "vapour_saturation"]

.. index:: equation of state (EOS); water / air / carbon dioxide ("wce")

Water, carbon dioxide and energy ("wce")
----------------------------------------

.. note::

   * **name**: "wce"
   * **component names**: ["water", "CO2", "energy"]
   * **phase conditions**: liquid, vapour, two-phase
   * **primary variable names**:

     * **liquid / vapour**: ["pressure", "temperature", "CO2_partial_pressure"]
     * **two-phase**: ["pressure", "vapour_saturation", "CO2_partial_pressure"]
   * **default primary variables**: [10\ :sup:`5` Pa, 20 :math:`^{\circ}`\ C, 0 Pa]
   * **default region**: 1 (liquid)
   * **default output fluid fields**: ["pressure", "temperature", "region", "CO2_partial_pressure", "vapour_saturation"]
       
