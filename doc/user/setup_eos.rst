.. index:: equation of state (EOS), simulation; EOS
.. _eos:

*****************
Equation of state
*****************

In different simulations there may be different combinations of mass and energy components present. In most subsurface simulations, water is present, but there may be others such as carbon dioxide or air. In non-isothermal simulations it is also necessary to model the energy component, to account for the effects of temperature on fluid properties. The set of equations describing the behaviour of the particular mixture of mass and energy components present under a range of thermodynamic conditions is referred to as the "equation of state" (EOS) module.

Given the primary thermodynamic variables, the EOS module calculates bulk fluid properties such as pressure and temperature, as well as the properties of the individual phases present, such as density, internal energy, viscosity etc. It also checks the primary variables to make sure they have not gone outside acceptable physical bounds, and handles primary variable switching when phase transitions occur.

In the Waiwera JSON input file, the **"eos"** value specifies the equation of state module to be used for the simulation. This can be either a string containing the abbreviated name of the EOS, or an object containing a **"name"** value. In general it is not necessary to specify the EOS as an object unless other EOS parameters besides the name need to be set.

.. note::
   **JSON object**: equation of state

   **JSON path**: eos

   +-------------+----------+-------------------+-----------------------+
   |**name**     |**type**  |**default**        |**value**              |
   +-------------+----------+-------------------+-----------------------+
   |"name"       |string    |"we"               |abbreviated EOS module |
   |             |          |                   |name                   |
   +-------------+----------+-------------------+-----------------------+
   |"primary"    |object    |{}                 |primary variable       |
   |             |          |                   |parameters             |
   |             |          |                   |                       |
   |             |          |                   |                       |
   +-------------+----------+-------------------+-----------------------+
   |"temperature"|number    |20\                |constant temperature ( |
   |             |          |:math:`^{\circ}`\ C|:math:`^{\circ}`\ C)   |
   |             |          |                   |for :ref:`water_eos`   |
   |             |          |                   |EOS                    |
   +-------------+----------+-------------------+-----------------------+

For example:

 .. code-block:: json

  {"eos": {"name": "wae"}}

selects the water / air / energy EOS. Since the only EOS parameter is the name, this can be specified more simply as:

.. code-block:: json

   {"eos": "wae"}

.. index:: primary variables
.. _primary_variable_parameters:

Primary variable parameters
===========================

Each EOS module has a particular set of primary thermodynamic variables which determine the fluid state in each cell (see :ref:`primary_variables`). Parameters related to the primary variables can be specified via the **"eos.primary"** value in the Waiwera input JSON file. This is an object with just one value, **"eos.primary.scale"**.

.. note::
   **JSON object**: primary variable parameters

   **JSON path**: eos.primary

   +-------------+----------+-------------------+-----------------------+
   |**name**     |**type**  |**default**        |**value**              |
   +-------------+----------+-------------------+-----------------------+
   |"scale"      |object    |{}                 |scaling parameters     |
   +-------------+----------+-------------------+-----------------------+

.. index:: primary variables; scaling

Scaling
-------

Waiwera internally non-dimensionalises the primary variables to improve numerical behaviour. In most cases this is carried out via a simple scaling by a fixed constant. These fixed constants have default values, but can be over-ridden via the **"eos.primary.scale"** value in the Waiwera input JSON file. This is also an object, with values specific to each EOS module (see below).

Note that all input and output thermodynamic variables are in their usual dimensional form (i.e. not scaled).

Water EOS modules
=================

.. index:: equation of state (EOS); water ("w")
.. _water_eos:

Water ("w")
-----------

+-------------------------------+--------------------------+
|**abbreviated name**:          |"w"                       |
+-------------------------------+--------------------------+
|**component names**:           |["water"]                 |
+-------------------------------+--------------------------+
|**phase names**:               |["liquid"]                |
+-------------------------------+--------------------------+
|**primary variable names**:    |["pressure"]              |
+-------------------------------+--------------------------+
|**default primary variables**: |[10\ :sup:`5` Pa]         |
|                               |                          |
+-------------------------------+--------------------------+
|**default region**:            |1 (liquid)                |
+-------------------------------+--------------------------+
|**default eos.primary.scale**: |{"pressure": 1e6}         |
|                               |                          |
+-------------------------------+--------------------------+
|**default output fluid         |["pressure", "region"]    |
|fields**:                      |                          |
+-------------------------------+--------------------------+

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

+-------------------------------+--------------------------------------------------+
|**abbreviated name**:          |"we"                                              |
+-------------------------------+--------------------------------------------------+
|**component names**:           |["water", "energy"]                               |
+-------------------------------+--------------------------------------------------+
|**phase names**:               |["liquid", "vapour"]                              |
+-------------------------------+--------------------------------------------------+
|**primary variable names**:    |**single-phase**: ["pressure", "temperature"]     |
|                               +--------------------------------------------------+
|                               |**two-phase**: ["pressure", "vapour_saturation"]  |
+-------------------------------+--------------------------------------------------+
|**default primary variables**: |[10\ :sup:`5` Pa, 20 :math:`^{\circ}`\ C]         |
|                               |                                                  |
+-------------------------------+--------------------------------------------------+
|**default region**:            |1 (liquid)                                        |
+-------------------------------+--------------------------------------------------+
|**default eos.primary.scale**: |{"pressure": 1e6, "temperature": 100}             |
|                               |                                                  |
+-------------------------------+--------------------------------------------------+
|**default output fluid         |["pressure", "temperature", "region",             |
|fields**:                      |"vapour_saturation"]                              |
+-------------------------------+--------------------------------------------------+

This is the simplest non-isothermal equation of state module, with only one mass component (water) but also including the energy component. Water may be in liquid, vapour or two-phase conditions, and may transition between these states. Primary variables are pressure and temperature for single-phase conditions but switch to pressure and vapour saturation under two-phase conditions.

Fluid properties are calculated directly from the thermodynamic formulation for water (see :ref:`water_thermodynamics`).

The **"eos.primary.scale"** object contains values for customising the non-dimensionalisation of pressure and temperature primary variables. (Vapour saturation is already non-dimensional.) For example:

 .. code-block:: json

  {"eos": {"name": "we", "primary": {"scale": {"temperature": 20}}}}

selects the water/energy equation of state and overrides the non-dimensionalisation of temperatures, so that they are scaled by a factor of 20.

.. _water_ncg_eos:

Water / NCG EOS modules
=======================

These EOS modules simulate mixtures of water and non-condensible gases (NCGs), together with energy. They work in much the same way as the water / energy EOS ("we") apart from modifications to the fluid properties resulting from the presence of the non-condensible gas.

The primary variables for these EOS modules are as for the water / energy EOS, but with an added third variable, the partial pressure of the non-condensible gas.

The **"eos.primary.scale"** contains values for customising the non-dimensionalisation of pressure, temperature and gas partial pressure primary variables. Gas partial pressures can be scaled either by a fixed constant, as for the pressure and temperature variables, or by the total pressure (the default). This can be selected by setting the **"eos.primary.scale.partial_pressure"** to **"pressure"**. For example:

 .. code-block:: json

  {"eos": {"name": "wce", "primary": {"scale": {"partial_pressure": "pressure"}}}}

selects the water/CO\ :sub:`2`/energy equation of state, and specifies that CO\ :sub:`2` partial pressures should be non-dimensionalised by scaling by the total pressure. Setting the **"eos.primary.scale.partial_pressure"** value to a number specifies scaling by a fixed constant, as for pressure and temperature variables. For example:

 .. code-block:: json

  {"eos": {"name": "wae", "primary": {"scale": {"partial_pressure": 1e5}}}}

selects the water/air/energy equation of state, and specifies that partial pressures of air should be non-dimensionalised by scaling by a fixed factor of 10\ :sup:`5`.

.. add detail on how NCG mixture EOS modules work? - using Henry's derivative to compute energy of solution etc.

.. index:: equation of state (EOS); water / air / energy ("wae")
.. _water_air_energy_eos:

Water, air and energy ("wae")
-----------------------------

+-------------------------------+-------------------------------------------------------------------------+
|**abbreviated name**:          |"wae"                                                                    |
+-------------------------------+-------------------------------------------------------------------------+
|**component names**:           |["water", "air", "energy"]                                               |
+-------------------------------+-------------------------------------------------------------------------+
|**phase names**:               |["liquid", "vapour"]                                                     |
+-------------------------------+-------------------------------------------------------------------------+
|**primary variable names**:    |**single-phase**: ["pressure", "temperature", "air_partial_pressure"]    |
|                               +-------------------------------------------------------------------------+
|                               |**two-phase**: ["pressure", "vapour_saturation", "air_partial_pressure"] |
+-------------------------------+-------------------------------------------------------------------------+
|**default primary variables**: |[10\ :sup:`5` Pa, 20 :math:`^{\circ}`\ C, 0 Pa]                          |
+-------------------------------+-------------------------------------------------------------------------+
|**default region**:            |1 (liquid)                                                               |
+-------------------------------+-------------------------------------------------------------------------+
|**default eos.primary.scale**: |{"pressure": 1e6, "temperature": 100, "partial_pressure": "pressure"}    |
|                               |                                                                         |
+-------------------------------+-------------------------------------------------------------------------+
|**default output fluid         |["pressure", "temperature", "region", "air_partial_pressure",            |
|fields**:                      |"vapour_saturation"]                                                     |
+-------------------------------+-------------------------------------------------------------------------+

.. index:: equation of state (EOS); water / air / carbon dioxide ("wce")

Water, carbon dioxide and energy ("wce")
----------------------------------------

+-------------------------------+-------------------------------------------------------------------------+
|**abbreviated name**:          |"wce"                                                                    |
+-------------------------------+-------------------------------------------------------------------------+
|**component names**:           |["water", "CO2", "energy"]                                               |
+-------------------------------+-------------------------------------------------------------------------+
|**phase names**:               |["liquid", "vapour"]                                                     |
+-------------------------------+-------------------------------------------------------------------------+
|**primary variable names**:    |**single-phase**: ["pressure", "temperature", "CO2_partial_pressure"]    |
|                               +-------------------------------------------------------------------------+
|                               |**two-phase**: ["pressure", "vapour_saturation", "CO2_partial_pressure"] |
+-------------------------------+-------------------------------------------------------------------------+
|**default primary variables**: |[10\ :sup:`5` Pa, 20 :math:`^{\circ}`\ C, 0 Pa]                          |
+-------------------------------+-------------------------------------------------------------------------+
|**default region**:            |1 (liquid)                                                               |
+-------------------------------+-------------------------------------------------------------------------+
|**default eos.primary.scale**: |{"pressure": 1e6, "temperature": 100, "partial_pressure": "pressure"}    |
|                               |                                                                         |
+-------------------------------+-------------------------------------------------------------------------+
|**default output fluid         |["pressure", "temperature", "region", "CO2_partial_pressure",            |
|fields**:                      |"vapour_saturation"]                                                     |
+-------------------------------+-------------------------------------------------------------------------+
       
