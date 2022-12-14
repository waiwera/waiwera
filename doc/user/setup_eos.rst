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

   +-----------------------+----------+-------------------+-----------------------+
   |**name**               |**type**  |**default**        |**value**              |
   +-----------------------+----------+-------------------+-----------------------+
   |"name"                 |string    |"we"               |abbreviated EOS module |
   |                       |          |                   |name                   |
   +-----------------------+----------+-------------------+-----------------------+
   |"primary"              |object    |{}                 |primary variable       |
   |                       |          |                   |parameters             |
   |                       |          |                   |                       |
   |                       |          |                   |                       |
   +-----------------------+----------+-------------------+-----------------------+
   |"temperature"          |number    |20\                |constant temperature ( |
   |                       |          |:math:`^{\circ}`\ C|:math:`^{\circ}`\ C)   |
   |                       |          |                   |for :ref:`water_eos`   |
   |                       |          |                   |EOS                    |
   +-----------------------+----------+-------------------+-----------------------+
   |"permeability_modifier"|object    |{}                 |parameters for effect  |
   |                       |          |                   |of fluid on            |
   |                       |          |                   |permeability           |
   |                       |          |                   |                       |
   +-----------------------+----------+-------------------+-----------------------+

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

.. index:: equation of state (EOS), permeability modification
.. _fluid_permeability_modification:

Permeability modification
=========================

For some equations of state, the fluid state can change the effective local permeability. For example, when the :ref:`water_salt_eos` is used, a solid halite phase can be present when salt precipitates out of solution.

This effect can be represented using the **"eos.permeability_modifier"** value. This object a **"type"** string value which determines how the permeability is reduced as the effective porosity decreases. Its possible values are "none", "power" and "verma-pruess".

.. note::
   **JSON object**: fluid permeability modifier

   **JSON path**: eos.permeability_modifier

   +-------------+----------+-------------------+-----------------------+
   |**name**     |**type**  |**default**        |**value**              |
   +-------------+----------+-------------------+-----------------------+
   |"type"       |string    |"none"             |permeability modifier  |
   |             |          |                   |type                   |
   +-------------+----------+-------------------+-----------------------+
   |"exponent"   |number    |3 for power law; 2 |exponent :math:`n` for |
   |             |          |for Verma-Pruess   |power law; for         |
   |             |          |                   |Verma-Pruess: 2 for    |
   |             |          |                   |tubes in series or 3   |
   |             |          |                   |for fractures in series|
   +-------------+----------+-------------------+-----------------------+
   |"phir"       |number    |0.1                |for Verma-Pruess,      |
   |             |          |                   |parameter              |
   |             |          |                   |:math:`\phi_r`         |
   +-------------+----------+-------------------+-----------------------+
   |"gamma"      |number    |0.7                |for Verma-Pruess,      |
   |             |          |                   |parameter              |
   |             |          |                   |:math:`\Gamma`         |
   +-------------+----------+-------------------+-----------------------+

If the type is "none", there is no permeability reduction (the default). If the type is "power", a power-law relationship is used to determine the permeability reduction :math:`k/k_0`:

.. math::

   \frac{k}{k_0} = \left(\frac{\phi}{\phi_0}\right)^n

Here :math:`k` and :math:`k_0` represent the permeabilities before and after modification respectively, :math:`\phi_0` is the rock porosity and :math:`\phi` is the effective porosity, reduced by e.g. the presence of solid halite. The exponent :math:`n` typically takes values between 2 and 3.

The Verma-Pruess model [Verma-Pruess]_ for permeability modification is more complex and is based on considering the possible geometries of the pores and how reduced porosity at the throats of the pores can have a larger effect on the effective permeability. Here the permeability reduction is given by:

.. math::

   \frac{k}{k_0} = \theta^n \frac{1 - \Gamma + \Gamma / \omega^n}{1 - \Gamma + \Gamma \left(\frac{\theta}{\theta + \omega - 1} \right)^n}

where

.. math::

   \theta = \frac{\frac{\phi}{\phi_0} - \phi_r}{1 - \phi_r}

and

.. math::

   \omega = 1 + \frac{1}{\Gamma (1/\phi_r - 1)}

When the parameter :math:`n` takes the value 2, the pores are represented by a series of one-dimensional tubes, whereas when it takes the value 3, the pores are represented by parallel-plate fracture segments. The parameter :math:`\phi_r` is the fraction of the original porosity at which the permeability is reduced to zero, and the parameter :math:`\Gamma` is the fractional length of the pore bodies.

.. [Verma-Pruess] Verma, A. and Pruess, K. (1988). "Thermohydrologic conditions and silica redistribution near high-level nuclear wastes emplaced in saturated geological formations", J. Geophysical Research, 93, B2, 1159 - 1173.

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
.. _water_CO2_energy_eos:

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
       
.. _water_salt_eos:

Water / salt EOS module
=======================

This EOS module simulates mixtures of water and salt (NaCl), i.e. brine, together with energy. It works in much the same way as the :ref:`water_energy_eos` EOS apart from modifications to the fluid properties resulting from the presence of salt. The main thermodynamic properties (e.g. density and internal energy) of brine are calculated using the formulation of [Driesner]_.

Salt can be present in dissolved form in the liquid phase, under either single-phase liquid or two-phase conditions. It is assumed there is no salt present in the vapour phase.

Salt in the liquid phase may be present in concentrations up to a limit defined by the solubility of salt in water. This is temperature-dependent but under typical conditions the maximum salt mass fraction is approximately 0.3. At higher concentrations the salt will precipitate out into solid-phase salt (halite). Hence, the salt EOS modules have a "solid" phase as well as the liquid and vapour phases. This solid phase is not considered mobile, and is omitted from flux calculations across mesh faces. However, when solid halite is present it does reduce the pore space available for brine. It can also optionally reduce the effective permeability (see :ref:`fluid_permeability_modification`).

The primary variables for this EOS module are as for the water / energy EOS, but with an added third variable for salt. This variable represents salt mass fraction (in the liquid phase), unless there is solid-phase halite present, in which case it switches to the solid-phase saturation, i.e. the volume fraction of halite.

For the water/salt EOS module, the thermodynamic "region" has an expanded meaning to differentiate between fluid with and without solid-phase halite present. When halite is not present, the region has the same meaning as for pure water (see :ref:`thermodynamic_regions`). However when halite is present the region is incremented by 4:

+------+----------+
|Region|Conditions|
+------+----------+
|1     |Liquid, no|
|      |halite    |
+------+----------+
|2     |Vapour, no|
|      |halite    |
+------+----------+
|4     |Two-phase,|
|      |no halite |
+------+----------+
|5     |Liquid,   |
|      |halite    |
+------+----------+
|6     |Vapour,   |
|      |halite    |
+------+----------+
|8     |Two-phase,|
|      |halite    |
+------+----------+

(Note that, as for the :ref:`water_energy_eos` EOS, supercritical fluid (region 3) is not supported.)

.. [Driesner] Driesner, T. (2007). "The system H2O - NaCl. Part II: Correlations for molar volume, enthalpy, and isobaric heat capacity from 0 1000 C, 1 to 5000 bar, and 0 to 1 XNaCl. Geochimica et Cosmochimica Acta, 71, 4902 -- 4919.

.. index:: equation of state (EOS); water / salt / energy ("wse")
.. _water_salt_energy_eos:

Water, salt and energy ("wse")
------------------------------

+-------------------------------+-------------------------------------------------------------------------+
|**abbreviated name**:          |"wse"                                                                    |
+-------------------------------+-------------------------------------------------------------------------+
|**component names**:           |["water", "salt", "energy"]                                              |
+-------------------------------+-------------------------------------------------------------------------+
|**phase names**:               |["liquid", "vapour", "solid"]                                            |
+-------------------------------+-------------------------------------------------------------------------+
|**primary variable names**:    |**single-phase, no halite**: ["pressure", "temperature",                 |
|                               |"salt_mass_fraction"]                                                    |
|                               +-------------------------------------------------------------------------+
|                               |**two-phase, no halite**: ["pressure", "vapour_saturation",              |
|                               |"salt_mass_fraction"]                                                    |
|                               +-------------------------------------------------------------------------+
|                               |**single-phase, halite**: ["pressure", "temperature", "solid_saturation"]|
|                               |                                                                         |
|                               +-------------------------------------------------------------------------+
|                               |**two-phase, halite**: ["pressure", "vapour_saturation",                 |
|                               |"solid_saturation"]                                                      |
+-------------------------------+-------------------------------------------------------------------------+
|**default primary variables**: |[10\ :sup:`5` Pa, 20 :math:`^{\circ}`\ C, 0]                             |
+-------------------------------+-------------------------------------------------------------------------+
|**default region**:            |1 (liquid)                                                               |
+-------------------------------+-------------------------------------------------------------------------+
|**default eos.primary.scale**: |{"pressure": 1e6, "temperature": 100,                                    |
|                               |"salt_mass_fraction/solid_saturation": 1}                                |
+-------------------------------+-------------------------------------------------------------------------+
|**default output fluid         |["pressure", "temperature", "region", "vapour_saturation",               |
|fields**:                      |"liquid_salt_mass_fraction", "solid_saturation"]                         |
+-------------------------------+-------------------------------------------------------------------------+

.. _water_salt_ncg_eos:

Water / salt / NCG EOS modules
==============================

These EOS modules simulate mixtures of water and salt (NaCl), i.e. brine, together with non-condensible gases (NCGs) and energy. They are essentially a cross between the :ref:`water_salt_eos` and the :ref:`water_ncg_eos`, using the same formulations for salt and NCG thermodynamics. In addition, the "salting out" effect of salt concentration on the dissolution of NCG into the liquid phase is simulated.

The primary variables for these EOS modules are as for the water / salt EOS, but with an added fourth variable for NCG partial pressure.

.. index:: equation of state (EOS); water / salt / air / energy ("wsae")
.. _water_salt_air_energy_eos:

Water, salt, air and energy ("wsae")
------------------------------------

+-------------------------------+-------------------------------------------------------------------------+
|**abbreviated name**:          |"wsae"                                                                   |
+-------------------------------+-------------------------------------------------------------------------+
|**component names**:           |["water", "salt", "air", "energy"]                                       |
+-------------------------------+-------------------------------------------------------------------------+
|**phase names**:               |["liquid", "vapour", "solid"]                                            |
+-------------------------------+-------------------------------------------------------------------------+
|**primary variable names**:    |**single-phase, no halite**: ["pressure", "temperature",                 |
|                               |"salt_mass_fraction", "air_partial_pressure"]                            |
|                               +-------------------------------------------------------------------------+
|                               |**two-phase, no halite**: ["pressure", "vapour_saturation",              |
|                               |"salt_mass_fraction", "air_partial_pressure"]                            |
|                               +-------------------------------------------------------------------------+
|                               |**single-phase, halite**: ["pressure", "temperature", "solid_saturation",|
|                               |"air_partial_pressure"]                                                  |
|                               +-------------------------------------------------------------------------+
|                               |**two-phase, halite**: ["pressure", "vapour_saturation",                 |
|                               |"solid_saturation", "air_partial_pressure"]                              |
+-------------------------------+-------------------------------------------------------------------------+
|**default primary variables**: |[10\ :sup:`5` Pa, 20 :math:`^{\circ}`\ C, 0, 0 Pa]                       |
+-------------------------------+-------------------------------------------------------------------------+
|**default region**:            |1 (liquid)                                                               |
+-------------------------------+-------------------------------------------------------------------------+
|**default eos.primary.scale**: |{"pressure": 1e6, "temperature": 100,                                    |
|                               |"salt_mass_fraction/solid_saturation": 1, "partial_pressure": "pressure"}|
+-------------------------------+-------------------------------------------------------------------------+
|**default output fluid         |["pressure", "temperature", "region", "air_partial_pressure",            |
|fields**:                      |"vapour_saturation", "liquid_salt_mass_fraction", "solid_saturation"]    |
+-------------------------------+-------------------------------------------------------------------------+

.. index:: equation of state (EOS); water / salt / carbon dioxide / energy ("wsce")
.. _water_salt_CO2_energy_eos:

Water, salt, carbon dioxide and energy ("wsce")
-----------------------------------------------

+-------------------------------+-------------------------------------------------------------------------+
|**abbreviated name**:          |"wsce"                                                                   |
+-------------------------------+-------------------------------------------------------------------------+
|**component names**:           |["water", "salt", "CO2", "energy"]                                       |
+-------------------------------+-------------------------------------------------------------------------+
|**phase names**:               |["liquid", "vapour", "solid"]                                            |
+-------------------------------+-------------------------------------------------------------------------+
|**primary variable names**:    |**single-phase, no halite**: ["pressure", "temperature",                 |
|                               |"salt_mass_fraction", "CO2_partial_pressure"]                            |
|                               +-------------------------------------------------------------------------+
|                               |**two-phase, no halite**: ["pressure", "vapour_saturation",              |
|                               |"salt_mass_fraction", "CO2_partial_pressure"]                            |
|                               +-------------------------------------------------------------------------+
|                               |**single-phase, halite**: ["pressure", "temperature", "solid_saturation",|
|                               |"CO2_partial_pressure"]                                                  |
|                               +-------------------------------------------------------------------------+
|                               |**two-phase, halite**: ["pressure", "vapour_saturation",                 |
|                               |"solid_saturation", "CO2_partial_pressure"]                              |
+-------------------------------+-------------------------------------------------------------------------+
|**default primary variables**: |[10\ :sup:`5` Pa, 20 :math:`^{\circ}`\ C, 0, 0 Pa]                       |
+-------------------------------+-------------------------------------------------------------------------+
|**default region**:            |1 (liquid)                                                               |
+-------------------------------+-------------------------------------------------------------------------+
|**default eos.primary.scale**: |{"pressure": 1e6, "temperature": 100,                                    |
|                               |"salt_mass_fraction/solid_saturation": 1, "partial_pressure": "pressure"}|
+-------------------------------+-------------------------------------------------------------------------+
|**default output fluid         |["pressure", "temperature", "region", "CO2_partial_pressure",            |
|fields**:                      |"vapour_saturation", "liquid_salt_mass_fraction", "solid_saturation"]    |
+-------------------------------+-------------------------------------------------------------------------+
