.. index:: equation of state (EOS), simulation; EOS
.. _eos:

*****************
Equation of state
*****************

In different simulations there may be different combinations of mass and energy components present. In most subsurface simulations, water is present, but there may be others such as carbon dioxide or air. In non-isothermal simulations it is also necessary to model the energy component, to account for the effects of temperature on fluid properties. The set of equations describing the behaviour of the particular mixture of mass and energy components present under a range of thermodynamic conditions is referred to as the "equation of state" (EOS) module.

Given the primary thermodynamic variables, the EOS module calculates bulk fluid properties such as pressure and temperature, as well as the properties of the individual phases present, such as density, internal energy, viscosity etc. It also checks the primary variables to make sure they have not gone outside acceptable physical bounds, and handles primary variable switching when phase transitions occur.

In the Waiwera JSON input file, the **"eos"** value specifies the equation of state module to be used for the simulation. This can be either a string containing the abbreviated name of the EOS, or an object containing a **"name"** value. In general it is not necessary to specify the EOS as an object unless other EOS parameters besides the name need to be set.

.. admonition:: JSON input

   **JSON object**: equation of state

   **JSON path**: eos

   +--------------------------------+----------+-------------------+-----------------------+
   |**name**                        |**type**  |**default**        |**value**              |
   +--------------------------------+----------+-------------------+-----------------------+
   |"name"                          |string    |"we"               |abbreviated EOS module |
   |                                |          |                   |name                   |
   +--------------------------------+----------+-------------------+-----------------------+
   |"primary"                       |object    |{}                 |primary variable       |
   |                                |          |                   |parameters             |
   |                                |          |                   |                       |
   |                                |          |                   |                       |
   +--------------------------------+----------+-------------------+-----------------------+
   |"temperature"                   |number    |20\                |constant temperature ( |
   |                                |          |:math:`^{\circ}`\ C|:math:`^{\circ}`\ C)   |
   |                                |          |                   |for :ref:`water_eos`   |
   |                                |          |                   |EOS                    |
   +--------------------------------+----------+-------------------+-----------------------+
   |"permeability_modifier"         |object |  |{}                 |parameters for effect  |
   |                                |``null``  |                   |of fluid on            |
   |                                |          |                   |permeability           |
   |                                |          |                   |                       |
   +--------------------------------+----------+-------------------+-----------------------+
   |"relative_permeability_modifier"|object |  |depends on EOS     |parameters for effect  |
   |                                |``null``  |                   |of fluid on relative   |
   |                                |          |                   |permeability           |
   |                                |          |                   |                       |
   +--------------------------------+----------+-------------------+-----------------------+
   |"conditions"                    |string    |depends on EOS     |alternative methods of |
   |                                |          |                   |specifying initial and |
   |                                |          |                   |boundary conditions    |
   |                                |          |                   |                       |
   +--------------------------------+----------+-------------------+-----------------------+

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

.. admonition:: JSON input

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

This effect can be represented using the **"eos.permeability_modifier"** value. This object has a **"type"** string value which determines how the permeability is reduced as the effective porosity decreases. Its possible values are "none", "power" and "verma-pruess".

.. admonition:: JSON input

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

If the type is "none", there is no permeability reduction (the default). This can also be achieved by setting the **"eos.permeability_modifier"** value to ``null``. If the type is "power", a power-law relationship is used to determine the permeability reduction :math:`k/k_0`:

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

Relative permeability modification
==================================

For some equations of state, the fluid state can change the effective local relative permeability. For example, when the :ref:`supercritical_water_energy_eos` EOS is used, relative permeabilities may be modified so that they approach simple saturation functions as the critical point is approached.

This effect can be represented using the **"eos.relative_permeability_modifier"** value. This object has a **"type"** string value which determines how the relative permeability is modified according to temperature. Its possible values are "none" and "linear". If the type is "none", there is no relative permeability modification (the default for most equations of state). This can also be achieved by setting the **"eos.relative_permeability_modifier"** value to ``null``.

.. admonition:: JSON input

   **JSON object**: fluid relative permeability modifier

   **JSON path**: eos.relative_permeability_modifier

   +---------------------+----------+-------------------+-----------------------+
   |**name**             |**type**  |**default**        |**value**              |
   +---------------------+----------+-------------------+-----------------------+
   |"type"               |string    |"none"             |relative permeability  |
   |                     |          |                   |modifier type          |
   +---------------------+----------+-------------------+-----------------------+
   |"minimum_temperature"|number    |350                |minimum temperature    |
   |                     |          |                   |:math:`T_0` for        |
   |                     |          |                   |relative permeability  |
   |                     |          |                   |modification           |
   +---------------------+----------+-------------------+-----------------------+

If the type is "linear" (the default for the :ref:`supercritical_water_energy_eos` EOS), the effective relative permeability :math:`K^r_p` is given by:

.. math::

   f = \begin{cases}
   k^r_p & T \le T_0 \\
   K^r_p = (1 - \xi) k^r_p + \xi S_p & T > T_0
   \end{cases}

where :math:`k^r_p` is the original unmodified relative permeability, :math:`T` is temperature (:math:`^{\circ}`\ C), :math:`T_0` is a specified minimum temperature (:math:`^{\circ}`\ C), :math:`\xi = (T - T_0) / (T_c - T_0)`, :math:`T_c` is critical temperature of water (:math:`^{\circ}`\ C) and :math:`S_p` is the saturation of phase :math:`p`. As the critical point is approached, the properties of the liquid and vapour phases converge, so that the relative permeability functions approach simple saturation functions (:math:`K^r_p = S_p`) and do not include any other interactions between phases.

.. index:: simulation; initial conditions, initial conditions
.. index:: simulation; boundary conditions, boundary conditions

Alternative initial and boundary conditions
===========================================

In general, :ref:`initial_conditions` and :ref:`boundary_conditions` are set by specifying the primary variables for the equation of state being used. However, for convenience, some equations of state allow the user to specify them using variables different from the primary variables. This can be done using the **"eos.conditions"** value. Its possible values depend on the EOS being used.

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
|**regions supported**:         |1                         |
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
|**regions supported**:         |1, 2, 4                                           |
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

The primary variables for these EOS modules are as for the water / energy EOS, but with an added third variable, the partial pressure of the non-condensible gas. (The first variable, pressure, now represents the total pressure, not the partial pressure of water.)

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
|**regions supported**:         |1, 2, 4                                                                  |
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
|**regions supported**:         |1, 2, 4                                                                  |
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

(Note that, as for the :ref:`water_energy_eos` EOS, water thermodynamic regions 3 and 5 are not supported.)

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
|**regions supported**:         |1, 2, 4 (regions 5, 6, 8 indicate halite)                                |
+-------------------------------+-------------------------------------------------------------------------+
|**default region**:            |1 (liquid, no halite)                                                    |
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

The primary variables for these EOS modules are as for the water / salt EOS, but with an added fourth variable for NCG partial pressure. (As for the :ref:`water_ncg_eos`, the first variable represents total pressure, not partial pressure of water.)

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
|**regions supported**:         |1, 2, 4 (regions 5, 6, 8 indicate halite)                                |
+-------------------------------+-------------------------------------------------------------------------+
|**default region**:            |1 (liquid, no halite)                                                    |
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
|**regions supported**:         |1, 2, 4 (regions 5, 6, 8 indicate halite)                                |
+-------------------------------+-------------------------------------------------------------------------+
|**default region**:            |1 (liquid, no halite)                                                    |
+-------------------------------+-------------------------------------------------------------------------+
|**default eos.primary.scale**: |{"pressure": 1e6, "temperature": 100,                                    |
|                               |"salt_mass_fraction/solid_saturation": 1, "partial_pressure": "pressure"}|
+-------------------------------+-------------------------------------------------------------------------+
|**default output fluid         |["pressure", "temperature", "region", "CO2_partial_pressure",            |
|fields**:                      |"vapour_saturation", "liquid_salt_mass_fraction", "solid_saturation"]    |
+-------------------------------+-------------------------------------------------------------------------+

.. _supercritical_eoses:

Supercritical water EOS modules
===============================

.. index:: equation of state (EOS); supercritical water / energy ("se")
.. _supercritical_water_energy_eos:

Supercritical water and energy ("se")
-------------------------------------

+-------------------------------+--------------------------------------------------+
|**abbreviated name**:          |"se"                                              |
+-------------------------------+--------------------------------------------------+
|**component names**:           |["water", "energy"]                               |
+-------------------------------+--------------------------------------------------+
|**phase names**:               |["liquid", "vapour", "supercritical"]             |
+-------------------------------+--------------------------------------------------+
|**primary variable names**:    |**regions 1, 2, 5**: ["pressure", "temperature"]  |
|                               +--------------------------------------------------+
|                               |**region 3**: ["density", "temperature"]          |
|                               +--------------------------------------------------+
|                               |**region 4**: ["pressure", "vapour_saturation"]   |
+-------------------------------+--------------------------------------------------+
|**default primary variables**: |[10\ :sup:`5` Pa, 20 :math:`^{\circ}`\ C]         |
|                               |                                                  |
+-------------------------------+--------------------------------------------------+
|**regions supported**:         |1, 2, 3, 4, 5                                     |
+-------------------------------+--------------------------------------------------+
|**default region**:            |1 (liquid)                                        |
+-------------------------------+--------------------------------------------------+
|**default eos.primary.scale**: |{"pressure": 1e6, "temperature": 100, "density":  |
|                               |322}                                              |
+-------------------------------+--------------------------------------------------+
|**default output fluid         |["pressure", "temperature", "region",             |
|fields**:                      |"vapour_saturation", "liquid_density",            |
|                               |"vapour_density", "supercritical_density",        |
|                               |"liquidlike_fraction", "supercritical_phases"]    |
+-------------------------------+--------------------------------------------------+

This is based on the :ref:`water_energy_eos` EOS, but extends its capabilities to supercritical water and high-temperature steam. It can only be used in conjunction with the IAPWS-97 thermodynamics module (see :ref:`water_thermodynamics`). Whereas the "we" EOS is limited to liquid water, dry steam and two-phase conditions, with liquid water and two-phase only simulated below temperatures of 350 :math:`^{\circ}`\ C, the "se" EOS module can also simulate IAPWS-97 region 3 (see :ref:`thermodynamic_regions`), which covers near-critical and supercritical fluids, and region 5. Hence all pressures and temperatures up to 100 MPa and 800 :math:`^{\circ}`\ C can be simulated, as well as temperatures up to 2000 :math:`^{\circ}`\ C for pressures below 50 MPa.

The primary variables for this EOS are the same as those for the "we" EOS in regions 1, 2 and 4 (and region 5 primary variables are the same as for region 2). For region 3, it is not possible to use pressure and temperature as primary variables, as the thermodynamic equations are poorly behaved near the critical point when expressed as functions of these variables. Instead, Waiwera follows the IAPWS-97 formulation and uses density and temperature as primary variables in region 3.

However, for convenience it is possible to specify region 3 initial and boundary conditions using the more familiar pressure and temperature variables, if desired. This can be done via the **"eos.conditions"** value in the Waiwera input JSON file. This is a string value, and setting its value to "pressure" means that all initial and boundary conditions are interpreted as pressures and temperatures (and are converted internally to densities and temperatures).

As supercritical fluid is single-phase but can behave in a liquid-like or vapour-like way, or somewhere in between the two (see :ref:`supercritical_thermodynamics`), EOS "se" results for supercritical cells are given not in terms of liquid or vapour phases but as a separate third "supercritical" phase. This has the usual phase properties such as density, viscosity etc. and it is possible to include results for them in the Waiwera HDF5 output in the usual way (see :ref:`output_fluid_fields`). By default the HDF5 output also includes cell values for the :ref:`liquidlike_fraction` of the supercritical fluid. (For convenience, results for liquid-like fraction are also given for sub-critical cells, with the value set to 1 for liquid water, 0 for dry steam and to the liquid phase saturation for two-phase fluid.)

.. index:: equation of state (EOS); supercritical water / air / energy ("sae")
.. _supercritical_water_air_energy_eos:

Supercritical water, air and energy ("sae")
-------------------------------------------

+-------------------------------+-------------------------------------------------------------------------+
|**abbreviated name**:          |"sae"                                                                    |
+-------------------------------+-------------------------------------------------------------------------+
|**component names**:           |["water", "air", "energy"]                                               |
+-------------------------------+-------------------------------------------------------------------------+
|**phase names**:               |["liquid", "vapour", "supercritical"]                                    |
+-------------------------------+-------------------------------------------------------------------------+
|**primary variable names**:    |**regions 1, 2, 5**: ["pressure", "temperature", "air_partial_pressure"] |
|                               +-------------------------------------------------------------------------+
|                               |**region 3**: ["density", "temperature", "air_partial_pressure"]         |
|                               |                                                                         |
|                               +-------------------------------------------------------------------------+
|                               |**region 4**: ["pressure", "vapour_saturation", "air_partial_pressure"]  |
+-------------------------------+-------------------------------------------------------------------------+
|**default primary variables**: |[10\ :sup:`5` Pa, 20 :math:`^{\circ}`\ C, 0 Pa]                          |
+-------------------------------+-------------------------------------------------------------------------+
|**regions supported**:         |1, 2, 3, 4, 5                                                            |
+-------------------------------+-------------------------------------------------------------------------+
|**default region**:            |1 (liquid)                                                               |
+-------------------------------+-------------------------------------------------------------------------+
|**default eos.primary.scale**: |{"pressure": 1e6, "temperature": 100, "partial_pressure": "pressure",    |
|                               |"density": 322}                                                          |
+-------------------------------+-------------------------------------------------------------------------+
|**default output fluid         |["pressure", "temperature", "region", "vapour_saturation",               |
|fields**:                      |"liquid_density", "vapour_density", "supercritical_density",             |
|                               |"liquidlike_fraction", "supercritical_phases", "air_partial_pressure"]   |
+-------------------------------+-------------------------------------------------------------------------+

This combines the :ref:`supercritical_water_energy_eos` and :ref:`water_air_energy_eos` EOS modules, so that mixtures of sub- or super-critical water and air can be simulated.

The primary variables are the same as those for the :ref:`supercritical_water_energy_eos` EOS, with a third variable added for the partial pressure of air. Note that for regions 1, 2, 4 and 5 the pressure variable represents the total pressure (not partial pressure of water), but for region 3 the density variable represents the water density (not total density). However, as for the :ref:`supercritical_water_energy_eos` EOS, it is possible to specify region 3 initial and boundary conditions using pressures instead of densities, by setting the **"eos.conditions"** value in the Waiwera input JSON file to "pressure". (In this case, the pressures for region 3 also represent total pressure.)

Note also that for this EOS, the computed :ref:`liquidlike_fraction` values represent the liquid-like fraction of the water component only, independently of the air component.

The behaviour of mixtures of water and air is not yet well understood at high temperatures and pressures, so Waiwera uses an approximate formulation for supercritical water in the presence of air:

- the density of air is computed under all conditions using the ideal gas law, which is valid only if the partial pressure of air is much lower than the air critical pressure (approximately 3.79 MPa)
- Henry's Law is used to model the dissolution of air into liquid, with the Henry's coefficient represented as a polynomial function of temperature. This polynomial is valid only for temperatures below the critical point of water. In Waiwera, for temperatures greater than 370 :math:`^{\circ}`\ C, the Henry's coefficient at 370 :math:`^{\circ}`\ C is used.
- similarly, the air energy of solution is computed from the temperature derivative of the Henry's coefficient, so for temperatures greater than 370 :math:`^{\circ}`\ C, the energy of solution at 370 :math:`^{\circ}`\ C is used
- Dalton's Law is assumed at all temperatures and pressures (i.e. total pressure is equal to the sum of the water and air partial pressures)
- for completely liquid-like or vapour-like supercritical water (see :ref:`supercritical_thermodynamics`), the effects of air on the fluid properties are modelled in the same way as for sub-critical liquid water or vapour
- in :ref:`widom_delta`, in which the supercritical water behaviour is in between liquid-like and vapour-like, the fluid properties are computed first as if the water component were completely liquid-like and again as if it were completely vapour-like. The actual fluid properties are then linearly interpolated between these two according to the :ref:`liquidlike_fraction`. 

For most applications, air enters the model through an atmospheric boundary and is present mainly near the surface, while concentrations of air in the deeper parts of the model, where supercritical water may be found, are very low. In such cases, the above approximations used to represent the properties of supercritical water containing air have little effect on the model results. However, they will not give accurate results if the model does contain significant concentrations of air in supercritical water. 
