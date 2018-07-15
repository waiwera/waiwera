*****************
How Waiwera works
*****************

.. focus on what user needs to know about to set up simulation and interpret results

.. two-point flux approximation, upstream weighting?
.. MINC for fractured media?

Mass and energy conservation equations
======================================

Waiwera works by solving time-dependent conservation equations over the simulation domain. In general there may be several mass 'components' present in the problem, for example water and CO\ :sub:`2`. A separate mass conservation equation is solved for each mass component present.

For non-isothermal problems, an additional energy conservation equation is also solved, simultaneously with the mass conservation equations.

The mass and energy conservation equations over an arbitrary volume :math:`V` are:

.. math::
   :label: conservation

   \frac{d}{dt} \int_{V} {M^c \, dV} = \int_{\partial V} {\mathbf{F^c} \cdot \hat{n} \, dA} + \int_{V} {q^c \, dV}

where :math:`c = 1,\ldots C+1` (and :math:`C` is the number of mass components, e.g. 1 for pure water). Component :math:`C+1` is the energy component.

For component :math:`c`, :math:`M^c` is the mass or energy density in :math:`V_n`, :math:`\mathbf{F^c}` is the flux and :math:`q^c` represents source or sink terms (per unit volume).

.. _finite_volume_discretisation:

Finite volume discretisation
============================

The flow domain is discretised using a finite volume mesh made up of :math:`N` cells. For each component we can define cell-averaged mass or energy densities :math:`M_n^c` for each component, and cell-averaged source terms :math:`q_n^c` as:

.. math::

   M_n^c = \frac{1}{V_n} \int_{V_n} {M^c \, dV}

   q_n^c = \frac{1}{V_n} \int_{V_n} {q^c \, dV}

We can also represent the flux integral in :eq:`conservation` as:

.. math::

   \int_{\partial V_n} {\mathbf{F^c} \cdot \hat{n} \, dA} = \sum_m {A_{nm} F_{nm}^c}

where :math:`A_{nm}` is the area of the face connecting cells :math:`n` and :math:`m`, and :math:`F_{nm}^c` is the normal component of the flux through it.

Then the discretised conservation equations for cell :math:`V_n` can be written:

.. math::
   :label: discretised_conservation

   \frac{d}{dt} M_n^c = \frac{1}{V_n} \sum_m {A_{nm} F_{nm}^c} + q_n^c

.. _primary_variables:

Primary variables
=================

Waiwera solves equation :eq:`discretised_conservation` for the thermodynamic state in each cell in the simulation mesh. The thermodynamic state in each cell is represented by a small set of 'primary variables', one for each conservation equation. The primary variables depend on the equation of state being used. The primary variables also depend on the phase conditions in the cell.

For example, for non-isothermal pure water simulations, there is just one mass component (water) and one energy component, so two primary variables are needed to describe the thermodynamic state. For single-phase conditions we may use pressure and temperature as the primary variables. All other fluid properties (e.g. density, viscosity etc.) can be calculated from the primary variables.

However, for two-phase conditions, the pressure and temperature are not independent, as they are related via the saturation curve. Hence, they cannot be used as primary variables to describe the thermodynamic state. For two-phase conditions, Waiwera uses pressure and vapour saturation as primary variables.

Because the choice of primary variables depends on the phase conditions, when the fluid in a cell changes phase, the primary variables must be changed.

Time evolution
==============

The discretised conservation equations :eq:`discretised_conservation` are of the form:

.. math::
   :label: RLeqn

   \frac{d}{dt} \mathit{L}(t, \mathbf{Y}) = \mathit{R}(t, \mathbf{Y})

where :math:`t` is time and :math:`\mathbf{Y}` is the vector of primary variables for all cells in the simulation mesh (of total length :math:`N(C+1)`). Here :math:`L` represents the cell-averaged mass and energy balances, as a function of time and the primary thermodynamic variables. Similarly, :math:`R` represents inflows into the cells (per unit volume) from flows through the cell faces, together with sources and sinks within the cell.

Solving the set of ordinary differential equations :eq:`RLeqn` with respect to time, we can compute the time evolution of :math:`\mathbf{Y}`, the thermodynamic state of the entire discretised simulation domain.

For solving the conservation equations, :math:`L` and :math:`R` are complex, non-linear functions of the primary variables :math:`\mathbf{Y}`. Hence equation :eq:`RLeqn` must be solved numerically, computing the solution :math:`\mathbf{Y}` at discrete times.

Waiwera contains a module for the numerical solution of ordinary differential equations of the form :eq:`RLeqn`, using different numerical methods. The simplest of these is the 'backwards Euler' method, which discretises equation :eq:`RLeqn` as follows:

.. math::
   :label: beuler

   \frac{L(t^{n+1}, \mathbf{Y}^{n+1}) - L(t^n, \mathbf{Y}^n)}{\Delta t} \approx R(t^{n+1}, \mathbf{Y}^{n+1})

where :math:`t^n` is the :math:`n^{th}` discretised time, and :math:`\Delta t` is the time step size, so that :math:`t^{n+1} = t^n + \Delta t`. For the backwards Euler method, at each time step we must solve equation :eq:`beuler` for the unknown new solution :math:`\mathbf{Y}^{n+1}`.

.. _function_evaluations:

Function evaluations
====================

Waiwera needs to evaluate the functions :math:`L` and :math:`R` for any given set of primary variables (and time). The function :math:`L`, representing the mass and energy densities :math:`M_n^c` in the cells, is relatively straightforward to evaluate, by summing the contributions of the different phases. Considering a particular cell:

.. math::

   M_n^c =
   \Biggl \lbrace
   {
   \phi_n \sum_p{S_p \rho_p X_p^c}, c \leq C
   \atop
   (1 - \phi_n) \rho_{r} c_{r} T + \phi_n \sum_p {S_p \rho_p u_p}, c = C + 1
   }

where the :math:`p` subscripts refer to phases, and the :math:`r` subscripts refer to rock properties. Here :math:`\phi_n` is the porosity in the cell, :math:`S` is phase saturation, :math:`\rho` is density, :math:`X` is mass fraction, :math:`u` is internal energy density, :math:`c_r` is the rock specific heat and :math:`T` is temperature.

The function :math:`R`, representing fluxes into the cells, has contributions from source and sink terms (which are easily evaluated), and from fluxes through faces. This latter contribution is computed by summing the component face fluxes in each phase:

.. math::

   F_{nm}^c = \sum_p{F_p^c}

where the phase fluxes are given by:

.. math::
   :label: flux

   F_p^c =
   \Biggl \lbrace
   {
   -k \frac{k_r^p}{\mu_p} \rho_p X_p^c (\frac{\partial P}{\partial n} - \bar{\rho}_p \mathbf{g}.\hat{n}), c \leq C
   \atop
   -K \frac{\partial T}{\partial n} + \sum_{i=1}^{C} {\sum_p{h_p^i F_p^i}} , c = C + 1
   }

Here :math:`k` is effective permeability normal to the face, :math:`k_r` is relative permeability, :math:`\mu` is viscosity, :math:`P` is pressure, :math:`\mathbf{g}` is the gravity vector, :math:`K` is rock heat conductivity and :math:`h` is enthalpy. :math:`\hat{n}` is the unit vector normal to the face, and :math:`\bar{\rho}_p` is the effective phase density on the face.

The normal gradients of pressure and temperature are evaluated by finite differencing across the phase, i.e. taking the difference between the values in the cells on either side of the face and dividing by the distance between the cell centres. This "two-point flux approximation" relies on the assumption that the mesh satisfies the "orthogonality criterion", i.e. that the line joining the cell centres is orthogonal to the face.

When evaluating the phase fluxes using equation :eq:`flux`, the flow quantities :math:`k_r`, :math:`\rho_p`, :math:`\mu`, :math:`X_c^p` and :math:`h_p` are "upstream weighted", i.e. their values are taken from the cell upstream from the face. This is needed for numerical stability. The rock permeability :math:`k` and heat conductivity :math:`K` on the face are evaluated using harmonic weighting of the values in the cells on either side of the face.

For the gravity term, Waiwera calculates the effective phase density on the face as a saturation-weighted average of the phase densities in the cells on either side:

.. math::

   \bar{\rho}_p = \frac{S_p^1 \rho_p^1 + S_p^2 \rho_p^2}{S_p^1 + S_p^2}

where :math:`S_p^1`, :math:`S_p^2` are the phase saturations in the two cells, and :math:`\rho_p^1`, :math:`\rho_p^2` are the corresponding phase densities. This formulation ensures a smooth variation in effective phase density on the face when the adjoining cells change phase. If both adjoining cells have the same saturation (e.g. in single-phase conditions) then this weighted average reduces to a simple arithmetic average.


Solution of equations at each time step
=======================================

Regardless of the time stepping method used, the discretised equations to be solved at each time step (e.g. :eq:`beuler`) are non-linear. If we write them in a generic form:

.. math::
   :label: fx0

   f(\mathbf{y}) = \mathbf{0}

then at each time step we must solve this for the solution :math:`\mathbf{y} = \mathbf{Y}^{n+1}`. Because of the non-linearity, it must be solved numerically using a non-linear solution technique such as Newton's method. This is an iterative method which starts from an initial estimate of the solution (here taken as :math:`\mathbf{y} = \mathbf{Y}^n`) and adjusts the provisional solution :math:`\mathbf{y}` at each iteration until equation :eq:`fx0` is satisfied, to within a pre-specified tolerance.

At each iteration, Newton's method adds an update :math:`\Delta \mathbf{y}` to the provisional solution :math:`\mathbf{y}` according to:

.. math::
   :label: newton

   \mathbf{J} \Delta \mathbf{y} = -f(\mathbf{y})

where :math:`\mathbf{J}` is the Jacobian matrix of the function :math:`f`, i.e. the matrix of partial derivatives of :math:`f` with respect to :math:`\mathbf{y}`.

At each iteration, the Newton update equation :eq:`newton` represents a large, sparse system of linear equations to be solved numerically. "Krylov subspace" iterative methods (e.g. conjugate gradient methods) are appropriate for solving such systems. For typical simulations of large problems, most of the computation time is spent in the solution of the linear equations.

Waiwera uses the "SNES" non-linear solver provided by the `PETSc <https://www.mcs.anl.gov/petsc/>`_ library to solve equation :eq:`fx0` at each time step. For problems in which the Jacobian matrix :math:`\mathbf{J}` is difficult to calculate, the SNES solver offers an option to calculate it automatically using finite differencing. In this case the Jacobian partial derivatives are evaluated approximately by adding small increments onto the primary variable vector :math:`\mathbf{y}` and re-evaluating the function :math:`f`. Waiwera makes use of this approach to calculate the Jacobian matrix.

The SNES solver in turn makes use of the "KSP" suite of linear solvers, also provided by PETSc, to solve the linear system :eq:`newton` at each Newton iteration.
