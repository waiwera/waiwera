.. index:: Waiwera; numerical methods, numerical methods

*****************
How Waiwera works
*****************

.. index:: conservation equations; fluid
.. _conservation_equations:

Mass and energy conservation equations
======================================

Waiwera works by solving time-dependent conservation equations for mass (and usually energy) over the simulation domain. In general there may be several mass 'components' present in the problem, for example water and CO\ :sub:`2`. A separate mass conservation equation is solved for each mass component present.

For non-isothermal problems, an additional energy conservation equation is also solved, simultaneously with the mass conservation equations.

The mass and energy conservation equations over an arbitrary volume :math:`V` are:

.. math::
   :label: conservation

   \frac{d}{dt} \int_{V} {M^\kappa \, dV} = \int_{\partial V} {\mathbf{F^\kappa} \cdot \hat{n} \, dA} + \int_{V} {q^\kappa \, dV}

where :math:`\kappa = 1,\ldots N+1` (and :math:`N` is the number of mass components, e.g. 1 for pure water). Component :math:`N+1` is the energy component.

For component :math:`\kappa`, :math:`M^\kappa` is the mass or energy density in :math:`V`, :math:`\mathbf{F^\kappa}` is the flux of that component through the boundary :math:`\partial V` and :math:`q^\kappa` represents source or sink terms (per unit volume).

.. index:: numerical methods; finite volume discretisation
.. _finite_volume_discretisation:

Finite volume discretisation
============================

The flow domain is discretised using a finite volume mesh. In the :math:`i`\ :sup:`th` cell with volume :math:`V_i`, for each component we can define cell-averaged mass or energy densities :math:`M_i^\kappa` for each component, and cell-averaged source terms :math:`q_i^\kappa` as:

.. math::

   M_i^\kappa = \frac{1}{V_i} \int_{V_i} {M^\kappa \, dV}

   q_i^\kappa = \frac{1}{V_i} \int_{V_i} {q^\kappa \, dV}

We can also represent the flux integral in :eq:`conservation` as:

.. math::

   \int_{\partial V_i} {\mathbf{F^\kappa} \cdot \hat{n} \, dA} = \sum_j {A_{ij} F_{ij}^\kappa}

where :math:`A_{ij}` is the area of the face connecting cells :math:`i` and :math:`j`, and :math:`F_{ij}^\kappa` is the normal component of the flux through it.

Then the discretised conservation equations for cell :math:`V_i` can be written:

.. math::
   :label: discretised_conservation

   \frac{d}{dt} M_i^\kappa = \frac{1}{V_i} \sum_j {A_{ij} F_{ij}^\kappa} + q_i^\kappa

.. index:: thermodynamics; primary variables, primary variables
.. _primary_variables:

Primary variables
=================

Waiwera solves equation :eq:`discretised_conservation` for the thermodynamic state in each cell in the simulation mesh. The thermodynamic state in each cell is represented by a small set of 'primary variables', one for each conservation equation. The primary variables depend on the equation of state being used. They also depend on the phase conditions in the cell.

For example, for non-isothermal pure water simulations, there is just one mass component (water) and one energy component, so two primary variables are needed to describe the thermodynamic state. For single-phase conditions we may use pressure and temperature as the primary variables. All other fluid properties (e.g. density, viscosity etc.) can be calculated from the primary variables.

However, for two-phase conditions, the pressure and temperature are not independent, as they are related via the saturation curve. Hence, they cannot be used as primary variables to describe the thermodynamic state. For two-phase conditions, Waiwera uses pressure and vapour saturation as primary variables.

Because the choice of primary variables depends on the phase conditions, when the fluid in a cell changes phase, the primary variables must be changed.

In practice, Waiwera does not solve for the raw primary variables in each cell, but for non-dimensionalised primary variables. This improves numerical behaviour, mostly because the raw primary variables typicallyhave very different magnitudes from each other. The non-dimensionalised variables are related to their corresponding raw primary variables in most cases via a simple scaling by a fixed constant. These scaling parameters have default values but can also be customised for particular problems (see :ref:`primary_variable_parameters`).

.. index:: numerical methods; time evolution
.. _time_evolution:

Time evolution
==============

The discretised conservation equations :eq:`discretised_conservation` are of the form:

.. math::
   :label: RLeqn

   \frac{d}{dt} \mathbf{L}(t, \mathbf{Y}) = \mathbf{R}(t, \mathbf{Y})

where :math:`t` is time and :math:`\mathbf{Y}` is the vector of primary variables for all cells in the simulation mesh. Here :math:`\mathbf{L}` represents the cell-averaged mass and energy balances, as a function of time and the primary thermodynamic variables. Similarly, :math:`\mathbf{R}` represents inflows (per unit volume) into the cells from flows through the cell faces, together with sources and sinks within the cell.

Solving the set of ordinary differential equations :eq:`RLeqn` with respect to time, we can compute the time evolution of :math:`\mathbf{Y}`, the thermodynamic state of the entire discretised simulation domain.

When solving the conservation equations :eq:`discretised_conservation`, :math:`\mathbf{L}` and :math:`\mathbf{R}` are complicated, non-linear functions of the primary variables :math:`\mathbf{Y}`. Hence equation :eq:`RLeqn` must be solved numerically, computing the solution :math:`\mathbf{Y}` at discrete times.

Waiwera contains a module for the numerical solution of ordinary differential equations of the form :eq:`RLeqn`, using different numerical methods. The simplest of these is the 'backwards Euler' method, which discretises equation :eq:`RLeqn` as follows:

.. math::
   :label: beuler

   \frac{\mathbf{L}(t^{n+1}, \mathbf{Y}^{n+1}) - \mathbf{L}(t^n, \mathbf{Y}^n)}{\Delta t} = \mathbf{R}(t^{n+1}, \mathbf{Y}^{n+1})

where :math:`t^n` is the :math:`n^{th}` discretised time, and :math:`\Delta t` is the time step size, so that :math:`t^{n+1} = t^n + \Delta t`. For the backwards Euler method, at each time step we must solve equation :eq:`beuler` for the unknown new solution :math:`\mathbf{Y}^{n+1}`.

.. index:: numerical methods; function evaluations
.. _function_evaluations:

Function evaluations
====================

Waiwera needs to evaluate the functions :math:`\mathbf{L}` and :math:`\mathbf{R}` for any given set of primary variables (and time). The function :math:`\mathbf{L}`, representing the mass and energy densities :math:`M_i^\kappa` in the cells, is relatively straightforward to evaluate, by summing the contributions of the different phases. Considering a particular cell:

.. math::

   M_i^\kappa =
   \begin{cases}
   \phi_i \sum_p{S_{i,p} \rho_{i,p} X_{i,p}^\kappa} & \kappa \leq N \\
   (1 - \phi_i) \rho_i^r c_i T_i + \phi_i \sum_p {S_{i,p} \rho_{i,p} u_{i,p}} & \kappa = N + 1
   \end{cases}

where for cell :math:`i`, :math:`\phi_i` is the rock porosity, :math:`c_i` is rock specific heat and :math:`T_i` is temperature. For phase :math:`p`, :math:`S_{i,p}` is the phase saturation, :math:`\rho_{i,p}` is phase density, :math:`X_{i,p}^\kappa` is phase mass fraction for component :math:`\kappa` and :math:`u_{i,p}` is phase internal energy density.

The function :math:`\mathbf{R}`, representing fluxes into the cells, has contributions from source and sink terms (which are easily evaluated), and from fluxes through faces. This latter contribution is computed by summing the component face fluxes in each phase:

.. math::

   F_{ij}^\kappa = \sum_p{F_{ij,p}^\kappa}

where the phase fluxes are given by:

.. math::
   :label: flux

   F_{ij,p}^\kappa =
   \begin{cases}
   -\frac{k \cdot k^r_p}{\nu_p} X_p^\kappa (\frac{\partial P}{\partial n} - \rho_{ij,p} \mathbf{g}.\hat{n}) & \kappa \leq N \\
   -K \frac{\partial T}{\partial n} + \sum_{m=1}^{N} {\sum_p{h_p^m F_{ij,p}^m}} & \kappa = N + 1
   \end{cases}

Here :math:`k` is effective rock permeability normal to the face, :math:`K` is rock heat conductivity, :math:`P` is pressure, :math:`\mathbf{g}` is the gravity vector and :math:`\hat{n}` is the unit vector normal to the face. For phase :math:`p`, :math:`k^r_p` is the phase relative permeability, :math:`\nu_p` is phase kinematic viscosity, :math:`h_p^m` is phase enthalpy for component :math:`m` and :math:`\rho_{ij, p}` is the effective phase density on the face.

The normal gradients of pressure and temperature are evaluated by finite differencing across the phase, i.e. taking the difference between the values in the cells on either side of the face and dividing by the distance between the cell centres. This "two-point flux approximation" relies on the assumption that the mesh satisfies the "orthogonality criterion", i.e. that the line joining the cell centres is orthogonal to the face (see :ref:`mesh_orthogonality`).

When evaluating the phase fluxes using equation :eq:`flux`, the flow quantities :math:`k^r_p`, :math:`\rho_p`, :math:`\nu_p`, :math:`X_p^\kappa` and :math:`h_p^m` are "upstream weighted", i.e. their values are taken from the cell upstream from the face. This is needed for numerical stability. The rock permeability :math:`k` and heat conductivity :math:`K` on the face are evaluated using harmonic weighting of the values in the cells on either side of the face.

For the gravity term, Waiwera calculates the effective phase density on the face as a saturation-weighted average of the phase densities in the cells on either side:

.. math::

   \rho_{ij,p} = \frac{S_{i,p} \rho_{i,p} + S_{j,p} \rho_{j,p}}{S_{i,p} + S_{j,p}}

This formulation can be derived by considering a force balance over the two cells, and ensures a smooth variation in effective phase density on the face when the adjoining cells change phase. If both adjoining cells have the same saturation (e.g. in single-phase conditions) then this weighted average reduces to a simple arithmetic average of :math:`\rho_{i,p}` and :math:`\rho_{j,p}`.

.. index:: solver
.. _nonlinear_equations:

Solution of equations at each time step
=======================================

Regardless of the time stepping method used, the discretised equations to be solved at each time step (e.g. :eq:`beuler`) are non-linear. If we write them in a generic form:

.. math::
   :label: fx0

   \mathbf{f}(\mathbf{Y}) = \mathbf{0}

then at each time step we must solve this for the solution :math:`\mathbf{Y} = \mathbf{Y}^{n+1}`. Because of the non-linearity, it must be solved numerically using a non-linear solution technique such as Newton's method. Newton's method is an iterative method which starts from an initial estimate of the solution (here taken as :math:`\mathbf{Y} = \mathbf{Y}^n`) and adjusts the provisional solution :math:`\mathbf{Y}` at each iteration until equation :eq:`fx0` is satisfied, to within a pre-specified tolerance.

At each iteration, Newton's method adds an update :math:`\Delta \mathbf{Y}` to the provisional solution :math:`\mathbf{Y}` according to:

.. math::
   :label: newton

   \mathbf{J} \Delta \mathbf{Y} = -\mathbf{f}

where :math:`\mathbf{J}` is the Jacobian matrix of the function :math:`\mathbf{f}`, i.e. the matrix of partial derivatives of :math:`\mathbf{f}` with respect to :math:`\mathbf{Y}`.

.. index:: solver; non-linear, numerical methods; non-linear equations

At each iteration, the Newton update equation :eq:`newton` represents a large, sparse system of linear equations to be solved numerically. "Krylov subspace" iterative methods (e.g. conjugate gradient methods) are appropriate for solving such systems. For typical simulations of large flow problems, most of the computation time is spent in the solution of the linear equations.

.. index:: PETSc; SNES

Waiwera uses the "SNES" non-linear solver provided by the `PETSc <https://www.mcs.anl.gov/petsc/>`_ library to solve equation :eq:`fx0` at each time step. For problems in which the Jacobian matrix :math:`\mathbf{J}` is difficult to calculate, the SNES solver offers an option to calculate it automatically using finite differencing. In this case the Jacobian partial derivatives are evaluated approximately by adding small increments onto the primary variable vector :math:`\mathbf{y}` and re-evaluating the function :math:`\mathbf{f}`. Waiwera makes use of this approach to calculate the Jacobian matrix (see :ref: `jacobian_matrix`).

.. index:: PETSc; KSP, solver; linear, numerical methods; linear equations

The SNES solver in turn makes use of the "KSP" suite of linear solvers, also provided by PETSc, to solve the linear system :eq:`newton` at each Newton iteration.

.. _simulating_tracers:

Simulating tracers
==================

.. index:: conservation equations; tracer, tracers; conservation equations
.. _tracer_eqns:

Tracer conservation equations
-----------------------------

Waiwera can also simulate the movement of tracers, by solving separate tracer mass conservation equations. These are very similar to the fluid mass conservation equation :eq:`conservation`, but with an extra term representing decay:

.. math::
   :label: tracer_conservation

   \frac{d}{dt} \int_{V} {M^T \, dV} = \int_{\partial V} {\mathbf{F}^T \cdot \hat{n} \, dA} + \int_{V} {q^T \, dV} - \int_{V} {\alpha M^T \, dV}

where :math:`M^T` is the mass density of tracer in :math:`V`, :math:`\mathbf{F}^T` is is the flux of tracer through the boundary :math:`\partial V`, :math:`q^T` represents tracer source or sink terms (per unit volume) and :math:`\alpha` is the decay rate of tracer in :math:`V`.

.. index:: numerical methods; finite volume discretisation

Finite volume discretisation of the tracer conservation equation is done in exactly the same way as for the fluid mass conservation equation, by introducing the cell-averaged quantities:

.. math::

   M_i^T = \frac{1}{V_i} \int_{V_i} {M^T \, dV}

   q_i^T = \frac{1}{V_i} \int_{V_i} {q^T \, dV}

   \alpha_i = \frac{1}{V_i} \int_{V_i} {\alpha \, dV}

We can represent the flux integral in :eq:`tracer_conservation` as

.. math::

   \int_{\partial V_i} {\mathbf{F^T} \cdot \hat{n} \, dA} = \sum_j {A_{ij} F_{ij}^T}

and approximate the decay term by

.. math::

   \int_{V} {\alpha M^T \, dV} \approx \alpha_i V_i M_i^T

Then the discretised tracer conservation equation can be written as

.. math::
   :label: discretised_tracer_conservation

   \frac{d}{dt} M_i^T = \frac{1}{V_i} \sum_j {A_{ij} F_{ij}^T} + q_i^T - \alpha_i M_i^T

Multiple tracers may be simulated simultaneously, in which case a discretised conservation equation of this form is solved for each one. Each equation is solved for the mass fraction of tracer :math:`X_i^T` in each cell.

.. _tracer_terms:

Evaluating the terms in the tracer equations
--------------------------------------------

The mass term :math:`M_i^T` represents the tracer mass per unit volume in each cell. Currently, tracers in Waiwera are assumed to be specific to a particular fluid phase (e.g. liquid or vapour), so that the mass term may be evaluated as:

.. math::

   M_i^T = \phi_i S_{i,p} \rho_{i,p} X_i^T

where :math:`p` is the tracer phase. The tracer flux term can be written as the sum of two terms, representing advection and diffusion respectively:

.. math::

   F_{ij}^T = X_{ij}^T F_{ij,p} - [\phi \rho_p \tau]_{ij} D^T \frac{\partial X^T}{\partial n}

In the advection term, :math:`F_{ij, p}` is the total mass flux (over all mass components) in phase :math:`p`. The quantity :math:`X_{ij}^T` represents the effective tracer mass fraction on the face between cells :math:`i` and :math:`j`, and is upstream weighted (e.g. :math:`X_{ij}^T = X_i^T` if flow is from cell :math:`i` to cell :math:`j`).

In the diffusion term, :math:`D^T` is the tracer diffusion coefficient and the normal derivative of tracer mass fraction is evaluated by finite differencing across the face, as for the pressure and temperature gradients (see :ref:`function_evaluations`). The factor pre-multiplying the diffusion coefficient, :math:`[\phi \rho_p \tau]_{ij}`, is a product of three cell quantities, rock porosity :math:`\phi`, fluid phase density :math:`\rho_p` and a dimensionless "tortuosity" :math:`\tau`. The value of this total factor at the face is found using harmonic weighting of the values at the cells on either side.

The tortuosity in cell :math:`i` is found from a product of rock and fluid contributions, :math:`\tau_{i,0}` and :math:`\tau_{i,p}` respectively:

.. math::

   \tau_{i} = \tau_{i,0} \tau_{i,p}

At present a simple "constant-diffusivity" formulation is used for tortuosity, in which the rock contribution :math:`\tau_{i,0}` is identically 1 (i.e. no rock contribution), and the fluid contribution is equal to the phase saturation:

.. math::

   \tau_{i,p} = S_{i,p}

The tracer source term :math:`q_i^T` is specified for injection, while for production it is evaluated from:

.. math::

   q_i^T = X_i^T f_{i,p} q_i

where :math:`f_{i,p}` is the fluid flow fraction for phase :math:`p` in cell :math:`i` (computed based on the phase mobilities) and :math:`q_i` is the mass production rate.

.. _tracer_temp_decay:

Temperature-dependent decay
---------------------------

.. index:: tracers; temperature-dependent decay

The tracer decay rate :math:`\alpha_i` defaults to zero, but can be given a non-zero value (specific to each tracer, but independent of the cell). Alternatively, each tracer can be assigned a temperature-dependent decay rate, evaluated according to the Arrhenius equation:

.. math::
   :label: arrhenius

   \alpha_i = \alpha^0 e^{-E_0 / (R T_i^k)}

where :math:`\alpha^0` is a constant decay rate, :math:`E_0` is the "activation energy" for the tracer, :math:`R` is the universal gas constant and :math:`T_i^k` is the cell temperature in Kelvin.

Solution of tracer equations at each time step
----------------------------------------------

Like the discretised fluid conservation equations :eq:`discretised_conservation`, the discretised tracer conservation equation :eq:`discretised_tracer_conservation` can be written in the form :eq:`RLeqn`. However, the terms in the tracer equation are all linear in the tracer mass fractions :math:`X_i^T`. Hence, for tracer, the functions :math:`\mathbf{L}` and :math:`\mathbf{R}` are also linear.

.. index:: PETSc; KSP, solver; linear, numerical methods; linear equations
.. index:: numerical methods; time evolution

In Waiwera, tracers are assumed to be passive, that is, their presence does not affect the properties of the fluid through which they move. This means that the tracer mass conservation equations can be solved independently from the fluid conservation equations. At each time step, the nonlinear flow equations are solved as described above (see :ref:`nonlinear_equations`), using a nonlinear solver, after which the tracer mass fractions can be found by solving a single auxiliary set of linear equations. Again, the "KSP" suite of linear equation solvers provided by PETSc are used for this.
