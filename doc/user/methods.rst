*****************
How Waiwera works
*****************

.. focus on what user needs to know about to set up simulation and interpret results

.. two-point flux approximation, upstream weighting?
.. MINC for fractured media?

Mass and energy conservation equations
**************************************

Waiwera works by solving time-dependent conservation equations over the simulation domain. In general there may be several mass 'components' present in the problem, for example water and CO\ :sub:`2`. A separate mass conservation equation is solved for each mass component present.

For non-isothermal problems, an additional energy conservation equation is also solved, simultaneously with the mass conservation equations.

The mass and energy conservation equations over an arbitrary volume :math:`V` are:

.. math::
   :label: conservation

   \frac{d}{dt} \int_{V} {M^c \, dV} = \int_{\partial V} {\vec{F^c} \cdot \hat{n} \, dA} + \int_{V} {q^c \, dV}

where :math:`c = 1,\ldots C+1` (and :math:`C` is the number of mass components, e.g. 1 for pure water). Component :math:`C+1` is the energy component.

For component :math:`c`, :math:`M^c` is the mass or energy density in :math:`V_n`, :math:`\vec{F^c}` is the flux and :math:`q^c` represents source or sink terms.

Finite volume discretisation
****************************

The flow domain is discretised using a finite volume mesh made up of :math:`N` cells. For each component we can define cell-averaged mass or energy densities :math:`M_n^c` for each component, and cell-averaged source terms :math:`q_n^c` as:

.. math::

   M_n^c = \frac{1}{V_n} \int_{V_n} {M^c \, dV}

   q_n^c = \frac{1}{V_n} \int_{V_n} {q^c \, dV}

We can also represent the flux integral in :eq:`conservation` as:

.. math::

   \int_{\partial V_n} {\vec{F^c} \cdot \hat{n} \, dA} = \sum_m {A_{nm} F_{nm}^c}

where :math:`A_{nm}` is the area of the face connecting cells :math:`n` and :math:`m`, and :math:`F_{nm}^c` is the normal component of the flux through it.

Then the discretised conservation equations for cell :math:`V_n` can be written:

.. math::
   :label: discretised_conservation

   \frac{d}{dt} M_n^c = \frac{1}{V_n} \sum_m {A_{nm} F_{nm}^c} + q_n^c

Primary variables
*****************

Waiwera solves equation :eq:`discretised_conservation` for the thermodynamic state in each cell in the simulation mesh. The thermodynamic state in each cell is represented by a small set of 'primary variables', one for each conservation equation. The primary variables depend on the equation of state being used. The primary variables also depend on the phase conditions in the cell.

For example, for non-isothermal pure water simulations, there is just one mass component (water) and one energy component, so two primary variables are needed to describe the thermodynamic state. For single-phase conditions we may use pressure and temperature as the primary variables. All other fluid properties (e.g. density, viscosity etc.) can be calculated from the primary variables.

However, for two-phase conditions, the pressure and temperature are not independent, as they are related via the saturation curve. Hence, they cannot be used as primary variables to describe the thermodynamic state. For two-phase conditions, Waiwera uses pressure and vapour saturation as primary variables.

Because the choice of primary variables depends on the phase conditions, when the fluid in a cell changes phase, the primary variables must be changed.

Time evolution
**************

The discretised conservation equations :eq:`discretised_conservation` are of the form:

.. math::
   :label: RLeqn

   \frac{d}{dt} \mathit{L}(t, \vec{Y}) = \mathit{R}(t, \vec{Y})

where :math:`t` is time and :math:`\vec{Y}` is the vector of primary variables for all cells in the simulation mesh (of total length :math:`N(C+1)`). Here :math:`L` represents the cell-averaged mass and energy balances, as a function of time and the primary thermodynamic variables. Similarly, :math:`R` represents inflows into the cells (per unit volume) from flows through the cell faces, together with sources and sinks within the cell.

Solving the set of ordinary differential equations :eq:`RLeqn` with respect to time, we can compute the time evolution of :math:`\vec{Y}`, the thermodynamic state of the entire discretised simulation domain.

For solving the conservation equations, :math:`L` and :math:`R` are complex, non-linear functions of the primary variables :math:`\vec{Y}`. Hence equation :eq:`RLeqn` must be solved numerically, computing the solution :math:`\vec{Y}` at discrete times.

Waiwera contains a module for the numerical solution of ordinary differential equations of the form :eq:`RLeqn`, using different numerical methods. The simplest of these is the 'backwards Euler' method, which discretises equation :eq:`RLeqn` as follows:

.. math::
   :label: beuler

   \frac{L(t^{n+1}, \vec{Y}^{n+1}) - L(t^n, \vec{Y}^n)}{\Delta t} \approx R(t^{n+1}, \vec{Y}^{n+1})

where :math:`t^n` is the :math:`n^{th}` discretised time, and :math:`\Delta t` is the time step size, so that :math:`t^{n+1} = t^n + \Delta t`. For the backwards Euler method, at each time step we must solve equation :eq:`beuler` for the unknown new solution :math:`\vec{Y}^{n+1}`.

Solution of equations at each time step
***************************************

Regardless of the time stepping method used, the discretised equations to be solved at each time step (e.g. :eq:`beuler`) are non-linear. If we write them in a generic form:

.. math::
   :label: fx0

   f(\vec{y}) = \vec{0}

then at each time step we must solve this for the solution :math:`\vec{y} = \vec{Y}^{n+1}`. Because of the non-linearity, it must be solved numerically using a non-linear solution technique such as Newton's method. This is an iterative method which starts from an initial estimate of the solution (here taken as :math:`\vec{y} = \vec{Y}^n`) and adjusts the provisional solution :math:`\vec{y}` at each iteration until equation :eq:`fx0` is satisfied, to within a pre-specified tolerance.

At each iteration, Newton's method adds an update :math:`\Delta \vec{y}` to the provisional solution :math:`\vec{y}` according to:

.. math::
   :label: newton

   \matrix{J} \Delta \vec{y} = -f(\vec{y})

where :math:`\matrix{J}` is the Jacobian matrix of the function :math:`f`, i.e. the matrix of partial derivatives of :math:`f` with respect to :math:`\vec{y}`.

At each iteration, the Newton update equation :eq:`newton` represents a large, sparse system of linear equations to be solved numerically. "Krylov subspace" iterative methods (e.g. conjugate gradient methods) are appropriate for solving such systems. For typical simulations of large problems, most of the computation time is spent in the solution of the linear equations.

Waiwera uses the "SNES" non-linear solver provided by the `PETSc <https://www.mcs.anl.gov/petsc/>`_ library to solve equation :eq:`fx0` at each time step. For problems in which the Jacobian matrix :math:`\matrix{J}` is difficult to calculate, the SNES solver offers an option to calculate it automatically using finite differencing. In this case the Jacobian partial derivatives are evaluated approximately by adding small increments onto the primary variable vector :math:`\vec{y}` and re-evaluating the function :math:`f`. Waiwera makes use of this approach to calculate the Jacobian matrix.

The SNES solver in turn makes use of the "KSP" suite of linear solvers, also provided by PETSc, to solve the linear system :eq:`newton` at each Newton iteration.
