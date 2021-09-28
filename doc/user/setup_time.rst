.. index:: time, time; stepping, time step
.. _time_stepping:

*************
Time stepping
*************

Waiwera solves for the :ref:`time_evolution` of the system being simulated using a numerical time-stepping method. Time is discretised, proceeding from the simulation start time to its end time in steps.

At each time step, the time-stepping method results in a system of non-linear equations to be solved, using the `PETSc <https://www.mcs.anl.gov/petsc/>`_ "SNES" iterative non-linear solver (see :ref:`nonlinear_equations`). This in turn leads to a sequence of systems of linear equations to be solved at each non-linear solver iteration, using the PETSc "KSP" linear solvers.

All time-related parameters are specified in the Waiwera JSON input file via the **"time"** value, which is an object.

.. note::
   **JSON object**: time parameters

   **JSON path**: time

   +------------+------------+------------+--------------------------+
   |**name**    |**type**    |**default** |**value**                 |
   +------------+------------+------------+--------------------------+
   |"start"     |number      |0           |simulation start time (s) |
   |            |            |            |                          |
   |            |            |            |                          |
   +------------+------------+------------+--------------------------+
   |"stop"      |number |    |``null``    |simulation stop time (s)  |
   |            |``null``    |            |                          |
   +------------+------------+------------+--------------------------+
   |"step"      |object      |see below   |time stepping parameters  |
   |            |            |            |                          |
   +------------+------------+------------+--------------------------+

.. index:: time; start, time; stop

The **"start"** and **"stop"** values in the "time" object are numbers specifying the start and stop times for the simulation. The "stop" value can also be given the ``null`` value (or simply not specified), in which case the simulation will not stop at a prescribed time -- other criteria will determine when the simulation stops (e.g. how many time steps have been taken).

For example:

.. code-block:: json

   {"time": {"start": 86400, "stop": 172800}}

specifies a simulation starting at time 86400 s and ending at 172800 s.

All parameters related to time-stepping are specified via the **"step"** value in the "time" object. The "step" value is itself an object.

.. note::
   **JSON object**: time stepping parameters

   **JSON path**: time.step

   +-----------------+------------+------------------+-----------------------+
   |**name**         |**type**    |**default**       |**value**              |
   +-----------------+------------+------------------+-----------------------+
   |"method"         |string      |"beuler"          |time stepping method   |
   |                 |            |                  |                       |
   +-----------------+------------+------------------+-----------------------+
   |"size"           |number |    |0.1 s             |time step sizes (s)    |
   |                 |array       |                  |                       |
   +-----------------+------------+------------------+-----------------------+
   |"adapt"          |object      |see below         |adaptive time-stepping |
   |                 |            |                  |parameters             |
   +-----------------+------------+------------------+-----------------------+
   |"maximum"        |object      |{"size": ``null``,|maximum time step size |
   |                 |            |"number": 100,    |and number, and number |
   |                 |            |"tries": 10}      |of tries per step      |
   |                 |            |                  |                       |
   +-----------------+------------+------------------+-----------------------+
   |"stop"           |object      |{"step":          |minimum and maximum    |
   |                 |            |{"minimum":       |time step sizes causing|
   |                 |            |``null``,         |simulation to stop     |
   |                 |            |"maximum":        |                       |
   |                 |            |``null``}}        |                       |
   +-----------------+------------+------------------+-----------------------+
   |"solver"         |object      |see below         |non-linear and linear  |
   |                 |            |                  |solver parameters      |
   +-----------------+------------+------------------+-----------------------+

.. index:: time step; methods, numerical methods; time evolution
.. _time_stepping_methods:

Time stepping methods
=====================

The time evolution of the solution vector of primary thermodynamic variables :math:`\mathbf{Y}` is found by solving the discretised mass and energy conservation equations, which can be written (see :ref:`time_evolution`) in the general form:

.. math::
   :label: RLeqn2

   \frac{d}{dt} \mathbf{L}(t, \mathbf{Y}) = \mathbf{R}(t, \mathbf{Y})

where :math:`t` is time and  the left- and right-hand side functions :math:`\mathbf{L}` and :math:`\mathbf{R}` represent the cell mass and energy balances, and the inflows into each cell from fluxes and source terms respectively.

Waiwera contains a module for solving ordinary differential equations of this form, using different numerical time-stepping methods. These solve for the solution :math:`\mathbf{Y}^n` at a sequence of discretised times :math:`t^n`. At the :math:`n^{th}` time step the new solution :math:`\mathbf{Y}^{n+1}` at time :math:`t^{n+1} = t^n + \Delta t^n` (where :math:`\Delta t^n` is the current time step size) is computed from :math:`\mathbf{Y}^n` (and possibly other previous solutions) by solving a set of non-linear equations, which depend on the time-stepping method being used.

At present, only a few relatively simple time-stepping methods are included. Particularly for geothermal flow models, the complex and highly non-linear nature of the equations being solved mean that the time-stepping methods need to be very stable.

The time-stepping method is specified via the **"method"** value in the "time.step" object. This is a simple string value which defaults to "beuler", selecting the :ref:`backwards_euler` method.

For example:

.. code-block:: json

   {"time": {"step": {"method": "beuler"}}}

specifies the backwards Euler time stepping method.

.. index:: numerical methods; backwards Euler
.. _backwards_euler:

Backwards Euler
---------------

The simplest time-stepping method included is the "backwards Euler" method, selected in the Waiwera JSON input file by setting the "method" value in the "time.step" object to **"beuler"**. This fully-implicit method is only first-order accurate, but is highly stable. At least for the present, it is recommended for most applications.

The backwards Euler method discretises equation :eq:`RLeqn2` as follows:

.. math::
   :label: beuler2

   \frac{1}{\Delta t^n} \big(\mathbf{L}^{n+1} - \mathbf{L}^n \big) = \mathbf{R}^{n+1}

where :math:`\mathbf{L}^n = \mathbf{L}(t^n, \mathbf{Y}^n)` and :math:`\mathbf{R}^n = \mathbf{R}(t^n, \mathbf{Y}^n)`.

.. index:: numerical methods; BDF2

BDF2
----

BDF2 (selected in the Waiwera JSON input file by setting the "method" value in the "time.step" object to **"bdf2"**) is one of a series of "backward differentiation formula" methods (also known as "Gear algorithms") designed for solving stiff differential equations. It is an implicit "linear multistep" method: the new solution :math:`\mathbf{Y}^{n+1}` is found not only from the previous solution :math:`\mathbf{Y}^n` but also from :math:`\mathbf{Y}^{n-1}`. BDF2 is second-order accurate but has a slightly smaller stability region than the backwards Euler method (which can be considered the lowest-order member of the family of BDF methods).

The variable-stepsize BDF2 method discretises equation :eq:`RLeqn2` as follows:

.. math::
   :label: BDF2

   \frac{1}{\Delta t^n} \Big(\frac{1 + 2r}{1 + r} \mathbf{L}^{n+1} - (1 + r) \mathbf{L}^n + \frac{r^2}{1 + r} \mathbf{L}^{n-1} \Big) = \mathbf{R}^{n+1})

where :math:`r = \Delta t^n / \Delta t^{n-1}` is the stepsize ratio.

Equation :eq:`BDF2` cannot be used for the first time step in the simulation (:math:`n = 1`), as there is no previous solution :math:`\mathbf{Y}^0`. Hence, backwards Euler is used as a startup method on the first time step. 

.. index:: time step; sizes
.. _specifying_time_step_sizes:

Specifying time step sizes
==========================

Time step sizes can be specified using the **"size"** value in the "time.step" object. This can be either a single number, to set a fixed time step size, or an array of numbers. For example:

.. code-block:: json

   {"time": {"step": {"size": 3600}}}

sets a fixed time step size of 3600 s. In the following example, an array of increasing time step sizes is specified:

.. code-block:: json

   {"time": {"step": {"size": [1e3, 2e3, 3e3, 4e3]}}}

If all specified time step sizes have been used, but the simulation has not yet finished, then the simulation will continue using the last specified time step size (unless :ref:`adaptive_time_stepping` has been selected, in which case the step size will be automatically adapted once the specified step sizes have been performed).

.. index:: time step; adaptive
.. _adaptive_time_stepping:

Adaptive time stepping
======================

When adaptive time stepping is used, the time stepper will automatically adjust the time step size as the simulation progresses. Only the initial time step size need be specified, via the "time.step.size" value (see :ref:`specifying_time_step_sizes`). If an array of step sizes is specified, these will be performed first before adaptive time stepping commences.

The time step adaption algorithm uses the concept of a "monitor value" :math:`\eta`, which is essentially a non-dimensional measure of how much the solution has changed over the course of the last time step. If :math:`\eta` is within a specified range :math:`[\eta_{min}, \eta_{max}]` then the time step size is left unchanged; otherwise it is decreased or increased. Specifically, the new time step size :math:`\Delta t^{n+1}` is given by:

.. math::

   \Delta t^{n+1} = \begin{cases}
     \alpha \Delta t^n & \eta < \eta_{min} \\
     \Delta t^n & \eta_{min} \le \eta \le \eta_{max} \\
     \beta \Delta t^n & \eta > \eta_{max}
   \end{cases}

where :math:`\alpha > 1` is an amplification factor for increasing the time step size, and :math:`\beta < 1` is a reduction factor for reducing it. These are specified in the JSON input file via the **"amplification"** and **"reduction"** values in the **"time.step.adapt"** object.

When the time step size is increased, it is limited by the maximum time step size (if any) specified using the **"time.step.maximum.size"** value. If this maximum is exceeded, the simulation will continue with the step size reduced to the specified maximum.

It is also possible to specify minimum and maximum time step sizes which, if exceeded, will cause the simulation to stop. These limits are specified via the **"time.step.stop.size.minimum"** and **"time.step.stop.size.maximum"** values. Either of these can take number values or ``null`` (the default), in which case no limit is enforced. If either limit is hit, a final time step is carried out at the limiting time step size before the simulation is stopped.

.. note::
   **JSON object**: time step stop sizes

   **JSON path**: time.step.stop.size

   +----------------+------------------+----------------+-----------------------+
   |**name**        |**type**          |**default**     |**value**              |
   +----------------+------------------+----------------+-----------------------+
   |"minimum"       |number | ``null`` |``null``        |Minimum time step size |
   |                |                  |                |causing simulation to  |
   |                |                  |                |stop                   |
   +----------------+------------------+----------------+-----------------------+
   |"maximum"       |number | ``null`` |``null``        |Maximum time step size |
   |                |                  |                |causing simulation to  |
   |                |                  |                |stop                   |
   +----------------+------------------+----------------+-----------------------+

The "time.step.adapt" object has a Boolean **"on"** value, which determines whether adaptive time stepping is to be used. Note, however, that it can still be useful to specify at least some of the other adaptor parameters even if the adaptor is switched off. This is because the adaptor is also used to handle :ref:`time_step_reductions`. If these parameters are not specified, default values will be used.

.. note::
   **JSON object**: time step adaptor

   **JSON path**: time.step.adapt

   +----------------+------------+------------+-------------------------+
   |**name**        |**type**    |**default** |**value**                |
   +----------------+------------+------------+-------------------------+
   |"on"            |boolean     |``false``   |whether adaptor is used  |
   +----------------+------------+------------+-------------------------+
   |"method"        |string      |"iteration" |method used for adapting |
   |                |            |            |time step size           |
   +----------------+------------+------------+-------------------------+
   |"minimum"       |number      |5           |minimum monitor value    |
   |                |            |            |:math:`\eta_{min}`       |
   |                |            |            |                         |
   +----------------+------------+------------+-------------------------+
   |"maximum"       |number      |8           |maximum monitor value    |
   |                |            |            |:math:`\eta_{max}`       |
   |                |            |            |                         |
   +----------------+------------+------------+-------------------------+
   |"amplification" |number      |2           |factor :math:`\alpha` for|
   |                |            |            |increasing time step size|
   |                |            |            |                         |
   +----------------+------------+------------+-------------------------+
   |"reduction"     |number      |0.2         |factor :math:`\beta` for |
   |                |            |            |reducing time step size  |
   |                |            |            |                         |
   +----------------+------------+------------+-------------------------+

Two different time step adaption methods are available, selected using the **"method"** string value in the "time.step.adapt" JSON object. They differ only in the way the monitor value :math:`\eta` is defined.

Non-linear iteration count method
---------------------------------

This method, selected by setting the "method" value to "iteration", uses the number of non-linear solver iterations in the latest time step as the monitor value :math:`\eta`. Because the non-linear solver starts from the previous solution :math:`\mathbf{Y}^n` as its initial estimate of the new solution :math:`\mathbf{Y}^{n+1}`, in general the difference between these values may be expected to be correlated with the number of iterations.

Relative change method
----------------------

This method, selected by setting the "method" value to "change", defines the monitor value as :math:`\eta = \|\mathbf{\Delta L}\|_{\infty}` where:

.. math::

   \Delta L_i = \frac{L_i^{n+1} - L_i^n} {\max{(|L_i^n|, \epsilon)}}

:math:`\epsilon` is a small constant (:math:`\epsilon = 10^{-3}` is used here), preventing division by zero. In this approach, the value of :math:`\eta` is essentially a measure of the relative change in the mass and energy balances in the cells.

Example
-------

In the example below, an initial time step size of 3600 s (1 hour) is used, after which time step sizes are chosen adaptively using the non-linear iteration count method. The time step size will be increased if the non-linear solver converges in fewer than 4 iterations, and decreased if it takes more than 8 iterations. The time step size is not allowed to exceed 86400 s (1 day), and the simulation should stop at time 2592000 s (30 days).

.. code-block:: json

   {"time": {"step": {"size": 3600,
                      "adapt": {"on": true,
                               "method": "iteration",
                               "minimum": 4, "maximum": 8},
                      "maximum": {"size": 86400}
                      },
             "stop": 2592000}}

.. index:: time step; reductions
.. _time_step_reductions:

Time step reductions
====================

If a time step cannot be completed with its original size, it is re-tried with a reduced step size. This may occur if, for example, the non-linear solver aborts or does not converge within the specified maximum allowed number of iterations.

The non-linear solver may abort if the linear solver does not converge, or if primary thermodynamic variables go outside the range of validity of the thermodynamic formulation (see :ref:`water_thermodynamics`). Slow convergence of the non-linear solver may be caused by a variety of factors, including large numbers of phase transitions within the time step.

The process of re-trying the time step with a reduced time step size may be carried out multiple times until the time step is successfully completed. There is, however, a limit on the number of allowable tries, specified by the **"time.step.maximum.tries"** value (which defaults to 10). If this limit is exceeded, the simulation will stop.

If specified-size time steps are being used (see :ref:`specifying_time_step_sizes`), the process of reducing the time step size is carried out by temporarily turning on the time step size adaptor (see :ref:`adaptive_time_stepping`). After a successful reduced-size time step has been completed, the adaptor will then try to increase the time step size again if possible. Once the original specified time step size has been attained the adaptor will be switched off, and the time stepper will resume using the specified time step sizes.

Time stepping termination
=========================

The time stepper can terminate in a number of ways:

* if the time reaches or exceeds the stop time, specified by the **"time.stop"** value in the JSON input (if the time exceeds the stop time, the time step size will be reduced to hit the stop time exactly)
* if the number of time steps reaches the limit specified by the **"time.step.maximum.number"** value
* if the time step size falls below a lower limit specified by **"time.step.stop.size.minimum"** or exceeds **"time.step.stop.size.maximum"**
* if a time step fails to complete, and the time step size reduction process is repeated more than the maximum allowable number of tries specified by the **"time.step.maximum.tries"** value

.. index:: time; steady state

.. _steady_state:

Steady-state simulations
========================

It is often necessary to solve for the steady-state behaviour of a system, for example to estimate the "natural state" of a geothermal reservoir before production. In this case the discretised conservation equations reduce to:

.. math::

   \mathbf{R}(\mathbf{Y}) = \mathbf{0}

Direct solution
---------------

In principle it is possible to solve these equations directly for the steady-state solution :math:`\mathbf{Y}`, without time stepping, in a single non-linear solution process. Waiwera does offer this as an option, by setting the "time.step.method" value to "directss", but this approach does not usually work well. The discretised conservation equations even in their transient form are numerically difficult to solve. Eliminating the time derivative term to give the steady-state form of the equations increases the difficulty further.

However, the main problem with this approach lies in the fact that the non-linear solver still needs a starting estimate of the solution, and in most cases will not converge unless this starting estimate is close to the steady-state solution. Hence, the direct solution approach is not usually recommended.

Using adaptive time-stepping
----------------------------

The usual approach to finding steady-state solutions is to solve the transient conservation equations using :ref:`adaptive_time_stepping`, using the :ref:`backwards_euler` time-stepping method, without limiting the time step size, and letting the time stepper run until a very large time step size has been achieved. As the time step size :math:`\Delta t^n` increases, it gradually reduces the left-hand side time derivative term in equation :eq:`beuler2`, until at very large time step sizes it is effectively zero.

This approach has the advantage that it usually still converges to the steady-state solution, even if it is started from an initial condition that is not close to the solution. The time-stepping process can be seen as effectively an outer iteration procedure that drives the problem from being transient to steady-state.

What constitutes a "very large" time step size is somewhat problem-dependent, and is determined mainly by numerical considerations rather than any physical time-scales of the transient problem. The main criterion is that the final time step size needs to be large enough to make the left-hand side derivative terms in equation :eq:`beuler2` negligibly small. For typical geothermal reservoir models a time step size of at least 10\ :sup:`15` s is usually needed for a reliable steady-state solution.

As the time step size increases and the left-hand side time derivative term in equation :eq:`beuler2` decreases in magnitude, the linear equations to be solved at each non-linear solver iteration generally become progressively more ill-conditioned. In the later stages of a steady-state simulation it is common for the linear solver to take more iterations to solve, or to fail. To obtain a properly converged steady state solution it may be necessary to experiment with different linear solvers and preconditioners (see :ref:`linear_equation_solution`).

Setting up a steady-state simulation using this approach can be done by specifying a large maximum stopping time step size (via "time.step.stop.size.maximum"), e.g. 10\ :sup:`15` s, and no stop time ("time.stop" = ``null``, the default). A limit on the total number of time steps ("time.step.maximum.number") is usually set, so that the simulation still stops even if a large time step size (and hence a true steady state) is never attained. After the simulation has finished, it is important to check that it has reached the specified stopping time step size rather than the maximum number of time steps.

For example:

.. code-block:: json

   {"time": {"step": {"size": 1e6,
                      "adapt": {"on": true,
                               "method": "iteration",
                               "minimum": 5, "maximum": 8},
                      "maximum": {"number": 500},
                      "method": "beuler",
                      "stop": {"size": {"maximum": 1e15}}
                      },
             "stop": null}}

sets up a steady-state simulation using adaptive time-stepping, with a starting time step size of 10\ :sup:`6` s and a large stopping time step size of 10\ :sup:`15` s, which must be attained within 500 time steps.

.. index:: numerical methods; non-linear equations, solver; non-linear, PETSc; SNES
.. _nonlinear_solution:

Solution of non-linear equations
================================

At each time step the `PETSc <https://www.mcs.anl.gov/petsc/>`_ "SNES" non-linear solver (with Newton-Raphson iteration by default) is used to solve the discretised mass and energy conservation equations, e.g. equation :eq:`beuler2` for the backwards Euler time-stepping method. The conservation equations are re-written as a function, known as the **residual** function, so that finding the root of this function corresponds to solving the original equation. For example, for the backwards Euler time-stepping method, the residual function :math:`\mathbf{f}` is:

.. math::

   \mathbf{f} = \mathbf{L}^{n+1} - \mathbf{L}^n - \Delta t \: \mathbf{R}^{n+1}

Convergence in the residual
---------------------------

The non-linear solution process is considered converged when all the elements of the residual :math:`\mathbf{f}` are sufficiently small. Note, however, that the left- and right-hand side vectors :math:`\mathbf{L}` and :math:`\mathbf{R}`, and hence also :math:`\mathbf{f}`, usually contain values of differing magnitudes, depending on whether they arise from mass or energy components. Hence, for the purpose of checking convergence, it is necessary to non-dimensionalise the residual :math:`\mathbf{f}` so that its elements are all of comparable sizes. The non-dimensionalised residual :math:`\mathbf{f'}` is defined as:

.. math::

   f'_i = \frac{f_i}{\max{(|L^n_i|, \epsilon_a)}}

and the non-linear solution process is then considered converged when :math:`\|\mathbf{f}'\|_{\infty} < \epsilon_r`. Here :math:`\epsilon_a` and :math:`\epsilon_r` are specified tolerances, set in the Waiwera JSON input file via the **"tolerance.function.absolute"** and **"tolerance.function.relative"** values respectively in the **"time.step.solver.nonlinear"** object.

.. note::
   **JSON object**: non-linear solver parameters

   **JSON path**: time.step.solver.nonlinear

   +------------+------------+------------------+--------------------+
   |**name**    |**type**    |**default**       |**value**           |
   +------------+------------+------------------+--------------------+
   |"maximum"   |object      |{"iterations": 8} |maximum number of   |
   |            |            |                  |iterations          |
   |            |            |                  |                    |
   +------------+------------+------------------+--------------------+
   |"minimum"   |object      |{"iterations": 0} |minimum number of   |
   |            |            |                  |iterations          |
   +------------+------------+------------------+--------------------+
   |"tolerance" |object      |{"function":      |relative and        |
   |            |            |{"relative": 10\  |absolute tolerances |
   |            |            |:sup:`-5`,        |on function value   |
   |            |            |"absolute": 1},   |and solution update |
   |            |            |"update":         |                    |
   |            |            |{"relative": 10\  |                    |
   |            |            |:sup:`-10`,       |                    |
   |            |            |"absolute": 1}}   |                    |
   +------------+------------+------------------+--------------------+

Convergence in the update
-------------------------

If the solution process has not converged in terms of the residual function value, then a second test is carried out to check the magnitude of the update :math:`\Delta \mathbf{Y}` to the solution :math:`\mathbf{Y}` during the latest iteration. If this update is sufficiently small then the process may also be considered converged, even if the residual function value check has not passed. Again, this solution update needs to be non-dimensionalised as different elements of the solution may be of different magnitudes. The non-dimensionalised solution update :math:`\Delta \mathbf{Y}'` is defined as:

.. math::

   \Delta Y'_i = \frac{\Delta Y_i}{\max{(|Y^n_i|, \delta_a)}}

and the solution process is considered converged in the update if :math:`\|\Delta \mathbf{Y}'\|_{\infty} < \delta_r`, where :math:`\delta_a` and :math:`\delta_r` are specified tolerances, set in the Waiwera JSON input file via the **"tolerance.update.absolute"** and **"tolerance.update.relative"** values respectively in the **"time.step.solver.nonlinear"** object.

Iteration limits
----------------

Limits on the number of non-linear solver iterations can be set via the **"maximum"** and **"minimum"** values in the "time.step.solver.nonlinear" JSON object. These values are both objects containing an **"iterations"** integer value, so for example the maximum number of non-linear solver iterations is set using "time.step.solver.nonlinear.maximum.iterations".

The "minimum.iterations" value defaults to zero, so that the non-linear solution process is allowed to converge without iterating, if it happens that either of the convergence checks are satisfied with the initial estimate of the solution. Under some conditions, it is useful to make sure the non-linear solver always takes at least one iteration. This can be done by setting the "minimum.iterations" value to 1.

.. index:: non-linear equations, Jacobian
.. _jacobian_matrix:

Computing the Jacobian matrix
-----------------------------

The nonlinear solution process requires the Jacobian matrix (i.e. the matrix of partial derivatives of the residual function :math:`\mathbf{f}` with respect to the primary variables) to be computed at each iteration (see :ref:`nonlinear_equations`). Waiwera makes use of PETSc's ability to compute the Jacobian automatically using finite differencing.

Specifically, the Jacobian matrix :math:`\mathbf{J}` is computed using the formula:

.. math::

   J_{i,j} = \frac{f_i(\mathbf{Y} + h_j \mathbf{e}_j) - f_i(\mathbf{Y})}{h_j}

where :math:`f_i` is the :math:`i^{th}` component of the residual function :math:`\mathbf{f}`, :math:`\mathbf{e}_j` is the :math:`j^{th}` standard basis vector and :math:`h_j` is the finite difference step size for the :math:`j^{th}` variable :math:`Y_j`.

The step size :math:`h_j` is calculated from:

.. math::

    h_j =
   \begin{cases}
   \epsilon Y_j & |Y_j| > \delta \\
   sgn(Y_j) \epsilon \delta & |Y_j| \leq \delta
   \end{cases}

where :math:`sgn` is the sign function,  :math:`\epsilon` is a fixed differencing increment parameter and :math:`\delta` is a differencing tolerance parameter (to prevent zero-division errors or loss of significance when the variable :math:`Y_j` is small). These two parameters may be set in the Waiwera JSON input file via the **"time.step.solver.nonlinear.jacobian.differencing"** object:

.. note::
   **JSON object**: Jacobian finite differencing parameters

   **JSON path**: time.step.solver.nonlinear.jacobian.differencing

   +------------+------------+------------------+--------------------------+
   |**name**    |**type**    |**default**       |**value**                 |
   +------------+------------+------------------+--------------------------+
   |"increment" |number      |10\ :sup:`-8`     |differencing increment    |
   |            |            |                  |:math:`\epsilon`          |
   |            |            |                  |                          |
   +------------+------------+------------------+--------------------------+
   |"tolerance" |number      |10\ :sup:`-2`     |differencing tolerance    |
   |            |            |                  |:math:`\delta`            |
   +------------+------------+------------------+--------------------------+

(These values can also be set using the :ref:`petsc_command_line_parameters` ``-mat_fd_coloring_err`` and ``-mat_fd_coloring_umin`` respectively.)

In the non-linear solution process, non-dimensionalised versions of the primary variables :math:`\mathbf{Y}` are used, rather than their raw values (see :ref:`primary_variable_parameters`). This is partly because raw variables representing different physical quantities (e.g. pressures and temperatures) may have very different magnitudes, which would make it difficult to choose differencing parameters appropriate for all of them.

In most cases the default differencing parameters should give satisfactory results. However, for some problems, adjusting them can improve performance. If, for example, the non-linear solver is often converging slowly, or not converging, it may be that the Jacobian is not being computed with sufficient accuracy. Trying a smaller differencing increment may help. Alternatively, if non-linear solver convergence seems to be held up by variables with near-zero values (e.g. partial pressures of non-condensible gases) it may be worth experimenting with the differencing tolerance.

Example
-------

In the following example, a steady-state simulation is specified with the maximum allowed number of non-linear solver iteration increased to 10, and the relative function tolerance reduced to 10\ :sup:`-6`. Additionally, the Jacobian finite differencing increment is set to :math:`10^{-9}`:

.. code-block:: json

   {"time": {"step": {"size": 1e6,
                      "adapt": {"on": true,
                               "method": "iteration",
                               "minimum": 5, "maximum": 8},
                      "maximum": {"number": 500},
                      "stop": {"size": {"maximum": 1e15}},
                      "method": "beuler",
                      "solver": {
                         "nonlinear": {"maximum": {"iterations": 10},
                                       "tolerance": {"function": {"relative": 1e-6}},
                                       "jacobian": {"differencing": {"increment": 1e-9}}}
                      }},
             "stop": null}}

.. index:: numerical methods; linear equations, solver; linear
.. _linear_equation_solution:

Solution of linear equations
============================

At each iteration of the non-linear solver (see :ref:`nonlinear_equations`), a large, sparse system of linear equations must be solved to find the latest Newton-Raphson update :math:`\Delta \mathbf{Y}` to the solution vector :math:`\mathbf{Y}`:

.. math::

   \mathbf{J} \Delta \mathbf{Y} = -\mathbf{f}

where :math:`\mathbf{J}` is the Jacobian matrix of the residual function :math:`\mathbf{f}`.

This system of linear equations is solved using the `PETSc <https://www.mcs.anl.gov/petsc/>`_ "KSP" suite of parallelised linear equation solvers. Linear solver parameters can be specified via the **"time.step.solver.linear"** value in the Waiwera JSON input file. This value is an object.

.. note::
   **JSON object**: linear solver parameters

   **JSON path**: time.step.solver.linear

   +-----------------+------------+---------------------+----------------+
   |**name**         |**type**    |**default**          |**value**       |
   +-----------------+------------+---------------------+----------------+
   |"type"           |string      |"bcgs"               |linear solver   |
   |                 |            |                     |type            |
   +-----------------+------------+---------------------+----------------+
   |"options"        |object      |{}                   |linear solver   |
   |                 |            |                     |options         |
   |                 |            |                     |                |
   +-----------------+------------+---------------------+----------------+
   |"tolerance"      |object      |{}                   |linear solver   |
   |                 |            |                     |tolerance       |
   +-----------------+------------+---------------------+----------------+
   |"maximum"        |object      |{}                   |maximum         |
   |                 |            |                     |iterations      |
   +-----------------+------------+---------------------+----------------+
   |"preconditioner" |object      |{"type": "asm",      |preconditioner  |
   |                 |            |"sub":               |options         |
   |                 |            |{"preconditioner":   |                |
   |                 |            |{"type": "ilu",      |                |
   |                 |            |"factor": {"levels": |                |
   |                 |            |0}}}}                |                |
   +-----------------+------------+---------------------+----------------+

Linear solver type
------------------

PETSc offers a range of different `KSP linear solver types <http://www.mcs.anl.gov/petsc/petsc-dev/docs/linearsolvertable.html>`_. For Waiwera, the most appropriate linear solvers are generally the "Krylov subspace" methods. The linear solver type can be specified in the Waiwera JSON input file via the **"type"** string value in the "time.step.solver.linear" object. The linear solver types that may be selected in this way are:

+------------+---------------+-----------------------+
|**name**    |**PETSc name** |**description**        |
|            |               |                       |
+------------+---------------+-----------------------+
|"gmres"     |KSPGMRES       |generalised minimum    |
|            |               |residual               |
+------------+---------------+-----------------------+
|"lgmres"    |KSPLGMRES      |augmented GMRES        |
+------------+---------------+-----------------------+
|"bcgs"      |KSPBCGS        |Bi-CGStab (stabilised  |
|            |               |bi-conjugate gradient) |
+------------+---------------+-----------------------+
|"bcgsl"     |KSPBCGSL       |Bi-CGStab(L)           |
+------------+---------------+-----------------------+

These represent the most commonly useful linear solver types for the linear equation systems solved by Waiwera. (Note that other PETSc linear solver types may be selected at run-time using :ref:`petsc_command_line_parameters`.) The GMRES and Bi-CGStab solvers generally perform adequately for many problems. For very ill-conditioned systems (e.g. near the end of steady-state simulations) the Bi-CGStab(L) solver may give better performance. If using the GMRES solver, increasing the restart parameter (see :ref:`solver_options`) may also help. If linear solver failures persist, it may be necessary to experiment with different :ref:`preconditioners`.

For example:

.. code-block:: json

   {"time": {"step": {"solver": {"linear": {"type": "gmres"}}}}}

selects the GMRES linear solver type.

Convergence parameters
----------------------

The above linear solver types are all iterative methods, so parameters may be set to control convergence criteria.

The convergence tolerance may be specified via the **"tolerance"** value in the "time.step.solver.linear" object. This value is itself an object, containing a **"relative"** number value for specifying the relative convergence tolerance. If not specified, then the PETSc default tolerance is used.

.. note::
   **JSON object**: linear solver tolerance

   **JSON path**: time.step.solver.linear.tolerance

   +------------+------------+--------------+---------------------+
   |**name**    |**type**    |**default**   |**value**            |
   +------------+------------+--------------+---------------------+
   |"relative"  |number      |PETSc default |relative convergence |
   |            |            |              |tolerance            |
   |            |            |              |                     |
   +------------+------------+--------------+---------------------+

The maximum allowed number of linear solver iterations can be specified using the **"maximum"** value in the "time.step.solver.linear" object, which again is itself an object, this time containing an **"iterations"** integer value for specifying the iteration limit. If not specified, then the PETSc default is used.

.. note::
   **JSON object**: linear solver iteration limit

   **JSON path**: time.step.solver.linear.maximum

   +-------------+------------+--------------+------------------+
   |**name**     |**type**    |**default**   |**value**         |
   +-------------+------------+--------------+------------------+
   |"iterations" |integer     |PETSc default |iteration limit   |
   |             |            |              |                  |
   +-------------+------------+--------------+------------------+

For example:

.. code-block:: json

   {"time": {"step": {"solver": {"linear": {"type": "gmres",
                                            "tolerance": {"relative": 1e-12},
                                            "maximum": {"iterations": 2000}
                                            }}}}}

selects a GMRES linear solver with relative tolerance 10\ :sup:`-12` and and iteration limit of 2000.

.. _solver_options:

Solver options
--------------

Some linear solvers may have options specific to the solver type, which can be specified via the **"options"** value in the "time.step.solver.linear" object.

.. note::
   **JSON object**: linear solver options

   **JSON path**: time.step.solver.options

   +------------+------------+------------+--------------+
   |**name**    |**type**    |**default** |**value**     |
   +------------+------------+------------+--------------+
   |"gmres"     |object      |see below   |GMRES options |
   |            |            |            |              |
   +------------+------------+------------+--------------+

Currently there is only one such linear solver (GMRES) with options available in this way. The GMRES solver in PETSc offers a "restarted GMRES" option, and the "options.gmres" object has a **"restart"** integer value to specify the number of Krylov search directions to orthogonalise against. For some Waiwera simulations the restarted GMRES linear solver performs well on difficult problems, particularly if the restart parameter is increased.

.. note::
   **JSON object**: GMRES linear solver options

   **JSON path**: time.step.solver.options.gmres

   +------------+------------+-------------------+--------------------+
   |**name**    |**type**    |**default**        |**value**           |
   +------------+------------+-------------------+--------------------+
   |"restart"   |integer     |PETSc default (30) |number of Krylov    |
   |            |            |                   |search directions   |
   +------------+------------+-------------------+--------------------+

For example:

.. code-block:: json

   {"time": {"step": {"solver": {"linear": {"type": "gmres",
                                            "options": {"gmres": {"restart": 200}}
                                            }}}}}

selects the restarted GMRES linear solver with a restart parameter of 200.

.. _preconditioners:

Preconditioners
---------------

Preconditioners are used to improve the convergence rate of iterative linear equation solvers. A preconditioner transforms the problem so that the resulting matrix has a lower condition number, allowing the linear solver to converge more rapidly. This is especially important when the original system of linear equations to be solved is ill-conditioned, as is often the case for the equations solved by Waiwera.

Preconditioning parameters can be specified using the **"preconditioner"** value in the "time.step.solver.linear" object.

.. note::
   **JSON object**: linear solver preconditioner

   **JSON path**: time.step.solver.linear.preconditioner

   +------------+------------+------------+---------------------+
   |**name**    |**type**    |**default** |**value**            |
   +------------+------------+------------+---------------------+
   |"type"      |string      |"asm"       |preconditioner type  |
   +------------+------------+------------+---------------------+
   |"sub"       |object      |see below   |sub-preconditioner   |
   |            |            |            |options              |
   +------------+------------+------------+---------------------+

PETSc offers a range of different preconditioners. The **"type"** string value in the "preconditioner" object can be used to specify the preconditioner type. The preconditioner types that may be selected in this way are:

+------------+----------------+-------------------+-------------------+
|**name**    |**PETSc name**  |**description**    |**parallel**       |
+------------+----------------+-------------------+-------------------+
|"bjacobi"   |PCBJACOBI       |Block Jacobi       |yes                |
|            |                |                   |                   |
+------------+----------------+-------------------+-------------------+
|"asm"       |PCASM           |Additive Schwarz   |yes                |
|            |                |method             |                   |
+------------+----------------+-------------------+-------------------+
|"ilu"       |PCILU           |incomplete LU      |no                 |
|            |                |factorisation      |                   |
+------------+----------------+-------------------+-------------------+
|"lu"        |PCLU            |LU factorisation   |no                 |
+------------+----------------+-------------------+-------------------+
|"none"      |PCNONE          |no preconditioning |no                 |
+------------+----------------+-------------------+-------------------+

Of these, only the "bjacobi" and "asm" preconditioners are suitable for parallel simulations. Which one works better depends on the problem and can only be determined by experiment.

The other preconditioner types are really only included here for testing purposes (e.g. if the linear solver is failing, and it is necessary to determine if the problem lies in the linear solver itself or the preconditioner). The "ilu" preconditioner could also be used for serial simulations.

Sub-preconditioners
-------------------

In a parallel simulation, the matrix is effectively treated as a block matrix, with one block on each processor by default (so the blocks are determined by the :ref:`mesh_partitioning`). The parallel preconditioner operates at the block level, and each block has its own sub-preconditioner, which operates in serial.

By default, the PETSc implementations of the Block Jacobi and Additive Schwarz parallel preconditioners use ILU(0) sub-preconditioning on each block by default (i.e. incomplete LU factorisation with no fill-in). Other sub-preconditioner types are available, but in general the ILU sub-preconditioner works adequately and there is little reason to use anything else.

For very demanding problems it may be necessary, however, to increase the level of fill-in in the ILU sub-preconditioner. The level of fill-in may be specified via the **"sub.preconditioner"** value in the "time.step.solver.linear.preconditioner" object. This is itself an object, which contains a **"factor.levels"** value specifying the level of fill-in. (There is also a **"type"** string value which can be used for changing the sub-preconditioner type.)

.. note::
   **JSON object**: linear solver sub-preconditioner

   **JSON path**: time.step.solver.linear.preconditioner.sub.preconditioner

   +------------+------------+--------------+-------------------------+
   |**name**    |**type**    |**default**   |**value**                |
   +------------+------------+--------------+-------------------------+
   |"type"      |string      |"ilu"         |sub-preconditioner type  |
   +------------+------------+--------------+-------------------------+
   |"factor"    |object      |{"levels": 0} |level of fill-in for     |
   |            |            |              |"ilu" sub-preconditioner |
   +------------+------------+--------------+-------------------------+

For example:

.. code-block:: json

   {"time": {"step": {"solver": {"linear": {"type": "bcgs",
                                            "preconditioner": {"sub":
                                              {"preconditioner":
                                                {"factor": {"levels": 3}}}}
                                            }}}}}

specifies a Bi-CGStab linear solver. The default ILU sub-preconditioner is used, but with the level of fill-in increased to 3.

Solution of tracer equations
============================

When tracers are being simulated (see :ref:`setup_tracers`), the tracer mass fractions are updated at each time step by solving a single auxiliary system of linear equations (separate from the main flow solution process). This is carried out using the `PETSc <https://www.mcs.anl.gov/petsc/>`_ "KSP" suite of parallelised linear equation solvers, as is done at each iteration of the non-linear flow solution process (see :ref:`linear_equation_solution`).

However, the linear solver parameters for the auxiliary tracer solution can be specified independently of those for the main flow solution, via the **"time.step.solver.auxiliary"** value in the Waiwera JSON input file.

This object is identical in structure to the **"time.step.solver.linear"** value, but has some slightly different default values. Specifically, the auxiliary linear solver and preconditioner types have different defaults (see table below), more suitable for solving the discretised tracer conservation equations.

All other values within the **"time.step.solver.auxiliary"** object can be specified in the same way as their counterparts in the **"time.step.solver.linear"** object, and have the same defaults.

.. note::
   **JSON object**: auxiliary linear solver parameters

   **JSON path**: time.step.solver.auxiliary

   +-----------------+------------+---------------------+----------------+
   |**name**         |**type**    |**default**          |**value**       |
   +-----------------+------------+---------------------+----------------+
   |"type"           |string      |"gmres"              |linear solver   |
   |                 |            |                     |type            |
   +-----------------+------------+---------------------+----------------+
   |"options"        |object      |{}                   |linear solver   |
   |                 |            |                     |options         |
   |                 |            |                     |                |
   +-----------------+------------+---------------------+----------------+
   |"tolerance"      |object      |{}                   |linear solver   |
   |                 |            |                     |tolerance       |
   +-----------------+------------+---------------------+----------------+
   |"maximum"        |object      |{}                   |maximum         |
   |                 |            |                     |iterations      |
   +-----------------+------------+---------------------+----------------+
   |"preconditioner" |object      |{"type": "bjacobi",  |preconditioner  |
   |                 |            |"sub":               |options         |
   |                 |            |{"preconditioner":   |                |
   |                 |            |{"type": "ilu",      |                |
   |                 |            |"factor": {"levels": |                |
   |                 |            |0}}}}                |                |
   +-----------------+------------+---------------------+----------------+

