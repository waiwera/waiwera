title: Code structure

## Main program

The [[waiwera(program)]] program is the main driver program for the code. Its main job is to create and use two objects:

- the [[waiwera:simulation]] object, representing the particular simulation to be run
- the [[waiwera:timestepper]] object

First the timestepper is initialized and associated with the simulation. Then its [[timestepper_type:run]] method is called, which steps the simulation through time until it completes.

## Timesteppers and ODEs

The [[waiwera:timestepper]] object is an instance of [[timestepper_type]]. This derived type is designed to step an ordinary differential equation (ODE) through time. It can do this with an arbitrary ODE of the form:

$$\frac{d}{dt} \mathit{LHS}(t, \vec{Y}) = \mathit{RHS}(t, \vec{Y})$$

where \(t\) is time and \(Y\) is a vector of unknowns (which can be distributed over multiple processors) over a discretized finite-volume mesh.

At each time step a set of non-linear residual equations is solved to update the solution vector \(Y\). A PETSc SNES parallel non-linear solver is used for this. The exact form of the non-linear equations depends on the time-stepping method (e.g. backwards Euler) being used.

The details of the ODE to be solved (including the LHS and RHS functions) are represented by an [[ode_type]] object. This is an abstract type (i.e. it does not implement any specific ODE). It contains a [[ode_type(type):solution]] component, which is a PETSc Vec object representing the vector \(Y\) above.

## Simulation object

The [[waiwera:simulation]] object is an instance of [[flow_simulation_type]], which extends the abstract [[ode_type]] to implement the particular ODE to be solved for simulating non-isothermal subsurface flow. This ODE is derived from discretized mass and energy balance equations to be solved over the finite-volume mesh.

In a [[flow_simulation_type]] object, the [[ode_type(type):solution]] component vector (inherited from its parent [[ode_type]]) represents the primary thermodynamic variables in each mesh cell.

The LHS function represents the mass and energy balances in each cell. The RHS function represents fluid inflows to each cell, via flows through the faces or source/sink terms.

A [[flow_simulation_type]] object also contains vector components for storing the [[flow_simulation_type(type):rock]] and [[flow_simulation_type(type):fluid]] properties in each mesh cell. Currently the [[flow_simulation_type(type):rock]] vector is simply initialized at the start of the simulation and then left unchanged during the run. However the option of including time-varying rock properties (e.g. from rock mechanics calculations) is left open.

The [[flow_simulation_type(type):fluid]] vector represents the full set of fluid properties (including bulk properties like pressure and temperature, and also properties of each phase like liquid density or vapour viscosity). The initial values of the fluid vector are calculated from the simulation initial conditions (i.e. initial solution vector) via the [[flow_simulation_type(type):fluid_init]] type-bound procedure. Subsequently, the fluid vector is updated from any given vector of primary variables via the [[flow_simulation_type(type):fluid_properties]] type-bound procedure. This is called before each function evaluation of the non-linear equation solution procedure.

## Non-linear solution process

The SNES non-linear solver uses Newton-Raphson iteration to solve the residual equations at each time step. At the beginning of each iteration, PETSc calculates the Jacobian matrix using finite differencing, incrementing each of the primary variables in turn and re-evaluating the residual function. Hence, the vector of fluid properties has to be updated for each of these function evaluations. However, the thermodynamic region and phase conditions are not updated for these small increments in primary variables (even if the incremented variables do in fact cross over into another thermodynamic region).

Once the Jacobian matrix has been calculated, the Newton iteration update is performed. The SNES solver includes capability for a subsequent "line search" after the default Newton update. However, this is disabled here, because changes to the primary variables during the line search are not necessarily small, and could result in thermodynamic region transitions. Because variable switching is used in these transitions, this would result in a discontinuous change in the solution vector during the line search, meaning it would not be able to converge.

After the line search (if any) is applied, and the final solution vector for the iteration has been determined, the SNES solver allows a user-defined function to be called. Here this is used to carry out any thermodynamic region transitions needed (via the [[flow_simulation_type(type):fluid_transitions]] type-bound procedure). This in turn calls the [[eos_type(type):transition]] procedure of the simulation's [[flow_simulation_type(type):eos]] component to carry out the region transition in a particular cell.

## Local-level objects: cell, face, rock and fluid

The [[cell_type(type)]] and [[face_type(type)]] derived types are provided for manipulating mesh cells and faces between them, when carrying out local-level computations. For computations within each cell, there are also [[rock_type(type)]] and [[fluid_type(type)]] derived types, for accessing and manipulating the rock and fluid properties in the cell. These types provide an object-oriented way of manipulating the data in the rock and fluid vectors.

The main components of these types are pointers which point to the corresponding entries in the vectors. The [[rock_type(type)]] has a [[rock_type(type):assign]] procedure which assigns the pointer components, and the other types have a similar assignment procedure. Because the components of these types are pointers, changing their component values updates the entries in the corresponding vector.

These local types are also nested: a face contains two cell components, and each cell contains a rock and a fluid component.
