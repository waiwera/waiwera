*****************
How Waiwera works
*****************

.. focus on what user needs to know about to set up simulation and interpret results

.. solution scheme:
.. finite volume discretisation of mass and energy balance equations
.. two-point flux approximation, upstream weighting
.. unknowns at each time step are the (non-dimensionalised) primary variables
.. Non-linear solve at each time step
.. Linear solver at each non-linear solver iteration
.. phase transitions: variable switching   
.. time stepping: abstracted
.. EOS handles calculation of fluid (mixture) properties from primary variables, phase transitions etc.
.. MINC for fractured media

.. software design: PETSc library (DMPlex for meshes, Vec/Mat for linear algebra, SNES, KSP), Fortran 2003
.. parallelism: MPI
