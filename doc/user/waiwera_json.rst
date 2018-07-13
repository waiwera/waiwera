===========================
Waiwera JSON file structure
===========================

All input data for a Waiwera simulation are contained within one object in the input JSON file.

The various aspects of the simulation are specified by different named values within that object, as listed below. These values are of various types, and many are objects themselves, with their own internal structure.

==============  ===================== =========================
Name            Type                  Specifies
==============  ===================== =========================
title           string                simulation title
mesh            string / object       simulation mesh
rock            object                rock properties
boundaries      array                 boundary conditions
source          array                 source and sink terms
initial         object                initial conditions
gravity         number / array / null gravity
thermodynamics  string                thermodynamic formulation
eos             object                equation of state
time            object                time stepping
logfile         Boolean / object      output log file
output          Boolean / object      output results file
==============  ===================== =========================

The exact structure of a Waiwera input JSON file is defined by the JSON `schema <http://json-schema.org/>`_ included in the `utils/` directory of the Waiwera source code.
