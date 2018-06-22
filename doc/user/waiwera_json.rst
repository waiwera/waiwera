======================================
How are Waiwera JSON files structured?
======================================

All input data for a Waiwera simulation are contained within one object in the input JSON file.

The various aspects of the simulation are specified by different named values within that object. These values are of various types, and many are objects themselves. They are listed in the table below.

==============  =====================
Name            Type
==============  =====================
title           string
mesh            string / object
rock            object
boundaries      array
source          array
initial         object
gravity         number / array / null   
thermodynamics  string
eos             object
time            object
logfile         Boolean / object
output          Boolean / object
==============  =====================

