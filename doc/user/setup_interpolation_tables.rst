.. index:: interpolation tables
.. _interpolation_tables:

********************
Interpolation tables
********************

Some Waiwera input parameters are specified in the form of tables of values, tabulated against some independent variable. Most commonly the independent variable is time, so the tables represent the time variation of the value. For example:

.. code-block:: json

   {"source": [
     {"cell": 99, "rate": [[0, -10.5], [3600, -9.9], [7200, -7.7]]}
   ]}

specifies a source producing from a cell, with flow rate decreasing with time.

Typically, it is necessary to perform two kinds of operation on the tabulated values:

* **interpolation**, to find values corresponding to times in between those given in the table
* **averaging**, to find the average value over a given time interval (e.g. a time step)

Waiwera provides different options for both interpolation and averaging of tabulated data. Different options may be appropriate depending on the quantity represented by the data.

.. index:: interpolation tables; interpolation
.. _table_interpolation:

Interpolation
=============

Waiwera offers three types of interpolation (see :numref:`interpolation_fig`), depending on how the data are assumed to vary between their tabulated values:

#. **Step** interpolation, in which the data are assumed to have a simple piecewise constant variation. This may be appropriate for quantities which generally remain constant but are adjusted at particular times. The value in between two tabulated times is simply assumed equal to the last tabulated value.
#. **Linear** interpolation, in which the data are assumed to vary linearly between the tabulated values. This is often appropriate for quantities that change continuously, with the tabulated values representing time-series measurements of the quantity over time.
#. **PCHIP** (Piecewise Cubic Hermite Interpolation Polynomial), a method which gives a smooth, differentiable curve between data points, without introducing any local extrema (i.e. oscillations) between them. The PCHIP method was introduced by Fritsch and Carlson (1980) [1]_ and refined by Fritsch and Butland (1984) [2]_.

.. _interpolation_fig:
.. figure:: interpolation.*
           :scale: 67 %
           :align: center

           Interpolation types

For points outside the table (i.e. before the first tabulated time or after the last one), the data are extrapolated by assuming a piecewise constant variation, regardless of the interpolation type.
 
Usually a JSON value containing tabulated data in the Waiwera input file will also have an associated **"interpolation"** string value. This can be set to **"step"**, **"linear"** or **"pchip"** to control the interpolation type. If the interpolation type is not specified, linear interpolation is used by default.

For example:

.. code-block:: json

   {"source": [
     {"cell": 99,
      "rate": [[0, -10.5], [3600, -9.9], [7200, -7.7]],
      "interpolation": "step"}
   ]}

specifies a source with step interpolation on the flow rate table.

.. index:: interpolation tables; averaging
.. _table_averaging:

Averaging
=========

Waiwera also offers two methods for averaging an interpolated data table over an interval.

* The first method uses integration. The integral of the data (using the specified interpolation method) is divided by the time interval to calculate the average value. This method is usually preferred, especially in situations when it is necessary to preserve the correct area under the data curve.
* The second method is simpler, calculating an approximate result from the average of the interpolated values at the start- and end-points of the interval. This method can be less expensive if there are many data points within the time interval. However, it is not as accurate, particularly if there are many tabulated values within the integration interval.

Usually a JSON value containing tabulated data in the Waiwera input file will also have an associated **"averaging"** string value, which can be set to either **"integrate"** or **"endpoint"** to control the averaging method. If the averaging method is not specified, the integration method is used by default.

For example:

.. code-block:: json

   {"source": [
     {"cell": 99,
      "rate": [[0, -10.5], [3600, -9.9], [7200, -7.7]],
      "interpolation": "linear", "averaging": "integrate"}
   ]}

specifies a source with linear interpolation and the integration averaging method.

.. [1] Fritsch, F.N. and Carlson, R.E. (1980), *Monotone piecewise cubic interpolation*, SIAM Journal on Numerical Analysis 17(2), pp. 238--246.
.. [2] Fritsch, F.N. and Butland, J. (1984), *A method for constructing local monotone piecewise cubic interpolants*, SIAM Journal on Scientific and Statistical Computing 5(2), pp. 300--304.
