API Documentation
=================
===================================================
 Lipid order --- 
===================================================

:Author: Ricardo Ramirez
:Year: 2024



This code is intended to compute the order parameters of lipids from all atom MD simulations.
The functions used to compute order parameters are sn1 and sn2 for the first and the second chain correpondingly.
sni function calls to the correponding function order_sni and it calls to the function get_vectors. Thus, getting the order
parameters for a period of time using the following equation:

.. math:: SCD = \left|\frac{3}{2}cos^2(\theta_i)-\frac{1}{2}\right|





Module
------

.. autosummary::
   :toctree: autosummary
   :recursive:

.. automodule:: lipidorder.lipid_order
   :members: