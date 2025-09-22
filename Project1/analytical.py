""" Module that provides the analytic solution to the following problems, with
necessary constants/conditions taken from function inputs:

ode1 (ODE): Steady 1D cylindrical heat conduction with no energy generation,
    assuming constant thermal conductivity, with specified temperature boundary
    conditions.
     
int1 (definite integral): Lift from an elliptical lift distribution for a
    finite wing, assuming incompressible, inviscid, and steady flow, using 
    lifting-line theory in a uniform freestream with constant density.
"""

import numpy as np

def ode1(temp1, temp2, rad1, rad2, d_rad):
    n = round((rad2-rad1) / d_rad)
    rad = np.linspace(rad1, rad2, n)
    temp = np.log(rad/rad1) / np.log(rad2/rad1) * (temp2-temp1) + temp1
    return temp, rad

def int1(rho, vel, lamda, span):
    lift = rho*vel*lamda*span*np.pi/4
    return lift
