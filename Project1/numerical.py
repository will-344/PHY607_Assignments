""" Module that provides numerical solution to the following problems, with
necessary constants/conditions taken from function inputs:

ode1 (ODE): Steady 1D cylindrical heat conduction with no energy generation,
    assuming constant thermal conductivity, with specified temperature boundary
    conditions.
     
int1 (definite integral): Elliptical lift distribution for a finite wing,
    assuming incompressible, inviscid, and steady flow, using lifting-line
    theory in a uniform freestream with constant density.
"""

import numpy as np

def euler(temp1, temp2, rad1, rad2, d_rad):
    x1 = np.array([temp1])
    x2 = np.array([(temp2-temp1) / (rad1 * np.log(rad2/rad1))])
    iterations = round((rad2-rad1) / d_rad)
    rad = np.linspace(rad1, rad2, iterations)
    
    for i in range(iterations-1):
        x1new = x1[i] + x2[i]*d_rad
        x2new = x2[i] + -x2[i]/rad[i]*d_rad
        x1 = np.append(x1, x1new)
        x2 = np.append(x2, x2new)
    
    temp = x1
    return temp, rad
