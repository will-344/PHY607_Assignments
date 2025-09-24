""" Module that provides numerical solutions to the following problems, with
necessary constants/conditions taken from function inputs:

ODE: Steady 1D cylindrical heat conduction with no energy generation, assuming
    constant thermal conductivity, with specified temperature boundary
    conditions.
     
definite integral: Lift from an e lliptical lift distribution for a finite
    wing, assuming incompressible, inviscid, and steady flow, using lifting-
    line theory in a uniform freestream with constant density.
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

def rk4(temp1, temp2, rad1, rad2, d_rad):
    x1 = np.array([temp1])
    x2 = np.array([(temp2-temp1) / (rad1 * np.log(rad2/rad1))])
    iterations = round((rad2-rad1) / d_rad)
    rad = np.linspace(rad1, rad2, iterations)
    
    for i in range(iterations-1):
        k1x1 = x2[i]
        k1x2 = -x2[i] / rad[i]
        
        x2a = x2[i] + k1x2*d_rad/2
        k2x1 = x2a
        k2x2 = -x2a / (rad[i] + d_rad/2)
        
        x2b = x2[i] + k2x2*d_rad/2
        k3x1 = x2b
        k3x2 = -x2b / (rad[i] + d_rad/2)
        
        x2c = x2[i] + k3x2*d_rad
        k4x1 = x2c
        k4x2 = -x2c / (rad[i] + d_rad)
        
        x1new = x1[i] + (k1x1 + 2*k2x1 + 2*k3x1 + k4x1)*d_rad/6
        x2new = x2[i] + (k1x2 + 2*k2x2 + 2*k3x2 + k4x2)*d_rad/6
        
        x1 = np.append(x1, x1new)
        x2 = np.append(x2, x2new)
    
    temp = x1
    return temp, rad

def riemann(rho, vel, lamda, span, intervals):
    area = 0
    width = span/intervals
    intervals += 1
    height = np.array([])
    y = np.linspace(-span/2, span/2, intervals)
    for i in range(intervals):
        height = np.append(height, 
                                rho*vel*lamda * np.sqrt(1 - (2*y[i]/span)**2))
        if i == intervals-1:
            break
        area += width*height[i]
    return area, height, y

def riemann_plot(height, y):
    intervals = len(y)
    height_bins = np.array([])
    y_bins = np.array([])
    for i in range(intervals):
        height_bins = np.append(height_bins, height[i])
        y_bins = np.append(y_bins, y[i])
        if i == intervals-1:
            break
        height_bins = np.append(height_bins, height[i])
        y_bins = np.append(y_bins, y[i+1])
    return height_bins, y_bins
    
