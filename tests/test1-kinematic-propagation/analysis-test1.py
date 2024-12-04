#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script description: This script does takes the kinematics and action computed by AIS++ and compares it to analytical solutions.
Author: Noam Mouelle
Date: 08-10-2024

Note: ensure that the physical constants are the same as in the AIS++ code!!
"""

import unittest
import numpy as np
from math import isclose, sqrt, cos, sin, cosh, sinh

# define constants here
g = 9.81
R = 6.37e6
m = 1.44e-25
h = 6.62607015e-34
pi = 3.141592653589793
hbar = h / (2 * pi)

def analytical_propagate_free(x0, y0, z0, vx0, vy0, vz0, t):
    """
    Analytical solution for a free particle propagation.
    """
    x = x0 + vx0 * t
    y = y0 + vy0 * t
    z = z0 + vz0 * t
    
    vx = vx0
    vy = vy0
    vz = vz0

    return x, y, z, vx, vy, vz

def analytical_propagate_uniform_acceleration(x0, y0, z0, vx0, vy0, vz0, t):
    """
    Analytical solution for a particle with uniform acceleration.
    """
    x = x0 + vx0 * t
    y = y0 + vy0 * t
    z = z0 + vz0 * t - 0.5 * g * t**2

    vx = vx0
    vy = vy0
    vz = vz0 - g * t

    return x, y, z, vx, vy, vz

def analytical_propagate_linear(x0, y0, z0, vx0, vy0, vz0, t):
    """
    Analytical solution for a particle with linear acceleration.
    """
    x = x0 * cos(sqrt(g/R) * t) + vx0 * sqrt(R/g) * sin(sqrt(g/R) * t)
    y = y0 * cos(sqrt(g/R) * t) + vy0 * sqrt(R/g) * sin(sqrt(g/R) * t)
    z = 1/2 * R - 1/2 * (R - 2*z0) * cosh(sqrt(2*g/R) * t) + vz0 * sqrt(R/(2*g)) * sinh(sqrt(2*g/R) * t)

    vx = -x0 * sqrt(g/R) * sin(sqrt(g/R) * t) + vx0 * cos(sqrt(g/R) * t)
    vy = -y0 * sqrt(g/R) * sin(sqrt(g/R) * t) + vy0 * cos(sqrt(g/R) * t)
    vz = -1/2 * (R - 2*z0) * sqrt(2*g/R) * sinh(sqrt(2*g/R) * t) + vz0 * cosh(sqrt(2*g/R) * t)

    return x, y, z, vx, vy, vz

def analytical_propagate(x0,y0,z0,vx0,vy0,vz0,t,acceleration):
    if acceleration == 'free':
        return analytical_propagate_free(x0,y0,z0,vx0,vy0,vz0,t)
    elif acceleration == 'uniform':
        return analytical_propagate_uniform_acceleration(x0,y0,z0,vx0,vy0,vz0,t)
    elif acceleration == 'linear':
        return analytical_propagate_linear(x0,y0,z0,vx0,vy0,vz0,t)
    
def Sfree(x0,y0,z0,vx0,vy0,vz0,t):
    """
    Action for a free particle.
    """
    return 0.5 * m * t * (vx0**2 + vy0**2 + vz0**2)

def Suniform(x0,y0,z0,vx0,vy0,vz0,t):
    """
    Action for a particle with uniform acceleration.
    """
    return 0.5 * m * t * (vx0**2 + vy0**2 + vz0**2) - m * g * t**2 * vz0 + m * g ** 2 * t**3 / 3 - m * g * t * z0

def S_linear(x0, y0, z0, vx0, vy0, vz0, t):
    """
    Action for a particle with linear acceleration.
    """
    part_1 = -2*g*R*t + 3*R*vz0 - 4*vx0*x0 - 4*vy0*y0 - 6*vz0*z0
    part_2 = 4*(vx0*x0 + vy0*y0)*cos(2 * sqrt(g/R) * t)
    part_3 = -2 * vz0 * (R - 2*z0) * cosh(sqrt(2*g/R) * t) ** 2
    part_4 = - vz0 * (R - 2*z0) * cosh(2 * sqrt(2*g/R) * t)
    part_5 = 4 * (R * (vx0**2 + vy0**2) - g * (x0**2 + y0**2)) * sin(2 * sqrt(g/R) * t) / (2 * sqrt(g * R))
    part_6 = sqrt(2) * (2 * R * vz0**2 + g * (R - 2 * z0)**2) * sinh(2 * sqrt(2*g/R) * t) / (2 * sqrt(g * R))
    
    # Final expression
    result = (1 / 8) * m * (part_1 + part_2 + part_3 + part_4 + part_5 + part_6)
    
    return result

from mpmath import mp, sqrt, cos, cosh, sin, sinh

# Set the precision to quad precision (approximately 113 bits, equivalent to 34 decimal places)
mp.dps = 34  # Set decimal precision to 34 digits

# define constants using mpmath for quad precision
g = mp.mpf('9.81')
R = mp.mpf('6.37e6')
au = mp.mpf('1.660539066e-27')
m = mp.mpf('86.90888') * au
h = mp.mpf('6.62607015e-34')
pi = mp.pi
hbar = h / (2 * pi)

def analytical_propagate_free(x0, y0, z0, vx0, vy0, vz0, t):
    """
    Analytical solution for a free particle propagation using quad precision.
    """
    x = x0 + vx0 * t
    y = y0 + vy0 * t
    z = z0 + vz0 * t
    
    vx = vx0
    vy = vy0
    vz = vz0

    return x, y, z, vx, vy, vz

def analytical_propagate_uniform_acceleration(x0, y0, z0, vx0, vy0, vz0, t):
    """
    Analytical solution for a particle with uniform acceleration using quad precision.
    """
    x = x0 + vx0 * t
    y = y0 + vy0 * t
    z = z0 + vz0 * t - mp.mpf('0.5') * g * t**2

    vx = vx0
    vy = vy0
    vz = vz0 - g * t

    return x, y, z, vx, vy, vz

def analytical_propagate_linear(x0, y0, z0, vx0, vy0, vz0, t):
    """
    Analytical solution for a particle with linear acceleration using quad precision.
    """
    sqrt_g_R = sqrt(g / R)
    sqrt_2g_R = sqrt(2 * g / R)
    
    x = x0 * cos(sqrt_g_R * t) + vx0 * sqrt(R/g) * sin(sqrt_g_R * t)
    y = y0 * cos(sqrt_g_R * t) + vy0 * sqrt(R/g) * sin(sqrt_g_R * t)
    z = mp.mpf('0.5') * R - mp.mpf('0.5') * (R - 2*z0) * cosh(sqrt_2g_R * t) + vz0 * sqrt(R/(2*g)) * sinh(sqrt_2g_R * t)

    vx = -x0 * sqrt_g_R * sin(sqrt_g_R * t) + vx0 * cos(sqrt_g_R * t)
    vy = -y0 * sqrt_g_R * sin(sqrt_g_R * t) + vy0 * cos(sqrt_g_R * t)
    vz = -mp.mpf('0.5') * (R - 2*z0) * sqrt_2g_R * sinh(sqrt_2g_R * t) + vz0 * cosh(sqrt_2g_R * t)

    return x, y, z, vx, vy, vz

def analytical_propagate(x0, y0, z0, vx0, vy0, vz0, t, acceleration):
    if acceleration == 'free':
        return analytical_propagate_free(x0, y0, z0, vx0, vy0, vz0, t)
    elif acceleration == 'uniform':
        return analytical_propagate_uniform_acceleration(x0, y0, z0, vx0, vy0, vz0, t)
    elif acceleration == 'linear':
        return analytical_propagate_linear(x0, y0, z0, vx0, vy0, vz0, t)
    
def Sfree(x0, y0, z0, vx0, vy0, vz0, t):
    """
    Action for a free particle using quad precision.
    """
    return mp.mpf('0.5') * m * t * (vx0**2 + vy0**2 + vz0**2)

def Suniform(x0, y0, z0, vx0, vy0, vz0, t):
    """
    Action for a particle with uniform acceleration using quad precision.
    """
    return (mp.mpf('0.5') * m * t * (vx0**2 + vy0**2 + vz0**2) 
            - m * g * t**2 * vz0 
            + m * g**2 * t**3 / mp.mpf('3') 
            - m * g * t * z0)

def S_linear(x0, y0, z0, vx0, vy0, vz0, t):
    """
    Action for a particle with linear acceleration using quad precision.
    """
    sqrt_g_R = sqrt(g / R)
    sqrt_2g_R = sqrt(2 * g / R)
    
    part_1 = -mp.mpf('2') * g * R * t + mp.mpf('3') * R * vz0 - mp.mpf('4') * vx0 * x0 - mp.mpf('4') * vy0 * y0 - mp.mpf('6') * vz0 * z0
    part_2 = mp.mpf('4') * (vx0 * x0 + vy0 * y0) * cos(mp.mpf('2') * sqrt_g_R * t)
    part_3 = -mp.mpf('2') * vz0 * (R - mp.mpf('2') * z0) * cosh(sqrt_2g_R * t) ** 2
    part_4 = - vz0 * (R - mp.mpf('2') * z0) * cosh(mp.mpf('2') * sqrt_2g_R * t)
    part_5 = mp.mpf('4') * (R * (vx0**2 + vy0**2) - g * (x0**2 + y0**2)) * sin(mp.mpf('2') * sqrt_g_R * t) / (mp.mpf('2') * sqrt(g * R))
    part_6 = sqrt(mp.mpf('2')) * (mp.mpf('2') * R * vz0**2 + g * (R - mp.mpf('2') * z0)**2) * sinh(mp.mpf('2') * sqrt_2g_R * t) / (mp.mpf('2') * sqrt(g * R))
    
    # Final expression
    result = (mp.mpf('1') / mp.mpf('8')) * m * (part_1 + part_2 + part_3 + part_4 + part_5 + part_6)
    
    return result


def S(x0,y0,z0,vx0,vy0,vz0,t,acceleration):
    if acceleration == 'free':
        return Sfree(x0,y0,z0,vx0,vy0,vz0,t)
    elif acceleration == 'uniform':
        return Suniform(x0,y0,z0,vx0,vy0,vz0,t)
    elif acceleration == 'linear':
        return S_linear(x0,y0,z0,vx0,vy0,vz0,t)
    
def get_errors(a,b):
    abs_err = np.abs(a-b)
    if b == 0:
        return abs_err, 0
    else:
        return abs_err, abs_err/np.abs(b)

def check_propagator_result(filename_initial, filename_final, tol, dt, acceleration):
    # load the initial and final states
    _,_,_,x0list,y0list,z0list,vx0list,vy0list,vz0list,_,_,_,_=np.loadtxt(filename_initial, unpack=True,skiprows=1,delimiter=',')
    _,_,_,x1list,y1list,z1list,vx1list,vy1list,vz1list,phase_double_list,phase_err_list,_,_=np.loadtxt(filename_final, unpack=True,skiprows=1,delimiter=',')
    # loop over the data and compare to analytical calculations
    abs_error_x, abs_error_y, abs_error_z, abs_error_vx, abs_error_vy, abs_error_vz, abs_error_phase_double = [], [], [], [], [], [], []
    rel_error_x, rel_error_y, rel_error_z, rel_error_vx, rel_error_vy, rel_error_vz, rel_error_phase_double = [], [], [], [], [], [], []
    for x0,y0,z0,vx0,vy0,vz0,x1,y1,z1,vx1,vy1,vz1,phase_double,phase_err in zip(x0list,y0list,z0list,vx0list,vy0list,vz0list,x1list,y1list,z1list,vx1list,vy1list,vz1list,phase_double_list,phase_err_list):
        # propagate the initial state
        x1a, y1a, z1a, vx1a, vy1a, vz1a = analytical_propagate(x0, y0, z0, vx0, vy0, vz0, dt, acceleration)
        phase_double_a = S(x0,y0,z0,vx0,vy0,vz0,dt,acceleration) / hbar
        # assert that numerical and analytical results are close, and collect the errors
        abs_error_x_, rel_error_x_ = get_errors(x1, x1a)
        abs_error_y_, rel_error_y_ = get_errors(y1, y1a)
        abs_error_z_, rel_error_z_ = get_errors(z1, z1a)
        abs_error_vx_, rel_error_vx_ = get_errors(vx1, vx1a)
        abs_error_vy_, rel_error_vy_ = get_errors(vy1, vy1a)
        abs_error_vz_, rel_error_vz_ = get_errors(vz1, vz1a)
        abs_error_phase_double_, rel_error_phase_double_ = get_errors(phase_double, phase_double_a)

        abs_error_x.append(abs_error_x_)
        abs_error_y.append(abs_error_y_)
        abs_error_z.append(abs_error_z_)
        abs_error_vx.append(abs_error_vx_)
        abs_error_vy.append(abs_error_vy_)
        abs_error_vz.append(abs_error_vz_)
        abs_error_phase_double.append(abs_error_phase_double_)

        rel_error_x.append(rel_error_x_)
        rel_error_y.append(rel_error_y_)
        rel_error_z.append(rel_error_z_)
        rel_error_vx.append(rel_error_vx_)
        rel_error_vy.append(rel_error_vy_)
        rel_error_vz.append(rel_error_vz_)
        rel_error_phase_double.append(rel_error_phase_double_)

    # convert all the list elements to floats
    abs_error_x = [float(i) for i in abs_error_x]
    abs_error_y = [float(i) for i in abs_error_y]
    abs_error_z = [float(i) for i in abs_error_z]
    abs_error_vx = [float(i) for i in abs_error_vx]
    abs_error_vy = [float(i) for i in abs_error_vy]
    abs_error_vz = [float(i) for i in abs_error_vz]
    abs_error_phase_double = [float(i) for i in abs_error_phase_double]
    rel_error_x = [float(i) for i in rel_error_x]
    rel_error_y = [float(i) for i in rel_error_y]
    rel_error_z = [float(i) for i in rel_error_z]
    rel_error_vx = [float(i) for i in rel_error_vx]
    rel_error_vy = [float(i) for i in rel_error_vy]
    rel_error_vz = [float(i) for i in rel_error_vz]
    rel_error_phase_double = [float(i) for i in rel_error_phase_double]

    # print the max abs and rel errors, as well as the mean abs and rel errors
    print('x: max abs err =', f'{max(abs_error_x):.3g}', 'mean abs err =', f'{np.mean(abs_error_x):.3g}', 'max rel err =', f'{max(rel_error_x):.3g}', 'mean rel err =', f'{np.mean(rel_error_x):.3g}')
    print('y: max abs err =', f'{max(abs_error_y):.3g}', 'mean abs err =', f'{np.mean(abs_error_y):.3g}', 'max rel err =', f'{max(rel_error_y):.3g}', 'mean rel err =', f'{np.mean(rel_error_y):.3g}')
    print('z: max abs err =', f'{max(abs_error_z):.3g}', 'mean abs err =', f'{np.mean(abs_error_z):.3g}', 'max rel err =', f'{max(rel_error_z):.3g}', 'mean rel err =', f'{np.mean(rel_error_z):.3g}')
    print('vx: max abs err =', f'{max(abs_error_vx):.3g}', 'mean abs err =', f'{np.mean(abs_error_vx):.3g}', 'max rel err =', f'{max(rel_error_vx):.3g}', 'mean rel err =', f'{np.mean(rel_error_vx):.3g}')
    print('vy: max abs err =', f'{max(abs_error_vy):.3g}', 'mean abs err =', f'{np.mean(abs_error_vy):.3g}', 'max rel err =', f'{max(rel_error_vy):.3g}', 'mean rel err =', f'{np.mean(rel_error_vy):.3g}')
    print('vz: max abs err =', f'{max(abs_error_vz):.3g}', 'mean abs err =', f'{np.mean(abs_error_vz):.3g}', 'max rel err =', f'{max(rel_error_vz):.3g}', 'mean rel err =', f'{np.mean(rel_error_vz):.3g}')
    print('phase_double: max abs err =', f'{max(abs_error_phase_double):.3g}', 'mean abs err =', f'{np.mean(abs_error_phase_double):.3g}', 'max rel err =', f'{max(rel_error_phase_double):.3g}', 'mean rel err =', f'{np.mean(rel_error_phase_double):.3g}')   


    for x0,y0,z0,vx0,vy0,vz0,x1,y1,z1,vx1,vy1,vz1,phase_double,phase_err in zip(x0list,y0list,z0list,vx0list,vy0list,vz0list,x1list,y1list,z1list,vx1list,vy1list,vz1list,phase_double_list,phase_err_list):
        x1a, y1a, z1a, vx1a, vy1a, vz1a = analytical_propagate(x0, y0, z0, vx0, vy0, vz0, dt, acceleration)
        phase_double_a = S(x0,y0,z0,vx0,vy0,vz0,dt,acceleration) / hbar

        assert isclose(x1, x1a, abs_tol=tol)
        assert isclose(y1, y1a, abs_tol=tol)
        assert isclose(z1, z1a, abs_tol=tol)
        assert isclose(vx1, vx1a, abs_tol=tol)
        assert isclose(vy1, vy1a, abs_tol=tol)
        assert isclose(vz1, vz1a, abs_tol=tol)
        assert isclose(phase_double, phase_double_a, abs_tol=tol)

    print('\nAll tests passed!')

def main():
    print("TEST 1: KINEMATIC PROPAGATION")
    tol = 1e-9
    dt  = mp.mpf('1')
    # check the free particle propagator
    print('\nFree particle propagator')
    check_propagator_result('output-files/zero_pot_initial_WPK.txt', 'output-files/zero_pot_final_WPK.txt', tol, dt, 'free')
    # check the uniform acceleration propagator
    print('\nUniform acceleration propagator')
    check_propagator_result('output-files/linear_pot_initial_WPK.txt', 'output-files/linear_pot_final_WPK.txt', tol, dt, 'uniform')
    # check the linear acceleration propagator
    print('\nLinear acceleration propagator')
    check_propagator_result('output-files/quadratic_pot_initial_WPK.txt', 'output-files/quadratic_pot_final_WPK.txt', tol, dt, 'linear')

if __name__ == '__main__':
    main()
