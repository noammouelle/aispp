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
from mpmath import mp, sqrt, cos, cosh, sin, sinh
import pandas as pd

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
freq0 = mp.mpf('429228004229873')
omega0 = 2 * pi * freq0

def analytical_propagate_free(x0, y0, z0, vx0, vy0, vz0, t, state, k):
    """
    Analytical solution for a free particle propagation using quad precision.
    """

    if state == 0:
        x = x0 + vx0 * t
        y = y0 + vy0 * t
        z = z0 + vz0 * t
        vx = vx0
        vy = vy0
        vz = vz0
    else:
        x = x0 + vx0 * t
        y = y0 + vy0 * t
        z = z0 + vz0 * t
        vx = vx0
        vy = vy0
        vz = vz0 + hbar * k / m

    return x, y, z, vx, vy, vz
    
def Sfree(x0, y0, z0, vx0, vy0, vz0, t, state, k):
    """
    Action for a free particle using quad precision.
    """
    if state == 0:
        vx = vx0
        vy = vy0
        vz = vz0
    else:
        vx = vx0
        vy = vy0
        vz = vz0 + hbar * k / m

    return mp.mpf('0.5') * m * t * (vx**2 + vy**2 + vz**2)

def get_delta(vz0, k, omega):
    return k * vz0 - (omega - omega0) + hbar * k**2 / (2 * m)

def analytical_phase_double(x0, y0, z0, vx0, vy0, vz0, t, state, k, omega, rabifreq):
    delta = get_delta(vz0, k, omega)
    X     = sqrt(rabifreq**2 + delta**2)
    phi   = X/2 * t
    action_phase = Sfree(x0, y0, z0, vx0, vy0, vz0, t, state, k) / hbar

    phase = action_phase - delta/2*t

    if state == 0:
        phase += mp.arg(cos(phi) + 1j * delta / X * sin(phi))
    else:
        phase += k*(z0 - hbar*k/(2*m)*t) - pi/2

    return phase

def analytical_phase_quad(t, state, omega):
    if state == 0:
        return 0
    else:
        return -omega*t
    
def get_amplitude(vz0, k, t, omega, rabifreq,state):
    delta = get_delta(vz0, k, omega)
    X     = sqrt(rabifreq**2 + delta**2)
    phi   = X/2 * t
    if state==0:
        return mp.fabs(cos(phi) - 1j * delta / X * sin(phi))
    else:
        return mp.fabs(rabifreq*sin(phi)/X)
    
def get_errors(a,b):
    abs_err = np.abs(a-b)
    if b == 0:
        return abs_err, 0
    else:
        return abs_err, abs_err/np.abs(b)

def get_kz_value(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('kz'):
                _, kz_value = line.split()
                return float(kz_value)
    return None

def get_omega_value(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('omega'):
                _, omega_value = line.split()
                return float(omega_value)
    return None

def get_rabifreq_value(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('rabifreq'):
                _, rabifreq = line.split()
                return float(rabifreq)*2*pi
    return None

def check_results(initial_df, final_df, k, omega, rabifreq):
    # assume all waps in the final_df are in the same state
    state = final_df['State'].values[0]
    # create a dictionary to store the errors
    errors = {'X_abs_err': [], 'X_rel_err': [], 'Y_abs_err': [], 'Y_rel_err': [],
              'Z_abs_err': [], 'Z_rel_err': [], 'VX_abs_err': [], 'VX_rel_err': [],
              'VY_abs_err': [], 'VY_rel_err': [], 'VZ_abs_err': [], 'VZ_rel_err': [],
              'Phase_abs_err': [], 'Phase_rel_err': [], 'PhaseQuad_abs_err': [], 'PhaseQuad_rel_err': [],
              'Amplitude_abs_err': [], 'Amplitude_rel_err': []}
    
    # loop over the rows of the dataframe
    for index, row_initial in initial_df.iterrows():
        row_final = final_df.iloc[index]
        # get the analytical kinematics
        x, y, z, vx, vy, vz = analytical_propagate_free(row_initial['X'], row_initial['Y'], row_initial['Z'],
                                                         row_initial['VX'], row_initial['VY'], row_initial['VZ'],
                                                         row_final['Time'] - row_initial['Time'], state, k)
        # get the analytical phase
        phase = analytical_phase_double(row_initial['X'], row_initial['Y'], row_initial['Z'],
                                        row_initial['VX'], row_initial['VY'], row_initial['VZ'],
                                        row_final['Time'] - row_initial['Time'], state, k, omega, rabifreq)
        # get the analytical phasequad
        phase_quad = analytical_phase_quad(row_final['Time'] - row_initial['Time'], state, omega)
        # get the analytical amplitude
        amplitude = get_amplitude(row_initial['VZ'], k, row_final['Time'] - row_initial['Time'], omega, rabifreq, state)

        # get the errors
        x_abs_err, x_rel_err = get_errors(row_final['X'], x)
        y_abs_err, y_rel_err = get_errors(row_final['Y'], y)
        z_abs_err, z_rel_err = get_errors(row_final['Z'], z)
        vx_abs_err, vx_rel_err = get_errors(row_final['VX'], vx)
        vy_abs_err, vy_rel_err = get_errors(row_final['VY'], vy)
        vz_abs_err, vz_rel_err = get_errors(row_final['VZ'], vz)
        phase_abs_err, phase_rel_err = get_errors(row_final['Phase'], phase)
        phase_quad_abs_err, phase_quad_rel_err = get_errors(row_final['PhaseQuad'], phase_quad)
        amplitude_abs_err, amplitude_rel_err = get_errors(row_final['Amplitude'], amplitude)

        # append the errors to the dictionary
        errors['X_abs_err'].append(x_abs_err)
        errors['X_rel_err'].append(x_rel_err)
        errors['Y_abs_err'].append(y_abs_err)
        errors['Y_rel_err'].append(y_rel_err)
        errors['Z_abs_err'].append(z_abs_err)
        errors['Z_rel_err'].append(z_rel_err)
        errors['VX_abs_err'].append(vx_abs_err)
        errors['VX_rel_err'].append(vx_rel_err)
        errors['VY_abs_err'].append(vy_abs_err)
        errors['VY_rel_err'].append(vy_rel_err)
        errors['VZ_abs_err'].append(vz_abs_err)
        errors['VZ_rel_err'].append(vz_rel_err)
        errors['Phase_abs_err'].append(phase_abs_err)
        errors['Phase_rel_err'].append(phase_rel_err)
        errors['PhaseQuad_abs_err'].append(phase_quad_abs_err)
        errors['PhaseQuad_rel_err'].append(phase_quad_rel_err)
        errors['Amplitude_abs_err'].append(amplitude_abs_err)
        errors['Amplitude_rel_err'].append(amplitude_rel_err)        
        
    return errors

def check_sim_results(filename_initial, filename_final, filename_input, tol_kin, tol_phase, tol_amplitude):
    # read the initial and final dataframes
    initial_df = pd.read_csv(filename_initial,skipinitialspace=True)
    final_df = pd.read_csv(filename_final,skipinitialspace=True)

    # separate into a ground state and excited state dataframes
    final_df_ground = final_df[final_df['State'] == 0]
    final_df_excited = final_df[final_df['State'] == 1]

    # get the kz value
    kz = get_kz_value(filename_input)
    # get the omega value
    omega = get_omega_value(filename_input)
    # get the rabifreq value
    rabifreq = get_rabifreq_value(filename_input)

    # check the results
    errors_ground = check_results(initial_df, final_df_ground, kz, omega, rabifreq)
    errors_excited = check_results(initial_df, final_df_excited, kz, omega, rabifreq)

    # convert erros to floats
    for key in errors_ground.keys():
        for i in range(len(errors_ground[key])):
            errors_ground[key][i] = float(errors_ground[key][i])
            errors_excited[key][i] = float(errors_excited[key][i])
    
    # convert to pd dataframe
    errors_ground = pd.DataFrame(errors_ground)
    errors_excited = pd.DataFrame(errors_excited)

    # printout the mean and max absolute and relative errors for each variable (ground then excited)
    print('GROUND STATE ERRORS')
    print('mean error')
    print(errors_ground.mean())
    print('max error')
    print(errors_ground.max())
    print('EXCITED STATE ERRORS')
    print('mean error')
    print(errors_excited.mean())
    print('max error')
    print(errors_excited.max())      

    # assert that all the errors are below the tolerance
    for key in errors_ground.keys():
        # assert that the absolute errors are below the tolerance
        if 'abs_err' in key:
            if ('X' in key) or ('Y' in key) or ('Z' in key):
                assert all(errors_ground[key] < tol_kin)
                assert all(errors_excited[key] < tol_kin)
            elif 'Phase' in key:
                assert all(errors_ground[key] < tol_phase)
                assert all(errors_excited[key] < tol_phase)
            elif 'Amplitude' in key:
                assert all(errors_ground[key] < tol_amplitude)
                assert all(errors_excited[key] < tol_amplitude)
            else:
                raise ValueError('Unknown key')

def main():
    print("TEST 2: SINGLE PULSE PROPAGATION")
    tol_kin = 1e-9
    tol_phase = 1e-6
    tol_amplitude = 1e-6
    check_sim_results('output-files/initial_WPK.txt', 'output-files/final_WPK.txt',  'input-files/input_final.aisi', tol_kin, tol_phase, tol_amplitude)

    print('All tests passed!')

if __name__ == '__main__':
    main()
