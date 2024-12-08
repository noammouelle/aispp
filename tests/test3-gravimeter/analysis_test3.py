import numpy as np
import mpmath as mp
import pandas as pd

# define constants as in ais++
g = mp.mpf('9.81')
R = mp.mpf('6.37e6')
au = mp.mpf('1.660539066e-27')
m = mp.mpf('86.90888') * au
h = mp.mpf('6.62607015e-34')
pi = mp.pi
hbar = h / (2 * pi)

def get_kz_value(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('kz'):
                _, kz_value, _, _ = line.split()
                return float(kz_value)
    return None

def get_T(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('detectiontime'):
                _, T_value = line.split()
                #return float(T_value)/2
                return mp.mpf('1.0')

def get_rabi_freq(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('rabifreq'):
                rabi_freq_value = line.split()[1]
                return 2*pi*float(rabi_freq_value)

def get_phase_analytical(kz,T,rabifreq):
    return -kz*g*T*(T+2*pi/(2*rabifreq))

def get_phase_simulated(final_df):
    # we are only interested in the interfering path ending up in the ground state 
    # (paths 0010 and 0100)
    phase_upper = final_df[final_df['Path'] == '0100']['Phase'].values[0]
    phase_upper_quad = final_df[final_df['Path'] == '0010']['PhaseQuad'].values[0]
    phase_lower = final_df[final_df['Path'] == '0010']['Phase'].values[0]
    phase_lower_quad = final_df[final_df['Path'] == '0010']['PhaseQuad'].values[0]

    # get the final positions and velocities
    z_upper = final_df[final_df['Path'] == '0100']['Z'].values[0]
    vz_upper = final_df[final_df['Path'] == '0100']['VZ'].values[0]
    z_lower = final_df[final_df['Path'] == '0010']['Z'].values[0]
    vz_lower = final_df[final_df['Path'] == '0010']['VZ'].values[0]

    # compute average momenta
    pbar = 0.5 * (m * vz_upper + m * vz_lower)
    # compute delta z
    delta_z = z_upper - z_lower

    return phase_upper + phase_upper_quad - phase_lower - phase_lower_quad - pbar * delta_z / hbar

def check_sim_result(input_filepath, initial_filepath, final_filepath, tol):
    # get relevant params
    kz = get_kz_value(input_filepath)
    T = get_T(input_filepath)
    rabifreq = get_rabi_freq(input_filepath)
    # load initial and final data
    initial_data = pd.read_csv(initial_filepath, skipinitialspace=True, dtype={'Path':str})
    final_data = pd.read_csv(final_filepath, skipinitialspace=True, dtype={'Path':str})
    # get z0 and vz0
    z0 = initial_data['Z'][0]
    vz0 = initial_data['VZ'][0]
    # get the analytical phase
    phase_analytical = get_phase_analytical(kz,T,rabifreq)

    # get the simulated phase
    phase_simulated = get_phase_simulated(final_data)

    # print both phases
    print('Analytical phase: ', phase_analytical)
    print('Simulated phase: ', phase_simulated)

    # assert within tolerance (tol)
    assert np.abs(phase_analytical - phase_simulated) < tol

    print('Test passed!')

if __name__ == '__main__':
    initial_filepath = 'output-files/data_initial_WPK.txt'
    final_filepath = 'output-files/data_final_WPK.txt'
    input_filepath = 'input-files/input_final.aisi'
    tol = 1e-1
    check_sim_result(input_filepath, initial_filepath, final_filepath, tol)