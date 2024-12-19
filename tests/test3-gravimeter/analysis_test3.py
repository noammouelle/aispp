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
                return mp.mpf('1')

def get_rabi_freq(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('rabifreq'):
                rabi_freq_value = line.split()[1]
                return 2*pi*float(rabi_freq_value)
            
def get_frequencychirp(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('frequencychirp'):
                frequency_chirp_value = line.split()[1]
                return float(frequency_chirp_value)

def get_phase_analytical(kz,T,rabifreq,frequencychirp):
    tau = pi/(2*rabifreq)
    T = T + 2*tau
    phi2 = -(kz*g+frequencychirp)*T**2 * (1 - 2*tau/T + 4*tau/(pi*T))
    phi2T = (kz*g+frequencychirp)/rabifreq**2 * (rabifreq*T*mp.cos(rabifreq*T) - mp.sin(rabifreq*T))
    deltaT = - kz*g*T - frequencychirp*T
    theta = tau * deltaT / 2
    dphi = 4*theta**2*mp.sin(2*phi2T)
    print('phi2: ', phi2)
    print('dphi: ', dphi)
    return phi2 + dphi

def get_phase_simulated(final_df):
    # we are only interested in the interfering path ending up in the ground state 
    # (paths 0010 and 0100)
    phase_upper = final_df[final_df['Path'] == '0100']['Phase'].values[0]
    phase_upper_quad = final_df[final_df['Path'] == '0100']['PhaseQuad'].values[0]
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

    # compute the phase
    phase = phase_upper + phase_upper_quad - phase_lower - phase_lower_quad - pbar * delta_z / hbar 

    # get the phase err
    err_upper = final_df[final_df['Path'] == '0100']['PhaseErr'].values[0]
    err_lower = final_df[final_df['Path'] == '0010']['PhaseErr'].values[0]
    err = np.sqrt(err_upper**2 + err_lower**2)

    return phase, err

def check_sim_result(input_filepath, initial_filepath, final_filepath, tol):
    # get relevant params
    kz = get_kz_value(input_filepath)
    T = get_T(input_filepath)
    rabifreq = get_rabi_freq(input_filepath)
    frequencychirp = get_frequencychirp(input_filepath)
    # load initial and final data
    initial_data = pd.read_csv(initial_filepath, skipinitialspace=True, dtype={'Path':str})
    final_data = pd.read_csv(final_filepath, skipinitialspace=True, dtype={'Path':str})
    # get z0 and vz0
    z0 = initial_data['Z'][0]
    vz0 = initial_data['VZ'][0]
    # get the analytical phase
    phase_analytical = get_phase_analytical(kz,T,rabifreq,frequencychirp)

    # get the simulated phase
    phase_simulated, phase_err = get_phase_simulated(final_data)

    # print both phases
    print('Analytical phase: ', phase_analytical)
    print('Simulated phase: ', phase_simulated, ' +/- ', phase_err)
    print('Difference: ', np.abs(phase_analytical - phase_simulated))

    # assert within tolerance (tol)
    assert np.abs(phase_analytical - phase_simulated) < tol

    print('Test passed!')

if __name__ == '__main__':
    initial_filepath = 'output-files/data_initial_WPK.txt'
    final_filepath = 'output-files/data_final_WPK.txt'
    input_filepath = 'input-files/input_final.aisi'
    tol = 1e-4
    check_sim_result(input_filepath, initial_filepath, final_filepath, tol)