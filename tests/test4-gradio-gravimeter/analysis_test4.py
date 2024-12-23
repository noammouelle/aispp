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

def get_phase_analytical(kz,T,rabifreq,frequencychirp,z0,vz0):
    tau = pi/(2*rabifreq)
    T = T + 2*tau
    eta = tau/T
    vr = hbar * kz / m
    vm = vz0 + vr/2
    gamma = 2*g/R
    # leading order phase
    phi2 = -(kz*g+frequencychirp-kz*gamma*z0)*T**2 * (1 - 2*eta + 4*eta/pi) \
           +kz*gamma*T**3 * (vm*(1 - 2*eta + 4*eta/pi) -g*T*(7/12 - 4/3*eta + 8/(3*pi)*eta))
    # Bertoldi corrections
    phi2T_term1 = -kz*vz0*gamma/rabifreq**3
    phi2T_term2 = kz*vz0/rabifreq
    phi2T_term3 = -g*kz*T*gamma*mp.cos(rabifreq*T)/rabifreq**3
    phi2T_term4 = kz*vz0*gamma*mp.cos(rabifreq*T)/rabifreq**3
    phi2T_term5 = g*kz*T*mp.cos(rabifreq*T)/rabifreq
    phi2T_term6 = -kz*vz0*mp.cos(rabifreq*T)/rabifreq
    phi2T_term7 = T*frequencychirp*mp.cos(rabifreq*T)/rabifreq
    phi2T_term8 = g*kz*T**3*gamma*mp.cos(rabifreq*T)/(6*rabifreq)
    phi2T_term9 = -kz*vz0*T**2*gamma*mp.cos(rabifreq*T)/(2*rabifreq)
    phi2T_term10 = -kz*T*z0*gamma*mp.cos(rabifreq*T)/rabifreq
    phi2T_term11 = g*kz*gamma*mp.sin(rabifreq*T)/rabifreq**4
    phi2T_term12 = -g*kz*mp.sin(rabifreq*T)/rabifreq**2
    phi2T_term13 = -frequencychirp*mp.sin(rabifreq*T)/rabifreq**2
    phi2T_term14 = -g*kz*T**2*gamma*mp.sin(rabifreq*T)/(2*rabifreq**2)
    phi2T_term15 = kz*vz0*T*gamma*mp.sin(rabifreq*T)/rabifreq**2
    phi2T_term16 = kz*z0*gamma*mp.sin(rabifreq*T)/rabifreq**2

    phi2T = phi2T_term1 + phi2T_term2 + phi2T_term3 + phi2T_term4 + phi2T_term5 + phi2T_term6 + phi2T_term7 + phi2T_term8 + phi2T_term9 + phi2T_term10 + phi2T_term11 + phi2T_term12 + phi2T_term13 + phi2T_term14 + phi2T_term15 + phi2T_term16
    '''
    print("term1: ", phi2T_term1)
    print("term2: ", phi2T_term2)
    print("term3: ", phi2T_term3)
    print("term4: ", phi2T_term4)
    print("term5: ", phi2T_term5)
    print("term6: ", phi2T_term6)
    print("term7: ", phi2T_term7)
    print("term8: ", phi2T_term8)
    print("term9: ", phi2T_term9)
    print("term10: ", phi2T_term10)
    print("term11: ", phi2T_term11)
    print("term12: ", phi2T_term12)
    print("term13: ", phi2T_term13)
    print("term14: ", phi2T_term14)
    print("term15: ", phi2T_term15)
    print("term16: ", phi2T_term16)
    '''
    deltaT = - kz*g*T - frequencychirp*T + kz*vz0*(1+gamma*T**2/2) \
             + kz*gamma*T*(z0-T**2*g/6)
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
    phase_analytical = get_phase_analytical(kz,T,rabifreq,frequencychirp,z0,vz0)

    # get the simulated phase
    phase_simulated, phase_err = get_phase_simulated(final_data)

    # print both phases
    print('\n-----------------------------------')
    print('TEST RESULTS:\n')
    print('Analytical phase: ', phase_analytical)
    print('Simulated phase: ', phase_simulated, ' +/- ', phase_err)
    print('Difference: ', np.abs(phase_analytical - phase_simulated))

    # assert within tolerance (tol)
    assert np.abs(phase_analytical - phase_simulated) < tol

    print('\nTEST PASSED!')
    print('-----------------------------------\n')

def main(i):
    initial_filepath = 'output-files/data_initial_{}_WPK.txt'.format(i)
    final_filepath = 'output-files/data_final_{}_WPK.txt'.format(i)
    input_filepath = 'input-files/input_final_{}.aisi'.format(i)
    tol = 5e-4
    check_sim_result(input_filepath, initial_filepath, final_filepath, tol)

if __name__ == '__main__':
    for i in range(10):
        main(i)