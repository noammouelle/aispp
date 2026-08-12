"""
Shared configuration for the rotating-frame end-to-end tests.

Builds .aisi inputs for short Mach-Zehnder sequences with an arbitrary number of
loops, runs ais++ on them, and loads the resulting phases. Kept separate from
the test script so the visualisation script can reuse it.

Environment
-----------
AISPP_BIN   path to the ais++ binary   (default: ../build/ais++ next to this file)
AISPY_PATH  path to the aispy checkout (default: guesses a sibling of the repo)
"""

import os
import subprocess
import sys

import mpmath as mp
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)

AISPP_BIN = os.environ.get('AISPP_BIN', os.path.join(_REPO, 'build', 'ais++'))


def _find_aispy():
    env = os.environ.get('AISPY_PATH')
    if env:
        return env
    # a worktree of aispp usually has a matching worktree of aispy alongside it
    candidates = [
        os.path.join(os.path.dirname(_REPO), 'aispy'),
        os.path.abspath(os.path.join(_REPO, '..', '..', '..', '..', 'aispy',
                                     '.claude', 'worktrees',
                                     os.path.basename(_REPO))),
        os.path.abspath(os.path.join(_REPO, '..', 'aispy')),
    ]
    for c in candidates:
        if os.path.isdir(os.path.join(c, 'aispy')):
            return c
    raise RuntimeError('could not locate aispy; set AISPY_PATH')


AISPY_PATH = _find_aispy()
sys.path.insert(0, AISPY_PATH)

from aispy.utils import AISFlow, pi, hbar, kz   # noqa: E402

# ── physical scale of the test interferometer ─────────────────────────────────
#
# Short and small: a single atom, a 3-pulse MZ (lmt_order = 1), and an
# interrogation time of 0.1 s. Long enough that the Sagnac phase is ~0.1 rad at
# Earth rate, short enough that the whole matrix of runs finishes in seconds.

INTERROGATION_TIME = mp.mpf('0.1')            # s
RABI_FREQ          = 2 * pi * mp.mpf('1e4')   # rad/s
V0Z                = mp.mpf('19.62')          # m/s, launch velocity
V0X                = mp.mpf('0.01')           # m/s, transverse velocity
ZR                 = 450.085                  # m, Rayleigh range

OMEGA_EARTH = 7.292115e-5                     # rad/s

# effective wavevector, rad/m — used for the analytic Coriolis phase
K_EFF = float(kz)


def build_param_dict(utype='linear_pot', rotation=None, loopnumber=1,
                     interrogation_time=None, v0x=None, v0y=0.0,
                     printtrajectory=False,
                     rabi_freq=None, v0z=None, fountain=False,
                     ultrafast=0, qag_abs=1e-13, qag_rel=1e-13):
    """
    Parameter dictionary for one run. Mirrors examples/build_inputs.py.

    ``fountain=True`` picks the launch velocity that brings the atom back to its
    starting height exactly at the end of the sequence (v0z = g * L * T), which
    is what a real long-baseline fountain does.
    """
    T = INTERROGATION_TIME if interrogation_time is None else mp.mpf(str(interrogation_time))
    vx = V0X if v0x is None else mp.mpf(str(v0x))
    rabi = RABI_FREQ if rabi_freq is None else 2 * pi * mp.mpf(str(rabi_freq))
    if fountain:
        # apex at the midpoint of the whole sequence, which lasts 2*L*T
        vz = mp.mpf('9.81') * loopnumber * T
    else:
        vz = V0Z if v0z is None else mp.mpf(str(v0z))

    # a full L-loop sequence lasts 2*L*T; leave a little margin before detection
    detection_time = 2 * loopnumber * T + mp.mpf('0.001')

    # AISFlow mirrors the interrogation-time list into a palindrome to build the
    # diamond chain, so a D-diamond sequence needs ceil(D/2) entries. Passing a
    # single entry for D >= 3 silently builds the wrong number of diamonds while
    # the path list is still built for D, and nothing interferes.
    t_list = [T] * max(1, (loopnumber + 1) // 2)

    potential_params = {'utype': utype}
    if rotation is not None:
        potential_params['rotation'] = [float(r) for r in rotation]

    return {
        'cloud_params': {
            'natoms':       1,
            'initialstate': 0,
            'sigma':        0.0,
            'longtemp':     0,
            'transtemp':    0,
            'x0':           [0.0, 0.0, 0.0],
            'v0':           [float(vx), float(v0y), float(vz)],
        },
        'potential_params': potential_params,
        'sequence_params': {
            't_init':             mp.mpf('0.0'),
            'detectiontime':      detection_time,
            'interrogation_time': t_list,
            'lmt_order':          1,
            'dt_lmt':             0,
            'automaticdetuning':  1,
            'frequencychirp':     0,
            'kchirp':             0,
            'ultranarrow':        True,
            'sequencename':       'MZ',
            'loopnumber':         loopnumber,
        },
        'pulse_params': {
            'rabi_freq':      rabi,
            'wtype':          'gaussian',
            'phi0':           0,
            'kx_psr':         0,
            'ky_psr':         0,
            'ky':             0,
            'waist':          mp.sqrt(2 * ZR / kz),
            'focallength':    0,
            'zupwardlaser':   0,
            'zdownwardlaser': 0,
            'beam_radius':    0.02,
            'baseline':       10,
            'zernike_params': {},
        },
        'simulation_params': {
            'amplitudethreshold': 0,
            'coherencelength':    1e-3,
            'usemcbranching':     0,
            'ignoredetuning':     0,
            'usestaticapprox':    0,
            'seed':               1,
            'usedetvolselection': 0,
            'usepathselection':   1,
            'xdet':               [-1.0, 1.0],
            'ydet':               [-1.0, 1.0],
            'zdet':               [-1.0, 1.0],
            # The action integral grows with the flight time (~1e3 m^2/s^2 for a
            # multi-second fountain), so an absolute tolerance of 1e-13 is below
            # double precision there and GSL bails out with a roundoff error.
            # Long sequences need a looser setting, or ultrafast=1 if only the
            # trajectories matter.
            'gslqagabserr':       qag_abs,
            'gslqagrelerr':       qag_rel,
            'gslkinodeabserr':    1e-12,
            'gslkinoderelerr':    0,
            'gslpulseodeabserr':  1e-12,
            'gslpulseoderelerr':  0,
            # the action phase carries the Sagnac term, so ultrafast must be off
            'ultrafast':          ultrafast,
        },
        'io_params': {
            # the main output holds one randomly *sampled* port per atom, which
            # with a single atom is often a non-interfering one; the _PROB file
            # holds every port and is deterministic
            'printprobs':       1,
            'printwavepackets': 0,
            'printtrajectory':  1 if printtrajectory else 0,
        },
    }


def apply_tilt_compensation(aisi_path, omega):
    """
    Rewrite the kx/ky lines of an .aisi file so that the effective wavevector is
    held fixed in the *inertial* frame rather than in the rotating one.

    This is the simulation counterpart of counter-rotating the retroreflector
    ("tip-tilt compensation"), the standard way real experiments recover the
    contrast a rotating frame would otherwise destroy.

    Why it is needed: a Mach-Zehnder does not close in a rotating frame. The two
    arms differ by the recoil velocity u = hbar k / m along k, and the Coriolis
    acceleration -2 Omega x v acts on that difference, leaving the arms
    transversely separated by 2 |Omega_perp| u T^2 at recombination. That is a
    *differential* effect, so no common-mode launch velocity can undo it -- only
    changing the direction of the momentum kicks can.

    Holding k fixed in inertial space means its rotating-frame components obey
    k(t) = R(-Omega t) k(0), i.e. to first order in Omega t

        k(t) = k(0) - (Omega x k(0)) t.

    Each pulse gets the tilt appropriate to its own firing time, with the sign
    following that pulse's beam direction.
    """
    omega = np.asarray(omega, dtype=float)

    lines = open(aisi_path).read().splitlines()
    values = {}
    for ln in lines:
        parts = ln.split()
        if parts and parts[0] in ('t0', 'kx', 'ky', 'kz'):
            values[parts[0]] = [float(x) for x in parts[1:]]

    for key in ('t0', 'kx', 'ky', 'kz'):
        if key not in values:
            raise RuntimeError(f'{aisi_path}: no "{key}" line to tilt-compensate')

    t0, kx, ky, kz = values['t0'], values['kx'], values['ky'], values['kz']
    if not (len(t0) == len(kx) == len(ky) == len(kz)):
        raise RuntimeError(f'{aisi_path}: pulse arrays have inconsistent lengths')

    new_kx, new_ky = [], []
    for t, x, y, z in zip(t0, kx, ky, kz):
        k = np.array([x, y, z])
        tilted = k - np.cross(omega, k) * t
        # plain Python floats: repr() of a numpy scalar is "np.float64(...)" on
        # numpy >= 2, which the ais++ parser silently fails to read
        new_kx.append(float(tilted[0]))
        new_ky.append(float(tilted[1]))

    out = []
    for ln in lines:
        parts = ln.split()
        if parts and parts[0] == 'kx':
            out.append('kx ' + ' '.join(f'{v:.17g}' for v in new_kx) + ' ')
        elif parts and parts[0] == 'ky':
            out.append('ky ' + ' '.join(f'{v:.17g}' for v in new_ky) + ' ')
        else:
            out.append(ln)
    open(aisi_path, 'w').write('\n'.join(out) + '\n')


def run_case(name, workdir, tilt_compensation=None, **kwargs):
    """
    Build an input, run ais++, and return a dict with the interfering phase
    shift and the paths of the output files.
    """
    import h5py

    os.makedirs(workdir, exist_ok=True)
    AISFlow(build_param_dict(**kwargs), name, workdir)

    aisi = os.path.join(workdir, name + '.aisi')
    out  = os.path.join(workdir, name + '.h5')
    prob = os.path.join(workdir, name + '_PROB.h5')
    traj = os.path.join(workdir, name + '_TRAJ.h5')

    if tilt_compensation is not None:
        apply_tilt_compensation(aisi, tilt_compensation)

    proc = subprocess.run([AISPP_BIN, '-i', aisi, '-o', out],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f'ais++ failed for {name}:\n{proc.stdout}\n{proc.stderr}')

    with h5py.File(prob, 'r') as f:
        interfering = f['interferingFlag'][:].astype(bool)
        phases = f['phaseShifts'][:]
        probabilities = f['probabilities'][:]
        positions = f['positions'][:]
        rotation = f['rotation'][:] if 'rotation' in f else np.zeros(3)

    if not interfering.any():
        raise RuntimeError(f'no interfering port for {name}')

    # take the brightest interfering port; the others carry the same phase
    idx = np.where(interfering)[0]
    best = idx[np.argmax(probabilities[idx])]

    return {
        'name': name,
        'phase': float(phases[best]),
        'positions': positions,
        'interfering': interfering,
        'n_interfering': int(interfering.sum()),
        'rotation': np.asarray(rotation, float),
        'output': out,
        'prob': prob,
        'trajectory': traj if kwargs.get('printtrajectory') else None,
        'stdout': proc.stdout,
    }


def analytic_coriolis_phase(omega_y, v_x, T, loopnumber=1):
    """
    Leading-order rotation phase of an L-loop MZ with k along z and Omega along y.

    For a single loop the standard result is |dphi| = 2 k_eff Omega_y v_x T^2.
    Successive loops are traversed in the opposite sense, so the signed sum over
    L loops alternates: it vanishes for even L and equals the single-loop value
    for odd L.
    """
    single = 2.0 * K_EFF * omega_y * v_x * T * T
    return single if loopnumber % 2 == 1 else 0.0
