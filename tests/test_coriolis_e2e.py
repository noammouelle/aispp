#!/usr/bin/env python3
"""
End-to-end tests for rotating-frame (Coriolis) support in ais++.

Runs real Mach-Zehnder sequences through the ais++ binary and checks the
resulting interferometric phase against independent predictions.

    python tests/test_coriolis_e2e.py            # full suite
    python tests/test_coriolis_e2e.py --quick    # skip the slower scans

Environment: AISPP_BIN, AISPY_PATH (see coriolis_common.py).

What is checked
---------------
1. Null test          omega = 0 reproduces the inertial potentials exactly, and
                      asking for rotation with an inertial utype is an error.
2. Sagnac area law    the odd-in-omega phase equals 2 m Omega.A / hbar, with A
                      the area enclosed by the *non-rotating* arms.
3. Linearity          the odd part scales as Omega, the even part as Omega^2.
4. T^2 scaling        the odd part scales as the square of the interrogation time.
5. Centrifugal        the even part has the size expected of k_eff Omega^2 z T^2.
6. Multi-loop         successive loops reverse the circulation, so the enclosed
                      area and the rotation phase collapse for even loop counts;
                      the residual is shown to be a finite-pulse-duration effect.
7. Closure            the arms recombine at the same point.

Phases are decomposed as
    odd  = (phi(+Omega) - phi(-Omega)) / 2     -> exactly first order in Omega
    even = (phi(+Omega) + phi(-Omega)) / 2 - phi(0)  -> exactly second and higher
which isolates the Sagnac and centrifugal contributions without any fitting.
"""

import argparse
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coriolis_common as cc                                   # noqa: E402

sys.path.insert(0, cc.AISPY_PATH)
from aispy.trajectory import load_trajectory, enclosed_area, sagnac_phase  # noqa: E402

OMEGA = cc.OMEGA_EARTH
WORKDIR = os.environ.get('CORIOLIS_WORKDIR',
                         os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      '_coriolis_e2e'))

_checks = 0
_failures = 0


def check(ok, what, detail=''):
    global _checks, _failures
    _checks += 1
    if ok:
        print(f'  ok    {what}' + (f'   [{detail}]' if detail else ''))
    else:
        _failures += 1
        print(f'  FAIL  {what}' + (f'   [{detail}]' if detail else ''))


def check_close(got, want, rtol, what):
    err = abs(got - want) / max(abs(want), 1e-300)
    check(err <= rtol, what, f'got {got:.8g}, want {want:.8g}, rel {err:.3g} (tol {rtol:g})')


def section(name):
    print(f'\n{name}')


_cache = {}


def phase(tag, **kw):
    """Run a case (memoised) and return its interfering phase."""
    if tag not in _cache:
        _cache[tag] = cc.run_case(tag, WORKDIR, **kw)
    return _cache[tag]['phase']


def run(tag, **kw):
    if tag not in _cache:
        _cache[tag] = cc.run_case(tag, WORKDIR, **kw)
    return _cache[tag]


def odd_even(tag, omega=OMEGA, **kw):
    """Return (odd, even) parts of the phase in Omega, about Omega = 0."""
    # trajectories are always written: they are small for a single atom, and it
    # keeps the memoised runs usable by the area-law and closure checks
    p0 = phase(tag + '_z', utype='rotating_linear_pot', rotation=[0, 0, 0],
               printtrajectory=True, **kw)
    pp = phase(tag + '_p', utype='rotating_linear_pot', rotation=[0,  omega, 0],
               printtrajectory=True, **kw)
    pm = phase(tag + '_m', utype='rotating_linear_pot', rotation=[0, -omega, 0], **kw)
    return (pp - pm) / 2.0, (pp + pm) / 2.0 - p0


# ── 1. null test ─────────────────────────────────────────────────────────────

def test_null():
    section('1. Null test: omega = 0 must reproduce the inertial potentials')

    inertial = phase('null_inertial', utype='linear_pot')
    rot0 = phase('null_rot0', utype='rotating_linear_pot', rotation=[0, 0, 0])
    check(inertial == rot0, 'rotating_linear_pot at omega=0 == linear_pot',
          f'{inertial!r} vs {rot0!r}')

    # omitting the rotation key entirely must behave the same
    rot_absent = phase('null_rotabs', utype='rotating_linear_pot')
    check(inertial == rot_absent, 'rotating_linear_pot with no rotation key == linear_pot')

    # asking for rotation with an inertial potential must be rejected, not ignored
    try:
        phase('null_bad', utype='linear_pot', rotation=[0, OMEGA, 0])
        check(False, 'rotation with an inertial utype is rejected')
    except RuntimeError:
        check(True, 'rotation with an inertial utype is rejected')


# ── 2. Sagnac area law ───────────────────────────────────────────────────────

def test_area_law():
    section('2. Sagnac area law for a single loop')

    odd, even = odd_even('L1', loopnumber=1)
    ref = run('L1_z', utype='rotating_linear_pot', rotation=[0, 0, 0],
              loopnumber=1)

    traj = load_trajectory(ref['trajectory'])
    area = enclosed_area(traj, potential='rotating_linear_pot')
    law = sagnac_phase(area, [0, OMEGA, 0])

    print(f'    enclosed area A_y = {area[1]:.6e} m^2')
    print(f'    odd-in-omega phase = {odd:+.8e} rad')
    print(f'    2 m Omega.A / hbar = {law:+.8e} rad')

    # the overall sign depends on which arm the port frame calls "upper"; the
    # magnitude is the physics
    check_close(abs(odd), abs(law), 1e-3,
                'odd-in-omega phase equals the Sagnac area law')

    # and both agree with the textbook single-loop formula 2 k Omega v_x T^2
    analytic = cc.analytic_coriolis_phase(OMEGA, float(cc.V0X),
                                          float(cc.INTERROGATION_TIME))
    check_close(abs(odd), analytic, 3e-3,
                'odd-in-omega phase equals 2 k_eff Omega v_x T^2')


# ── 3. scaling in omega ──────────────────────────────────────────────────────

def test_omega_scaling():
    section('3. Scaling with the rotation rate')

    odd1, even1 = odd_even('L1', loopnumber=1)
    odd2, even2 = odd_even('L1x2', omega=2 * OMEGA, loopnumber=1)

    check_close(odd2 / odd1, 2.0, 1e-4, 'odd part is linear in omega')

    # The even part is a ~3e-3 rad residual extracted from a total phase of
    # ~9e5 rad, so it carries about a percent of numerical scatter: sweeping the
    # rate over 0.5x-4x Earth gives even/Omega^2 constant to +-1%, but
    # non-monotonically, which is noise rather than a higher-order term. The
    # tolerance below reflects that floor; the check still separates a quadratic
    # (ratio 4) from a linear (ratio 2) dependence unambiguously.
    check_close(even2 / even1, 4.0, 5e-2, 'even part is quadratic in omega')


# ── 4. scaling in T ──────────────────────────────────────────────────────────

def test_time_scaling(quick):
    section('4. Scaling with the interrogation time')
    if quick:
        print('    skipped (--quick)')
        return

    T = float(cc.INTERROGATION_TIME)
    odd1, _ = odd_even('L1', loopnumber=1)
    odd2, _ = odd_even('L1T2', loopnumber=1, interrogation_time=2 * T)

    # finite pulse durations make the effective T slightly longer than the
    # nominal one, which biases the ratio by ~0.1%
    check_close(odd2 / odd1, 4.0, 5e-3, 'odd part scales as T^2')


# ── 5. centrifugal term ──────────────────────────────────────────────────────

def test_centrifugal():
    section('5. Even-in-omega part is the centrifugal acceleration phase')

    _, even = odd_even('L1', loopnumber=1)
    ref = run('L1_z', utype='rotating_linear_pot', rotation=[0, 0, 0],
              loopnumber=1)
    z_mean = float(np.mean(load_trajectory(ref['trajectory'])['positions'][:, 2]))
    T = float(cc.INTERROGATION_TIME)

    # a centrifugal acceleration Omega^2 z acts like a small extra gravity, so
    # the phase is of order k_eff Omega^2 z T^2. This is an order-of-magnitude
    # anchor, not an identity: z varies over the flight and there is an Omega^2 x
    # term too.
    estimate = cc.K_EFF * OMEGA**2 * z_mean * T * T
    print(f'    even part  = {even:+.4e} rad')
    print(f'    k Omega^2 z T^2 ~ {estimate:.4e} rad  (z_mean = {z_mean:.3f} m)')
    check(0.2 < abs(even) / estimate < 5.0,
          'even part matches the centrifugal scale',
          f'ratio {abs(even)/estimate:.2f}')


# ── 6. multi-loop cancellation ───────────────────────────────────────────────

def test_multiloop(quick):
    section('6. Multi-loop: Coriolis cancels for even loop counts')

    rows = []
    for L in (1, 2, 4):
        odd, even = odd_even(f'ML{L}', loopnumber=L)
        ref = run(f'ML{L}_z', utype='rotating_linear_pot', rotation=[0, 0, 0],
                  loopnumber=L)
        traj = load_trajectory(ref['trajectory'])
        area = enclosed_area(traj, potential='rotating_linear_pot')
        rows.append({'L': L, 'odd': odd, 'even': even, 'A_y': area[1],
                     'traj': traj})
        print(f'    L={L}:  A_y = {area[1]:+.4e} m^2   odd = {odd:+.4e} rad   '
              f'even = {even:+.4e} rad')

    a1 = abs(rows[0]['A_y'])
    for r in rows[1:]:
        check(abs(r['A_y']) / a1 < 1e-4,
              f'enclosed area cancels for L={r["L"]}',
              f'|A_y|/|A_y(L=1)| = {abs(r["A_y"])/a1:.2e}')

    o1 = abs(rows[0]['odd'])
    for r in rows[1:]:
        check(abs(r['odd']) / o1 < 1e-2,
              f'rotation phase suppressed for L={r["L"]}',
              f'|odd|/|odd(L=1)| = {abs(r["odd"])/o1:.2e}')

    # The residual for even L is not physics: it is the finite pulse duration
    # breaking the symmetry between successive loops. Shortening the pulses
    # by 10x must shrink it while leaving the single-loop phase alone.
    if quick:
        print('    (pulse-duration refinement skipped: --quick)')
        return

    odd_L2_fast, _ = odd_even('MLfast2', loopnumber=2, rabi_freq=1e5)
    odd_L1_fast, _ = odd_even('MLfast1', loopnumber=1, rabi_freq=1e5)
    print(f'    10x shorter pulses:  odd(L=1) = {odd_L1_fast:+.4e}   '
          f'odd(L=2) = {odd_L2_fast:+.4e}')

    check(abs(odd_L2_fast) < 0.5 * abs(rows[1]['odd']),
          'even-L residual shrinks with shorter pulses (finite-pulse artefact)',
          f'{abs(rows[1]["odd"]):.3e} -> {abs(odd_L2_fast):.3e}')
    check_close(abs(odd_L1_fast), abs(rows[0]['odd']), 1e-2,
                'single-loop phase is unaffected by pulse duration')


# ── 7. path closure ──────────────────────────────────────────────────────────

def test_closure():
    section('7. Interferometer closure')

    from aispy.trajectory import arm_loop

    def gap(tag, **kw):
        ref = run(tag, utype='rotating_linear_pot', **kw)
        loop = arm_loop(load_trajectory(ref['trajectory']),
                        potential='rotating_linear_pot')
        return float(np.linalg.norm(loop['upper'][-1] - loop['lower'][-1]))

    gaps = {}
    for L in (1, 2, 4):
        gaps[L] = gap(f'ML{L}_p', rotation=[0, OMEGA, 0], loopnumber=L)
        print(f'    L={L}: residual arm separation |dr| = {gaps[L]:.3e} m')

    # A single loop does NOT close in a rotating frame: the Coriolis force acts
    # differently on the two arms because their vertical velocities differ by the
    # recoil, and nothing undoes it. This is the effect that real gyro-compensated
    # gravimeters cancel with a tip-tilt mirror. The residual should be of order
    # Omega (hbar k / m) T^2.
    recoil_v = 1.0546e-34 * cc.K_EFF / (86.90888 * 1.660539066e-27)
    T = float(cc.INTERROGATION_TIME)
    scale = OMEGA * recoil_v * T * T
    print(f'    Omega (hbar k / m) T^2 = {scale:.3e} m')
    check(0.3 < gaps[1] / scale < 10.0,
          'single-loop residual has the Coriolis-opening scale',
          f'ratio {gaps[1]/scale:.2f}')

    # ...and it must grow linearly with Omega, confirming it is the Coriolis
    # opening rather than a sequence-timing artefact
    g1 = gap('CL1x1', rotation=[0, OMEGA, 0], loopnumber=1, printtrajectory=True)
    g2 = gap('CL1x2', rotation=[0, 2 * OMEGA, 0], loopnumber=1, printtrajectory=True)
    check_close(g2 / g1, 2.0, 1e-2, 'single-loop residual is linear in omega')

    # Even loop counts reverse the deflection and do close.
    for L in (2, 4):
        check(gaps[L] < 0.1 * gaps[1],
              f'even loop count closes far better (L={L})',
              f'{gaps[L]:.3e} m vs {gaps[1]:.3e} m at L=1')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true',
                    help='skip the slower scans')
    args = ap.parse_args()

    print('ais++ rotating-frame end-to-end tests')
    print(f'  binary : {cc.AISPP_BIN}')
    print(f'  aispy  : {cc.AISPY_PATH}')
    print(f'  workdir: {WORKDIR}')

    t0 = time.time()
    test_null()
    test_area_law()
    test_omega_scaling()
    test_time_scaling(args.quick)
    test_centrifugal()
    test_multiloop(args.quick)
    test_closure()

    print(f'\n{_checks} checks, {_failures} failures  '
          f'({len(_cache)} simulations, {time.time()-t0:.1f} s)')
    if _failures == 0:
        print('ALL TESTS PASSED')
    return 1 if _failures else 0


if __name__ == '__main__':
    sys.exit(main())
