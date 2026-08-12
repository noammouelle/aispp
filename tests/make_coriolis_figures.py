#!/usr/bin/env python3
"""
Generate the rotating-frame figures.

    python tests/make_coriolis_figures.py [--outdir examples/coriolis]
    python tests/make_coriolis_figures.py --only strategy

Everything is drawn at TRUE SCALE. With a single photon recoil the two arms end
up ~1e6 times closer together than the trajectory is long, which is unplottable;
the honest fix is LMT rather than a drawing trick, so these examples use
lmt_order = 101, giving an arm separation of n * (hbar k/m) * T = tens of cm.

A stronger pulse (rabi_freq = 100 kHz) is used with LMT because the sequence
builder's LMT blocks leave an uncorrected vertical separation of order
(hbar k/m) * dt_pi * n(n-1)/2, which at n = 101 and a 10 kHz Rabi frequency is
3.3 mm -- larger than the coherence length, so nothing interferes. This is a
property of the pulse sequence, not of the rotating frame: it is identical at
Omega = 0. It scales with the pi-pulse duration, so a 100 kHz Rabi frequency
brings it to 0.34 mm and the interferometer closes.

Geometry: fountain sequences, launched so the atom returns to its starting
height at recombination. Loop counts are compared at FIXED total free-fall time
(T = T_tot / 2L), which is the real design question, rather than at fixed T.

Figures
-------
coriolis_effect.png          what rotation does to the trajectories
coriolis_fountain.png        uniform gravity + Coriolis over a full fountain
coriolis_separation.png      arm separations vs time for 1, 2 and 4 loops
coriolis_planes.png          arm traces on the x-y, x-z and y-z planes
coriolis_loops.png           1/2/4-loop enclosed areas at fixed total time
coriolis_closure.png         launch velocity vs tilting k, for closure
coriolis_longbaseline.png    the long-baseline strategy: why tip-tilt does not
                             scale, and what to do instead

Rotation is Earth's, at latitude 45 deg, with z vertical and x east.
"""

import argparse
import os
import sys
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coriolis_common as cc                                     # noqa: E402

sys.path.insert(0, cc.AISPY_PATH)
from aispy.trajectory import (load_trajectory,                    # noqa: E402
                              plot_trajectory_planes, plot_loop_comparison,
                              reconstruct_trajectories)

LAT = np.deg2rad(45.0)
OMEGA_VEC = cc.OMEGA_EARTH * np.array([0.0, np.cos(LAT), np.sin(LAT)])
OMEGA_PERP = OMEGA_VEC[1]
POT = 'rotating_linear_pot'

NLMT = 101          # LMT order: sets the arm separation to n*(hbar k/m)*T
RABI = 1e5          # Hz; short pi pulses keep the LMT blocks closing
T_TOT = 2.5         # s of total free fall for the main figures

G = 9.81
HBAR = 1.054571817e-34
MASS = 86.90888 * 1.660539066e-27
RECOIL_V = HBAR * cc.K_EFF / MASS
KB = 1.380649e-23
COH_LEN = HBAR / (MASS * np.sqrt(KB * 1e-9 / MASS))   # 1 nK cloud

BASE = dict(fountain=True, ultrafast=1, lmt_order=NLMT, rabi_freq=RABI)


def run(tag, workdir, loops=1, total_time=T_TOT, rotation=None, **kw):
    return cc.run_case(tag, workdir, utype=POT,
                       rotation=OMEGA_VEC if rotation is None else rotation,
                       loopnumber=loops, interrogation_time=total_time / (2 * loops),
                       printtrajectory=True, **BASE, **kw)


def arms(r):
    sm = reconstruct_trajectories(load_trajectory(r['trajectory']),
                                  potential=POT, arm_grouping=True)
    return sm['upper'], sm['lower']


def grid(up, lo, n=3000):
    t = np.linspace(max(up['t'][0], lo['t'][0]), min(up['t'][-1], lo['t'][-1]), n)
    return t, (lambda d, k: np.interp(t, d['t'], d[k]))


def save(fig, outdir, name):
    p = os.path.join(outdir, name)
    fig.savefig(p, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('  wrote', p)
    return p


# ── what rotation does ───────────────────────────────────────────────────────

def fig_effect(runs, outdir):
    up_r, lo_r = arms(runs[('rot', 1)])
    up_i, lo_i = arms(runs[('inert', 1)])
    t, ip = grid(up_r, lo_r)
    com = lambda u, l, k: 0.5 * (ip(u, k) + ip(l, k))

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))

    axs[0].plot(t, (ip(up_r, 'z') - ip(lo_r, 'z')) * 100, color='#1f77b4', lw=1.4)
    axs[0].set_ylabel('$\\Delta z$ between arms [cm]')
    axs[0].set_title('(a) the interferometer diamond\n(unchanged by rotation)',
                     fontsize=9)

    axs[1].plot(t, (com(up_r, lo_r, 'x') - com(up_i, lo_i, 'x')) * 1e3,
                color='#d62728', lw=1.4, label='$x$ (east)')
    axs[1].plot(t, (com(up_r, lo_r, 'y') - com(up_i, lo_i, 'y')) * 1e3,
                color='#2ca02c', lw=1.4, label='$y$ (north)')
    axs[1].set_ylabel('rotating $-$ inertial [mm]')
    axs[1].set_title('(b) common-mode deflection of the cloud\n'
                     '(a launch velocity can cancel this)', fontsize=9)
    axs[1].legend(fontsize=8)

    axs[2].plot(t, (ip(up_r, 'x') - ip(lo_r, 'x')) * 1e6, color='#d62728',
                lw=1.4, label='rotating')
    axs[2].plot(t, (ip(up_i, 'x') - ip(lo_i, 'x')) * 1e6, color='#888888',
                lw=1.4, ls='--', label='inertial')
    axs[2].set_ylabel('$\\Delta x$ between arms [$\\mu$m]')
    axs[2].set_title('(c) differential deflection: encloses area,\n'
                     'and leaves the arms open (a launch velocity cannot)',
                     fontsize=9)
    axs[2].legend(fontsize=8)

    for ax in axs:
        ax.set_xlabel('$t$ [s]')
        ax.axhline(0, color='k', lw=0.6, alpha=0.4)
        ax.grid(alpha=0.25, lw=0.5)
    fig.suptitle(f'Single-loop LMT-{NLMT} fountain, {T_TOT} s free fall, at Earth '
                 'rotation rate. All quantities are true scale.', fontsize=11)
    fig.tight_layout()
    return save(fig, outdir, 'coriolis_effect.png')


# ── the fountain itself ──────────────────────────────────────────────────────

def fig_fountain(runs, outdir):
    up_r, lo_r = arms(runs[('rot', 1)])
    up_i, lo_i = arms(runs[('inert', 1)])
    t, ip = grid(up_r, lo_r)
    z_com = 0.5 * (ip(up_r, 'z') + ip(lo_r, 'z'))
    dx = 0.5 * (ip(up_r, 'x') + ip(lo_r, 'x')) - 0.5 * (ip(up_i, 'x') + ip(lo_i, 'x'))
    dy = 0.5 * (ip(up_r, 'y') + ip(lo_r, 'y')) - 0.5 * (ip(up_i, 'y') + ip(lo_i, 'y'))

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.3))
    axs[0].plot(t, z_com, color='#1f77b4', lw=1.5)
    axs[0].fill_between(t, ip(lo_r, 'z'), ip(up_r, 'z'), color='#1f77b4',
                        alpha=0.25, lw=0)
    axs[0].set_ylabel('$z$ [m]')
    axs[0].set_title(f'(a) the fountain, arms shaded\napex {z_com.max():.2f} m, '
                     f'free fall {t[-1]-t[0]:.2f} s', fontsize=9)

    axs[1].plot(t, dx * 1e3, color='#d62728', lw=1.5, label='$x$ (east)')
    axs[1].plot(t, dy * 1e3, color='#2ca02c', lw=1.5, label='$y$ (north)')
    axs[1].set_ylabel('rotating $-$ inertial [mm]')
    axs[1].set_title('(b) Coriolis deflection of the cloud\n'
                     '$\\ddot{r} = -2\\,\\Omega \\times \\dot{r}$', fontsize=9)
    axs[1].legend(fontsize=8)

    axs[2].plot(t, (ip(up_r, 'z') - ip(lo_r, 'z')) * 100, color='#1f77b4', lw=1.5)
    axs[2].set_ylabel('$\\Delta z$ between arms [cm]')
    axs[2].set_title(f'(c) arm separation from {NLMT} photon recoils\n'
                     f'peak {np.abs(ip(up_r,"z")-ip(lo_r,"z")).max()*100:.1f} cm',
                     fontsize=9)

    apex = t[np.argmax(z_com)]
    for ax in axs:
        ax.axvline(apex, color='k', lw=0.7, ls=':', alpha=0.6)
        ax.set_xlabel('$t$ [s]')
        ax.grid(alpha=0.25, lw=0.5)
    fig.suptitle(f'LMT-{NLMT} fountain with uniform gravity and Coriolis: '
                 f'{t[-1]-t[0]:.1f} s free fall, {z_com.max():.1f} m apex. The '
                 'eastward deflection grows fastest at the apex and then levels '
                 'off, because $\\dot z$ reverses there.', fontsize=11)
    fig.tight_layout()
    return save(fig, outdir, 'coriolis_fountain.png')


# ── arm separations vs loop count ────────────────────────────────────────────

def fig_separation(runs, outdir, loops=(1, 2, 4)):
    fig, axs = plt.subplots(2, len(loops), figsize=(13, 6), sharex='col')
    for col, L in enumerate(loops):
        up, lo = arms(runs[('rot', L)])
        t, ip = grid(up, lo)
        d = {k: ip(up, k) - ip(lo, k) for k in 'xyz'}
        axs[0, col].plot(t, d['z'] * 100, color='#1f77b4', lw=1.3)
        axs[0, col].set_title(f'{L} loop' + ('s' if L > 1 else '') +
                              f'  ($T$ = {T_TOT/(2*L):.3f} s)', fontsize=10)
        axs[0, col].set_ylabel('$\\Delta z$ [cm]' if col == 0 else '')
        axs[1, col].plot(t, d['x'] * 1e6, color='#d62728', lw=1.3, label='$\\Delta x$')
        axs[1, col].plot(t, d['y'] * 1e6, color='#2ca02c', lw=1.3, label='$\\Delta y$')
        axs[1, col].set_ylabel('transverse [$\\mu$m]' if col == 0 else '')
        axs[1, col].set_xlabel('$t$ [s]')
        for ax in (axs[0, col], axs[1, col]):
            ax.axhline(0, color='k', lw=0.6, alpha=0.4)
            ax.grid(alpha=0.25, lw=0.5)
    axs[1, 0].legend(fontsize=8)
    fig.suptitle(f'Arm separation vs time, LMT-{NLMT}, fixed {T_TOT} s total free '
                 'fall. $\\Delta z$ reverses sign once per loop, so both the '
                 'enclosed area and the differential Coriolis kick cancel for an '
                 'even number of loops.', fontsize=11)
    fig.tight_layout()
    return save(fig, outdir, 'coriolis_separation.png')


# ── plane traces and loop comparison ─────────────────────────────────────────

def fig_planes(runs, outdir):
    traj = load_trajectory(runs[('rot', 1)]['trajectory'])
    fig, axes, info = plot_trajectory_planes(traj, potential=POT,
                                             planes=('xy', 'xz', 'yz'))
    fig.suptitle(fig._suptitle.get_text() +
                 f'\nLMT-{NLMT} fountain, {T_TOT} s free fall — true scale',
                 fontsize=10)
    fig.tight_layout()
    print(f'      A = {info["area"]} m^2, dphi = {info["sagnac_phase"]:.4e} rad, '
          f'gap = {info["gap"]:.3e} m')
    return save(fig, outdir, 'coriolis_planes.png')


def fig_loops(runs, outdir, loops=(1, 2, 4)):
    entries = [(f'{L} loop' + ('s' if L > 1 else '') + f' — $T$ = {T_TOT/(2*L):.3f} s',
                load_trajectory(runs[('inert', L)]['trajectory'])) for L in loops]
    fig, axes, rows = plot_loop_comparison(entries, potential=POT, plane='xz',
                                           sagnac_rotation=OMEGA_VEC)
    fig.suptitle('Enclosed area vs number of loops at fixed total free-fall time '
                 f'({T_TOT} s, LMT-{NLMT}, true scale) — successive loops reverse '
                 'the circulation, so even loop counts cancel', fontsize=10)
    fig.tight_layout()
    for r in rows:
        print(f"      {r['label'][:14]:16s} A = {r['area'][1]:+.4e} m^2   "
              f"dphi = {r['sagnac_phase']:+.4e} rad")
    return save(fig, outdir, 'coriolis_loops.png')


# ── what closes the interferometer ───────────────────────────────────────────

def fig_closure(outdir, workdir):
    predicted = 2 * OMEGA_PERP * NLMT * RECOIL_V * (T_TOT / 2) ** 2

    def dx(tag, **kw):
        up, lo = arms(run(tag, workdir, loops=1, **kw))
        t, ip = grid(up, lo)
        return t, ip(up, 'x') - ip(lo, 'x')

    v0x_list = [0.0, 0.001, 0.005, 0.02]
    tilt_scan = [0.0, 0.5, 0.9, 1.0, 1.1, 1.5, 2.0]
    vel_runs = [(v, dx(f'cl_vx{v}', v0x=v)) for v in v0x_list]
    tilt_runs = [(f, dx(f'cl_t{f}', v0x=0.001, tilt_compensation=OMEGA_VEC * f))
                 for f in tilt_scan]

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.4))
    cmap = plt.get_cmap('viridis')

    for i, (v, (t, d)) in enumerate(vel_runs):
        axs[0].plot(t, d * 1e6, lw=1.4, color=cmap(i / max(1, len(vel_runs) - 1)),
                    label=f'$v_{{0x}}$ = {v*1e3:g} mm/s')
    axs[0].axhline(-predicted * 1e6, color='k', ls=':', lw=1.0)
    axs[0].set_ylabel('$\\Delta x$ between arms [$\\mu$m]')
    axs[0].set_title('(a) launch velocity cannot close it —\n'
                     'all four curves coincide', fontsize=9)
    axs[0].legend(fontsize=7)

    for i, (f, (t, d)) in enumerate(tilt_runs):
        axs[1].plot(t, d * 1e6, lw=2.0 if f == 1.0 else 1.1,
                    color='#d62728' if f == 1.0 else cmap(i / (len(tilt_runs) - 1)),
                    label=f'{f:g}$\\,\\Omega$' + (' (matched)' if f == 1.0 else ''))
    axs[1].set_ylabel('$\\Delta x$ between arms [$\\mu$m]')
    axs[1].set_title('(b) counter-rotating $\\vec{k}$ does:\n'
                     'tilt rate as a fraction of $\\Omega$', fontsize=9)
    axs[1].legend(fontsize=7, ncol=2)

    resid = [abs(d[-1]) for _, (t, d) in tilt_runs]
    axs[2].plot(tilt_scan, np.array(resid) * 1e6, 'o-', color='#1f77b4', lw=1.4)
    axs[2].axhline(COH_LEN * 1e6, color='#2ca02c', ls='--', lw=1.2)
    axs[2].annotate('coherence length (1 nK)', xy=(0.35, COH_LEN * 1e6),
                    fontsize=8, color='#2ca02c', va='bottom')
    axs[2].set_yscale('log')
    axs[2].set_xlabel('tilt rate / $\\Omega$')
    axs[2].set_ylabel('$|\\Delta x|$ at recombination [$\\mu$m]')
    axs[2].set_title('(c) residual opening vs tilt rate', fontsize=9)

    for ax in axs[:2]:
        ax.set_xlabel('$t$ [s]')
    for ax in axs:
        ax.grid(alpha=0.25, lw=0.5)
    fig.suptitle('Closing a rotating-frame interferometer (single loop, '
                 f'LMT-{NLMT}). The arms differ by $n\\hbar k/m$ along '
                 '$\\vec{k}$ and the Coriolis force acts on that difference — a '
                 '*differential* effect, so no common-mode launch velocity '
                 'removes it.', fontsize=11)
    fig.tight_layout()
    print(f'      predicted opening = {predicted:.4e} m')
    for f, (t, d) in tilt_runs:
        print(f'      tilt {f:4.2f} x Omega -> |dx| = {abs(d[-1]):.4e} m')
    return save(fig, outdir, 'coriolis_closure.png')


# ── the long-baseline strategy ───────────────────────────────────────────────

def fig_longbaseline(outdir, workdir, beam_radius=0.02):
    """
    Why holding k in the inertial frame does not scale to a long baseline, and
    what to use instead.
    """
    totals = [2.5, 5.0, 10.0]
    baselines = np.logspace(0, 3.5, 200)          # 1 m .. 3 km

    # (b) even loop counts, no tilt at all
    loop_gap = {}
    for L in (1, 2, 4):
        up, lo = arms(run(f'lb_L{L}', workdir, loops=L))
        t, ip = grid(up, lo)
        loop_gap[L] = abs((ip(up, 'x') - ip(lo, 'x'))[-1])

    # (c) cloud excursion with and without the closing launch velocity
    excursion = {}
    for T_tot in totals:
        guess = cc.analytic_closing_v0x(OMEGA_PERP, T_tot)
        v0x, resid, info = cc.solve_closing_v0x(
            f'lb_v{T_tot}', workdir, guess=guess, utype=POT, rotation=OMEGA_VEC,
            loopnumber=2, interrogation_time=T_tot / 4, **BASE)
        excursion[T_tot] = (info['max_excursion_uncompensated'],
                            info['max_excursion'], v0x, resid)

    fig, axs = plt.subplots(1, 3, figsize=(14.5, 4.5))

    # (a) beam walk-off
    for T_tot, c in zip(totals, ('#1f77b4', '#ff7f0e', '#d62728')):
        axs[0].loglog(baselines, OMEGA_PERP * T_tot * baselines * 1e3, lw=1.6,
                      color=c, label=f'$T_{{tot}}$ = {T_tot:g} s')
    axs[0].axhline(beam_radius * 1e3, color='k', ls='--', lw=1.2)
    axs[0].annotate(f'beam radius {beam_radius*1e3:.0f} mm', xy=(1.5, beam_radius*1e3),
                    fontsize=8, va='bottom')
    for Lb, name in ((100, 'MAGIS-100'), (2000, '2 km')):
        axs[0].axvline(Lb, color='gray', ls=':', lw=1.0)
        axs[0].annotate(name, xy=(Lb * 0.85, 0.2), rotation=90, fontsize=7,
                        color='gray', ha='right', va='bottom')
    axs[0].set_xlabel('baseline $L$ [m]')
    axs[0].set_ylabel('beam walk-off at far station [mm]')
    axs[0].set_title('(a) tip-tilt does not scale:\n'
                     'walk-off $= \\Omega_\\perp T_{tot} L$', fontsize=9)
    axs[0].legend(fontsize=7)
    axs[0].grid(alpha=0.25, lw=0.5, which='both')

    # (b) even loops close it for free
    Ls = sorted(loop_gap)
    axs[1].bar([str(L) for L in Ls], [loop_gap[L] * 1e6 for L in Ls],
               color=['#d62728', '#2ca02c', '#2ca02c'], width=0.55)
    axs[1].axhline(COH_LEN * 1e6, color='k', ls='--', lw=1.2)
    axs[1].annotate('coherence length (1 nK)', xy=(-0.4, COH_LEN * 1e6 * 1.3),
                    fontsize=8)
    axs[1].set_yscale('log')
    axs[1].set_xlabel('number of loops')
    axs[1].set_ylabel('$|\\Delta x|$ at recombination [$\\mu$m]')
    axs[1].set_title('(b) even loop counts close it with\n'
                     'NO tilt (LMT-101, fixed $T_{tot}$)', fontsize=9)
    axs[1].grid(alpha=0.25, lw=0.5, axis='y')

    # (c) cloud excursion
    w = 0.35
    idx = np.arange(len(totals))
    axs[2].bar(idx - w/2, [excursion[T][0] * 1e3 for T in totals], w,
               color='#d62728', label='$v_{0x}$ = 0')
    axs[2].bar(idx + w/2, [excursion[T][1] * 1e3 for T in totals], w,
               color='#2ca02c', label='closing $v_{0x}$')
    axs[2].axhspan(1e-2, 5, color='#2ca02c', alpha=0.10, lw=0)
    axs[2].annotate('"a few mm" — acceptable', xy=(-0.45, 2.2), fontsize=8,
                    color='#2ca02c', va='bottom')
    axs[2].set_yscale('log')
    axs[2].set_xticks(idx)
    axs[2].set_xticklabels([f'{T:g} s' for T in totals])
    axs[2].set_xlabel('total free-fall time')
    axs[2].set_ylabel('max transverse excursion [mm]')
    axs[2].set_title('(c) a launch velocity keeps the cloud\n'
                     'on axis (it scales as $T_{tot}^3$)', fontsize=9)
    axs[2].legend(fontsize=8)
    axs[2].grid(alpha=0.25, lw=0.5, axis='y')

    fig.suptitle('Long-baseline strategy: use an even number of loops for '
                 'closure (not tip-tilt, which walks the beam off the far '
                 'station), and a small launch velocity to keep the cloud on '
                 'the beam axis.', fontsize=11)
    fig.tight_layout()

    print('      loop gaps:', {L: f'{loop_gap[L]:.3e} m' for L in Ls})
    for T_tot in totals:
        e0, e1, v0x, resid = excursion[T_tot]
        print(f'      T_tot={T_tot:4.1f}s: v0x={v0x*1e3:7.4f} mm/s  '
              f'excursion {e0*1e3:7.3f} -> {e1*1e3:6.3f} mm  '
              f'(final x residual {resid*1e6:+.2f} um)')
    return save(fig, outdir, 'coriolis_longbaseline.png')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', default=None)
    ap.add_argument('--workdir', default=None)
    ap.add_argument('--only', choices=('traj', 'closure', 'strategy'), default=None)
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    outdir = args.outdir or os.path.join(os.path.dirname(here), 'examples', 'coriolis')
    workdir = args.workdir or os.path.join(here, '_coriolis_figs')
    os.makedirs(outdir, exist_ok=True)

    print(f'rotation: {OMEGA_VEC} rad/s (Earth, latitude 45 deg)')
    print(f'LMT-{NLMT}, Rabi {RABI:.0e} Hz, {T_TOT} s total free fall')
    print(f'arm separation n*(hbar k/m)*T = '
          f'{NLMT*RECOIL_V*T_TOT/2*100:.1f} cm at 1 loop')
    print(f'output: {outdir}\n')

    if args.only in (None, 'traj'):
        runs = {}
        for L in (1, 2, 4):
            for tag, rot in (('rot', OMEGA_VEC), ('inert', [0, 0, 0])):
                t0 = time.time()
                runs[(tag, L)] = run(f'fig_{tag}_L{L}', workdir, loops=L, rotation=rot)
                print(f'  {tag}/L={L}: {time.time()-t0:.1f} s')
        fig_effect(runs, outdir)
        fig_fountain(runs, outdir)
        fig_separation(runs, outdir)
        fig_planes(runs, outdir)
        fig_loops(runs, outdir)

    if args.only in (None, 'closure'):
        fig_closure(outdir, workdir)

    if args.only in (None, 'strategy'):
        fig_longbaseline(outdir, workdir)


if __name__ == '__main__':
    main()
