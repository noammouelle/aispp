#!/usr/bin/env python3
"""
Generate the rotating-frame figures.

    python tests/make_coriolis_figures.py [--outdir examples/coriolis]
    python tests/make_coriolis_figures.py --only short     # skip the slow ones

Two interferometer scales are covered.

``short``  T = 0.1 s, launched at 19.62 m/s and still rising at detection.
           Cheap, and the one the end-to-end tests use.

``long``   A long-baseline fountain: T = 1.25 s with the launch velocity chosen
           so the atom returns to its starting height exactly at recombination.
           Total free-fall time is 2*L*T, i.e. 2.5 s, 5 s and 10 s for 1, 2 and
           4 loops, reaching 7.7 m, 31 m and 123 m. The 1- and 2-loop cases sit
           in the MAGIS-100 / AION regime; the 4-loop one is past it and is
           included to show the scaling.

Figures
-------
coriolis_effect.png                what rotation does to a short sequence
coriolis_planes.png                short: arm traces on the x-y, x-z, y-z planes
coriolis_loops.png                 short: 1/2/4-loop enclosed areas
coriolis_separation.png            short: arm separations vs time
coriolis_fountain.png              long: uniform gravity + Coriolis, the full
                                   parabolic fountain and its deflections
coriolis_longbaseline_planes.png   long: arm traces on the three planes
coriolis_longbaseline_loops.png    long: 1/2/4-loop enclosed areas, 2.5-10 s
coriolis_closure.png               long: why no launch velocity closes a
                                   rotating-frame interferometer, and how
                                   counter-rotating k does

The rotation vector is Earth's, oriented for latitude 45 deg with z vertical and
x pointing east, so all three coordinate planes carry some signal.
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
POT = 'rotating_linear_pot'

# name -> kwargs handed to cc.run_case
CONFIGS = {
    'short': dict(interrogation_time=0.1),
    # The figures only need trajectories, and ultrafast skips just the action
    # phase, leaving the kinematics untouched. That matters here: over several
    # seconds of flight the action reaches ~1e3 m^2/s^2, so the quadrature
    # cannot hit the tight absolute tolerance the short runs use. The Sagnac
    # phases quoted on the loop figures come from the enclosed area, not from
    # the simulated phase, so nothing is lost.
    'long':  dict(interrogation_time=1.25, fountain=True, ultrafast=1),
}


def get_runs(cfg_name, loops=(1, 2, 4), workdir=None):
    """Run (and cache on disk) the rotating and inertial cases for one config."""
    kw = CONFIGS[cfg_name]
    runs = {}
    for L in loops:
        for tag, rot in (('rot', OMEGA_VEC), ('inert', [0, 0, 0])):
            t0 = time.time()
            runs[(tag, L)] = cc.run_case(f'{cfg_name}_{tag}_L{L}', workdir,
                                         utype=POT, rotation=rot, loopnumber=L,
                                         printtrajectory=True, **kw)
            print(f'    {cfg_name}/{tag}/L={L}: {time.time()-t0:.1f} s')
    return runs


def arms(run, potential=POT):
    sm = reconstruct_trajectories(load_trajectory(run['trajectory']),
                                  potential=potential, arm_grouping=True)
    return sm['upper'], sm['lower']


def common_grid(up, lo, n=3000):
    t = np.linspace(max(up['t'][0], lo['t'][0]),
                    min(up['t'][-1], lo['t'][-1]), n)
    return t, (lambda d, k: np.interp(t, d['t'], d[k]))


# ── figure: what rotation does (short) ───────────────────────────────────────

def fig_effect(runs, outdir):
    up_r, lo_r = arms(runs[('rot', 1)])
    up_i, lo_i = arms(runs[('inert', 1)])
    t, ip = common_grid(up_r, lo_r)

    com_r = {k: 0.5 * (ip(up_r, k) + ip(lo_r, k)) for k in 'xy'}
    com_i = {k: 0.5 * (ip(up_i, k) + ip(lo_i, k)) for k in 'xy'}

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))

    axs[0].plot(t, (ip(up_r, 'z') - ip(lo_r, 'z')) * 1e6, color='#1f77b4', lw=1.4)
    axs[0].set_ylabel('$\\Delta z$ between arms [$\\mu$m]')
    axs[0].set_title('(a) the interferometer diamond\n(unchanged by rotation)',
                     fontsize=9)

    axs[1].plot(t, (com_r['x'] - com_i['x']) * 1e6, color='#d62728', lw=1.4, label='$x$')
    axs[1].plot(t, (com_r['y'] - com_i['y']) * 1e6, color='#2ca02c', lw=1.4, label='$y$')
    axs[1].set_ylabel('rotating $-$ inertial [$\\mu$m]')
    axs[1].set_title('(b) common-mode Coriolis deflection\nof the cloud centre',
                     fontsize=9)
    axs[1].legend(fontsize=8)

    axs[2].plot(t, (ip(up_r, 'x') - ip(lo_r, 'x')) * 1e9, color='#d62728',
                lw=1.4, label='rotating')
    axs[2].plot(t, (ip(up_i, 'x') - ip(lo_i, 'x')) * 1e9, color='#888888',
                lw=1.4, ls='--', label='inertial')
    axs[2].set_ylabel('$\\Delta x$ between arms [nm]')
    axs[2].set_title('(c) differential deflection: this is what\n'
                     'encloses area — and leaves the arms open', fontsize=9)
    axs[2].legend(fontsize=8)

    for ax in axs:
        ax.set_xlabel('$t$ [s]')
        ax.axhline(0, color='k', lw=0.6, alpha=0.4)
        ax.grid(alpha=0.25, lw=0.5)

    fig.suptitle('Single-loop Mach-Zehnder at Earth rotation rate ($T$ = 0.1 s): '
                 'the arms separate in $z$ by design, and the Coriolis force '
                 'then deflects them differently in $x$', fontsize=11)
    fig.tight_layout()
    return save(fig, outdir, 'coriolis_effect.png')


# ── figure: uniform gravity + Coriolis, the full fountain (long) ─────────────

def fig_fountain(runs, outdir):
    up_r, lo_r = arms(runs[('rot', 1)])
    up_i, lo_i = arms(runs[('inert', 1)])
    t, ip = common_grid(up_r, lo_r)

    z_com = 0.5 * (ip(up_r, 'z') + ip(lo_r, 'z'))
    dx_com = 0.5 * (ip(up_r, 'x') + ip(lo_r, 'x')) - 0.5 * (ip(up_i, 'x') + ip(lo_i, 'x'))
    dy_com = 0.5 * (ip(up_r, 'y') + ip(lo_r, 'y')) - 0.5 * (ip(up_i, 'y') + ip(lo_i, 'y'))
    vz = np.gradient(z_com, t)

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.3))

    axs[0].plot(t, z_com, color='#1f77b4', lw=1.5)
    axs[0].set_ylabel('$z$ [m]')
    axs[0].set_title(f'(a) uniform gravity: the fountain\n'
                     f'apex {z_com.max():.2f} m, free fall {t[-1]-t[0]:.2f} s',
                     fontsize=9)

    axs[1].plot(t, dx_com * 1e3, color='#d62728', lw=1.5, label='$x$ (east)')
    axs[1].plot(t, dy_com * 1e3, color='#2ca02c', lw=1.5, label='$y$ (north)')
    axs[1].set_ylabel('rotating $-$ inertial [mm]')
    axs[1].set_title('(b) Coriolis deflection of the cloud\n'
                     '$\\ddot{r} = -2\\,\\Omega \\times \\dot{r}$', fontsize=9)
    axs[1].legend(fontsize=8)

    axs[2].plot(t, (ip(up_r, 'z') - ip(lo_r, 'z')) * 1e3, color='#1f77b4', lw=1.5,
                label='$\\Delta z$ [mm]')
    ax2 = axs[2].twinx()
    ax2.plot(t, (ip(up_r, 'x') - ip(lo_r, 'x')) * 1e9, color='#d62728', lw=1.5,
             label='$\\Delta x$ [nm]')
    ax2.set_ylabel('$\\Delta x$ between arms [nm]', color='#d62728')
    ax2.tick_params(axis='y', labelcolor='#d62728')
    axs[2].set_ylabel('$\\Delta z$ between arms [mm]', color='#1f77b4')
    axs[2].tick_params(axis='y', labelcolor='#1f77b4')
    axs[2].set_title('(c) arm separation: $\\Delta z$ from the recoil,\n'
                     '$\\Delta x$ from Coriolis', fontsize=9)

    # The vertical velocity reverses at the apex, and with it the sign of the
    # Coriolis acceleration. The deflection itself keeps growing -- what changes
    # is its curvature, so the apex shows up as an inflection in panel (b).
    apex = t[np.argmax(z_com)]
    for ax in axs:
        ax.axvline(apex, color='k', lw=0.7, ls=':', alpha=0.6)
        ax.set_xlabel('$t$ [s]')
        ax.grid(alpha=0.25, lw=0.5)
    axs[0].annotate('apex ($\\dot z = 0$)', xy=(apex, z_com.max()),
                    xytext=(0.55, 0.35), textcoords='axes fraction', fontsize=8,
                    arrowprops=dict(arrowstyle='->', lw=0.8))

    fig.suptitle('Long-baseline fountain with uniform gravity and Coriolis: '
                 f'$T$ = 1.25 s, {t[-1]-t[0]:.1f} s total free fall, '
                 f'{z_com.max():.1f} m apex. The eastward deflection grows '
                 'fastest at the apex and then levels off — $\\dot z$ reverses '
                 'there, so the Coriolis acceleration does too.', fontsize=11)
    fig.tight_layout()
    return save(fig, outdir, 'coriolis_fountain.png')


# ── figure: arm separations vs time (short) ──────────────────────────────────

def fig_separation(runs, outdir, loops=(1, 2, 4)):
    fig, axs = plt.subplots(2, len(loops), figsize=(13, 6), sharex='col')
    for col, L in enumerate(loops):
        up, lo = arms(runs[('rot', L)])
        t, ip = common_grid(up, lo)
        d = {k: ip(up, k) - ip(lo, k) for k in 'xyz'}

        axs[0, col].plot(t, d['z'] * 1e6, color='#1f77b4', lw=1.3)
        axs[0, col].set_title(f'{L} loop' + ('s' if L > 1 else ''), fontsize=10)
        axs[0, col].set_ylabel('$\\Delta z$ [$\\mu$m]' if col == 0 else '')
        axs[1, col].plot(t, d['x'] * 1e9, color='#d62728', lw=1.3, label='$\\Delta x$')
        axs[1, col].plot(t, d['y'] * 1e9, color='#2ca02c', lw=1.3, label='$\\Delta y$')
        axs[1, col].set_ylabel('transverse [nm]' if col == 0 else '')
        axs[1, col].set_xlabel('$t$ [s]')
        for ax in (axs[0, col], axs[1, col]):
            ax.axhline(0, color='k', lw=0.6, alpha=0.4)
            ax.grid(alpha=0.25, lw=0.5)
    axs[1, 0].legend(fontsize=8)
    fig.suptitle('Arm separation vs time. $\\Delta z$ reverses sign once per '
                 'loop, so the signed area — and with it the rotation phase — '
                 'cancels for an even number of loops', fontsize=11)
    fig.tight_layout()
    return save(fig, outdir, 'coriolis_separation.png')


# ── figures: plane traces and loop comparison ────────────────────────────────

def fig_planes(runs, outdir, name, subtitle):
    traj = load_trajectory(runs[('rot', 1)]['trajectory'])
    fig, axes, info = plot_trajectory_planes(traj, potential=POT,
                                             planes=('xy', 'xz', 'yz'))
    fig.suptitle(fig._suptitle.get_text() + '\n' + subtitle, fontsize=10)
    fig.tight_layout()
    p = save(fig, outdir, name)
    print(f'      A = {info["area"]} m^2, dphi = {info["sagnac_phase"]:.4e} rad')
    return p


def fig_loops(runs, outdir, name, loops=(1, 2, 4), labels=None):
    # The areas come from the NON-rotating runs: rotation deflects the two arms
    # differently, and that O(Omega) piece of the area shows up as an O(Omega^2)
    # term in 2 m Omega.A / hbar, which for even loop counts swamps the
    # first-order phase it is meant to predict.
    entries = [(labels[i] if labels else f'{L} loop' + ('s' if L > 1 else ''),
                load_trajectory(runs[('inert', L)]['trajectory']))
               for i, L in enumerate(loops)]
    fig, axes, rows = plot_loop_comparison(entries, potential=POT, plane='xz',
                                           sagnac_rotation=OMEGA_VEC)
    p = save(fig, outdir, name)
    for r in rows:
        print(f"      {r['label']:22s} A = {r['area'][1]:+.4e} m^2   "
              f"dphi = {r['sagnac_phase']:+.4e} rad")
    return p


# ── figure: what actually closes the interferometer ──────────────────────────

def fig_closure(outdir, workdir):
    """
    A Mach-Zehnder does not close in a rotating frame. Show why no launch
    velocity fixes it, and what does.
    """
    T, L = 1.25, 1
    hbar, mass = 1.054571817e-34, 86.90888 * 1.660539066e-27
    u = hbar * cc.K_EFF / mass                       # single-photon recoil
    predicted = 2 * OMEGA_VEC[1] * u * T * T

    def case(tag, **kw):
        r = cc.run_case(tag, workdir, utype=POT, rotation=OMEGA_VEC,
                        loopnumber=L, interrogation_time=T, fountain=True,
                        ultrafast=1, printtrajectory=True, **kw)
        up, lo = arms(r)
        t, ip = common_grid(up, lo)
        return t, np.array([ip(up, k) - ip(lo, k) for k in 'xyz'])

    v0x_list = [0.0, 0.01, 0.05, 0.2]
    vel_runs = [(v, case(f'close_vx{v}', v0x=v)) for v in v0x_list]
    tilt_scan = [0.0, 0.5, 0.9, 1.0, 1.1, 1.5, 2.0]
    tilt_runs = [(f, case(f'close_tilt{f}', v0x=0.01,
                          tilt_compensation=OMEGA_VEC * f)) for f in tilt_scan]

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.4))

    # (a) launch velocity does nothing to the closure
    cmap = plt.get_cmap('viridis')
    for i, (v, (t, d)) in enumerate(vel_runs):
        axs[0].plot(t, d[0] * 1e9, lw=1.4, color=cmap(i / max(1, len(vel_runs) - 1)),
                    label=f'$v_{{0x}}$ = {v:g} m/s')
    axs[0].axhline(-predicted * 1e9, color='k', ls=':', lw=1.0)
    axs[0].annotate('$-2\\,\\Omega_\\perp\\,(\\hbar k/m)\\,T^2$',
                    xy=(t[len(t)//2], -predicted * 1e9), fontsize=8,
                    xytext=(0.08, 0.12), textcoords='axes fraction',
                    arrowprops=dict(arrowstyle='->', lw=0.8))
    axs[0].set_ylabel('$\\Delta x$ between arms [nm]')
    axs[0].set_title('(a) launch velocity cannot close it —\n'
                     'all four curves coincide', fontsize=9)
    axs[0].legend(fontsize=7)

    # (b) tilting the wavevector does
    for i, (f, (t, d)) in enumerate(tilt_runs):
        lw = 2.0 if f == 1.0 else 1.1
        col = '#d62728' if f == 1.0 else cmap(i / max(1, len(tilt_runs) - 1))
        axs[1].plot(t, d[0] * 1e9, lw=lw, color=col,
                    label=f'{f:g}$\\,\\Omega$' + (' (matched)' if f == 1.0 else ''))
    axs[1].set_ylabel('$\\Delta x$ between arms [nm]')
    axs[1].set_title('(b) counter-rotating $\\vec{k}$ does:\n'
                     'tilt rate as a fraction of $\\Omega$', fontsize=9)
    axs[1].legend(fontsize=7, ncol=2)

    # (c) residual vs tilt rate
    resid = [abs(d[0][-1]) for _, (t, d) in tilt_runs]
    axs[2].plot(tilt_scan, np.array(resid) * 1e9, 'o-', color='#1f77b4', lw=1.4)
    axs[2].set_yscale('log')
    axs[2].set_xlabel('tilt rate / $\\Omega$')
    axs[2].set_ylabel('$|\\Delta x|$ at recombination [nm]')
    axs[2].set_title('(c) residual opening vs tilt rate\n'
                     f'{resid[0]*1e9:.0f} nm $\\rightarrow$ '
                     f'{resid[tilt_scan.index(1.0)]*1e12:.2f} pm at matched rate',
                     fontsize=9)

    for ax in axs[:2]:
        ax.set_xlabel('$t$ [s]')
    for ax in axs:
        ax.grid(alpha=0.25, lw=0.5)
        ax.axhline(0, color='k', lw=0.6, alpha=0.4)

    fig.suptitle('Closing a rotating-frame interferometer. The arms differ by the '
                 'recoil velocity $\\hbar k/m$ along $\\vec{k}$, and the Coriolis '
                 'force acts on that difference — a *differential* effect, so no '
                 'common-mode launch velocity can undo it.', fontsize=11)
    fig.tight_layout()
    p = save(fig, outdir, 'coriolis_closure.png')
    print(f'      predicted opening 2*Om*(hbar k/m)*T^2 = {predicted:.4e} m')
    for f, (t, d) in tilt_runs:
        print(f'      tilt {f:4.2f} x Omega -> |dx| = {abs(d[0][-1]):.4e} m')
    return p


def save(fig, outdir, name):
    p = os.path.join(outdir, name)
    fig.savefig(p, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('  wrote', p)
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', default=None)
    ap.add_argument('--workdir', default=None)
    ap.add_argument('--only', choices=('short', 'long'), default=None)
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    outdir = args.outdir or os.path.join(os.path.dirname(here), 'examples', 'coriolis')
    workdir = args.workdir or os.path.join(here, '_coriolis_figs')
    os.makedirs(outdir, exist_ok=True)

    print(f'rotation: {OMEGA_VEC} rad/s  (Earth rate, latitude 45 deg)')
    print(f'output:   {outdir}')

    if args.only != 'long':
        print('\nshort sequence (T = 0.1 s)')
        runs = get_runs('short', workdir=workdir)
        fig_effect(runs, outdir)
        fig_planes(runs, outdir, 'coriolis_planes.png',
                   '$T$ = 0.1 s, single loop')
        fig_loops(runs, outdir, 'coriolis_loops.png')
        fig_separation(runs, outdir)

    if args.only != 'short':
        print('\nlong-baseline fountain (T = 1.25 s)')
        runs = get_runs('long', workdir=workdir)
        fig_fountain(runs, outdir)
        fig_planes(runs, outdir, 'coriolis_longbaseline_planes.png',
                   'long-baseline fountain: $T$ = 1.25 s, 2.5 s free fall, '
                   '7.7 m apex')
        fig_loops(runs, outdir, 'coriolis_longbaseline_loops.png',
                  labels=['1 loop — 2.5 s, 7.7 m',
                          '2 loops — 5.0 s, 31 m',
                          '4 loops — 10.0 s, 123 m'])
        fig_closure(outdir, workdir)


if __name__ == '__main__':
    main()
