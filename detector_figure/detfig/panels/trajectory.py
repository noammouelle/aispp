"""Panel (d): what the atoms do during one shot, as height z versus time t.

Each cloud is a ballistic parabola z(t) = z0 + v0 t - g t^2 / 2, released at
t = 0 (a drop has v0 = 0).  The clock pulses are vertical lines: at these
scales light crosses the baseline instantaneously, and the *same* pulses
address every cloud, which is the gradiometer's common-mode rejection.
Between the first and last pulse each cloud is split into two arms
(Mach-Zehnder, separation exaggerated by ``split``).  Brackets under the
time axis mark the cycle steps (1)-(4) used in panel (a).

By default the clouds are the fountains/drops of ``instrument.atoms``, so
panels (a) and (d) cannot disagree.

With ``trajectory.engine: mzplots`` the panel is drawn by the LMT sequence
generator of mz-plots (github.com/noammouelle/mz-plots) instead: order-n LMT
pulse ladders with alternating directions, finite light travel time and a
gradiometer pair, embedded here in a scaled scope.  Only the cycle-step
brackets are added on top, aligned with its t0 / t_close / t_det.
"""

import importlib
import math
import sys
from pathlib import Path

from .. import glyphs as g
from ..tikz import Canvas, pt


def clouds_from_instrument(cfg):
    """(label, z0, v0, z_end) for each trajectory drawn in panel (a)."""
    grav = cfg.trajectory.g
    out = []
    for t in cfg.instrument.atoms.trajectories:
        rise = t["apex"] - t["start"]
        out.append((t.get("label", t.get("id", "")), t["start"], math.sqrt(2 * grav * rise) if rise > 0 else 0.0, t["end"]))
    return out


def landing_time(z0, v0, z_end, grav):
    """Time at which a cloud, on its way down, reaches z_end (None: never)."""
    disc = v0 * v0 + 2 * grav * (z0 - z_end)
    return (v0 + math.sqrt(disc)) / grav if disc >= 0 else None


def load_mzplots(cfg):
    """Import mz-plots' lmt_sequence_diagram (installed, or from trajectory.mzplots.path)."""
    path = cfg.trajectory.mzplots.path
    if path:  # relative paths are relative to the detector_figure directory
        base = Path(__file__).resolve().parents[2]
        path = str((base / Path(path).expanduser()).resolve())
        if path not in sys.path:
            sys.path.insert(0, path)
    try:
        return importlib.import_module("lmt_sequence_diagram")
    except ImportError:
        return None


def latex_packages(cfg):
    """Extra LaTeX packages this panel needs (mz-plots uses amsmath, braket)."""
    if "trajectory" not in cfg.layout.panels or cfg.trajectory.engine != "mzplots":
        return ()
    lmt = load_mzplots(cfg)
    return tuple(lmt.LATEX_PACKAGES) if lmt else ()


def draw(cfg):
    if cfg.trajectory.engine == "mzplots":
        lmt = load_mzplots(cfg)
        if lmt is not None:
            return _draw_mzplots(cfg, lmt)
        if not cfg.trajectory.mzplots.fallback:
            raise ValueError("trajectory.engine is mzplots but lmt_sequence_diagram cannot be imported: "
                             "pip install -e <mz-plots checkout>, or set trajectory.mzplots.path")
        print("detfig: mz-plots not found, drawing panel (d) with the built-in engine", file=sys.stderr)
    return _draw_builtin(cfg)


def _draw_mzplots(cfg, lmt):
    tr = cfg.trajectory
    c = Canvas()
    mcfg = lmt.config_from_mapping(dict(tr.mzplots.config), "trajectory.mzplots.config")
    body, info = lmt.picture_body(mcfg)
    sx = tr.width / (info["tmax"] - info["tmin"])
    sy = tr.height / (info["zmax"] - info["zmin"])

    def X(t):
        return (t - info["tmin"]) * sx

    with c.on("content"):
        c.raw(f"\\begin{{scope}}[shift={{({-info['tmin'] * sx:.4f},{-info['zmin'] * sy:.4f})}}, "
              f"xscale={sx:.5f}, yscale={sy:.5f}, font={{{tr.mzplots.font}}}]")
        for line in body:
            c.raw(line)
        c.raw("\\end{scope}")
    # cycle steps under the time axis, below mz-plots' own t labels (at z = -0.3)
    yb = min(0.0, (-0.3 - info["zmin"]) * sy) - 0.32
    t_prep_end = info["tmin"] + 0.45 * (info["t0"] - info["tmin"])
    spans = [(info["tmin"], t_prep_end), (t_prep_end, info["t0"]),
             (info["t0"], info["t_close"]), (info["t_close"], info["t_det"])]
    for n, (a, b) in enumerate(spans, 1):
        xa, xb = X(a) + 0.03, X(b) - 0.03
        with c.on("annotations"):
            c.draw(f"{pt(xa, yb + 0.08)} -- {pt(xa, yb)} -- {pt(xb, yb)} -- {pt(xb, yb + 0.08)}",
                   "dfleader, line width=0.45pt")
        g.marker(c, (xa + xb) / 2, yb - 0.2, n)
    c.anchor("_letter", -0.75, tr.height + 0.65)
    return c


def _draw_builtin(cfg):
    tr = cfg.trajectory
    c = Canvas()
    grav, W, Hp = tr.g, tr.width, tr.height
    clouds = ([(cl.get("label", ""), cl["z0"], cl.get("v0", 0.0), cl.get("z_end")) for cl in tr.clouds]
              if tr.clouds else clouds_from_instrument(cfg))
    if not clouds:
        raise ValueError("trajectory panel: no clouds (set trajectory.clouds or instrument.atoms.trajectories)")
    if tr.clouds and any("z0" not in cl for cl in tr.clouds):
        raise ValueError("trajectory.clouds: each cloud needs z0 (and optionally v0, z_end, label)")
    ifo = tr.interferometer
    t_pulses = [ifo.t0 + k * ifo.T for k in ifo.pulses]
    t_end = tr.t_max if tr.t_max else max(t_pulses) + 0.6 * ifo.T
    t_start = -tr.t_prepare

    def z_at(z0, v0, t):
        return z0 + v0 * t - grav * t * t / 2

    def t_stop(z0, v0, z_end):
        t_land = landing_time(z0, v0, z_end, grav) if z_end is not None else None
        return min(t_end, t_land) if t_land else t_end

    t_split0, t_split1 = min(t_pulses), max(t_pulses)
    t_mid = (t_split0 + t_split1) / 2

    def kick(t):
        """Extra height of the second arm (exaggerated by ``split``)."""
        if not ifo.show_arms or not t_split0 <= t <= t_split1:
            return 0.0
        return ifo.split * (t - t_split0 if t <= t_mid else t_split1 - t)

    zs = [z_at(z0, v0, t) + k for _, z0, v0, z_end in clouds
          for t in _samples(0, t_stop(z0, v0, z_end), 40) for k in (0.0, kick(t))]
    z_lo, z_hi = tr.z_range if tr.z_range else (min(zs), max(zs))
    pad = 0.06 * (z_hi - z_lo or 1)
    z_lo, z_hi = z_lo - pad, z_hi + pad

    def X(t):
        return (t - t_start) / (t_end - t_start) * W

    def Y(z):
        return (z - z_lo) / (z_hi - z_lo) * Hp

    font = cfg.document.font
    with c.on("back"):
        c.rect(0, 0, W, Hp, "fill=dfinterior, draw=none", cmd="filldraw")
    # clock pulses: same light for every cloud
    for k, tp in enumerate(t_pulses):
        with c.on("content"):
            c.draw(f"{pt(X(tp), 0)} -- {pt(X(tp), Hp)}", "dfbeam, line width=1.1pt, opacity=0.75")
        if k < len(ifo.labels):
            c.text(X(tp), Hp, ifo.labels[k], f"anchor=south, inner sep=1.5pt, text=dfbeam!70!black, font={{{font}}}")
        if ifo.pulse_arrows:  # each pulse goes down and, retro-reflected, comes back up
            style = "dfbeam!70!black, line width=0.6pt, -{Stealth[length=1.1mm]}"
            with c.on("front"):
                c.draw(f"{pt(X(tp) - 0.07, Hp - 0.06)} -- {pt(X(tp) - 0.07, Hp - 0.36)}", style)
                c.draw(f"{pt(X(tp) + 0.07, Hp - 0.36)} -- {pt(X(tp) + 0.07, Hp - 0.06)}", style)
    if ifo.show_T and len(t_pulses) > 1:  # mark the pulse separation T
        ya = Hp + 0.48
        with c.on("annotations"):
            c.draw(f"{pt(X(t_pulses[0]), ya)} -- {pt(X(t_pulses[1]), ya)}",
                   "dfleader, line width=0.45pt, {Stealth[length=1.1mm]}-{Stealth[length=1.1mm]}")
            c.text((X(t_pulses[0]) + X(t_pulses[1])) / 2, ya, "$T$",
                   f"anchor=south, inner sep=1pt, text=dftext, font={{{font}}}")
    # clouds: waiting in the source (t < 0), then ballistic, split between pulses
    for name, z0, v0, z_end in clouds:
        t1 = t_stop(z0, v0, z_end)
        with c.on("front"):
            c.draw(f"{pt(X(t_start), Y(z0))} -- {pt(X(0), Y(z0))}",
                   "dfatoms, line width=0.9pt, dash pattern=on 1.5pt off 1.2pt, opacity=0.7")
            path = " -- ".join(pt(X(t), Y(z_at(z0, v0, t))) for t in _samples(0, t1, 60))
            c.draw(path, "dfatoms, line width=0.9pt")
            if ifo.show_arms and t_split1 <= t1:
                path = " -- ".join(pt(X(t), Y(z_at(z0, v0, t) + kick(t))) for t in _samples(t_split0, t_split1, 40))
                c.draw(path, "dfatoms!55, line width=0.9pt")
        if ifo.arm_labels and t_split1 <= t1 and name == clouds[-1][0]:
            # label the two arms of one cloud in the first half of the sequence
            tq = t_split0 + 0.55 * (t_mid - t_split0)
            lo, hi = sorted([(Y(z_at(z0, v0, tq)), 0), (Y(z_at(z0, v0, tq) + kick(tq)), 1)])
            for (yv, k), anchor in ((lo, "north west"), (hi, "south east")):
                c.text(X(tq), yv, ifo.arm_labels[k],
                       f"anchor={anchor}, inner sep=1pt, text=dftext, font=\\sffamily\\scriptsize")
        if t1 < t_end:  # detected where it lands
            g.cloud(c, X(t1), Y(z_at(z0, v0, t1)), 0.05, 0.9)
        g.cloud(c, X(0), Y(z0), 0.06)
        if tr.label_clouds and name:
            c.text(X(t_start), Y(z0), name, f"anchor=east, inner sep=2pt, text=dftext, font={{{font}}}")
        c.anchor(f"cloud:{name}" if name else "cloud", X(0), Y(z0))

    # axes
    with c.on("annotations"):
        c.draw(f"{pt(0, 0)} -- {pt(W + 0.15, 0)}", "dfleader, line width=0.5pt, -{Stealth[length=1.3mm]}")
        c.draw(f"{pt(0, 0)} -- {pt(0, Hp + 0.15)}", "dfleader, line width=0.5pt, -{Stealth[length=1.3mm]}")
        c.text(W + 0.15, 0, tr.t_label, f"anchor=west, inner sep=1.5pt, text=dftext, font={{{font}}}")
        c.text(0, Hp + 0.15, tr.z_label, f"anchor=south, inner sep=1.5pt, text=dftext, font={{{font}}}")
        c.draw(f"{pt(X(0), 0)} -- {pt(X(0), -0.08)}", "dfleader, line width=0.5pt")
        c.text(X(0), -0.08, "0", "anchor=north, inner sep=1pt, text=dftext, font=\\sffamily\\scriptsize")

    # cycle steps under the time axis
    steps = {"prepare": (t_start, -tr.t_prepare * 0.25), "launch": (-tr.t_prepare * 0.25, t_split0),
             "interrogate": (t_split0, t_split1), "detect": (t_split1, t_end)}
    yb = -0.42
    for n, key in enumerate(("prepare", "launch", "interrogate", "detect"), 1):
        a, b = steps[key]
        if b - a <= 1e-9:
            continue
        xa, xb = X(a) + 0.04, X(b) - 0.04
        with c.on("annotations"):
            c.draw(f"{pt(xa, yb + 0.08)} -- {pt(xa, yb)} -- {pt(xb, yb)} -- {pt(xb, yb + 0.08)}",
                   "dfleader, line width=0.45pt")
        g.marker(c, (xa + xb) / 2, yb - 0.2, n)
    c.anchor("_letter", -0.75, Hp + 0.65)
    return c


def _samples(a, b, n):
    return [a + (b - a) * i / n for i in range(n + 1)]
