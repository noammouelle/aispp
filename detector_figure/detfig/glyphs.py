"""Reusable drawing components (vacuum hardware, optics, atoms, markers).

Every function draws onto a :class:`~detfig.tikz.Canvas` in absolute page
coordinates (cm).  Colours are referred to by their style name (``beam``,
``metal``, ...), which :func:`color` turns into the TikZ colour defined in
the document preamble.  The panels compose these glyphs; nothing in here
knows about YAML.
"""

import itertools
import math

from .tikz import num, pt

PALETTE = set()  # filled by build.py with the names under style.colors
_MARKER_IDS = itertools.count(1)  # unique TikZ node names for markers


def color(name):
    """Map a style colour name to its TikZ name; pass other expressions through."""
    head = name.split("!")[0]
    return ("df" + name) if head in PALETTE else name


# -- vacuum hardware ---------------------------------------------------------

def interior_shading(depth, vertical=True):
    """Cylinder-like shading for the inside of a tube or chamber."""
    edge = f"dfmetal!{num(40 * depth)}!dfinterior"
    if vertical:
        return f"left color={edge}, right color={edge}, middle color=dfinterior"
    return f"top color={edge}, bottom color={edge}, middle color=dfinterior"


def _chamfered(x0, y0, x1, y1, ch):
    return [(x0 + ch, y0), (x1 - ch, y0), (x1, y0 + ch), (x1, y1 - ch),
            (x1 - ch, y1), (x0 + ch, y1), (x0, y1 - ch), (x0, y0 + ch)]


def chamber(c, x0, y0, x1, y1, wall, depth, chamfer=None, vertical=True):
    """A closed vacuum chamber drawn in section: metal wall + shaded interior."""
    ch = min(x1 - x0, y1 - y0) * 0.18 if chamfer is None else chamfer
    with c.on("walls"):
        c.polyline(_chamfered(x0, y0, x1, y1, ch), "fill=dfmetal, draw=dfmetal!70!black",
                   closed=True, cmd="filldraw")
    with c.on("interiors"):
        inner = _chamfered(x0 + wall, y0 + wall, x1 - wall, y1 - wall, max(ch - wall * 0.4, 0))
        c.polyline(inner, interior_shading(depth, vertical), closed=True, cmd="shade")


def vtube(c, x, y0, y1, r, wall, depth):
    """Vertical tube in longitudinal section between heights y0 < y1."""
    with c.on("walls"):
        c.rect(x - r - wall, y0, x + r + wall, y1, "fill=dfmetal, draw=dfmetal!70!black")
    with c.on("interiors"):
        c.rect(x - r, y0, x + r, y1, interior_shading(depth, True), cmd="shade")


def htube(c, x0, x1, y, r, wall, depth):
    """Horizontal tube in longitudinal section between x0 < x1."""
    with c.on("walls"):
        c.rect(x0, y - r - wall, x1, y + r + wall, "fill=dfmetal, draw=dfmetal!70!black")
    with c.on("interiors"):
        c.rect(x0, y - r, x1, y + r, interior_shading(depth, False), cmd="shade")


def end_cap(c, x, y, r, wall, up):
    """Blank flange closing a vertical tube."""
    s = 1 if up else -1
    with c.on("front"):
        c.rect(x - r - 1.6 * wall, y, x + r + 1.6 * wall, y + s * 1.6 * wall,
               "fill=dfmetal, draw=dfmetal!70!black")


def shield(c, x, y0, y1, r_in, thickness, gap):
    """Magnetic shield (both walls of the section) around a vertical tube."""
    with c.on("back"):
        for side in (-1, 1):
            xa = x + side * (r_in + gap)
            xb = xa + side * thickness
            c.rect(min(xa, xb), y0, max(xa, xb), y1, "fill=dfshield, draw=dfshield!70!black")


def coils(c, x, y0, y1, r, pitch, size, style="bars"):
    """Bias/guide coils seen in a longitudinal section.

    ``bars``: conductors running along the tube (MAGIS-style bars giving a
    transverse bias field) drawn as lines; ``turns``: solenoid windings drawn
    as rows of wire cross-sections."""
    with c.on("back"):
        if style == "bars":
            for side in (-1, 1):
                c.draw(f"{pt(x + side * r, y0)} -- {pt(x + side * r, y1)}",
                       f"dfcoil, line width={num(size * 2 * 28.45)}pt")
            return
        n = max(int((y1 - y0) / pitch), 1)
        step = (y1 - y0) / n
        for i in range(n):
            y = y0 + (i + 0.5) * step
            for side in (-1, 1):
                c.circle(x + side * r, y, size, "fill=dfcoil, draw=dfcoil!60!black, line width=0.2pt")


def break_mark(c, x, y, half_w):
    """Slanted 'tube continues' cut across a vertical tube end."""
    h = half_w * 0.25
    with c.on("front"):
        c.polyline([(x - half_w, y - h), (x + half_w, y + h), (x + half_w, y + h + 0.09),
                    (x - half_w, y - h + 0.09)], "fill=white, draw=none", closed=True, cmd="filldraw")
        c.draw(f"{pt(x - half_w, y - h)} -- {pt(x + half_w, y + h)}", "dfmetal!60!black")
        c.draw(f"{pt(x - half_w, y - h + 0.09)} -- {pt(x + half_w, y + h + 0.09)}", "dfmetal!60!black")


def baseline_gap(c, x0, x1, y0, y1, half_w, slope=0.12, dots=True):
    """Cut out the band y0..y1 of the drawing (for x0 < x < x1): the baseline
    continues, but is not drawn to scale.  Slanted cut lines across the tube
    (|x| < half_w) and a vertical ellipsis mark the break."""
    def edge(y, x):
        return y + slope * x
    with c.on("annotations"):
        c.polyline([(x0, edge(y0, x0)), (x1, edge(y0, x1)), (x1, edge(y1, x1)), (x0, edge(y1, x0))],
                   "fill=white, draw=none", closed=True, cmd="filldraw")
        for y in (y0, y1):
            c.draw(f"{pt(-half_w, edge(y, -half_w))} -- {pt(half_w, edge(y, half_w))}",
                   "dfmetal!60!black, line width=0.6pt")
        if dots:
            ym = (y0 + y1) / 2
            for k in (-1, 0, 1):
                c.circle(0, ym + k * (y1 - y0) * 0.22, 0.035, "fill=dfleader, draw=none", cmd="fill")


# -- light -------------------------------------------------------------------

def beam_shading(name="beam"):
    col = color(name)
    return f"left color={col}!15, right color={col}!15, middle color={col}!85, fill opacity=0.9"


def vbeam(c, x, y0, y1, w, name="beam"):
    """Interferometry beam of finite width (Gaussian-looking profile)."""
    with c.on("content"):
        c.rect(x - w / 2, y0, x + w / 2, y1, beam_shading(name), cmd="shade")


def hbeam(c, x0, x1, y, w, name="beam"):
    col = color(name)
    with c.on("content"):
        c.rect(x0, y - w / 2, x1, y + w / 2,
               f"top color={col}!15, bottom color={col}!15, middle color={col}!85, fill opacity=0.9",
               cmd="shade")


def beam_cone(c, x, y_narrow, y_wide, w_narrow, w_wide, name="beam"):
    """Beam expanding (or contracting) between two heights, e.g. inside a telescope."""
    with c.on("content"):
        c.polyline([(x - w_narrow / 2, y_narrow), (x + w_narrow / 2, y_narrow),
                    (x + w_wide / 2, y_wide), (x - w_wide / 2, y_wide)],
                   beam_shading(name), closed=True, cmd="shade")


def beam_arrows(c, x, y, w, retro, size=0.5, gap=0.05):
    """Propagation arrows beside the beam edges: down on the left (incoming) and,
    with ``retro``, up on the right (retro-reflected)."""
    style = "dfbeam!65!black, line width=0.8pt, -{Stealth[length=1.6mm, width=1.4mm]}"
    xl, xr = x - w / 2 - gap, x + w / 2 + gap
    with c.on("front"):
        c.draw(f"{pt(xl, y + size / 2)} -- {pt(xl, y - size / 2)}", style)
        if retro:
            c.draw(f"{pt(xr, y - size / 2)} -- {pt(xr, y + size / 2)}", style)


def lattice_line(c, x0, y0, x1, y1, width=0.6):
    """Optical lattice drawn as a dashed line (transport or launch)."""
    with c.on("content"):
        c.draw(f"{pt(x0, y0)} -- {pt(x1, y1)}",
               f"dflattice, line width={num(width)}pt, dash pattern=on 1.2pt off 0.8pt")


def lattice_fringes(c, x, y0, y1, half_w, spacing):
    """Standing-wave lattice drawn as short horizontal fringes along an axis."""
    n = max(int((y1 - y0) / spacing), 1)
    with c.on("content"):
        for i in range(n + 1):
            y = y0 + i * (y1 - y0) / n
            c.draw(f"{pt(x - half_w, y)} -- {pt(x + half_w, y)}", "dflattice, line width=0.7pt, opacity=0.85")


def lens(c, x, y, half_w, thick):
    """Biconvex lens for a vertical beam."""
    with c.on("front"):
        c.ellipse(x, y, half_w, thick, "fill=dfoptics!25, draw=dfoptics, fill opacity=0.85")


def mirror(c, x, y, half_w, facing_up, tip_tilt, t=0.06):
    """Mirror for a vertical beam, reflective face towards the beam."""
    s = 1 if facing_up else -1
    back = y - s * t
    with c.on("front"):
        c.rect(x - half_w, min(y, back), x + half_w, max(y, back), "fill=dfoptics!35, draw=dfoptics")
        c.draw(f"{pt(x - half_w, y)} -- {pt(x + half_w, y)}", "dfoptics!50!black, line width=0.9pt")
        if tip_tilt:
            for side in (-1, 1):
                xa = x + side * (half_w + 0.12)
                c.draw(f"{pt(xa, y - 0.13)} to[bend {'right' if side > 0 else 'left'}=35] {pt(xa, y + 0.13)}",
                       "{Stealth[length=1.1mm]}-{Stealth[length=1.1mm]}, dfoptics!80!black, line width=0.5pt")


def fold(c, x, y, size, d_in, d_out, tip_tilt=False):
    """45-degree mirror turning a beam travelling along ``d_in`` into ``d_out``.

    Directions are unit vectors such as (-1, 0) (travelling left) or (0, -1)
    (travelling down).  The mirror surface is perpendicular to d_out - d_in."""
    nx, ny = d_out[0] - d_in[0], d_out[1] - d_in[1]
    norm = math.hypot(nx, ny)
    tx, ty = -ny / norm, nx / norm  # along the mirror surface
    a, b = (x - tx * size, y - ty * size), (x + tx * size, y + ty * size)
    with c.on("front"):
        c.draw(f"{pt(*a)} -- {pt(*b)}", "dfoptics!50!black, line width=1.6pt")
        if tip_tilt:  # rocking arrow behind the mirror (opposite to its normal)
            bx, by = x - nx / norm * 0.2, y - ny / norm * 0.2
            c.draw(f"{pt(bx - tx * 0.15, by - ty * 0.15)} to[bend left=40] {pt(bx + tx * 0.15, by + ty * 0.15)}",
                   "{Stealth[length=1.2mm]}-{Stealth[length=1.2mm]}, dfoptics!80!black, line width=0.6pt")


# -- atoms -------------------------------------------------------------------

def cloud(c, x, y, r, opacity=1.0):
    with c.on("front"):
        c.shade(f"{pt(x, y)} circle[radius={num(r)}]",
                f"inner color=dfatoms, outer color=dfatoms!25, opacity={num(opacity)}")


def trajectory(c, x, y_start, y_apex, y_end, dx, clouds, r, leg_gap=0.05):
    """Fountain trajectory: up on one side, turn at the apex, down on the other.

    The two legs are ``leg_gap`` apart.  When start == apex the atoms are
    simply dropped (one downward leg)."""
    w = leg_gap / 2
    style = "dfatoms, line width=0.9pt, opacity=0.75, -{Stealth[length=1.5mm]}"
    fountain = y_apex - y_start > 1e-6
    with c.on("front"):
        if fountain:
            c.draw(f"{pt(x + dx - w, y_start)} -- {pt(x + dx - w, y_apex - w)} "
                   f"arc[start angle=180, end angle=0, radius={num(w)}] -- {pt(x + dx + w, y_end)}", style)
        else:
            c.draw(f"{pt(x + dx, y_start)} -- {pt(x + dx, y_end)}", style)
    for i in range(clouds):
        f = (i + 1) / (clouds + 1)
        if fountain:  # clouds on the way up, slowing down towards the apex
            y = y_start + (1 - (1 - f) ** 2) * (y_apex - y_start)
            cloud(c, x + dx - w, y, r, 0.45 + 0.5 * f)
        else:
            cloud(c, x + dx, y_start + f * f * (y_end - y_start), r, 0.45 + 0.5 * f)


def zeeman_slower(c, x0, x1, y, r, s):
    """Zeeman slower: coil windings around the atomic-beam tube, denser towards
    the MOT end (``s``: direction from the MOT towards the oven, +1 or -1)."""
    n = 9
    with c.on("front"):
        for i in range(n):
            f = (i / (n - 1)) ** 0.6  # windings bunch up towards the MOT end
            xw = x1 - f * (x1 - x0) if s < 0 else x0 + f * (x1 - x0)
            c.draw(f"{pt(xw, y - r - 0.03)} -- {pt(xw, y + r + 0.03)}", "dfcoil, line width=1.1pt")


def mot2d(c, x, y, size):
    """2D MOT: transverse cooling beams (one pair in the plane, one pair
    perpendicular to it, drawn as dot/cross) around the atomic beam."""
    with c.on("content"):
        for sy in (-1, 1):
            c.draw(f"{pt(x, y + sy * size)} -- {pt(x, y + sy * size * 0.3)}",
                   "dfcooling, line width=0.7pt, -{Stealth[length=1mm]}")
        for dx, sym in ((-0.11, "dot"), (0.11, "cross")):
            c.circle(x + dx, y, 0.04, "draw=dfcooling, fill=white, line width=0.5pt")
            if sym == "dot":
                c.circle(x + dx, y, 0.012, "fill=dfcooling, draw=none", cmd="fill")
            else:
                c.draw(f"{pt(x + dx - 0.028, y - 0.028)} -- {pt(x + dx + 0.028, y + 0.028)} "
                       f"{pt(x + dx - 0.028, y + 0.028)} -- {pt(x + dx + 0.028, y - 0.028)}",
                       "dfcooling, line width=0.4pt")


def mot(c, x, y, size):
    """Magneto-optical trap: cooling beams converging on an atom cloud."""
    with c.on("content"):
        for ang in (0, 60, 120):
            a = math.radians(ang)
            dx, dy = math.cos(a) * size, math.sin(a) * size
            for s in (-1, 1):
                c.draw(f"{pt(x + s * dx, y + s * dy)} -- {pt(x + s * dx * 0.35, y + s * dy * 0.35)}",
                       "dfcooling, line width=0.7pt, -{Stealth[length=1.1mm]}")
    cloud(c, x, y, size * 0.26)


def camera(c, x, y, facing, size):
    """Camera looking in direction ``facing`` ('left' or 'right'); (x, y) is the lens tip."""
    s = -1 if facing == "left" else 1
    body_x = x - s * size * 0.5
    with c.on("front"):
        c.polyline([(x, y - size * 0.18), (x, y + size * 0.18),
                    (body_x, y + size * 0.3), (body_x, y - size * 0.3)],
                   "fill=dflaser!70, draw=dflaser", closed=True, cmd="filldraw")
        c.rect(body_x - s * size, y - size * 0.42, body_x, y + size * 0.42, "fill=dflaser, draw=dflaser")


def laser_box(c, x0, y0, x1, y1, text, font):
    with c.on("front"):
        c.rect(x0, y0, x1, y1, "fill=dflaser!8, draw=dflaser, rounded corners=1pt")
        c.text((x0 + x1) / 2, (y0 + y1) / 2, text,
               f"align=center, text=dftext, font={{{font}}}, inner sep=1pt")


def marker(c, x, y, n, target=None):
    """Circled step number of the experimental cycle, optionally with a short
    leader to ``target`` (x, y)."""
    with c.on("annotations"):
        name = f"marker{next(_MARKER_IDS)}"
        c.raw(f"\\node[circle, draw=dfmarker, fill=white, line width=0.5pt, inner sep=0.6pt, "
              f"minimum size=3.4mm, font=\\sffamily\\bfseries\\scriptsize, text=dfmarker] ({name}) "
              f"at {pt(x, y)} {{{n}}};")
        if target is not None:
            c.raw(f"\\draw[dfmarker, line width=0.4pt] ({name}) -- {pt(*target)};")


def dot(c, x, y, r=0.022):
    with c.on("annotations"):
        c.circle(x, y, r, "fill=dfleader, draw=none", cmd="fill")


def zoom_box(c, x0, y0, x1, y1, letter):
    with c.on("annotations"):
        c.rect(x0, y0, x1, y1, "draw=dfleader, dashed, rounded corners=2pt, line width=0.5pt", cmd="draw")
        c.text(x0, y1, f"({letter})", "anchor=south west, inner sep=1pt, font=\\sffamily\\scriptsize, text=dfleader")


def cut_line(c, x, y, half_w, letter):
    """Section line (engineering style) marking where a cross-section is taken."""
    style = "dfleader, line width=0.6pt"
    with c.on("annotations"):
        c.draw(f"{pt(x - half_w, y)} -- {pt(x + half_w, y)}", style + ", dash pattern=on 3pt off 1pt on 1pt off 1pt")
        for s in (-1, 1):
            xe = x + s * half_w
            c.draw(f"{pt(xe, y)} -- {pt(xe, y - 0.18)}", style + ", -{Stealth[length=1.2mm]}")
        c.text(x - half_w, y + 0.02, f"({letter})",
               "anchor=south east, inner sep=1pt, font=\\sffamily\\scriptsize, text=dfleader")
