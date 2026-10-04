"""Reusable drawing components (vacuum hardware, optics, atoms, markers).

Every function draws onto a :class:`~detfig.tikz.Canvas` in absolute page
coordinates (cm).  Colours are referred to by their style name (``beam``,
``metal``, ...), which :func:`color` turns into the TikZ colour defined in
the document preamble.  The panels compose these glyphs; nothing in here
knows about YAML.
"""

import math

from .tikz import num, opts, pt

PALETTE = set()  # filled by build.py with the names under style.colors


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


def coils(c, x, y0, y1, r, pitch, size):
    """Guide/bias coil windings seen in section: rows of wire cross-sections."""
    n = max(int((y1 - y0) / pitch), 1)
    step = (y1 - y0) / n
    with c.on("back"):
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


def beam_arrows(c, x, y, w, retro, size=0.28):
    """White arrows on the beam: down (incoming) and, with ``retro``, up (reflected)."""
    style = f"-{{Stealth[length={num(size * 0.45)}cm, width={num(max(w * 0.7, 0.05))}cm]}}, white, line width=0.6pt"
    with c.on("front"):
        if retro:
            c.draw(f"{pt(x - w * 0.22, y + size / 2)} -- {pt(x - w * 0.22, y - size / 2)}", style)
            c.draw(f"{pt(x + w * 0.22, y - size / 2)} -- {pt(x + w * 0.22, y + size / 2)}", style)
        else:
            c.draw(f"{pt(x, y + size / 2)} -- {pt(x, y - size / 2)}", style)


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


def fold_mirror(c, x, y, size, beam_from_right, tip_tilt):
    """45-degree mirror turning a horizontal beam downward."""
    s = 1 if beam_from_right else -1
    a, b = (x - s * size, y - size), (x + s * size, y + size)
    with c.on("front"):
        c.draw(f"{pt(*a)} -- {pt(*b)}", "dfoptics!50!black, line width=1.6pt")
        if tip_tilt:
            cx, cy = x - s * 0.24, y + 0.22
            c.draw(f"{pt(cx - s * 0.14, cy - 0.12)} to[bend {'left' if s > 0 else 'right'}=45] {pt(cx + s * 0.12, cy + 0.12)}",
                   "{Stealth[length=1.1mm]}-{Stealth[length=1.1mm]}, dfoptics!80!black, line width=0.5pt")


# -- atoms -------------------------------------------------------------------

def cloud(c, x, y, r, opacity=1.0):
    with c.on("front"):
        c.shade(f"{pt(x, y)} circle[radius={num(r)}]",
                f"inner color=dfatoms, outer color=dfatoms!25, opacity={num(opacity)}")


def trajectory(c, x, y_start, y_apex, y_end, dx, clouds, r):
    """Fountain trajectory: up on one side, turn at the apex, down on the other.

    When start == apex the atoms are simply dropped (no upward leg)."""
    w = 0.025
    with c.on("front"):
        style = "dfatoms, line width=0.5pt, opacity=0.55"
        if y_apex - y_start > 1e-6:
            c.draw(f"{pt(x + dx - w, y_start)} -- {pt(x + dx - w, y_apex - w)} "
                   f"arc[start angle=180, end angle=0, radius={num(w)}] -- {pt(x + dx + w, y_end)}",
                   style + ", -{Stealth[length=1.1mm]}")
        else:
            c.draw(f"{pt(x + dx, y_start)} -- {pt(x + dx, y_end)}", style + ", -{Stealth[length=1.1mm]}")
    for i in range(clouds):
        f = (i + 1) / (clouds + 1)
        cloud(c, x + dx, y_start + f * (y_apex - y_start) if y_apex > y_start + 1e-6
              else y_start + f * (y_end - y_start), r, 0.35 + 0.55 * f)


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


def marker(c, x, y, n):
    """Circled step number of the experimental cycle."""
    with c.on("annotations"):
        c.raw(f"\\node[circle, draw=dfmarker, fill=white, line width=0.5pt, inner sep=0.6pt, "
              f"minimum size=3.4mm, font=\\sffamily\\bfseries\\scriptsize, text=dfmarker] at {pt(x, y)} {{{n}}};")


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
        c.text(x + half_w, y, f"({letter})", "anchor=west, inner sep=1.5pt, font=\\sffamily\\scriptsize, text=dfleader")
