"""Panel (a): the whole instrument as a longitudinal cutaway.

Vertical positions are physical (metres mapped linearly onto
``overview.height`` cm); horizontal sizes are schematic.  Top to bottom:

    laser -- input/conditioning pipe(s) -- steering mirror(s)
                                                |
                                            telescope
                                                |
    source -- transport arm -- [node] -- camera |
                                                |   tube + coils + shield
                                               ...
    source -- transport arm -- [node] -- camera
                                                |
                                          retro mirror

The laser route is either ``direct`` (one horizontal pipe and a single
fold mirror above the telescope, MAGIS-like) or ``periscope`` (input arm low
down -> fixed steering mirror -> vertical conditioning pipe beside the tower
-> two-mirror transfer pipe across the top, AION-like).
"""

from .. import glyphs as g
from ..labels import callouts, markers
from ..tikz import Canvas, pt

DOWN, UP, LEFT, RIGHT = (0, -1), (0, 1), (-1, 0), (1, 0)


def draw(cfg):
    ins, ov, st = cfg.instrument, cfg.overview, cfg.style
    c = Canvas()
    H, L = ov.height, float(ins.baseline)

    brk = ov.baseline_break
    z_brk = brk.z if brk.show else None

    def y_of(z):
        """Page height of physical height z (m); above a baseline break the
        drawing is shifted up by the break's gap."""
        return z / L * H + (brk.gap if z_brk is not None and z > z_brk else 0.0)

    r, wall, depth = ins.tube.radius, st.wall, st.depth
    hw, hh = ov.node.half_width, ov.node.half_height
    bw = ins.beam.width
    sh, co = ins.shield, ins.coils
    r_coil = r + wall + sh.gap / 2
    r_shield = r + wall + sh.gap + (sh.layers - 1) * sh.spacing  # outermost layer
    nodes = sorted(ins.nodes, key=lambda n: n["z"])
    if not nodes:
        raise ValueError("instrument.nodes must contain at least one node")

    # The vacuum tube runs between the outermost nodes, or further if
    # instrument.tube_extent says so (e.g. AION-10: sources at 0 and 5 m of 10).
    z_lo, z_hi = nodes[0]["z"], nodes[-1]["z"]
    if ins.tube_extent:
        z_lo, z_hi = min(z_lo, ins.tube_extent[0]), max(z_hi, ins.tube_extent[1])
    tube_bot = y_of(z_lo) if z_lo < nodes[0]["z"] else y_of(z_lo) - hh
    tube_top = y_of(z_hi) if z_hi > nodes[-1]["z"] else y_of(z_hi) + hh

    # -- tube segments between/beyond the nodes, with coils and shield --------
    edges = [tube_bot] + [v for n in nodes for v in (y_of(n["z"]) - hh, y_of(n["z"]) + hh)] + [tube_top]
    for y0, y1 in zip(edges[::2], edges[1::2]):
        if y1 - y0 < 0.05:
            continue
        g.vtube(c, 0, y0 - 0.01, y1 + 0.01, r, wall, depth)
        if sh.show:
            for j in range(sh.layers):
                g.shield(c, 0, y0 + 0.08, y1 - 0.08, r + wall, sh.thickness, sh.gap + j * sh.spacing)
        if co.show:
            g.coils(c, 0, y0 + 0.1, y1 - 0.1, r_coil, co.pitch, min(sh.gap * 0.3, 0.035), co.style)
        if ins.support.show:
            xs = r_shield + sh.thickness + ins.support.offset
            with c.on("back"):
                for s in (-1, 1):
                    c.draw(f"{pt(s * xs, y0)} -- {pt(s * xs, y1)}", "dfsupport, line width=0.8pt, dashed")

    # -- nodes, sources, cameras -----------------------------------------------
    sw = ov.source
    for n in nodes:
        _node(c, n, y_of(n["z"]), ov, wall, depth, bw)

    # -- bottom: retro-reflection mirror or blank flange -----------------------
    stub = ov.stub
    y_beam_bot = tube_bot
    if ins.retro_mirror.show:
        ym = tube_bot - stub
        g.vtube(c, 0, ym - 0.01, tube_bot + 0.01, r, wall, depth)
        g.chamber(c, -hw * 0.8, ym - 0.28, hw * 0.8, ym + 0.01, wall * 1.4, depth)
        g.mirror(c, 0, ym - 0.12, r * 1.1, True, ins.retro_mirror.tip_tilt)
        c.anchor("mirror", r * 1.1, ym - 0.15)
        y_beam_bot = ym - 0.12
    else:
        g.end_cap(c, 0, tube_bot, r, wall, up=False)
        c.anchor("flange", r, tube_bot)

    # -- top: telescope, then the laser route ---------------------------------
    y_cur = tube_top
    tel = ins.telescope
    cone = None
    if tel.show:
        r_tel = r * 1.35
        y_t0 = y_cur + stub
        y_t1 = y_t0 + tel.length
        g.vtube(c, 0, y_cur - 0.01, y_t0 + 0.01, r, wall, depth)
        g.chamber(c, -r_tel - wall, y_t0, r_tel + wall, y_t1, wall, depth, chamfer=0.04)
        g.lens(c, 0, y_t0 + 0.12, r * 1.05, 0.05)
        g.lens(c, 0, y_t1 - 0.12, ins.beam.input_width * 1.3, 0.035)
        cone = (y_t0 + 0.12, y_t1 - 0.12)
        c.anchor("telescope", r_tel + wall, (y_t0 + y_t1) / 2)
        y_cur = y_t1
    lz = ins.laser
    if lz.show:
        y_beam_top = _laser(c, cfg, y_cur, y_of, tube_bot)
    else:
        g.end_cap(c, 0, y_cur, r, wall, up=True)
        y_beam_top = y_cur

    # -- beam --------------------------------------------------------------------
    narrow = ins.beam.input_width
    if cone:
        g.vbeam(c, 0, y_beam_bot, cone[0], bw)
        g.beam_cone(c, 0, cone[1], cone[0], narrow, bw)
        g.vbeam(c, 0, cone[1], y_beam_top, narrow)
    else:
        g.vbeam(c, 0, y_beam_bot, y_beam_top, bw)
    arrows_at = ins.beam.arrows_at
    for k, f in enumerate(arrows_at if isinstance(arrows_at, list) else [arrows_at]):
        y_arrows = tube_bot + (tube_top - tube_bot) * f
        g.beam_arrows(c, 0, y_arrows, bw, ins.retro_mirror.show, size=0.5, gap=ins.beam.arrow_gap)
        if k == 0:
            c.anchor("beam", bw * 0.3, y_arrows)

    # -- atoms ---------------------------------------------------------------------
    if ins.atoms.show:
        for t in ins.atoms.trajectories:
            g.trajectory(c, 0, y_of(t["start"]), y_of(t["apex"]), y_of(t["end"]), t.get("dx", 0.0),
                         ins.atoms.clouds, ins.atoms.cloud_radius, ins.atoms.leg_gap)
            if t.get("id"):
                c.anchor(f"apex:{t['id']}", t.get("dx", 0.0), y_of(t["apex"]))

    # -- pivot point of the Coriolis-compensating beam rotation -----------------
    if ins.pivot.show:
        yp = y_of(ins.pivot.z)
        with c.on("front"):
            c.circle(0, yp, 0.07, "draw=dfoptics!50!black, fill=white, line width=0.6pt")
            c.draw(f"{pt(-0.1, yp)} -- {pt(0.1, yp)} {pt(0, yp - 0.1)} -- {pt(0, yp + 0.1)}",
                   "dfoptics!50!black, line width=0.5pt")
        c.anchor("pivot", 0.07, yp)

    # -- scale: a dimension line along the tube, or a ticked axis -------------
    sb = ov.scale_bar
    if sb.show:
        _scale_bar(c, sb, y_of, z_lo, z_hi, L, nodes, hw, sw, wall)

    # -- baseline break: the tube continues, not drawn to scale -----------------
    if z_brk is not None:
        y0 = z_brk / L * H
        half = r_shield + sh.thickness + 0.06
        x_left = min(-half, sb.x - 0.15) if sb.show else -half
        g.baseline_gap(c, x_left, half + 0.02, y0, y0 + brk.gap, half)
        c.anchor("break", half, y0 + brk.gap / 2)

    # -- zoom indicators --------------------------------------------------------
    zm = ov.zoom
    letters = {name: chr(ord("a") + i) for i, name in enumerate(cfg.layout.panels)}
    if cfg.layout.zoom_indicators:
        if zm.node and "node" in letters:
            n = next((n for n in nodes if n["id"] == zm.node), None)
            if n is None:
                raise ValueError(f"overview.zoom.node '{zm.node}' is not one of the node ids")
            y = y_of(n["z"])
            side = n.get("source", "none")
            ext = hw + sw.arm + sw.width + 0.5
            x0 = -ext if side == "left" else -hw - 0.35
            x1 = ext if side == "right" else hw + 0.35
            half = max(hh, sw.height / 2) + 0.12
            g.zoom_box(c, x0, y - half, x1, y + half, letters["node"])
        if "trajectory" in letters and cfg.trajectory.zoom and ins.atoms.trajectories:
            zs = [z for t in ins.atoms.trajectories for z in (t["start"], t["apex"], t["end"])]
            yb0, yb1 = y_of(min(zs)), y_of(max(zs))
            xb = -(r_shield + sh.thickness + 0.1)
            with c.on("annotations"):
                c.draw(f"{pt(xb + 0.08, yb0)} -- {pt(xb, yb0)} -- {pt(xb, yb1)} -- {pt(xb + 0.08, yb1)}",
                       "dfleader, line width=0.5pt, dashed")
                y_tag = yb0 + cfg.trajectory.zoom_at * (yb1 - yb0)
                c.text(xb, y_tag, f"({letters['trajectory']})",
                       "anchor=east, fill=white, inner sep=1pt, font=\\sffamily\\scriptsize, text=dfleader")
        if zm.section_z is not None and "section" in letters:
            g.cut_line(c, 0, y_of(zm.section_z), r_shield + sh.thickness + 0.18, letters["section"])

    # -- labels, markers, legend -------------------------------------------------
    def along(base, z):
        y = y_of(z)
        x = {"tube": (bw / 2 + r) / 2 + 0.02, "beam": bw * 0.3, "wall": r + wall / 2, "coil": r_coil,
             "shield": r_shield + sh.thickness / 2,
             "support": r_shield + sh.thickness + ins.support.offset}.get(base)
        return None if x is None else (x, y)

    callouts(c, ov.labels, ov.label_column, ov.label_gap, cfg.document.font, along)
    markers(c, ov.markers, along)
    if ov.cycle_legend.show:
        _legend(c, ov.cycle_legend, y_beam_bot + ov.cycle_legend.y, cfg.document.font)

    c.anchor("_letter", ov.label_column.left - 1.2, y_beam_top + 0.4)
    return c


def _node(c, n, y, ov, wall, depth, bw):
    """One launch/detection node with its camera and (optional) source sidearm."""
    hw, hh, sw = ov.node.half_width, ov.node.half_height, ov.source
    nid = n["id"]
    g.chamber(c, -hw, y - hh, hw, y + hh, wall * 1.4, depth)
    c.anchor(f"node:{nid}", hw, y + hh * 0.5)
    c.anchor(f"node:{nid}.center", 0, y)
    # launch lattice: two dashed lines hugging the beam inside the node
    for s in (-1, 1):
        g.lattice_line(c, s * (bw / 2 + 0.035), y - hh + wall * 1.4, s * (bw / 2 + 0.035), y + hh - wall * 1.4)
    c.anchor(f"node:{nid}.launch", bw / 2 + 0.035, y - hh * 0.4)
    side = n.get("source", "none")
    cam_side = 1 if side == "left" else -1
    if n.get("camera", True):
        g.camera(c, cam_side * (hw + 0.02), y, "left" if cam_side > 0 else "right", 0.2)
        c.anchor(f"node:{nid}.camera", cam_side * (hw + 0.25), y + 0.08)
    if side not in ("left", "right"):
        return
    s = -1 if side == "left" else 1
    xa, xb = s * hw, s * (hw + sw.arm)
    g.htube(c, min(xa, xb) - 0.01, max(xa, xb) + 0.01, y, sw.arm_radius, wall, depth)
    xc0, xc1 = xb, xb + s * sw.width
    g.chamber(c, min(xc0, xc1), y - sw.height / 2, max(xc0, xc1), y + sw.height / 2,
              wall * 1.4, depth, vertical=False)
    xm = (xc0 + xc1) / 2
    # oven stub feeding the chamber from the far side
    xo = xc1 + s * 0.25
    g.htube(c, min(xc1, xo) - 0.01, max(xc1, xo) + 0.01, y, sw.arm_radius * 0.8, wall, depth)
    with c.on("front"):
        c.rect(min(xo, xo + s * 0.18), y - 0.13, max(xo, xo + s * 0.18), y + 0.13,
               "fill=dfmetal!80!black, draw=dfmetal!50!black")
    g.mot(c, xm, y, min(sw.width, sw.height) * 0.33)
    g.lattice_line(c, xm - s * 0.12, y, 0, y)
    c.anchor(f"source:{nid}", xc1 + s * 0.05, y + sw.height * 0.3)
    c.anchor(f"source:{nid}.mot", xm, y)
    c.anchor(f"source:{nid}.oven", xo + s * 0.09, y)
    c.anchor(f"arm:{nid}", (xa + xb) / 2, y)


def _laser(c, cfg, y_cur, y_of, tube_bot):
    """Pipes, mirrors and laser box feeding the telescope from above.

    Returns the height at which the beam enters the vertical axis."""
    ins, ov, st = cfg.instrument, cfg.overview, cfg.style
    lz, wall, depth = ins.laser, st.wall, st.depth
    s = 1 if lz.side == "right" else -1
    w_in = ins.beam.input_width
    pipe_r = w_in * 1.6
    stub = ov.stub
    y_f = y_cur + stub + 0.2  # height of the top transfer level
    neck = pipe_r if ins.telescope.show else ins.tube.radius  # no telescope: full-width beam
    g.vtube(c, 0, y_cur - 0.01, y_f - 0.2 + 0.01, neck, wall, depth)
    g.chamber(c, -0.26, y_f - 0.22, 0.26, y_f + 0.22, wall * 1.4, depth, chamfer=0.06)
    tip_on_axis = ins.steering_mirror.tip_tilt and (lz.route == "direct" or lz.tip_tilt_on == "axis")
    g.fold(c, 0, y_f, 0.12, (-s, 0), DOWN, tip_on_axis)
    c.anchor("fold", -s * 0.12, y_f + 0.12)

    if lz.route == "direct":
        x_end = s * (0.26 + lz.pipe_length)
        g.htube(c, min(s * 0.25, x_end), max(s * 0.25, x_end), y_f, pipe_r, wall, depth)
        g.hbeam(c, min(0, x_end), max(0, x_end), y_f, w_in)
        c.anchor("pipe", x_end - s * lz.pipe_length / 2, y_f + pipe_r)
        _laser_box(c, cfg, x_end, y_f, s)
        return y_f

    if lz.route != "periscope":
        raise ValueError(f"instrument.laser.route must be 'direct' or 'periscope', not '{lz.route}'")
    # periscope: transfer pipe across the top, vertical conditioning pipe
    # beside the tower, steering mirror in a cross where the input arm enters
    xc = s * lz.offset
    y_in = y_of(lz.input_z) if lz.input_z is not None else tube_bot + 0.5
    g.htube(c, min(s * 0.25, xc), max(s * 0.25, xc), y_f, pipe_r, wall, depth)
    g.hbeam(c, min(0, xc), max(0, xc), y_f, w_in)
    g.chamber(c, xc - 0.22, y_f - 0.22, xc + 0.22, y_f + 0.22, wall * 1.4, depth, chamfer=0.06)
    g.fold(c, xc, y_f, 0.12, UP, (-s, 0), ins.steering_mirror.tip_tilt and lz.tip_tilt_on == "transfer")
    c.anchor("transfer", (xc) / 2, y_f + pipe_r)
    g.vtube(c, xc, y_in + 0.2, y_f - 0.2, pipe_r, wall, depth)
    g.vbeam(c, xc, y_in, y_f, w_in)
    c.anchor("pipe", xc + s * pipe_r, (y_in + y_f) / 2)
    g.chamber(c, xc - 0.2, y_in - 0.2, xc + 0.2, y_in + 0.2, wall * 1.4, depth, chamfer=0.05)
    g.fold(c, xc, y_in, 0.11, (-s, 0), UP, False)
    c.anchor("cross", xc - s * 0.2, y_in - 0.15)
    x_end = xc + s * (0.2 + lz.pipe_length)
    g.htube(c, min(xc + s * 0.19, x_end), max(xc + s * 0.19, x_end), y_in, pipe_r, wall, depth)
    g.hbeam(c, min(xc, x_end), max(xc, x_end), y_in, w_in)
    c.anchor("input", x_end - s * lz.pipe_length / 2, y_in - pipe_r)
    _laser_box(c, cfg, x_end, y_in, s)
    return y_f


def _laser_box(c, cfg, x_end, y, s):
    """Viewport at the end of the pipe, then the (in-air) laser box."""
    lz = cfg.instrument.laser
    pipe_r = cfg.instrument.beam.input_width * 1.6
    with c.on("front"):
        c.rect(x_end - 0.03, y - pipe_r - 0.09, x_end + 0.03, y + pipe_r + 0.09,
               "fill=dfoptics!30, draw=dfoptics")
    bx0, bx1 = x_end + s * 0.12, x_end + s * (0.12 + lz.box_width)
    g.hbeam(c, min(x_end, bx0), max(x_end, bx0), y, cfg.instrument.beam.input_width)
    hb = lz.box_height / 2
    g.laser_box(c, min(bx0, bx1), y - hb, max(bx0, bx1), y + hb, lz.text, cfg.document.font)
    c.anchor("laser", (bx0 + bx1) / 2, y + hb)


def _scale_bar(c, sb, y_of, z_lo, z_hi, L, nodes, hw, sw, wall):
    x = sb.x
    with c.on("annotations"):
        if sb.style == "dimension":
            # stop the arrowheads on a sidearm rather than inside it
            def inset(z):
                crosses = any(n["z"] == z and n.get("source") == "left" for n in nodes) and -hw - sw.arm < x < -hw
                return sw.arm_radius + wall * 1.5 if crosses else 0.0
            c.draw(f"{pt(x, y_of(z_lo) + inset(z_lo))} -- {pt(x, y_of(z_hi) - inset(z_hi))}",
                   "dfleader, line width=0.5pt, {Stealth[length=1.4mm]}-{Stealth[length=1.4mm]}")
            z = z_lo + sb.step
            while sb.step and z < z_hi - 1e-9:
                c.draw(f"{pt(x - 0.05, y_of(z))} -- {pt(x + 0.05, y_of(z))}", "dfleader, line width=0.4pt")
                z += sb.step
            c.text(x, y_of(z_lo + sb.at * (z_hi - z_lo)), _fill_length(sb.label, z_hi - z_lo),
                   "rotate=90, fill=white, inner sep=1.5pt, font=\\sffamily\\scriptsize, text=dftext")
        else:
            c.draw(f"{pt(x, y_of(0))} -- {pt(x, y_of(L))}", "dfleader, line width=0.5pt")
            z = 0.0
            while z <= L + 1e-9:
                c.draw(f"{pt(x, y_of(z))} -- {pt(x - 0.08, y_of(z))}", "dfleader, line width=0.5pt")
                c.text(x - 0.1, y_of(z), f"{z:g}", "anchor=east, inner sep=1pt, font=\\sffamily\\scriptsize, text=dftext")
                z += sb.step or L
            c.text(x - 0.55, y_of(L / 2), _fill_length(sb.label, L),
                   "rotate=90, anchor=south, inner sep=1pt, font=\\sffamily\\scriptsize, text=dftext")


def _fill_length(text, length):
    """Substitute {length} / {length:g} (not str.format: labels contain LaTeX braces)."""
    return text.replace("{length:g}", f"{length:g}").replace("{length}", f"{length:g}")


def _legend(c, lg, y, font):
    """One-line legend of the experimental cycle: (1) prepare -> (2) ... ."""
    with c.on("annotations"):
        prev = None
        for i, step in enumerate(lg.steps, 1):
            where = f"at {pt(lg.x, y)}" if prev is None else f"[right=7pt of {prev}]"
            c.raw(f"\\node[circle, draw=dfmarker, fill=white, line width=0.5pt, inner sep=0.6pt, minimum size=3.4mm, "
                  f"font=\\sffamily\\bfseries\\scriptsize, text=dfmarker, anchor=west] (lgs{i}) {where} {{{i}}};")
            if prev:
                c.raw(f"\\draw[dfleader, line width=0.6pt, -{{Stealth[length=1.5mm]}}] "
                      f"({prev}.east) ++(1pt,0) -- (lgs{i}.west);")
            c.raw(f"\\node[right=1pt of lgs{i}, inner sep=1pt, text height=1.6ex, text depth=0.4ex, "
                  f"text=dftext, font={{{font}}}] (lgt{i}) {{{step}}};")
            prev = f"lgt{i}"
