"""Panel (a): the whole instrument as a longitudinal cutaway.

Vertical positions are physical (metres mapped linearly onto
``overview.height`` cm); horizontal sizes are schematic.  Top to bottom:

    laser box -- conditioning pipe -- fold/steering mirror
                                          |
                                      telescope
                                          |
    source -- transport arm -- [node] -- camera
                                          |   tube + coils + shield
                                         ...
    source -- transport arm -- [node] -- camera
                                          |
                                    retro mirror (optional)
"""

from .. import glyphs as g
from ..labels import callouts, markers
from ..tikz import Canvas, num, pt


def draw(cfg):
    ins, ov, st = cfg.instrument, cfg.overview, cfg.style
    c = Canvas()
    H, L = ov.height, float(ins.baseline)

    def y_of(z):
        return z / L * H

    r, wall, depth = ins.tube.radius, st.wall, st.depth
    hw, hh = ov.node.half_width, ov.node.half_height
    bw = ins.beam.width
    sh, co = ins.shield, ins.coils
    r_coil = r + wall + sh.gap / 2
    r_shield = r + wall + sh.gap
    nodes = sorted(ins.nodes, key=lambda n: n["z"])
    if not nodes:
        raise ValueError("instrument.nodes must contain at least one node")
    y_bot, y_top = y_of(nodes[0]["z"]), y_of(nodes[-1]["z"])

    # -- main tube between nodes, with coils and shield ------------------------
    for a, b in zip(nodes, nodes[1:]):
        y0, y1 = y_of(a["z"]) + hh, y_of(b["z"]) - hh
        g.vtube(c, 0, y0 - 0.01, y1 + 0.01, r, wall, depth)
        if sh.show:
            g.shield(c, 0, y0 + 0.08, y1 - 0.08, r + wall, sh.thickness, sh.gap)
        if co.show:
            g.coils(c, 0, y0 + 0.1, y1 - 0.1, r_coil, co.pitch, min(sh.gap * 0.3, 0.035))
        if ins.support.show:
            xs = r_shield + sh.thickness + ins.support.offset
            with c.on("back"):
                for s in (-1, 1):
                    c.draw(f"{pt(s * xs, y0)} -- {pt(s * xs, y1)}", "dfsupport, line width=0.8pt, dashed")

    # -- nodes, sources, cameras -----------------------------------------------
    sw = ov.source
    for n in nodes:
        y, nid = y_of(n["z"]), n["id"]
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
        if side in ("left", "right"):
            s = -1 if side == "left" else 1
            xa, xb = s * hw, s * (hw + sw.arm)
            g.htube(c, min(xa, xb) - 0.01, max(xa, xb) + 0.01, y, sw.arm_radius, wall, depth)
            xc0, xc1 = xb, xb + s * sw.width
            g.chamber(c, min(xc0, xc1), y - sw.height / 2, max(xc0, xc1), y + sw.height / 2,
                      wall * 1.4, depth, vertical=False)
            xm = (xc0 + xc1) / 2
            # oven / 2D-MOT stub feeding the chamber from the far side
            xo = xc1 + s * 0.25
            g.htube(c, min(xc1, xo) - 0.01, max(xc1, xo) + 0.01, y, sw.arm_radius * 0.8, wall, depth)
            with c.on("front"):
                c.rect(min(xo, xo + s * 0.18), y - 0.13, max(xo, xo + s * 0.18), y + 0.13,
                       "fill=dfmetal!80!black, draw=dfmetal!50!black")
            g.mot(c, xm, y, min(sw.width, sw.height) * 0.33)
            g.lattice_line(c, xm + s * -0.12, y, -s * 0.0, y)
            c.anchor(f"source:{nid}", xc1 + s * 0.05, y + sw.height * 0.3)
            c.anchor(f"source:{nid}.mot", xm, y)
            c.anchor(f"source:{nid}.oven", xo + s * 0.09, y)
            c.anchor(f"arm:{nid}", (xa + xb) / 2, y)

    # -- bottom: retro-reflection mirror or blank flange -----------------------
    stub = ov.stub
    y_beam_bot = y_bot - hh
    if ins.retro_mirror.show:
        ym = y_bot - hh - stub
        g.vtube(c, 0, ym - 0.01, y_bot - hh + 0.01, r, wall, depth)
        g.chamber(c, -hw * 0.8, ym - 0.28, hw * 0.8, ym + 0.01, wall * 1.4, depth)
        g.mirror(c, 0, ym - 0.12, r * 1.1, True, ins.retro_mirror.tip_tilt)
        c.anchor("mirror", r * 1.1, ym - 0.15)
        y_beam_bot = ym - 0.12
    else:
        g.end_cap(c, 0, y_bot - hh, r, wall, up=False)
        c.anchor("mirror", r, y_bot - hh)

    # -- top: telescope, fold/steering mirror, conditioning pipe, laser --------
    y_beam_top = y_top + hh
    beam_top_parts = []
    y_cur = y_top + hh
    tel = ins.telescope
    r_tel = r * 1.35
    if tel.show:
        y_t0 = y_cur + stub
        y_t1 = y_t0 + tel.length
        g.vtube(c, 0, y_cur - 0.01, y_t0 + 0.01, r, wall, depth)
        g.chamber(c, -r_tel - wall, y_t0, r_tel + wall, y_t1, wall, depth, chamfer=0.04)
        g.lens(c, 0, y_t0 + 0.12, r * 1.05, 0.05)
        g.lens(c, 0, y_t1 - 0.12, ins.beam.input_width * 1.3, 0.035)
        beam_top_parts.append(("cone", y_t0 + 0.12, y_t1 - 0.12))
        c.anchor("telescope", r_tel + wall, (y_t0 + y_t1) / 2)
        y_cur = y_t1
    lz = ins.laser
    if lz.show:
        s = 1 if lz.side == "right" else -1
        y_f = y_cur + stub + 0.2
        g.vtube(c, 0, y_cur - 0.01, y_f - 0.2 + 0.01, ins.beam.input_width * 1.6, wall, depth)
        g.chamber(c, -0.26, y_f - 0.22, 0.26, y_f + 0.22, wall * 1.4, depth, chamfer=0.06)
        pipe_r = ins.beam.input_width * 1.6
        x_end = s * (0.26 + lz.pipe_length)
        g.htube(c, min(s * 0.25, x_end), max(s * 0.25, x_end), y_f, pipe_r, wall, depth)
        g.hbeam(c, min(0, x_end), max(0, x_end), y_f, ins.beam.input_width)
        g.fold_mirror(c, 0, y_f, 0.12, s > 0, ins.steering_mirror.tip_tilt)
        bx0, bx1 = (x_end, x_end + s * 1.5)
        g.laser_box(c, min(bx0, bx1), y_f - 0.32, max(bx0, bx1), y_f + 0.32, lz.text, cfg.document.font)
        c.anchor("fold", -s * 0.12, y_f + 0.12)
        c.anchor("pipe", x_end - s * lz.pipe_length / 2, y_f + pipe_r)
        c.anchor("laser", (bx0 + bx1) / 2, y_f + 0.32)
        y_beam_top = y_f
    else:
        g.end_cap(c, 0, y_cur, r, wall, up=True)
        y_beam_top = y_cur

    # -- beam --------------------------------------------------------------------
    narrow = ins.beam.input_width
    if beam_top_parts:
        _, y_wide, y_narrow = beam_top_parts[0]
        g.vbeam(c, 0, y_beam_bot, y_wide, bw)
        g.beam_cone(c, 0, y_narrow, y_wide, narrow, bw)
        g.vbeam(c, 0, y_narrow, y_beam_top, narrow)
    else:
        g.vbeam(c, 0, y_beam_bot, y_beam_top, bw)
    y_arrows = y_bot + (y_top - y_bot) * ins.beam.arrows_at
    g.beam_arrows(c, 0, y_arrows, bw, ins.retro_mirror.show, size=0.45)
    c.anchor("beam", bw * 0.3, y_arrows)

    # -- atoms ---------------------------------------------------------------------
    if ins.atoms.show:
        for t in ins.atoms.trajectories:
            g.trajectory(c, 0, y_of(t["start"]), y_of(t["apex"]), y_of(t["end"]), t.get("dx", 0.0),
                         ins.atoms.clouds, 0.045)

    # -- scale bar: a dimension line between the end nodes, or a ticked axis ------
    sb = ov.scale_bar
    if sb.show:
        x = sb.x
        with c.on("annotations"):
            if sb.style == "dimension":
                c.draw(f"{pt(x, y_bot)} -- {pt(x, y_top)}",
                       "dfleader, line width=0.5pt, {Stealth[length=1.4mm]}-{Stealth[length=1.4mm]}")
                z = sb.step
                while z < L - 1e-9:
                    y = y_of(z)
                    c.draw(f"{pt(x - 0.05, y)} -- {pt(x + 0.05, y)}", "dfleader, line width=0.4pt")
                    z += sb.step
                c.text(x, y_of(sb.at * L), sb.label,
                       "rotate=90, fill=white, inner sep=1.5pt, font=\\sffamily\\scriptsize, text=dftext")
            else:
                c.draw(f"{pt(x, y_of(0))} -- {pt(x, y_of(L))}", "dfleader, line width=0.5pt")
                z = 0.0
                while z <= L + 1e-9:
                    y = y_of(z)
                    c.draw(f"{pt(x, y)} -- {pt(x - 0.08, y)}", "dfleader, line width=0.5pt")
                    c.text(x - 0.1, y, f"{z:g}", "anchor=east, inner sep=1pt, font=\\sffamily\\scriptsize, text=dftext")
                    z += sb.step
                c.text(x - 0.55, y_of(L / 2), sb.label,
                       "rotate=90, anchor=south, inner sep=1pt, font=\\sffamily\\scriptsize, text=dftext")

    # -- zoom indicators --------------------------------------------------------
    zm = ov.zoom
    letters = _panel_letters(cfg)
    if cfg.layout.get("zoom_indicators", True):
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
        _legend(c, ov.cycle_legend, y_bot - hh + ov.cycle_legend.y, cfg.document.font)

    c.anchor("_letter", ov.label_column.left - 1.2, max(y_beam_top, y_top + hh) + 0.4)
    return c


def _panel_letters(cfg):
    return {name: chr(ord("a") + i) for i, name in enumerate(cfg.layout.panels)}


def _legend(c, lg, y, font):
    """One-line legend of the experimental cycle: (1) prepare -> (2) ... ."""
    steps = lg.steps
    with c.on("annotations"):
        prev = None
        for i, step in enumerate(steps, 1):
            where = f"at {pt(lg.get('x', -2.6), y)}" if prev is None else f"[right=5pt of {prev}]"
            arrow = "" if prev is None else f"\\draw[dfleader, -{{Stealth[length=1.1mm]}}] ({prev}.east) ++(1pt,0) -- ++(3.5pt,0);"
            c.raw(f"\\node[circle, draw=dfmarker, fill=white, line width=0.5pt, inner sep=0.6pt, minimum size=3.4mm, "
                  f"font=\\sffamily\\bfseries\\scriptsize, text=dfmarker, anchor=west] (lgs{i}) {where} {{{i}}};")
            if arrow:
                c.raw(arrow)
            c.raw(f"\\node[right=1pt of lgs{i}, inner sep=1pt, text=dftext, font={{{font}}}] (lgt{i}) {{{step}}};")
            prev = f"lgt{i}"
