"""Panel (b): one launch/detection node with its atom-source sidearm.

Shows the local part of a shot: atoms are cooled in the source chamber (MOT),
carried along the sidearm by a horizontal transport lattice, handed to the
vertical launch lattice inside the node, launched into the tube, and imaged
by the node cameras when they come back.  The common interferometry beam
passes straight through the node.  The node centre is at (0, 0).
"""

from .. import glyphs as g
from ..labels import callouts, markers
from ..tikz import Canvas, pt


def draw(cfg):
    nd, st, ins = cfg.node, cfg.style, cfg.instrument
    c = Canvas()
    r, wall, depth = nd.tube_radius, nd.wall, st.depth
    hw, hh = nd.half_width, nd.half_height
    s = -1 if nd.source_side == "left" else 1

    # tube above and below, with coils and shield, broken at the far ends
    sh_gap, sh_t = 0.2, 0.07
    for y0, y1 in ((hh, hh + nd.tube_length), (-hh - nd.tube_length, -hh)):
        g.vtube(c, 0, y0 - 0.01, y1 + 0.01, r, wall, depth)
        inner0, inner1 = (y0 + 0.15, y1) if y0 > 0 else (y0, y1 - 0.15)
        if ins.shield.show:
            g.shield(c, 0, inner0, inner1, r + wall, sh_t, sh_gap)
        if ins.coils.show:
            g.coils(c, 0, inner0 + 0.05, inner1 - 0.05, r + wall + sh_gap / 2, 0.22, 0.045)
    y_top, y_bot = hh + nd.tube_length, -hh - nd.tube_length
    edge = r + wall + sh_gap + sh_t + 0.05
    g.break_mark(c, 0, y_top - 0.12, edge)
    g.break_mark(c, 0, y_bot + 0.04, edge)

    # node chamber
    g.chamber(c, -hw, -hh, hw, hh, wall, depth)
    c.anchor("node", hw * 0.75, hh * 0.85)

    # source sidearm: arm, chamber with MOT, oven/2D-MOT stub
    xa, xb = s * hw, s * (hw + nd.arm)
    g.htube(c, min(xa, xb) - 0.01, max(xa, xb) + 0.01, 0, nd.arm_radius, wall, depth)
    sw, sh = nd.source.width, nd.source.height
    xc0, xc1 = xb, xb + s * sw
    g.chamber(c, min(xc0, xc1), -sh / 2, max(xc0, xc1), sh / 2, wall, depth, vertical=False)
    xm = (xc0 + xc1) / 2
    g.mot(c, xm, 0, min(sw, sh) * 0.33)
    c.anchor("mot", xm, -min(sw, sh) * 0.12)
    c.anchor("source", xm - s * sw * 0.3, sh / 2)
    if nd.oven:
        xo = xc1 + s * 0.45
        g.htube(c, min(xc1, xo) - 0.01, max(xc1, xo) + 0.01, 0, nd.arm_radius * 0.7, wall, depth)
        with c.on("front"):
            c.rect(min(xo, xo + s * 0.35), -0.24, max(xo, xo + s * 0.35), 0.24,
                   "fill=dfmetal!80!black, draw=dfmetal!50!black, rounded corners=1pt")
        c.anchor("oven", xo + s * 0.17, -0.24)

    # interferometry beam straight through
    g.vbeam(c, 0, y_bot, y_top, nd.beam_width)
    g.beam_arrows(c, 0, hh + nd.tube_length * 0.55, nd.beam_width, ins.retro_mirror.show, size=0.4)
    c.anchor("beam", nd.beam_width * 0.3, -hh - nd.tube_length * 0.45)

    # transport lattice along the arm, with a cloud on its way in
    g.lattice_line(c, xm, 0, 0, 0, width=0.9)
    xt = (xa + xb) / 2
    g.cloud(c, xt, 0, 0.06, 0.85)
    with c.on("front"):
        c.arrow(xt - s * 0.05, 0.16, xt - s * 0.38, 0.16, "dfatoms, line width=0.5pt")
    c.anchor("transport", xt + s * 0.2, -0.0)

    # launch lattice: fringes along the axis inside the node, cloud launched up
    g.lattice_fringes(c, 0, -hh + wall + 0.08, hh - wall - 0.08, nd.beam_width * 0.75, 0.07)
    g.cloud(c, 0, 0, 0.07)
    with c.on("front"):
        c.arrow(-nd.beam_width * 0.9, -0.05, -nd.beam_width * 0.9, hh * 0.85, "dfatoms, line width=0.6pt")
    c.anchor("launch", nd.beam_width * 0.75, -hh * 0.45)
    # a returning cloud in the tube above (the one that will be imaged)
    g.cloud(c, 0.0, hh + nd.tube_length * 0.25, 0.06, 0.7)
    c.anchor("atoms", 0.06, hh + nd.tube_length * 0.25)

    # cameras on viewports, imaging the node centre
    for side in nd.cameras:
        k = 1 if side == "right" else -1
        if k == s:
            raise ValueError(f"node.cameras: no room for a camera on the source side ({side})")
        xv = k * (hw + 0.02)
        with c.on("front"):
            c.rect(min(xv, xv + k * 0.1), -0.22, max(xv, xv + k * 0.1), 0.22,
                   "fill=dfoptics!25, draw=dfoptics")
        g.camera(c, xv + k * 0.3, 0, "left" if k > 0 else "right", 0.34)
        with c.on("content"):
            c.polyline([(0, 0), (xv, 0.18), (xv, -0.18)], "fill=dfoptics, fill opacity=0.12, draw=none",
                       closed=True, cmd="filldraw")
        c.anchor(f"camera:{side}", xv + k * 0.7, 0.18)
    if nd.cameras:
        c.anchor("camera", *c.anchors[f"camera:{nd.cameras[0]}"])

    def along(base, z):
        # here z is a page offset (cm) above the node centre, not metres
        x = {"tube": (nd.beam_width / 2 + r) / 2, "wall": r + wall / 2, "coil": r + wall + sh_gap / 2,
             "shield": r + wall + sh_gap + sh_t / 2, "beam": nd.beam_width * 0.3}.get(base)
        return None if x is None else (x, z)

    callouts(c, nd.labels, nd.label_column, nd.label_gap, cfg.document.font, along)
    markers(c, nd.markers, along)
    c.anchor("_letter", nd.label_column.left, y_top + 0.3)
    return c
