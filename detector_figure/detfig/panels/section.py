"""Panel (c): cross-section of the interferometry region, layer by layer.

The layers come straight from ``section.layers`` (outside -> in), so the
same code draws a round MAGIS-like shield, an octagonal one, two shield
layers, or no support frame at all.  Each layer with a ``label`` gets a
callout; its leader starts on the layer itself at angle ``angle`` (degrees,
default spread automatically from upper right to lower right).
"""

import math

from .. import glyphs as g
from ..labels import callouts
from ..tikz import Canvas, num, pt


def _rot(sides):
    """Polygon rotation that puts a flat edge at the top (and bottom)."""
    return 90 - 180 / sides if sides else 0


def _ring_path(radius, sides):
    rot = _rot(sides)
    if not sides:
        return f"(0,0) circle[radius={num(radius)}]"
    pts = [(radius * math.cos(math.radians(rot + 360 * k / sides)),
            radius * math.sin(math.radians(rot + 360 * k / sides))) for k in range(sides)]
    return " -- ".join(pt(x, y) for x, y in pts) + " -- cycle"


def _edge_radius(radius, sides, angle):
    """Distance from the centre to a polygon's edge in direction ``angle``."""
    if not sides:
        return radius
    rot = _rot(sides)
    seg = 360 / sides
    rel = (angle - rot) % seg - seg / 2
    return radius * math.cos(math.radians(seg / 2)) / math.cos(math.radians(rel))


def draw(cfg):
    sec = cfg.section
    c = Canvas()
    layers = sec.layers
    labelled = [lay for lay in layers if lay.get("label")]
    # labels are evenly spaced top-to-bottom in layer order on each side; each
    # leader starts on its layer at the angle where the layer is level with it
    r_out = max(lay.get("radius", lay.get("size", 1.0)) for lay in layers) if layers else 1.0
    wanted, auto = {}, {}
    for side in ("right", "left"):
        group = [lay for lay in labelled if lay.get("side", "right") == side]
        n_lab = max(len(group) - 1, 1)
        for i, lay in enumerate(group):
            y = r_out * 0.9 * (1 - 2 * i / n_lab) if len(group) > 1 else 0.0
            wanted[id(lay)] = y
            R0 = lay.get("radius", lay.get("size", 1.0))
            a = math.degrees(math.asin(max(-0.9, min(0.9, y / R0))))
            auto[id(lay)] = a if side == "right" else 180 - a
    labels = []
    outer = 0.0
    for lay in layers:
        kind, col = lay["kind"], g.color(lay.get("color", _default_color(lay["kind"])))
        R = lay.get("radius", lay.get("size", 1.0))
        outer = max(outer, R)
        t = lay.get("thickness", 0.06)
        sides = lay.get("sides", 0)
        ang = lay.get("angle", auto.get(id(lay), 0))
        if kind == "square":
            with c.on("walls"):
                c.raw(f"\\filldraw[fill={col}, draw={col}!60!black, even odd rule] "
                      f"{pt(-R, -R)} rectangle {pt(R, R)} {pt(-R + t, -R + t)} rectangle {pt(R - t, R - t)};")
                for sx in (-1, 1):
                    for sy in (-1, 1):
                        c.polyline([(sx * (R - t), sy * (R - t)), (sx * (R - t - 0.25), sy * (R - t)),
                                    (sx * (R - t), sy * (R - t - 0.25))],
                                   f"fill={col}, draw={col}!60!black", closed=True, cmd="filldraw")
            rr = R - t / 2
            a = math.radians(ang)
            k = min(1 / abs(math.cos(a)) if math.cos(a) else 9, 1 / abs(math.sin(a)) if math.sin(a) else 9)
            anchor = (rr * math.cos(a) * k, rr * math.sin(a) * k)
        elif kind == "ring":
            n_layers = lay.get("layers", 1)
            spacing = lay.get("spacing", t * 1.8)
            with c.on("walls"):
                for j in range(n_layers):
                    ro = R - j * spacing
                    c.raw(f"\\filldraw[fill={col}, draw={col}!60!black, even odd rule] "
                          f"{_ring_path(ro, sides)} {_ring_path(ro - t, sides)};")
            rr = _edge_radius(R - t / 2, sides, ang)
            anchor = (rr * math.cos(math.radians(ang)), rr * math.sin(math.radians(ang)))
        elif kind == "coils":
            n = lay.get("count", 16)
            size = lay.get("wire", 0.05)
            with c.on("walls"):
                for k in range(n):
                    a = 2 * math.pi * k / n
                    c.circle(R * math.cos(a), R * math.sin(a), size,
                             f"fill={col}, draw={col}!60!black, line width=0.2pt")
            a = 2 * math.pi * round(math.radians(ang) / (2 * math.pi / n)) / n
            anchor = (R * math.cos(a), R * math.sin(a))
        elif kind == "vacuum":
            with c.on("interiors"):
                c.raw(f"\\shade[inner color=dfinterior, outer color=dfmetal!{num(35 * cfg.style.depth)}!dfinterior] "
                      f"{_ring_path(R, sides)};")
            anchor = (R * 0.8 * math.cos(math.radians(ang)), R * 0.8 * math.sin(math.radians(ang)))
        elif kind == "beam":
            with c.on("content"):
                c.shade(f"(0,0) circle[radius={num(R)}]", f"inner color={col}!90, outer color={col}!0")
            anchor = (R * 0.55 * math.cos(math.radians(ang)), R * 0.55 * math.sin(math.radians(ang)))
        elif kind == "atoms":
            g.cloud(c, 0, 0, R)
            anchor = (R * 0.5 * math.cos(math.radians(ang)), R * 0.5 * math.sin(math.radians(ang)))
        else:
            raise ValueError(f"section.layers: unknown kind '{kind}' "
                             "(use square, ring, coils, vacuum, beam or atoms)")
        if lay.get("label"):
            name = f"layer{len(labels)}"
            c.anchor(name, *anchor)
            labels.append({"text": lay["label"], "anchor": name, "side": lay.get("side", "right"),
                           "shift": wanted[id(lay)] - anchor[1]})

    callouts(c, labels, {"right": sec.label_column, "left": -sec.label_column}, sec.label_gap,
             cfg.document.font)
    if sec.title:
        c.text(0, -outer - 0.15, sec.title,
               f"anchor=north, align=center, inner sep=1pt, text=dftext, font={{{cfg.document.font}}}")
    c.anchor("_letter", -outer - 0.2, outer + 0.35)
    return c


def _default_color(kind):
    return {"square": "support", "ring": "metal", "coils": "coil", "beam": "beam",
            "atoms": "atoms", "vacuum": "interior"}.get(kind, "metal")
