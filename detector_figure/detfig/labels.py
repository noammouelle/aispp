"""Callout labels and step markers attached to named anchors.

Panels register anchors (``c.anchor("node:lower", x, y)``); labels in the
YAML refer to them by name.  Labels are stacked in a column at a fixed x on
each side of the panel and pushed apart vertically so they never overlap,
then connected to their targets with thin leader lines.

Anchor names may carry a height suffix for anything that runs along the
baseline, e.g. ``shield@50`` means "the shield at z = 50 m"; the panel
supplies a resolver for those.
"""

import re

from . import glyphs
from .tikz import pt


class AnchorError(KeyError):
    pass


def resolve(c, name, along=None):
    """Return (x, y) of anchor ``name``; ``along(base, z)`` handles ``base@z``."""
    if name in c.anchors:
        return c.anchors[name]
    if "@" in name and along is not None:
        base, _, z = name.partition("@")
        point = along(base, float(z))
        if point is not None:
            return point
    known = ", ".join(sorted(c.anchors))
    raise AnchorError(f"unknown anchor '{name}'. Known anchors: {known}"
                      + ("; plus <tube|beam|wall|coil|shield|support>@<z in m>" if along else ""))


def spread(desired, heights, gap):
    """1-D label repulsion.

    ``desired`` are target y values, ``heights`` the number of text lines.
    Returns positions in the same order, at least ``gap * lines`` apart and
    as close as possible (in the mean) to where they want to be."""
    order = sorted(range(len(desired)), key=lambda i: desired[i])
    clusters = []  # each: [indices], centred on the mean of their desired y
    for i in order:
        clusters.append([i])
        while len(clusters) > 1:
            a, b = clusters[-2], clusters[-1]
            ya, yb = _layout(a, desired, heights, gap), _layout(b, desired, heights, gap)
            top_a = ya[-1] + gap * heights[a[-1]] / 2
            bottom_b = yb[0] - gap * heights[b[0]] / 2
            if bottom_b >= top_a:
                break
            clusters[-2:] = [a + b]
    out = [0.0] * len(desired)
    for cl in clusters:
        for i, y in zip(cl, _layout(cl, desired, heights, gap)):
            out[i] = y
    return out


def _layout(cluster, desired, heights, gap):
    offs, y = [], 0.0
    for k, i in enumerate(cluster):
        if k:
            y += gap * (heights[cluster[k - 1]] + heights[i]) / 2
        offs.append(y)
    shift = sum(desired[i] - o for i, o in zip(cluster, offs)) / len(cluster)
    return [o + shift for o in offs]


def text_width(text, char=0.155):
    """Rough width (cm) of a label: longest line, at ~footnotesize sans."""
    lines = text.split("\\\\")
    plain = [re.sub(r"\\[a-zA-Z]+|[{}$^_\\]", "", ln) for ln in lines]
    return max(len(p) for p in plain) * char


def callouts(c, labels, columns, gap, font, along=None):
    """Draw YAML label entries ``{text, anchor, side, dx, dy, shift}``.

    ``side`` is left/right (stacked in a column at ``columns[side]``) or
    top/bottom (spread in a row at height ``columns[side]``).  ``shift``
    moves the text along its column/row without moving the target."""
    for side in ("left", "right", "top", "bottom"):
        items = [lab for lab in labels if lab.get("side", "right") == side]
        if not items:
            continue
        if side not in columns:
            raise KeyError(f"labels with side '{side}' need a '{side}' entry in label_column")
        targets = []
        for lab in items:
            if str(lab.get("text", "")).endswith("\\") and not str(lab["text"]).endswith("\\\\"):
                raise ValueError(f"label text {lab['text']!r} ends in a backslash: in YAML flow style "
                                 "({text: ..., ...}) a comma ends the text, so quote it: text: 'a\\,b'")
            x, y = resolve(c, lab["anchor"], along)
            targets.append((x + lab.get("dx", 0.0), y + lab.get("dy", 0.0)))
        vertical = side in ("left", "right")
        if vertical:
            sizes = [lab["text"].count("\\\\") + 1 for lab in items]
            pos = spread([t[1] + lab.get("shift", 0.0) for t, lab in zip(targets, items)], sizes, gap)
        else:
            sizes = [text_width(lab["text"]) + 0.25 for lab in items]
            pos = spread([t[0] + lab.get("shift", 0.0) for t, lab in zip(targets, items)], sizes, 1.0)
        line = columns[side]
        for lab, (tx, ty), p in zip(items, targets, pos):
            glyphs.dot(c, tx, ty)
            if vertical:
                knee = line + (-0.12 if side == "right" else 0.12)
                path = f"{pt(tx, ty)} -- {pt(knee, p)} -- {pt(line, p)}"
                where, anchor = (line, p), ("west" if side == "right" else "east")
                align = "left" if side == "right" else "right"
            else:
                knee = line + (-0.1 if side == "top" else 0.1)
                path = f"{pt(tx, ty)} -- {pt(p, knee)} -- {pt(p, line)}"
                where, anchor, align = (p, line), ("south" if side == "top" else "north"), "center"
                align += ", text height=1.6ex, text depth=0.4ex"  # common baseline along the row
            with c.on("annotations"):
                c.draw(path, "dfleader, line width=0.35pt")
                c.text(*where, lab["text"],
                       f"anchor={anchor}, align={align}, inner sep=1.5pt, text=dftext, font={{{font}}}")


def markers(c, items, along=None):
    for m in items:
        x, y = resolve(c, m["anchor"], along)
        glyphs.marker(c, x + m.get("dx", 0.0), y + m.get("dy", 0.0), m["n"])
