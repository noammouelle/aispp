"""Minimal TikZ writer.

Everything the panels draw goes through a :class:`Canvas`.  A canvas keeps
its commands in named *layers* that are emitted in a fixed order, so a
panel can draw a chamber wall after the atoms and still have the wall end
up underneath them.  The layer order is what makes the "cutaway" look work:
all metal walls are drawn first, then all vacuum interiors (which therefore
punch open every port where two chambers meet), then the physics content
(beam, atoms, lattices), then annotations.
"""

from contextlib import contextmanager

LAYERS = ("back", "walls", "interiors", "content", "front", "annotations")


def num(x):
    """Format a length in cm compactly: 1.250 -> '1.25', 2.0 -> '2'."""
    s = f"{x:.3f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


def pt(x, y):
    return f"({num(x)},{num(y)})"


def opts(*parts):
    """Join TikZ option fragments, skipping empty ones."""
    return ", ".join(p for p in parts if p)


class Canvas:
    def __init__(self):
        self._layers = {name: [] for name in LAYERS}
        self.layer = "content"
        self.anchors = {}  # name -> (x, y), used by labels and markers

    # -- bookkeeping ---------------------------------------------------------
    @contextmanager
    def on(self, layer):
        """Temporarily switch the active layer: ``with c.on("walls"): ...``."""
        old, self.layer = self.layer, layer
        try:
            yield self
        finally:
            self.layer = old

    def anchor(self, name, x, y):
        self.anchors[name] = (x, y)

    def raw(self, line):
        self._layers[self.layer].append(line)

    def comment(self, text):
        self.raw(f"% {text}")

    def tikz(self):
        out = []
        for name in LAYERS:
            if self._layers[name]:
                out.append(f"% --- layer: {name}")
                out.extend(self._layers[name])
        return "\n".join(out)

    # -- primitives ----------------------------------------------------------
    def draw(self, path, style=""):
        self.raw(f"\\draw[{style}] {path};")

    def fill(self, path, style=""):
        self.raw(f"\\fill[{style}] {path};")

    def filldraw(self, path, style=""):
        self.raw(f"\\filldraw[{style}] {path};")

    def shade(self, path, style=""):
        self.raw(f"\\shade[{style}] {path};")

    def text(self, x, y, text, style=""):
        self.raw(f"\\node[{style}] at {pt(x, y)} {{{text}}};")

    def rect(self, x0, y0, x1, y1, style="", cmd="filldraw"):
        getattr(self, cmd)(f"{pt(x0, y0)} rectangle {pt(x1, y1)}", style)

    def polyline(self, points, style="", closed=False, cmd="draw"):
        path = " -- ".join(pt(x, y) for x, y in points)
        if closed:
            path += " -- cycle"
        getattr(self, cmd)(path, style)

    def circle(self, x, y, r, style="", cmd="filldraw"):
        getattr(self, cmd)(f"{pt(x, y)} circle[radius={num(r)}]", style)

    def ellipse(self, x, y, rx, ry, style="", cmd="filldraw"):
        getattr(self, cmd)(f"{pt(x, y)} ellipse[x radius={num(rx)}, y radius={num(ry)}]", style)

    def arrow(self, x0, y0, x1, y1, style=""):
        self.draw(f"{pt(x0, y0)} -- {pt(x1, y1)}", opts("-{Stealth[length=1.6mm]}", style))
