from . import node, overview, section, trajectory

# name used in layout.panels -> function(cfg) returning a Canvas
PANELS = {"overview": overview.draw, "node": node.draw, "section": section.draw,
          "trajectory": trajectory.draw}
