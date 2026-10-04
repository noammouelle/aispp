# detfig: YAML → TikZ schematics of long-baseline atom interferometers

Draws a three-panel overview figure of a vertical atom-interferometric
detector (MAGIS-100, AION-10, or a generic one) as a functionally annotated
cutaway. The panels are:

* **(a)** the whole instrument, with the vertical axis drawn to scale;
* **(b)** one source, launch and detection node;
* **(c)** a cross-section of the interferometry region.

The output is plain TikZ, so fonts match the thesis and everything stays
vector. See [ANALYSIS.md](ANALYSIS.md) for the reasoning behind the design
and for which numbers are confirmed and which are not.

```
pip install pyyaml            # plus a TeX installation with TikZ (pdflatex)
cd detector_figure
python -m detfig configs/generic.yaml -o out
python -m detfig configs/magis100.yaml -o out --set overview.height=12
python -m detfig configs/aion10.yaml -o out
```

Each run writes these files to `out/`:

| file | use |
|---|---|
| `<name>.tex` | standalone document |
| `<name>.pdf` / `.png` | compiled figure / preview (needs `pdflatex`, `pdftoppm`) |
| `<name>-body.tex` | only the picture, for `\input{}` in the thesis (needs `\usepackage{tikz}` and `\usetikzlibrary{arrows.meta,positioning}`) |

In the thesis, `\includegraphics[width=\textwidth]{out/generic_detector.pdf}`
is usually the simplest option.

## Configuration

`configs/defaults.yaml` lists **every** option, with a comment on each. A
figure config holds only what differs from the defaults, and it can build
on another config:

```yaml
extends: generic.yaml            # start from the thesis figure
output: {name: my_variant}
instrument:
  baseline: 10
  nodes:                         # lists replace, dicts merge
    - {id: lower, z: 0, source: left}
    - {id: upper, z: 10, source: right}
  retro_mirror: {tip_tilt: false}
style:
  colors: {beam: "B2182B"}       # one colour changed, the rest kept
```

Mistakes stop the run with a message rather than being silently ignored.
This covers misspelt keys (including keys inside list items such as labels
and nodes), invalid choices (`laser.side: rigth`), duplicate node ids, node
heights outside the baseline, and unquoted `#` colours:
`unknown key 'config.instrument.shiled' (did you mean 'shield'?)`.

### Instrument layouts

| option | what it changes |
|---|---|
| `instrument.nodes` | number and height of launch/detection nodes; `source: left/right/none` per node |
| `instrument.tube_extent: [0, 10]` | tube continues beyond the outermost nodes (AION-10: nodes at 0 and 5 m of a 10 m tower) |
| `instrument.laser.route` | `direct`: one pipe and fold mirror above the telescope (MAGIS-like). `periscope`: input arm → steering mirror in a cross → vertical conditioning pipe → two-mirror transfer pipe (AION-like) |
| `steering_mirror.tip_tilt`, `retro_mirror.{show,tip_tilt}` | the two actuated mirrors |
| `coils.style` | `bars` (conductors along the tube, transverse bias) or `turns` (solenoid) |
| `shield.layers`, `section.layers[].sides` | number of shield layers; circular or polygonal (octagonal for MAGIS) cross-section |
| `atoms.trajectories` | fountains (`start < apex`) or drops (`start == apex`), in metres |
| `node.below` | `tube` (shielded tube continues below the node in (b)) or `stub` (bare stub to the retro mirror) |

**Units.** Heights along the baseline (`z`, `baseline`, `apex`,
`section_z`, `anchor@z`) are in metres. All other sizes are page
centimetres. The horizontal direction is schematic: a 15 cm tube on a
100 m tower would be invisible if drawn to scale.

**Text** is LaTeX. Use `\\` for a line break. In YAML's `{...}` flow style
a comma ends the value, so quote any text that contains one:
`{text: 'launch\,/\,detection', ...}`. The tool raises an error if you
forget.

### Labels and markers

Labels point at named **anchors**, not at coordinates, so moving a part
moves its label with it:

```yaml
overview:
  labels:
    - {text: magnetic shield, anchor: shield@55, side: right}
    - {text: atom source, anchor: source:upper.oven, side: top, dy: 0.1}
  markers:                       # circled cycle-step numbers
    - {n: 3, anchor: tube@30, dx: 0.4}
```

* `side`: `left`/`right` stacks labels in a column at
  `label_column.left/right`; `top`/`bottom` puts them in a row at
  `label_column.top/bottom`. Labels push each other apart automatically.
* `dx`, `dy`: move the point the leader line touches.
* `shift`: move only the text along its column or row.
* `line`: put this label in its own column or row, at this x (left/right)
  or y (top/bottom) instead of the shared `label_column` value.
* Markers take `leader: true` to draw a short line from the circled
  number to its anchor.

Placement is per config. After moving sources to the other side, or
changing the laser route, re-check the labels and `layout` offsets in the
rendered PNG; the geometry follows automatically, but text placement
does not.

| panel | anchors |
|---|---|
| (a) overview | `node:<id>`, `node:<id>.center`, `.launch`, `.camera`, `source:<id>`, `source:<id>.mot`, `.oven`, `arm:<id>`, `apex:<trajectory id>`, `telescope`, `fold` (mirror above the telescope), `laser`, `pipe`, `mirror` (retro) or `flange`, `beam`; periscope route adds `transfer`, `cross`, `input`; plus `tube@z`, `beam@z`, `wall@z`, `coil@z`, `shield@z`, `support@z` |
| (b) node | `node`, `mot`, `source`, `oven`, `transport`, `launch`, `atoms`, `beam`, `camera`, `camera:right`, plus `tube@y` etc., where *y* is a page offset (cm) above the node centre |

An unknown anchor name raises an error that lists all valid anchors.

### Cross-section layers

Panel (c) is built entirely from YAML, listed from the outside in:

```yaml
section:
  layers:
    - {kind: square, size: 1.45, thickness: 0.07, label: structural support}
    - {kind: ring, radius: 1.15, thickness: 0.06, sides: 8, color: shield, layers: 2, label: magnetic shield}
    - {kind: coils, radius: 0.78, count: 18, label: guide-field coils}
    - {kind: ring, radius: 0.62, color: metal, label: UHV tube}
    - {kind: vacuum, radius: 0.55}
    - {kind: beam, radius: 0.36, label: clock beam}
    - {kind: atoms, radius: 0.08, label: atom clouds}
```

`sides: 0` gives a circle; `sides: 8` gives an octagon, as in the MAGIS shield.

## Code map

```
detfig/
  config.py      YAML loading, `extends`, deep merge, typo check
  tikz.py        Canvas: TikZ commands in ordered layers
  glyphs.py      tubes, chambers, beam, lenses, mirrors, MOT, camera, ...
  labels.py      anchors, collision-free callouts, cycle markers
  panels/        overview.py (a), node.py (b), section.py (c)
  build.py       assembles panels, writes .tex, runs pdflatex
```

The cutaway look comes from the layer order in `tikz.py`: all metal walls
are drawn first, then all vacuum interiors. Wherever a sidearm or tube
meets a chamber, the interior is therefore drawn open, without any
port-cutting code. To add a new element, write a glyph and call it from a
panel. To add a new panel, write a `draw(cfg) -> Canvas` function, register
it in `panels/__init__.py`, and add its defaults to `defaults.yaml`.

Rendered previews of all three configs are in `previews/`.

Tests: `python -m pytest tests` (the compile tests are skipped without
`pdflatex`).
