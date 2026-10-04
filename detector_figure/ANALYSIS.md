# Figure analysis: MAGIS-100 and AION-10 → a generic detector schematic

All facts below were taken from the arXiv LaTeX sources and figure files of
the three papers:

| short | paper | where facts come from |
|---|---|---|
| **M** | Abe et al., *MAGIS-100*, QST 6, 044003 (2021), arXiv:2104.02835 | `main.tex` line numbers, Figs. 1–18 |
| **G** | Glick et al., *Coriolis force compensation and laser beam delivery for 100-m baseline atom interferometry*, arXiv:2311.05714 | `main.tex` line numbers, beam-delivery and pivot-point figures |
| **A** | AION Collaboration, *AION-10 technical design report*, arXiv:2508.03491 | `CDR.tex` line numbers; dimensions read off Fig. 1 (drawing PDR-D03) |

Values that are a calculation of ours, or schematic licence in the
drawing, are marked as such.

## 1. What the original figures show

| figure | content | used for |
|---|---|---|
| M Fig. 4 | spacetime diagram of one clock interferometer, π/2–π–π/2 at 0, T, 2T | panel (d) |
| M Fig. 5 | clock gradiometer: two clouds at z₁ and z₂, light from both ends, extra opposite-direction π pulses for LMT | panel (d), beam arrows |
| M Fig. 6 | operating modes A–D in the 100 m tube (drops; short launches; tall dual-isotope launch) | trajectories |
| M Fig. 7 | (a) shaft, three sources on the **same side**, 100 m arrow source to source; (b) surface optics; (c) one module | panel (a) |
| M Fig. 8 | shield cross-section: tube, square aluminium scaffolding, 8 bias-coil bars, octagonal mu-metal | panel (c) |
| M Figs. 10–12 | atom source CAD; connection node with crossed launch-lattice beams on scaffolding; detection optics | panel (b) |
| M Fig. 14 | interferometry laser chain ending in top tip-tilt → telescope → tube → bottom tip-tilt | top of panel (a) |
| M Fig. 17 | bottom retro-reflection tip-tilt chamber | bottom of panel (a) |
| G beam-delivery figures | laser room → 4f relay → fibre → up 3.5 m → across 1 m → 45° tip-tilt mirror → ×30 telescope | MAGIS top route |
| G pivot-point figures | single-cloud dual-isotope mode; three-cloud gradiometer mode, x–z plots | panel (d) |
| A Fig. 1 | engineering drawing with dimensions (heights, BCP offset, telescope length) | AION panel (a) |
| A beam-pipe layers, shielding, coil figures | octagonal double shield, saddle coils, square frame | AION panel (c) |
| A interconnect and launch-lattice figures | octagonal chamber, sidearm, camera and lattice ports, scaffold ±0.5 m | panel (b) |

## 2. Architecture, per experiment

| element | MAGIS-100 | AION-10 |
|---|---|---|
| baseline | "100 m vertical vacuum pipe", MINOS shaft (M L386); 17 modules of 5.3 m (M L396) | "aligned pair of 5 m baseline atom interferometers in a 10 m vertical vacuum tube" (A L262) |
| sources | top, middle and bottom, 50 m apart (M L367, Fig. 6), all on one side (M Fig. 7a) | interconnect chambers at 0.5 m and 5.5 m above the retro mirror (A Fig. 1) |
| source chain | AOSense beam source: "Sr oven, Zeeman slower, and 2D MOT" → blue 3D MOT → red MOT → crossed 1064 nm dipole trap with evaporation → matter-wave lensing (M L446–463) | "strontium ovens, 2D and 3D magneto-optical traps" (A L264); 2D MOT on an inclined axis, 3D MOT 0.5 m from the beam axis (A Fig. 1, Detail A); no Zeeman slower mentioned |
| into the tube | horizontal lattice shuttle, then a fibre switch sends **the same beams** to the vertical launch lattice (M L476–489) | "transported through XHV connections into the interconnect chambers" (A L264) |
| launch lattice | two beams crossing at a shallow angle, overlapping for ~10 cm, mirrors on scaffolding 50 cm above and below the node (M L486–489, Fig. 11b) | "two interfering laser beams form a moving optical lattice", mirrors on a scaffold ~0.5 m above and below each chamber (A L467) |
| detection | 461 nm fluorescence, cameras on two perpendicular axes through each node centre, in-vacuum lenses (M L508) | camera assemblies at the base of each pipe; horizontal viewports plus two at 45° above (A L262, L448) |
| magnetic | single-layer octagonal shield of four overlapping mu-metal sheets around a square aluminium scaffolding truss; 8 coil bars, 1 G horizontal bias (M L401–409) | double octagonal shield, 3 × 0.8 mm each, 15 mm apart, 4.6 m per pipe (A L1542–1565); saddle coils on a 272 mm circle, horizontal field (A L1630–1652) |
| tube | ~15 cm ID (scaled from M Fig. 8; not stated in the text) | 153 mm (vacuum model, A L1735) |
| beam route | laser room ~10 m away; 4f relay in rough vacuum (~1 Torr) at ~3 m height; 1 m fibre in air; then in UHV up 3.5 m, across 1 m, and down via the 45° tip-tilt mirror (G L372–446) | input arm (~1.6 m, at the lower-chamber level) → steering mirror in a six-way cross → vertical beam-conditioning pipe 0.72 m from the axis → two mirrors in the beam-transfer pipe turn it through 180° → telescope (A L262, L281, L325, Fig. 1) |
| telescope | ×30, f = 150 mm and 4500 mm (4.65 m long), 300 µm → ~1 cm waist; lens 1 is 50 mm below the mirror (G L249, L276) | f = 55 mm and 1947 mm, lenses ~2 m apart, ~1 cm waist (A L339) |
| tip-tilt mirrors | two, "on each end of the baseline" (G L157): the 45° top mirror (G L253) and the bottom retro mirror in its own chamber (M L599–609) | one of the two transfer-pipe mirrors is dynamic, plus the bottom retro mirror on the phase-shear detection platform (A L262, L325) |
| pulse directions | LMT π pulses come "from alternating directions"; the upward ones are the retro-reflected light, selected by Doppler shift (M L327, L553) | — |
| operating modes | A: two clouds dropped together, each falling 50 m in >3 s; B/C: launches of a few metres, detected at the source; D: dual-isotope launch of over 50 m (M L371) | both clouds launched upwards (A L467); heights and T not given |

**Pivot point (G Eq. 1, L246).** The top mirror rotates the beam about an
actuation point in front of lens 1. The telescope images that point to a
pivot point a distance d_p = f₂ + M²f₁ − M²d_a beyond lens 2. The paper's
examples put the pivot at the bottom source (G L205, L313). Our own
thin-lens estimate: with d_a = 50 mm, d_p ≈ 94.5 m. A single mirror at
one end fails at 100 m: the beam walks ~2.2 cm across the baseline, more
than its 1 cm waist (G L142).

## 3. How the brief and the first draft were corrected

1. **Two tip-tilt mirrors**, one at each end (G L157), rather than the single one in the brief's sketch.
2. **The bottom retro mirror is required.** LMT π pulses alternate direction (M L327), and the upward pulses are the retro-reflected light (M L553).
3. **The top steering mirror is a 45° fold before the telescope.** In both experiments the beam reaches it by a periscope: up beside the telescope, then across. This is `laser.route: periscope`.
4. **The source chain is longer than "oven → 3D MOT".** It is oven → (Zeeman slower in MAGIS) → 2D MOT → 3D MOT → narrow-line MOT → dipole trap. Panel (b) draws the stages listed in `node.stages`; MAGIS's label adds "red MOT, dipole trap".
5. **The launch lattice is two crossed beams folded by scaffold mirrors**, not a standing wave along the axis (M Fig. 11b, A L467).
6. **Cross-section nesting.** In MAGIS the square scaffolding sits *inside* the octagonal shield, around the tube, with the 8 coil bars between them. AION has a double octagonal shield inside a square frame, with saddle coils and a heater/insulation layer.
7. **AION-10 geometry.** The sources are at the bottom and middle and both launch upward. The input arm is at the lower-chamber level, the telescope is 2 m, and the beam-transfer pipe sits on top.
8. **MAGIS trajectories show mode A** (two simultaneous 50 m drops); panel (d) uses the same clouds.

## 4. Design of the figure

* **(a) Longitudinal cutaway**, with real vertical positions and a schematic horizontal scale. It shows:
  * sources on sidearms, launch/detection nodes, the shielded tube with bias coils;
  * a finite-width beam with down/up arrow pairs, the two tip-tilt mirrors, and the pivot point (optional);
  * fountain and drop trajectories, and cycle markers ①–④.
  * The generic version breaks the baseline (`overview.baseline_break`), so the same drawing stands for 10 m, 100 m or 1 km.
* **(b) Node inset**: oven → 2D MOT → 3D MOT → transport lattice → crossed launch lattice → camera, with the common beam passing straight through. Source preparation and transport are *local* to each node; only the clock beam spans the baseline.
* **(c) Cross-section**, built layer by layer from YAML.
* **(d) z(t) panel**. MAGIS and AION use the built-in engine: the clouds of panel (a), released at t = 0, with common π/2–π–π/2 pulses and the arms split (exaggerated). The generic figure uses [mz-plots](https://github.com/noammouelle/mz-plots), which draws:
  * order-n LMT pulse ladders with alternating directions (M L327, L553);
  * finite light travel time, so a pulse reaches the upper cloud L/c later;
  * |e⟩/|g⟩ output ports.

  In both engines, the ①–④ brackets under the time axis link "where" in (a) to "when".
* **Left out on purpose:** pumps, bellows, bakeout hardware, flanges, truss members, the laser-locking chain, and tip-tilt stage internals.

## 5. Schematic licence and remaining uncertainty

* AION's sidearms and conditioning pipe are about 90° apart in plan; the drawing puts them on opposite sides.
* AION launch heights and T are not given in the TDR, so its trajectories are illustrative.
* The exact MAGIS source heights are not stated. They are drawn 50 m apart (M Fig. 6), at 0, 50 and 100 m.
* In panel (d), the arm separation is exaggerated, and LMT pulse trains are drawn as single π/2–π–π/2 lines.
