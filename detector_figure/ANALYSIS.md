# Figure analysis: MAGIS-100 and AION-10 → a generic detector schematic

## 0. What could and could not be checked

This environment's network policy blocks `arxiv.org` and every other
academic host tried: ar5iv, export.arxiv, OSTI, IOP, Fermilab's lss archive,
the Stanford lab page, Liverpool's repository and ResearchGate. **The figure
images themselves could not be downloaded or inspected.** What follows
combines:

* facts **confirmed** from web-search result text that quotes the papers
  (marked ✔, with the source),
* the figure-by-figure reading in the original brief (the user's notes), and
* recollection, marked *unverified*.

To redo the extraction properly, allow `arxiv.org` in the environment's
network settings and re-run the analysis (Figs. 5–7, 10–12, 14, 17 of
2104.02835; Figs. 1, 5–8, 15–16 of 2508.03491).

## 1. Confirmed architecture

| Element | MAGIS-100 | AION-10 |
|---|---|---|
| Baseline | ~100 m vertical, MINOS access shaft, Fermilab ✔ | ~10 m tower in the Beecroft building, Oxford: "a gradiometer made of two 5-metre interferometers stacked in one tower" ✔ |
| Atom sources | 3, at top / middle / bottom ✔ | 2, at the **bottom and middle**; "these clouds are launched upwards" ✔ (AION project site) |
| Source → baseline | the **same pair of lattice beams** shuttles atoms horizontally out of the source and launches them vertically, switched by a fibre switch ✔ | atoms "transported through XHV connections into the interconnect chambers" at the base of each 5 m pipe ✔ |
| Detection | 461 nm fluorescence imaging at the nodes ✔ | camera assemblies at the base of each pipe ✔ |
| Species / clock line | Sr, 698 nm ¹S₀–³P₀ ✔; 689 nm narrow-line cooling ✔; source = oven + Zeeman slower + 2D MOT (AOSense) ✔ | ⁸⁷Sr, 698 nm ✔; ovens, 2D and 3D MOTs ✔ |
| Vacuum tube | "support frame containing a 6" diameter vacuum tube" ✔. 17 modules of ~5.3 m, from a 2022 MAGIS shield abstract, not 2104.02835 | modular beam pipe; pipe → heater tape + insulation → coils → shield ✔ |
| Magnetic environment | **single-layer octagonal** mu-metal shield ("four overlapping mu-metal sheets") ✔; **coil bars** along the length giving a 1 G **horizontal** bias ✔ | 10 m magnetic shield ✔ (shape and layering *unverified*) |
| Laser entry | from the top; laser room at the surface; in-vacuum 4f relay lenses down to the top of the shaft ✔ | input arm → fixed steering mirror in a six-way cross → **vertical** beam-conditioning pipe → **pair of mirrors** in the beam-transfer pipe across the top → telescope ✔ |
| Top steering | in-vacuum piezo tip-tilt mirror **before** the telescope (Mad City Labs Nano-MTA2, vacuum version) ✔ | "at least one of the two [transfer] mirrors" dynamically adjusted, together with the bottom mirror, for a tunable pivot point ✔ |
| Telescope | f = 150 mm and 4500 mm (×30), ~300 µm waist at the mirror → ~1 cm waist ✔ | Keplerian, lenses ~2 m apart ✔ |
| Bottom | **tip-tilt retro-reflection mirror** in a vacuum chamber at the bottom; tripod of 3 piezo stages (2 angles + height) ✔ | piezo retro-reflecting mirror on the phase-shear detection platform ✔ |
| Operating modes | (A) top and middle clouds dropped ~50 m; (B) short upward launches detected at the origin ✔ | both clouds launched upwards ✔ |

### Two tip-tilt mirrors, not one

Glick et al. (arXiv:2311.05714) ✔: *"Traditionally … done with a single
piezo-actuated tip-tilt mirror at one end of the baseline, but this approach
breaks down for larger interferometer baselines."* The MAGIS scheme instead
has *"piezo-actuated mirrors on each end of the baseline and tun[es] the path
length of the interferometer beam before a telescope to set the position of
the pivot point."* In more detail:

* **Top mirror.** In vacuum, about 50 mm before the first (f = 150 mm) lens.
  Stepper motors outside the vacuum set the mirror–lens distance, and so the
  pivot point, which is the mirror's image through the ×30 telescope. During
  each shot the piezos rotate the beam for Coriolis compensation.
* **Bottom mirror.** The retro-reflector, *"angled to keep the reflected beam
  aligned with the downgoing beam."*

The bottom mirror is required, not optional. LMT clock sequences need
pulses travelling in both directions, so a single downward beam is not
enough. The figure therefore draws the retro mirror by default, with
tip-tilt arrows on **both** mirrors, and a down/up arrow pair on the beam.

## 2. Corrections to the original brief, and to the first draft of this figure

1. The brief's ASCII layout has a single "tip–tilt / retro-reflection
   mirror" at the bottom. Both experiments have **two** actuated mirrors,
   one at each end. The figure labels both.
2. In MAGIS the top steering mirror sits **before** the telescope, inside the
   vacuum system. Its position relative to the telescope is physics, not
   decoration, because it sets the pivot point. So the figure orders the top
   as: laser → conditioning pipe → steering mirror → telescope → baseline.
3. The brief's eight labels are kept, but the three source-scale ones
   (source, cooling/trapping, transport lattice) move mainly into inset (b),
   where they can be drawn legibly. Panel (a) keeps one "atom source" label.
4. **AION-10 is not "sources at both ends".** The sources are at the bottom
   and the middle, both launch upwards, and the tube continues 5 m above the
   middle node to the telescope (`instrument.tube_extent`). The beam does not
   come straight across the top either: it comes up a vertical conditioning
   pipe beside the tower and is turned through 180° by two mirrors
   (`laser.route: periscope`).
5. **Bias coils are bars, not a solenoid** (MAGIS). They are drawn as
   conductors running along the tube (`coils.style: bars`), not as rows of
   windings, and labelled "bias(-field) coils" rather than "guide-field
   coils".
6. **The MAGIS shield is single-layer.** The first draft drew two layers in
   the cross-section.
7. **Transport and launch lattice** are drawn as two geometries, which is
   right. The MAGIS config's labels add that the two use the same beams, and
   no lattice wavelength is given because none was confirmed.
8. **MAGIS trajectories** now show one documented mode, the ~50 m drops,
   instead of a mixture of modes.

## 3. Design decisions carried into the code

These come from the brief (target realism 5–6/10) and were checked against
the confirmed facts above:

* **(a) Longitudinal cutaway**, physical vertical scale and schematic
  horizontal scale. It shows: sources on sidearms, launch/detection nodes,
  the shielded UHV tube with guide coils, a finite-width beam (Gaussian
  shading) running the full length, two faint fountain trajectories,
  numbered cycle markers ①–④ and a scale dimension.
* **(b) Node inset** (an abstraction of MAGIS Fig. 11 and the AION
  interconnect figures): oven/2D MOT → 3D MOT → horizontal transport
  lattice → vertical launch lattice in the node → camera on a viewport. The
  common beam passes straight through. This makes the point that source and
  transport are *local*, and only the interferometry beam spans the baseline.
* **(c) Cross-section** (the idea of AION Fig. 16): support ⊃ shield ⊃
  coils ⊃ tube ⊃ beam + atoms.
* **Left out on purpose:** pumps, bellows, bakeout, flanges, truss members,
  the laser-locking chain, and tip-tilt stage internals.

## 4. Open items to verify once arXiv is reachable

* Exact source/node heights for MAGIS (the middle one is assumed to be at 50 m).
* Which side each sidearm attaches on, and whether MAGIS sidearms all point
  the same way.
* The AION input-arm height, which side of the tower the conditioning pipe is
  on, and the true telescope length (TDR Fig. 1).
* Shield layer count (both), and the AION tube diameter and module length.
* Camera positions on the nodes (MAGIS Fig. 12; AION Figs. 5–8).
