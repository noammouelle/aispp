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
| Baseline | ~100 m vertical, MINOS access shaft, Fermilab ✔ | ~10 m vertical tower, Beecroft building, Oxford ✔ |
| Atom sources | 3, at top / middle / bottom ✔ | 2 (lower and upper sidearm) ✔ |
| Source → baseline | horizontal lattice shuttle from source into the main tube, then vertical lattice launch ✔ | sidearms → interconnect chambers (from the brief) |
| Species / clock line | Sr, 698 nm ¹S₀–³P₀ ✔; 689 nm narrow-line cooling ✔ | ⁸⁷Sr, 698 nm ✔ |
| Vacuum tube | 6-inch tube ✔; 17 modules of ~5.3 m (*search snippet, possibly a later document*) | modular beam pipe (lengths *unverified*) |
| Magnetic environment | bias coils inside an **octagonal** mu-metal shield ✔ | 10 m magnetic shield ✔ (layering *unverified*) |
| Laser entry | from the top; laser room at the surface ✔ | beam-conditioning pipe → transfer pipe across the top → telescope (from the brief) |
| Top steering | **in-vacuum piezo tip-tilt mirror before the telescope** (Mad City Labs Nano-MTA2, vacuum version) ✔ [2311.05714] | *unverified* |
| Telescope | f = 150 mm and 4500 mm (×30), ~300 µm waist at the mirror → ~1 cm waist ✔ [2311.05714] | *unverified* |
| Bottom | **tip-tilt retro-reflection mirror** in a vacuum chamber at the bottom; tripod of 3 piezo stages (2 angles + height), coarse actuators through bellows ✔ [2104.02835] | bottom mirror on a platform (from the brief) |

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

## 2. Corrections to the original brief

1. The brief's ASCII layout has a single "tip–tilt / retro-reflection
   mirror" at the bottom. MAGIS has **two** actuated mirrors (above).
   The figure labels both.
2. In MAGIS the top steering mirror sits **before** the telescope, inside the
   vacuum system. Its position relative to the telescope is physics, not
   decoration, because it sets the pivot point. So the figure orders the top
   as: laser → conditioning pipe → steering mirror → telescope → baseline.
3. The brief's eight labels are kept, but the three source-scale ones
   (source, cooling/trapping, transport lattice) move mainly into inset (b),
   where they can be drawn legibly. Panel (a) keeps one "atom source" label.

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

* Exact source/node heights for MAGIS (the middle one is assumed to be at 50 m)
  and for AION-10 (both drawn at the ends of the baseline).
* Which side each sidearm attaches on, and whether MAGIS sidearms all point
  the same way.
* The AION input-arm position and the true geometry of the conditioning and
  transfer pipes (TDR Fig. 1).
* Shield layer count (both), and the AION tube diameter and module length.
* Camera positions on the nodes (MAGIS Fig. 12; AION Figs. 5–8).
