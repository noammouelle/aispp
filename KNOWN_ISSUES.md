# KNOWN_ISSUES.md — ais++ v0.0.2

Open problems that are understood but not fixed in this release. Each one is
pinned by a test asserting the *current* behaviour, so the suite stays green
while the issue is open; fixing an issue will make its test fail, which is the
signal to update the expectation.

---

## 1. `GetDelPhi` returns zero for every beam type

**Severity: high — this is a silent physics gap, not a missing convenience.**

`AISLaserBeam::GetDelPhi` is supposed to return the local correction to the
plane-wave wavevector from the beam's curved wavefront, i.e. the gradient of the
spatially varying phase with the plane-wave `k·z` term removed. It is what gives
the atom the wavefront-gradient contribution to the photon recoil.

It currently returns `{0, 0, 0}` in all cases:

- **`gaussian` and `confocal`** dispatch to `gaussianGradientWavefront`
  (`03_aisbeams/src/AISWavefronts.cc`), which computes the gradient into a local
  variable and then unconditionally returns `{0.0, 0.0, 0.0}`, with the comment
  *"wavefront curvature accounted for in the zernike coeffs in this branch"* —
  a leftover from the `zernike-aberration-only` line of work, where the curvature
  was meant to come from Zernike coefficients instead of the analytic beam.

  The consequence is that commit `d89bbbc`, *"Wire GetDelPhi to the existing
  gaussianGradientWavefront for Gaussian beams"*, is a no-op in practice: it
  wired the call site to a function that always returns zero.

- **`interpolated`** falls through to the same zero return, because there is no
  gradient of the sampled field. This is the piece deliberately deferred in
  v0.0.2.

- **Zernike aberrations** have the same gap in the other direction: `GetPhi`
  applies them (restored in v0.0.2) but `GetDelPhi` never did, so an aberrated
  beam contributes its phase but not its gradient.

**What this means for results.** Any simulation whose signal depends on
wavefront-gradient recoil is currently computing that term as zero. Phase
observables that depend only on `GetPhi` — which is most of the wavefront-
aberration work — are unaffected.

**Fixing it.** Three separate pieces of work:

1. *Analytic beams.* Return `gradient` instead of the zero literal in
   `gaussianGradientWavefront`. Before doing so, resolve the double-counting the
   original comment was worried about: if a user supplies Zernike defocus *and*
   the analytic beam supplies curvature, the curvature is counted twice. The
   likely answer is that they are alternative descriptions and the input should
   pick one, but that needs a decision rather than a patch.
2. *Zernike gradient.* `GradZmn` exists in `03_aisbeams` and implements Noll's
   recursive formula for the Cartesian derivatives. It is currently uncalled.
   An earlier attempt to use it (`zernike-aberration-only-analytical-diff-recoil`)
   was abandoned at "Added modifications to avoid division by zero errors.
   Unsuccessful for now", so expect the ρ → 0 limit to need care.
3. *Interpolated beams.* `AISLaserBeam::GetPhiWrapper` and the `gsl_deriv.h`
   include are retained for exactly this: numerical differentiation of `GetPhi`.
   Note that trilinear interpolation is only C⁰, so the differentiated gradient
   is discontinuous across cell boundaries. If the gradient feeds the dynamics,
   that discontinuity lands directly in the equations of motion, and a smoother
   interpolant (tricubic, or a fitted local polynomial) is probably needed before
   the result is trustworthy.

**Pinned by:** `tests/test_interpolation.cc`, section *"Analytic beam types are
unaffected"*.

---

## 2. Sampled beams extrapolate as a constant outside the grid

`trilinearInterpolation` clamps both the cell index and the interpolation
weights, so a point outside the sampled volume takes the nearest edge value
rather than diverging. This is a deliberate improvement over the previous
behaviour, where a coordinate below `grid[0]` selected the *top* cell and
produced a large negative weight — but constant extrapolation is still a
modelling choice, not physics.

If a meaningful fraction of the ensemble leaves the grid, the beam is not being
modelled where those atoms are. There is currently no diagnostic for how many
atoms sample outside the grid; adding a counter would make this visible instead
of silent.

`aisoptics` has a `bounds_policy` on the Python side (`raise` by default) with no
C++ counterpart. Aligning the two would be the natural fix.

**Pinned by:** `tests/test_interpolation.cc`, section *"Out-of-grid sampling
clamps at both ends"*.

---

## 3. Trilinear interpolation is C⁰

Independently of issue 1, the interpolant itself is only continuous, not
continuously differentiable. Phases are fine; anything derived from a gradient
will show cell-boundary artefacts. Relevant if and when issue 1 is addressed for
sampled beams.

---

## 4. `usestaticapprox` is obsolete but still live

`AISPulsePropagator.cc` still branches on it and carries its own
`// TODO: Remove. This is now obsolete.` It is documented in
`INPUT_REFERENCE.md` as obsolete-but-accepted. Removing it is a breaking input
change, so it is deferred rather than done quietly.

---

## 5. `zlaser` has no effect

Parsed, stored on `AISLaserBeam`, never read. Beam geometry comes entirely from
`focallength` and the sign of `kz`. A remote branch is literally named after the
discovery (*"Ignores laser position now (useless)"*). It is retained for input
compatibility and documented as inert in `INPUT_REFERENCE.md`.

Unlike `beamradius` and `baseline` — which were inert for the same reason and are
live again in v0.0.2 now that the Zernike block is restored — `zlaser` has no
remaining consumer, so the choice is between wiring it to something meaningful or
removing it.
