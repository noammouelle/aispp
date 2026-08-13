# INPUT_REFERENCE.md — AIS++ input file format

AIS++ reads plain-text `.aisi` files. The recommended way to generate them is via
the **aispy** Python package, which handles the pulse schedule automatically (see
`aispy/utils.py`, class `AISFlow`). This document describes every parameter for
reference or for writing custom input files.

Lines beginning with `#` are comments. Parameter order within each section does
not matter. Blank lines are ignored.

---

## How to generate input files

The intended workflow is:

```python
from aispy.utils import AISFlow

param_dict = {
    'cloud_params':      {...},
    'potential_params':  {...},
    'sequence_params':   {...},
    'pulse_params':      {...},
    'simulation_params': {...},
    'io_params':         {...},
}

flow = AISFlow(param_dict, flowdir='path/to/output', workdir='path/to/workdir')
```

`AISFlow` writes the complete `.aisi` file including the full pulse schedule
(start/end times, wavevectors, laser frequencies, detunings). The per-pulse arrays
in the **Pulse parameters** section below are therefore rarely written by hand.

---

## 1. Cloud parameters

These define the initial atomic ensemble.

| Key           | Type    | Unit  | Description |
|---------------|---------|-------|-------------|
| `natoms`      | int     | —     | Number of Monte Carlo atoms to simulate. |
| `initialstate`| int     | —     | Initial internal state. `0` = ground state \|g⟩, `1` = excited state \|e⟩. |
| `sigma`       | float   | m     | RMS radius of the initial position distribution (Gaussian). Set to `0` for a point cloud. |
| `transtemp`   | float   | K     | Transverse temperature of the atom cloud (sets transverse velocity spread). |
| `longtemp`    | float   | K     | Longitudinal temperature (sets velocity spread along the launch direction). |
| `x0`          | 3 floats| m     | Initial centre-of-mass position `x y z`. |
| `v0`          | 3 floats| m/s   | Initial centre-of-mass velocity `vx vy vz`. For a vertically launched cloud, `vz` is the launch velocity. |

**Example:**
```
natoms 100000
initialstate 0
sigma 0.001
transtemp 1e-9
longtemp 0.0
x0 0.0 0.0 0.0
v0 0.0 0.0 19.62
```

### Phase-space grid mode

| Key        | Type            | Description |
|------------|-----------------|-------------|
| `initmode` | string          | `gaussian` (default) samples the thermal cloud described above. `psgrid` instead places atoms on a deterministic 6-D grid. |
| `xgrid`, `ygrid`, `zgrid` | `min max n` | Position grid limits (m) and point count per axis. |
| `vxgrid`, `vygrid`, `vzgrid` | `min max n` | Velocity grid limits (m/s) and point count per axis. These are **frame** velocities; see the rotating-frame note above. |

All six axis keys are required when `initmode psgrid` is set, and each `n` must
be at least 1. The Gaussian cloud keys must still be present — the parser
requires them — but they are ignored in this mode. The total atom count is the
product of the six point counts, not `natoms`.

**Example** — a transverse position/velocity map at fixed $z$:
```
initmode psgrid
xgrid  -0.002 0.002 41
ygrid  -0.002 0.002 41
zgrid   0.0   0.0    1
vxgrid -0.01  0.01  21
vygrid -0.01  0.01  21
vzgrid  0.0   0.0    1
```

---

## 2. Potential parameters

| Key     | Type   | Description |
|---------|--------|-------------|
| `utype` | string | External potential type. Inertial: `zero_pot` (no force), `linear_pot` (uniform gravity, $g = 9.81$ m/s²), `quadratic_pot` (gravity gradient). Rotating-frame counterparts: `rotating_pot`, `rotating_linear_pot`, `rotating_quadratic_pot`. |
| `rotation` | 3 × double | *Optional, default `0 0 0`.* Angular velocity $\vec\Omega$ of the frame in rad/s, in simulation-frame components. Only accepted with a `rotating_*` `utype`; supplying it with an inertial `utype` is an error. |

**Example:**
```
utype linear_pot
```

### Rotating frame

The `rotating_*` potentials add the Coriolis and centrifugal forces of a frame
rotating at constant $\vec\Omega$. The equation of motion integrated is

$$m\ddot{\vec r} = -\nabla V - 2m\,\vec\Omega\times\dot{\vec r} - m\,\vec\Omega\times(\vec\Omega\times\vec r).$$

Internally this is done by promoting the Hamiltonian to
$H = p^2/2m - \vec\Omega\cdot(\vec r\times\vec p) + V(\vec r)$, so the centrifugal
term is generated automatically from $p^2/2m$ and is *not* added separately.

Because the canonical momentum is $\vec p = m(\dot{\vec r} + \vec\Omega\times\vec r)$,
it is not $m\dot{\vec r}$ when the frame rotates. This is invisible at the input
and output boundaries: `v0`, the `vgrid` limits and the thermal spreads are all
**frame** velocities $\dot{\vec r}$, and every velocity written to file is a frame
velocity too. The conversion is the identity when $\vec\Omega = 0$, so inertial
runs are bit-for-bit unchanged.

Every output file gains a `rotation` dataset holding the $\vec\Omega$ actually used.

**Example** — Earth rotation at latitude 45° N, with $z$ vertical and $x$ pointing east:
```
utype rotating_linear_pot
rotation 0.0 5.156e-5 5.156e-5
```

---

## 3. Simulation parameters

Control numerical behaviour and performance trade-offs.

| Key                  | Type  | Default | Description |
|----------------------|-------|---------|-------------|
| `amplitudethreshold` | float | —       | Wavepackets with amplitude below this value are discarded. Typical value `0.001`. Reduce to keep more paths at the cost of memory and runtime. |
| `coherencelength`    | float | —       | Maximum spatial separation (m) for two wavepackets to be considered interfering at detection. |
| `usemcbranching`     | 0/1   | —       | Enable stochastic path pruning (Monte Carlo branching). When `1`, a single-atom path-finder run is performed first to identify the interfering paths; those are propagated deterministically, and all other branches are stochastically pruned. |
| `usepathselection`   | 0/1   | —       | Only simulate the paths listed in `pathstosimulate`. Faster than MC branching for known interferometer geometries. |
| `pathstosimulate`    | strings| —      | Space-separated list of binary path strings (e.g. `0101...0`). Each character encodes the internal state at successive beam-splitter interactions (0 = ground, 1 = excited). Generated automatically by `AISFlow`. |
| `ignoredetuning`     | 0/1   | —       | Set $\Delta = 0$ in the beam-splitter ODE. Useful for isolating wavefront effects from Doppler shifts. |
| `usestaticapprox`    | 0/1   | —       | Use the static approximation (atom does not move during a pulse). Obsolete; kept for backward compatibility. |
| `ultrafast`          | 0/1   | —       | Skip the classical action phase (kinematic phase). Only wavefront/laser phases are accumulated. Significantly faster; valid when the interferometric signal is dominated by wavefront effects. |
| `usedetvolselection` | 0/1   | —       | Only propagate wavepackets whose final position falls within the detection volume defined by `xdet`, `ydet`, `zdet`. |
| `xdet`               | 2 floats| — | Detection volume x-range: `xmin xmax` (m). |
| `ydet`               | 2 floats| — | Detection volume y-range: `ymin ymax` (m). |
| `zdet`               | 2 floats| — | Detection volume z-range: `zmin zmax` (m). |
| `seed`               | float | —       | RNG seed for atom ensemble initialisation. Set to `-1` for a random seed. |
| `gslqagabserr`       | float | 1e-12   | Absolute error tolerance for GSL QAGS numerical integration (classical action). |
| `gslqagrelerr`       | float | 1e-12   | Relative error tolerance for GSL QAGS integration. |
| `gslkinodeabserr`    | float | 1e-9    | Absolute error tolerance for the kinematic ODE solver (Hamilton's equations). |
| `gslkinoderelerr`    | float | 0       | Relative error tolerance for the kinematic ODE solver. |
| `gslpulseodeabserr`  | float | 1e-9    | Absolute error tolerance for the beam-splitter ODE solver (two-level Schrödinger equation). |
| `gslpulseoderelerr`  | float | 0       | Relative error tolerance for the beam-splitter ODE solver. |

**Notes on tolerances.** The default combination of `absErr = 1e-9` and `relErr = 0` for
the ODE solvers means only the absolute tolerance is active. For simulations where phases
accumulate to large values, tightening `gslpulseodeabserr` to `1e-12` improves accuracy at
the cost of runtime. The QAGS tolerances affect the action-phase integral and are only
relevant when `ultrafast = 0`.

---

## 4. Sequence parameters

| Key             | Type    | Unit | Description |
|-----------------|---------|------|-------------|
| `detectiontime` | float   | s    | Time at which atoms are detected. Must be ≥ the end time of the last pulse. Atoms are propagated ballistically from the last pulse to this time. |

---

## 5. IO parameters

| Key              | Type | Description |
|------------------|------|-------------|
| `printprobs`     | 0/1  | Write detection probabilities to the HDF5 output file. |
| `printwavepackets`| 0/1 | Write full wavepacket data (position, velocity, amplitude, phase) to the HDF5 output file. Produces large files for big ensembles. |
| `printtrajectory`| 0/1  | Write a snapshot of every wavepacket at every pulse boundary to a companion `*_TRAJ.h5` file. Consumed by `aispy.trajectory` for arm reconstruction, enclosed-area and Sagnac-phase analysis. Optional; default `0`. |

---

## 6. Pulse parameters

Each laser pulse is described by a set of per-pulse arrays. All arrays must have
the same length $N_{\text{pulses}}$. **These are normally generated by `AISFlow`;
hand-editing is error-prone for large LMT sequences.**

Each row below gives the keyword, its meaning, and the units. Values are written
as space-separated lists on a single line.

| Key              | Unit     | Description |
|------------------|----------|-------------|
| `t0`             | s        | Start time of each pulse. |
| `t1`             | s        | End time of each pulse. Pulse duration = `t1 - t0 = π/Ω₀` for a π-pulse. |
| `kx`, `ky`, `kz` | rad/m    | Wavevector components of the laser field for each pulse. For the clock transition at 698 nm, `|k| ≈ 9 × 10⁶` rad/m. The sign of `kz` encodes the propagation direction (upward vs downward beam). |
| `omega`          | rad/s    | Angular frequency of the laser field for each pulse. For each pulse, `omega` is computed to put the transition on resonance accounting for the Doppler shift of the atom at that time. Computed automatically by `AISFlow`. |
| `rabifreq`       | rad/s    | On-axis Rabi frequency $\Omega_0$ for each pulse. |
| `phi0`           | rad      | Constant phase offset added to the laser phase $\varphi(\boldsymbol{r})$ for each pulse. Usually `0`. |
| `wtype`          | string   | Beam profile type. Options: `gaussian` (analytic Gaussian intensity and wavefront), `confocal` (concave mirror with coincident in/out foci), `flat_square` (uniform intensity, flat wavefront), `interpolated` (sampled beam read from an HDF5 grid — see below). |
| `waist`          | m        | Beam waist $w_0$ at focus. Used by `gaussian` and `confocal`; ignored by `interpolated` and `flat_square`. |
| `focallength`    | m        | Focus position $f$ along the beam axis (distance from the mirror at $z = 0$). Set to `0` for focus at the mirror (HEHN configuration). |
| `zlaser`         | m        | **Inert.** Parsed and stored but never read; beam geometry comes from `focallength` and the sign of `kz`. Retained for input compatibility — see `KNOWN_ISSUES.md` #5. |
| `beamradius`     | m        | Normalisation radius $r_\text{beam}$ for the Zernike expansion below. Required (and must be positive) when any `zernikecoeff_N` is non-zero. It does **not** clip the beam: atoms outside this radius still receive intensity from the envelope. |
| `baseline`       | m        | Vertical separation between the two interferometers in a gradiometer configuration. |
| `beaminterpolationparamsfilenames` | path | Per-pulse path to the HDF5 beam file, required when `wtype` is `interpolated` and ignored otherwise. |
| `tiptiltx`, `tiptilty` | deg | Tip/tilt of the retroreflector about the $x$ and $y$ axes. Applied to `interpolated` beams as a shear of the sample point, $x \to x - z\tan\theta_x$. Optional; default `0`. |
| `kxchirp`, `kychirp`, `kzchirp` | rad/m/s | Linear time-derivative of the wavevector (wavevector chirp). Set to `0` for no chirp. |
| `frequencychirp` | rad/s²   | Linear time-derivative of the laser frequency. Set to `0` for no chirp. |
| `zernikecoeff_N` | rad/m²   | Zernike coefficient $N$ of the wavefront aberration (see below). Multiple coefficients can be included with different indices. |

### Zernike wavefront coefficients

The wavefront phase deviation from a plane wave is expanded in Zernike polynomials
over the beam aperture:

$$\varphi(\boldsymbol{r}) = \sum_N c_N Z_N(\boldsymbol{r}/r_{\text{beam}})$$

Each coefficient is specified as a per-pulse array with key `zernikecoeff_N` where
`N` is the Zernike Noll index. Set all coefficients to `0` to use the analytic Gaussian wavefront profile alone
(wavefront curvature from the beam geometry is still included via `focallength` and
`waist`). `beamradius` sets $r_\text{beam}$ and must be positive whenever any
coefficient is non-zero.

> **Restored in v0.0.2.** Between commit `adbebdf` and this release the Zernike
> block was absent from `GetPhi`, so `zernikecoeff_N` was parsed and silently
> ignored. Results produced in that window contain no aberration. Note also that
> the aberration contributes to the phase but **not** to the wavefront gradient —
> see `KNOWN_ISSUES.md` #1.

---

## 6b. Interpolated (sampled) beams

Setting `wtype interpolated` replaces the analytic profile with a beam sampled on
a regular 3-D grid, as written by the `AISPPExporter` in the companion
**aisoptics** package:

```python
from aisoptics import AISPPExporter
AISPPExporter().export_total_field("beam.h5", field_grid)
```

The file must contain 1-D coordinate datasets `x`, `y`, `z` and flattened
`phase` and `amplitude` datasets in C order, indexed `i*Ny*Nz + j*Nz + k`. It
should carry the attribute `export_format = "aispp_current_interpolation_hdf5"`;
ais++ rejects any other value and warns if the attribute is absent.

Conventions, matching the exporter:

- `phase` is the slowly varying envelope phase in radians and **excludes** the
  $\pm k z$ carrier, which ais++ handles separately.
- `amplitude` is field amplitude relative to the configured `rabifreq`, so the
  effective Rabi frequency is `rabifreq × amplitude`.

Points outside the sampled volume take the nearest edge value (constant
extrapolation) rather than extrapolating linearly. `waist`, and the Zernike
coefficients are not applied to interpolated beams — the grid is taken as the
complete description of the beam.

---

## 7. Complete minimal example

A simple Mach–Zehnder (π/2 – T – π – T – π/2) sequence for a single atom,
no LMT, Gaussian beam, no wavefront aberrations:

```
# Cloud
natoms 1000
initialstate 0
sigma 0.0
transtemp 0.0
longtemp 0.0
x0 0.0 0.0 0.0
v0 0.0 0.0 0.0

# Potential
utype zero_pot

# Simulation
amplitudethreshold 0.001
coherencelength 1.0
usemcbranching 0
usepathselection 1
pathstosimulate 000 001 010 011
ignoredetuning 0
usestaticapprox 0
ultrafast 0
usedetvolselection 0
xdet -1.0 1.0
ydet -1.0 1.0
zdet -1.0 1.0
seed -1
gslqagabserr 1e-12
gslqagrelerr 1e-12
gslkinodeabserr 1e-9
gslkinoderelerr 0
gslpulseodeabserr 1e-9
gslpulseoderelerr 0

# Sequence
detectiontime 0.0016

# IO
printprobs 1
printwavepackets 0

# Pulses  (3 pulses: pi/2, pi, pi/2)
# Generated by AISFlow for T = 0.5 ms, Omega_0 = 2pi*1000 rad/s
t0 0.0 0.0005 0.001
t1 0.00025 0.00075 0.00125
kx 0.0 0.0 0.0
ky 0.0 0.0 0.0
kz 9000000.0 9000000.0 9000000.0
omega 2696919267000000.0 2696919267000000.0 2696919267000000.0
rabifreq 6283.185307 6283.185307 6283.185307
phi0 0.0 0.0 0.0
wtype gaussian gaussian gaussian
waist 0.03 0.03 0.03
focallength 0.0 0.0 0.0
zlaser 0.0 0.0 0.0
beamradius 0.05 0.05 0.05
baseline 0.0 0.0 0.0
kxchirp 0.0 0.0 0.0
kychirp 0.0 0.0 0.0
kzchirp 0.0 0.0 0.0
frequencychirp 0.0 0.0 0.0
zernikecoeff_1 0.0 0.0 0.0
zernikecoeff_4 0.0 0.0 0.0
```

> **In practice**, use `AISFlow` to generate files like this. For an LMT sequence of
> order $n = 101$, the pulse arrays each contain ~400 entries with individually
> Doppler-tuned frequencies; constructing these by hand is infeasible.