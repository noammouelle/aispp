# THEORY.md — Semi-classical beam-splitter model in AIS++

This document derives the equations solved by AIS++ during a laser pulse. The derivation follows
the unitary transformation approach of Antoine [1], adapted to the specific form implemented in
the code. The reader is assumed familiar with basic quantum mechanics and atom interferometry.

---

## 1. Setup and notation

We consider a two-level atom with ground state $|g\rangle$ and excited state $|e\rangle$,
separated by frequency $\omega_0$, in an external potential $U(\hat{\boldsymbol{r}}, \hat{\boldsymbol{p}})$.
The full Hamiltonian is

$$\hat{H} = \hat{H}_{\text{ext}}(\hat{\boldsymbol{r}}, \hat{\boldsymbol{p}}) \otimes \hat{\mathbf{1}} + \hbar\omega_0 |e\rangle\langle e| + \hbar\Omega(\hat{\boldsymbol{r}}, t)\cos\!\bigl(\boldsymbol{k}\cdot\hat{\boldsymbol{r}} - \omega t + \varphi(\hat{\boldsymbol{r}})\bigr)\bigl(|e\rangle\langle g| + |g\rangle\langle e|\bigr),$$

where $\hat{H}_{\text{ext}} = \hat{\boldsymbol{p}}^2/2m + U(\hat{\boldsymbol{r}}, \hat{\boldsymbol{p}})$
is the centre-of-mass Hamiltonian, $\boldsymbol{k}$ is the laser wavevector, $\omega$ the laser
frequency, $\Omega(\hat{\boldsymbol{r}}, t)$ the (position- and time-dependent) Rabi frequency,
and $\varphi(\hat{\boldsymbol{r}})$ captures any deviation of the laser phase from an ideal plane
wave (e.g. wavefront curvature).

The state of the atom is a spinor wavepacket,

$$|\psi\rangle = \begin{pmatrix} \psi_e(\boldsymbol{r}, t) \\ \psi_g(\boldsymbol{r}, t) \end{pmatrix}.$$

In the semi-classical approximation we track a finite ensemble of Gaussian wavepackets, each
characterised by a centre $(\boldsymbol{r}_0, \boldsymbol{p}_0)$ in phase space, a complex
amplitude $A$, and a phase $\phi$. Operators evaluated on such a wavepacket reduce to their
classical values at the wavepacket centre.

---

## 2. Step 1 — Interaction picture: removing $\hat{H}_{\text{ext}}$ and $\hbar\omega_0$

Define the first unitary transformation

$$\hat{U}_1(t, t_0) = e^{\frac{i}{\hbar}\hat{H}_{\text{ext}}(t - t_0)} \otimes \Bigl(e^{i\omega_0(t-t_0)}|e\rangle\langle e| + |g\rangle\langle g|\Bigr),$$

and the transformed state $|\phi_1\rangle = \hat{U}_1|\psi\rangle$. A standard calculation gives
the transformed Hamiltonian

$$\hat{H}_1 = \hbar\hat{U}_{\text{ext}}^\dagger\,\Omega(\hat{\boldsymbol{r}}, t)\cos\!\bigl(\boldsymbol{k}\cdot\hat{\boldsymbol{r}} - \omega t + \varphi(\hat{\boldsymbol{r}})\bigr)\hat{U}_{\text{ext}} \otimes \bigl(e^{i\omega_0(t-t_0)}|e\rangle\langle g| + \text{h.c.}\bigr),$$

where $\hat{U}_{\text{ext}}(t, t_0) = \exp\!\bigl(-\frac{i}{\hbar}\hat{H}_{\text{ext}}(t-t_0)\bigr)$.

Applying the rotating-wave approximation (dropping terms oscillating at $2\omega$) gives

$$\hat{H}_1 = \frac{\hbar}{2}\hat{U}_{\text{ext}}^\dagger\,\Omega(\hat{\boldsymbol{r}}, t) \Bigl(e^{i\Phi(\hat{\boldsymbol{R}}(t,t_0),\,t)}|e\rangle\langle g| + \text{h.c.}\Bigr)\hat{U}_{\text{ext}},$$

where the full laser phase is

$$\Phi(\hat{\boldsymbol{r}}, t) = \boldsymbol{k}\cdot\hat{\boldsymbol{r}} - (\omega - \omega_0)t - \omega_0 t_0 + \varphi(\hat{\boldsymbol{r}}),$$

and $\hat{\boldsymbol{R}}(t, t_0) \equiv \hat{U}_{\text{ext}}^\dagger\hat{\boldsymbol{r}}\hat{U}_{\text{ext}}$
is the position operator in the Heisenberg picture under $\hat{H}_{\text{ext}}$ — i.e., the
quantum analogue of the classical trajectory. In the semi-classical limit,

$$\hat{\boldsymbol{R}}(t, t_0)\,\zeta(\boldsymbol{r}_0, \boldsymbol{p}_0, t_0) = \boldsymbol{r}_{\text{cl}}(\boldsymbol{r}_0, \boldsymbol{p}_0, t_0, t)\,\zeta(\boldsymbol{r}_0, \boldsymbol{p}_0, t_0),$$

so $\hat{\boldsymbol{R}}$ simply evaluates the classical trajectory at time $t$.

---

## 3. Step 2 — Linearised rotating frame: removing the laser phase

The laser phase $\Phi$ is a non-linear function of $\hat{\boldsymbol{r}}$, which prevents a clean
separation of internal and external dynamics. We remove it approximately by a second unitary
transformation

$$\hat{U}_2(t, t_0) = e^{-i\Theta(\hat{\boldsymbol{R}}'(t,t_0),\,t,\,t_0)}|e\rangle\langle e| + |g\rangle\langle g|,$$

where the phase is

$$\Theta(\hat{\boldsymbol{R}}'(t,t_0), t, t_0) = \boldsymbol{k}'\cdot\hat{\boldsymbol{R}}'(t,t_0) - (\omega - \omega_0)t - \omega_0 t_0 + \varphi',$$

with

$$\boldsymbol{k}' = \boldsymbol{k} + \nabla\varphi(\boldsymbol{r}^*), \qquad \varphi' = \varphi(\boldsymbol{r}^*).$$

Here $\boldsymbol{r}^*$ is the phase-space point about which we linearise (the wavepacket centre at
$t_0$), and $\hat{\boldsymbol{R}}'(t, t_0)$ is the linearised propagator — the solution to
Hamilton's equations expanded to first order around $(\boldsymbol{r}^*, \boldsymbol{p}^*)$:

$$\hat{\boldsymbol{R}}'(t, t_0) = A(t,t_0)\cdot\hat{\boldsymbol{r}} + B(t,t_0)\cdot\frac{\hat{\boldsymbol{p}}}{m} + \boldsymbol{\xi}(t,t_0).$$

The matrices $A$, $B$ and the vector $\boldsymbol{\xi}$ are obtained by solving the linearised
Hamilton equations; for AIS++'s current potentials (uniform gravity) $A = I$ and $B = 0$,
simplifying the calculation considerably.

**Effect of $\hat{U}_2$ on a wavepacket.** At $t = t_0$, $\hat{U}_2$ acts on the excited
component as

$$\hat{U}_2(t_0, t_0)|\zeta_e(\boldsymbol{r}_0^e, \boldsymbol{p}_0^e)\rangle = e^{-i\Phi(\boldsymbol{r}_0^e, t_0)}|\zeta_e(\boldsymbol{r}_0^e,\; \boldsymbol{p}_0^e - \hbar\boldsymbol{k}')\rangle,$$

imprinting the local laser phase and shifting the momentum of the excited wavepacket by
$-\hbar\boldsymbol{k}'$. The ground-state component is unchanged. This is precisely what
`AISPulsePropagator::ApplyU2` implements.

**Transformed Hamiltonian.** Under $\hat{U}_2$, the Hamiltonian $\hat{H}_1$ becomes (to leading
order in $\delta \equiv \Phi - \Theta$, which is small by construction)

$$\hat{H}_2^{(0)} = \hbar\Delta(\dot{\hat{\boldsymbol{R}}}', t, t_0)|e\rangle\langle e| + \frac{\hbar\Omega(\hat{\boldsymbol{R}}, t)}{2}\bigl(|e\rangle\langle g| + |g\rangle\langle e|\bigr),$$

plus an additive contribution from $i\hbar\dot{\hat{U}}_2\hat{U}_2^\dagger$,

$$i\hbar\dot{\hat{U}}_2\hat{U}_2^\dagger = \hbar\left(\boldsymbol{k}'\cdot\frac{d\hat{\boldsymbol{R}}'}{dt} - (\omega - \omega_0) + \frac{\hbar|\boldsymbol{k}'|^2}{2m}\right)|e\rangle\langle e|,$$

where the last term is a recoil correction. Combining, the effective detuning is

$$\Delta(\dot{\hat{\boldsymbol{R}}}', t, t_0) = \boldsymbol{k}'\cdot\frac{d\hat{\boldsymbol{R}}'}{dt} - (\omega - \omega_0) + \frac{\hbar|\boldsymbol{k}'|^2}{2m}.$$

In the semi-classical limit, $d\hat{\boldsymbol{R}}'/dt$ is replaced by the classical velocity
along the linearised trajectory, so

$$\Delta = \boldsymbol{k}'\cdot\boldsymbol{v}'_{\text{cl}}(t) - (\omega - \omega_0) + \frac{\hbar|\boldsymbol{k}'|^2}{2m}.$$

This is the Doppler shift (first term), the bare detuning (second term), and the two-photon
recoil (third term). Because $\boldsymbol{k}' = \boldsymbol{k} + \nabla\varphi(\boldsymbol{r}^*)$,
the wavefront gradient $\nabla\varphi$ contributes directly to $\Delta$ — this is the primary
mechanism by which beam imperfections affect the interferometric phase.

---

## 4. Step 3 — Solving the two-level ODE

In the frame defined by $\hat{U}_1\hat{U}_2$, the Schrödinger equation reduces to a
position-dependent two-level system. Writing the spinor as
$(c_e(t),\, c_g(t))^T$, the ODE is

$$\frac{d}{dt}\begin{pmatrix}c_e \\ c_g\end{pmatrix} = -i\begin{pmatrix}\Delta(\boldsymbol{r}'_{\text{cl}}, \boldsymbol{v}'_{\text{cl}}, t) & \Omega(\boldsymbol{r}_{\text{cl}}, t)/2 \\ \Omega(\boldsymbol{r}_{\text{cl}}, t)/2 & 0\end{pmatrix}\begin{pmatrix}c_e \\ c_g\end{pmatrix}.$$

Splitting into real and imaginary parts (as required by GSL's real-valued ODE interface),

$$\frac{d}{dt}\begin{pmatrix}\text{Re}\,c_e \\ \text{Re}\,c_g \\ \text{Im}\,c_e \\ \text{Im}\,c_g\end{pmatrix} = \begin{pmatrix}0 & 0 & \Delta & \Omega/2 \\ 0 & 0 & \Omega/2 & 0 \\ -\Delta & -\Omega/2 & 0 & 0 \\ -\Omega/2 & 0 & 0 & 0\end{pmatrix}\begin{pmatrix}\text{Re}\,c_e \\ \text{Re}\,c_g \\ \text{Im}\,c_e \\ \text{Im}\,c_g\end{pmatrix}.$$

This 4-dimensional real system is what `AISPulsePropagator::funcU3` evaluates, and it is
integrated by `AISPulsePropagator::ApplyU3` using GSL's `gsl_odeiv2_step_rk8pd` (8th-order
Runge–Kutta–Prince–Dormand) with user-specified absolute and relative tolerances.

The classical quantities $\boldsymbol{r}_{\text{cl}}$, $\boldsymbol{r}'_{\text{cl}}$,
$\boldsymbol{v}'_{\text{cl}}$ entering $\Delta$ and $\Omega$ are evaluated at each ODE step using
`AISKinematicPropagator::CalculateNewPhaseSpaceCoordsLinearized`.

---

## 5. Step 4 — Inverse transformations: $\hat{U}_2^\dagger$ and $\hat{U}_1$

After the ODE integration, the wavepacket is in the $(\hat{U}_1\hat{U}_2)$ frame. The inverse
transformations recover the wavepacket in the original frame.

**$\hat{U}_2^\dagger$ (`ApplyU2Dagger`).** The effect on the excited wavepacket centre is

$$\boldsymbol{r}_0^e \;\to\; \boldsymbol{r}_0^e + \delta\boldsymbol{r}, \qquad \boldsymbol{p}_0^e \;\to\; \boldsymbol{p}_0^e + \hbar A^T\boldsymbol{k}',$$

and the phase acquires the contribution

$$\delta\phi = \varphi(\boldsymbol{r}^*) - \boldsymbol{r}^*\cdot\nabla\varphi(\boldsymbol{r}^*) + \boldsymbol{k}'\cdot\bigl(A\boldsymbol{r}_0^e + B\,\delta\boldsymbol{k}' + \boldsymbol{\xi}\bigr) - (\omega_{\text{laser}} - \omega_{\text{Sr}})\,t_1 - \omega_{\text{Sr}}\,t_0,$$

where $\delta\boldsymbol{k}' = -\frac{\hbar}{2m}B A^T\boldsymbol{k}'$ is a small recoil-induced
position shift arising from the BCH commutator $\frac{1}{2}[\boldsymbol{k}'\cdot A\hat{\boldsymbol{r}},\, \boldsymbol{k}'\cdot B\hat{\boldsymbol{p}}/m] = -\frac{i\hbar}{2m}\boldsymbol{k}'^T(AB^T)\boldsymbol{k}'$.
The optical-frequency part $(\omega_{\text{laser}} - \omega_{\text{Sr}})t_1 + \omega_{\text{Sr}}t_0$
is accumulated in `__float128`.

**$\hat{U}_1$ (`ApplyU1`).** The kinematic propagator is run backwards from $t_1$ to $t_0$
(note the reversed time argument in the code), advancing the wavepacket centre and accumulating
the kinematic (action) phase for the return journey. The `fAddEnergyPhase` flag controls whether
the energy-phase term $-H_{\text{ext}} \cdot t / \hbar$ is included; it is temporarily set to
`true` during this step.

---

## 6. Wavepacket splitting and the interferometer

After applying $\hat{U}_1\hat{U}_2\hat{U}_3\hat{U}_2^\dagger\hat{U}_1$ to each incoming
wavepacket, the result is two daughter wavepackets: one that absorbed a photon (path "1") and one
that did not (path "0"). Each carries the appropriate amplitude and phase from the ODE solution.

Over a full $n$-pulse LMT Mach–Zehnder sequence, the wavepacket tree has up to $2^{2n+2}$
branches. AIS++ handles this with a stochastic pruning scheme (`ApplyMCBranching`): at each
beam-splitter, non-interfering branches are collapsed onto one randomly chosen daughter with
probability proportional to its squared amplitude, and the survivor is reweighted. Branches on
paths that must interfere at the final beam-splitter are always propagated deterministically.

At detection, pairs of wavepackets that share an internal state and overlap in phase space are
identified. Their probability of detection in each output port is

$$P = |A_1|^2 + |A_2|^2 + 2|A_1||A_2|\cos(\Delta\varphi),$$

where $\Delta\varphi = \phi_1 - \phi_2$ is the accumulated phase difference.

---

## 7. Phase precision

The total phase of each wavepacket is maintained as a sum of two parts:

- **`phaseDouble` (64-bit `double`):** kinematic contributions (classical action, wavefront
  terms, $\hat{U}_2$ and $\hat{U}_2^\dagger$ real-valued phase corrections). These are of order
  $kv_{\text{rec}} T \sim \mathcal{O}(10^3)$ rad.

- **`phaseQuad` (`__float128`, ~33 significant digits):** optical-frequency contributions,
  specifically terms of the form $\omega_{\text{laser}} t$ and $\omega_0 t_0$. These are of
  order $10^{15}$ rad. The interferometric phase difference is the small difference between two
  such large numbers; `double` would lose all significant figures.

The final phase difference $\Delta\varphi = \Delta\phi_{\text{double}} + \Delta\phi_{\text{quad}}$
is computed by summing both parts at detection time.

---

## References

[1] C. Antoine, *Contribution à la théorie des interféromètres atomiques*, PhD thesis,
Université Pierre et Marie Curie (2004). Sections 2.4, 3.1–3.2, Appendix A.

[2] N. Mouelle, J. Mitchell, V. Gibson, U. Schneider, *Wavefront Curvature and Transverse Atomic
Motion in Time-Resolved Atom Interferometry: Impact and Mitigation*, arXiv:2510.26739 (2025).

[3] N. Mouelle, *Numerical treatment of the semi-classical beam-splitter* (internal notes, 2024).