# Benscoter restrained warping (`LinearBenscoterBeamElement3D`)

`LinearBenscoterBeamElement3D` is the Timoshenko element with one more field: the
warping amplitude \(\psi(x)\), carried on the seventh nodal DOF of a warping mesh
(`[warping] = 1`). It needs a thin-walled section from `element_section.txt`.

## Section level

`ThinWalledSectionAnalysis.benscoter_stiffness` (8 × 8, right-handed) and
`SectionAnalysis.benscoter_stiffness` (element signs) act on

\[
\boldsymbol\varepsilon_8 = [\varepsilon_x,\ \kappa_y,\ \kappa_z,\ \gamma_{xy},\ \gamma_{xz},\ \phi_x,\ d,\ \psi']^{\mathsf T},
\qquad d = \psi - \phi_x .
\]

The section warps with the Saint-Venant warping of the six beam strains plus
\(\omega(s)\,d\), where \(\omega\) is the torsional warping about the shear centre,
orthogonalised against rigid axial motion. The wall strains gain

\[
\varepsilon_{xx} \mathrel{+}= \omega\,\psi', \qquad \gamma_{xs} \mathrel{+}= \omega_{,s}\, d .
\]

So with \(\psi = \phi_x\) the section is exactly the Saint-Venant one. The stiffness is
\(\int \mathbf G^{\mathsf T}\mathbf K_w\mathbf G\,\mathrm ds\) over the walls, with the
leading 6 × 6 block set to the Saint-Venant `stiffness`. Then:

- `K[6, 6]` is the warping shear stiffness \(S\). For an isotropic I-beam it is
  \(G\,2 b t_f (h/2)^2\).
- `K[7, 7]` is the primary warping stiffness \(E\Gamma\). For an I-beam it is
  \(E t_f b^3 h^2/24\). Secondary (through-thickness) warping is not included here.
- The off-diagonal terms couple warping to extension, bending and shear for anisotropic walls.

`warping_inertia` is \(\int \rho_A\,\omega^2\,\mathrm ds\).

## Element

- Element DOFs: the 12 standard ones (node-major), then \(\psi_1, \psi_2\).
- \(\psi\) and \(\theta_x\) are linear, so \(d\) is linear like the shear strains.
- Stiffness: a full rule at `max_full_order`, then the \(\{\gamma_{xy}, \gamma_{xz}, d\}\)
  block is replaced by its value at the shear-block order (default one point). This avoids
  warping-shear locking the same way it avoids shear locking.
- Mass: the rigid-section mass of the `SectionAnalysis` plus a consistent
  \(m_\omega\) mass on \(\psi\).
- Boundary conditions: `W` (or `CHI`) = 0 restrains warping. A free end needs nothing
  (the bimoment \(B = \partial W/\partial\psi'\) vanishes naturally).

The exact restrained-torsion solution is in
[proofs/benscoter/restrained_torsion.md](../proofs/benscoter/restrained_torsion.md); the
tests compare against it in `tests/test_benscoter_beam_element.py`.

## Limits

- Linear static and modal only; `linear_geometric_stiffness_matrix` raises (second order
  is a later step).
- One warping amplitude per section (torsional warping). Distortion and other
  cross-section modes come with the generalised beam theory step.
