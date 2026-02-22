# Euler–Bernoulli: Shear force from constitutive law

This note proves that **Euler–Bernoulli (EB) elements yield zero shear resultants** \(V_y\), \(V_z\) from the constitutive relation \(\sigma = D\,\varepsilon\). Shear force in EB is therefore only defined by equilibrium (\(V = \mathrm{d}M/\mathrm{d}x\)), not by the stress output at Gauss points.

For the contrasting Timoshenko case (non-zero shear from constitutive law), see [../timoshenko/shear_force_from_constitutive.md](../timoshenko/shear_force_from_constitutive.md).

---

## 1. Common setup (beam strain and stress resultants)

We use a 6-component beam strain and stress resultant ordering:

- **Strain** (conjugate to stress resultants):
  \[
  \varepsilon = \bigl[\, \varepsilon_x \;\; \kappa_y \;\; \kappa_z \;\; \gamma_{xy} \;\; \gamma_{xz} \;\; \phi_x \,\bigr]^\top.
  \]

- **Stress resultants** (from constitutive law \(\sigma = D\,\varepsilon\)):
  \[
  \sigma = \bigl[\, N \;\; M_y \;\; M_z \;\; V_y \;\; V_z \;\; T \,\bigr]^\top.
  \]

The constitutive matrix \(D\) is **diagonal**, so the fourth and fifth stress resultants are

\[
V_y = D_{44}\,\gamma_{xy}, \qquad V_z = D_{55}\,\gamma_{xz}.
\]

Hence: if \(\gamma_{xy} = \gamma_{xz} = 0\) **or** \(D_{44} = D_{55} = 0\), then \(V_y = V_z = 0\).

---

## 2. Euler–Bernoulli element

**Kinematic assumption:** plane sections remain plane and normal to the neutral axis. Transverse shear deformation is neglected, so the shear strains are **defined to be zero**:

\[
\gamma_{xy} = 0, \qquad \gamma_{xz} = 0 \quad \text{(everywhere in the element)}.
\]

**Strain–displacement matrix \(B\):**  
The EB \(B\) matrix has no rows that depend on \((\partial u_y/\partial x - \theta_z)\) or \((\partial u_z/\partial x - \theta_y)\). The fourth and fifth rows of \(B\) are zero, so for any element displacement vector \(\mathbf{u}_e\),

\[
\varepsilon = B\,\mathbf{u}_e \quad \Rightarrow \quad \gamma_{xy} = 0,\;\; \gamma_{xz} = 0.
\]

**Constitutive matrix \(D\):**  
The EB material matrix sets the shear stiffness entries to zero (no constitutive law for shear):

\[
D_{44} = 0, \qquad D_{55} = 0.
\]

So

\[
V_y = D_{44}\,\gamma_{xy} = 0, \qquad V_z = D_{55}\,\gamma_{xz} = 0
\]

at every Gauss point, for any load.

**Conclusion:**  
The EB element **does not produce shear force from \(\sigma = D\,\varepsilon\)**. In EB theory, shear force is defined by **equilibrium**: \(V = \mathrm{d}M/\mathrm{d}x\). If needed, \(V\) must be obtained by differentiating the computed bending moment \(M\) along the beam, not from the Gauss-point stress resultants.

**References in code:**  
`pre_processing/element_library/euler_bernoulli/utilities/B_matrix.py` (γ_xy, γ_xz = 0), `D_matrix.py` (shear rows zero).
