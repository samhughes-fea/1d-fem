# Timoshenko: Shear force from constitutive law

This note proves that **Timoshenko (T) elements can yield non-zero shear resultants** \(V_y\), \(V_z\) from the constitutive relation \(\sigma = D\,\varepsilon\) at Gauss points, because they include shear deformation and non-zero shear stiffness in \(D\).

For the contrasting Euler–Bernoulli case (zero shear from constitutive law), see [../euler_bernoulli/shear_force_from_constitutive.md](../euler_bernoulli/shear_force_from_constitutive.md).

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

Hence: if \(\gamma_{xy},\,\gamma_{xz}\) are not identically zero **and** \(D_{44},\,D_{55} > 0\), then \(V_y,\,V_z\) can be non-zero.

---

## 2. Timoshenko element

**Kinematic assumption:** plane sections remain plane but **not** necessarily normal to the neutral axis. Transverse shear deformation is included:

\[
\gamma_{xy} = \frac{\partial u_y}{\partial x} - \theta_z, \qquad \gamma_{xz} = \frac{\partial u_z}{\partial x} - \theta_y.
\]

**Strain–displacement matrix \(B\):**  
The Timoshenko \(B\) matrix has non-zero fourth and fifth rows that compute \(\gamma_{xy}\) and \(\gamma_{xz}\) from the nodal displacements and rotations. So in general, for a non-trivial \(\mathbf{u}_e\),

\[
\varepsilon = B\,\mathbf{u}_e \quad \Rightarrow \quad \gamma_{xy},\;\gamma_{xz} \;\text{ are not necessarily zero}.
\]

**Constitutive matrix \(D\):**  
The Timoshenko material matrix includes shear stiffness (shear correction factor \(\kappa\), area \(A\), shear modulus \(G\)):

\[
D_{44} = \kappa GA > 0, \qquad D_{55} = \kappa GA > 0.
\]

So

\[
V_y = D_{44}\,\gamma_{xy} = \kappa GA\,\gamma_{xy}, \qquad V_z = D_{55}\,\gamma_{xz} = \kappa GA\,\gamma_{xz}.
\]

Whenever \(\gamma_{xy} \neq 0\) or \(\gamma_{xz} \neq 0\) (e.g. under transverse load), \(V_y\) and/or \(V_z\) are non-zero.

**Conclusion:**  
The Timoshenko element **does produce shear force from \(\sigma = D\,\varepsilon\)** at Gauss points, because it has non-zero shear strains and non-zero shear stiffness in \(D\).

**References in code:**  
`pre_processing/element_library/timoshenko/utilities/B_matrix.py` (γ_xy, γ_xz from displacement/rotation), `D_matrix.py` (κGA on diagonal entries 3, 4).
