# Saint-Venant problem for thin-walled sections

This page derives the equations solved by `ThinWalledSectionAnalysis`. Notation is as in [../../element_library/thin_walled_section_analysis.md](../../element_library/thin_walled_section_analysis.md).

## 1. Displacement field

The section is rigid in its own plane. Points on the mid-surface of a wall, at \((y, z)\), move by

\[
u_x = u - y\theta_z + z\theta_y + w(x, s),\qquad u_y = v - z\theta_x,\qquad u_z = w_z + y\theta_x,
\]

where \(w(x,s)\) is the warping. A point a distance \(n\) off the mid-surface, along \(\mathbf n\), moves with the same rigid part at its own position, plus \(n\beta\), the wall's through-thickness slope. Kirchhoff walls have \(\gamma_{xn} = 0\), which fixes \(\beta = -\partial_x u_n\) with \(u_n = \mathbf u\cdot\mathbf n\). Removing the rigid rotation part leaves

\[
\beta_K = -(n_y\gamma_{xy} + n_z\gamma_{xz} + r_n\phi).
\]

## 2. Wall strains

Differentiating gives the four wall strains of the element-library page.

- \(\varepsilon_{xx}\) and \(\kappa_{xx}\) are the mid-surface value and the through-thickness gradient of \(\partial_x u_x\).
- \(\gamma_{xs}\) and \(\kappa_{xs}\) are those of \(\partial_s u_x + \partial_x u_s\).

On a straight wall \(\mathbf n\) is constant along \(s\), so \(\partial_s r_n = t_y n_z - t_z n_y = 1\). The twisting curvature is then \(\kappa_{xs} = -\phi + \partial_s\beta_K = -2\phi\).

In matrix form, with nodal warping \(\hat{\mathbf w}\),

\[
\mathbf e_w = \mathbf Z\boldsymbol\psi + \mathbf Y\boldsymbol\psi' + \mathbf B\hat{\mathbf w} + \mathbf N\hat{\mathbf w}'.
\]

- \(\mathbf B\) puts \(\partial_s w\) into \(\gamma_{xs}\).
- \(\mathbf N\) puts \(w\) into \(\varepsilon_{xx}\), multiplying \(\hat{\mathbf w}'\).
- \(\mathbf Y\) holds the \(-(n_y, n_z, r_n)\) entries of \(\kappa_{xx}\).

## 3. Rigid modes

The warping fields \(w = 1\), \(w = y\) and \(w = z\) duplicate rigid motions. Take \(\hat{\mathbf w} = y\,a(x)\) together with \(\Delta\gamma_{xy} = -a\) and \(\Delta\kappa_z = a'\).

- \(\gamma_{xs}\): the warping adds \(t_y a\) and \(\Delta\gamma_{xy}\) adds \(-t_y a\), so the change is zero.
- \(\varepsilon_{xx}\): the warping adds \(ya'\) and \(\Delta\kappa_z\) adds \(-ya'\), so the change is zero.
- \(\kappa_{xx}\): \(\Delta\kappa_z\) adds \(-n_y a'\) and \(\mathbf Y\) adds \(+n_y a'\), so the change is zero.

So the strain is unchanged. These three modes must be removed from the warping. The code requires \(\int t\, w\, m\, ds = 0\) for each mode \(m\in\{1, y, z\}\), written \(\mathbf D^\top\hat{\mathbf w} = 0\) with multipliers \(\boldsymbol\lambda\). This defines the section's rotation as the thickness-weighted average rotation.

## 4. Equilibrium

The strain energy per length is \(\tfrac12\int\mathbf e_w^\top\mathbf K\mathbf e_w\,ds\).

**Warping equation.** Vary \(\hat{\mathbf w}(x)\) and integrate by parts along \(x\):

\[
\int\mathbf B^\top\mathbf K\mathbf e_w\,ds - \frac{d}{dx}\int\mathbf N^\top\mathbf K\mathbf e_w\,ds + \mathbf D\boldsymbol\lambda = 0 .
\]

**Section forces.** Vary \(\boldsymbol\psi(x)\) in the internal virtual work. The result is \(\int\mathbf F^\top\delta\boldsymbol\psi\,dx\), with

\[
\mathbf F = \int\mathbf Z^\top\mathbf K\mathbf e_w\,ds - \frac{d}{dx}\int\mathbf Y^\top\mathbf K\mathbf e_w\,ds .
\]

The second term is the moment of the walls' Kirchhoff transverse shear \(Q_x = \partial_x M_{xx}\). These \(\mathbf F\) are the stress resultants of the table in the element-library page.

**Beam equilibrium** with no distributed load gives \(\mathbf F' = \mathbf T\mathbf F\), where the only non-zero entries of \(\mathbf T\) are \(M_y' = V_z\) and \(M_z' = -V_y\).

## 5. Linear-in-\(x\) solution

Seek \(\hat{\mathbf w} = \hat{\mathbf w}_0 + x\hat{\mathbf w}_1\) and \(\boldsymbol\psi = \boldsymbol\psi_0 + x\boldsymbol\psi_1\). Write

- \(\mathbf E = \int\mathbf B^\top\mathbf K\mathbf B\), \(\mathbf C = \int\mathbf B^\top\mathbf K\mathbf N\), \(\mathbf R = \int\mathbf B^\top\mathbf K\mathbf Z\), \(\mathbf L = \int\mathbf N^\top\mathbf K\mathbf Z\),
- \(\mathbf L_y = \int\mathbf B^\top\mathbf K\mathbf Y\), \(\mathbf A = \int\mathbf Z^\top\mathbf K\mathbf Z\) and \(\mathbf Q = \int\mathbf Z^\top\mathbf K\mathbf Y\).

The strain is \(\mathbf e_w = \mathbf Z(\boldsymbol\psi_0 + x\boldsymbol\psi_1) + \mathbf Y\boldsymbol\psi_1 + \mathbf B(\hat{\mathbf w}_0 + x\hat{\mathbf w}_1) + \mathbf N\hat{\mathbf w}_1\), and \(\partial_x\mathbf e_w = \mathbf Z\boldsymbol\psi_1 + \mathbf B\hat{\mathbf w}_1\).

**Terms in \(x\).** These give the warping equation \(\mathbf E\hat{\mathbf w}_1 + \mathbf R\boldsymbol\psi_1 + \mathbf D\boldsymbol\lambda_1 = 0\), together with \(\mathbf R^\top\hat{\mathbf w}_1 + \mathbf A\boldsymbol\psi_1 = \mathbf F_1 = \mathbf T\mathbf F_0\).

**Constant terms.** In the warping equation, \(\frac{d}{dx}\int\mathbf N^\top\mathbf K\mathbf e_w = \mathbf C^\top\hat{\mathbf w}_1 + \mathbf L\boldsymbol\psi_1\), so

\[
\mathbf E\hat{\mathbf w}_0 + \mathbf R\boldsymbol\psi_0 + \mathbf D\boldsymbol\lambda_0 = -(\mathbf C - \mathbf C^\top)\hat{\mathbf w}_1 + (\mathbf L - \mathbf L_y)\boldsymbol\psi_1 .
\]

In the forces, \(\frac{d}{dx}\int\mathbf Y^\top\mathbf K\mathbf e_w = \mathbf Q^\top\boldsymbol\psi_1 + \mathbf L_y^\top\hat{\mathbf w}_1\), so

\[
\mathbf R^\top\hat{\mathbf w}_0 + \mathbf A\boldsymbol\psi_0 = \mathbf F_0 - (\mathbf L - \mathbf L_y)^\top\hat{\mathbf w}_1 - (\mathbf Q - \mathbf Q^\top)\boldsymbol\psi_1 .
\]

Both levels use the same symmetric matrix \(\begin{bmatrix}\mathbf E & \mathbf R & \mathbf D\\ \mathbf R^\top & \mathbf A & 0\\ \mathbf D^\top & 0 & 0\end{bmatrix}\). It is non-singular because, once the rigid modes are removed, the only null vectors of \(\begin{bmatrix}\mathbf E & \mathbf R\\ \mathbf R^\top & \mathbf A\end{bmatrix}\) are gone. These are the pairs of §3, together with \(w=1\), which lies in the null space of \(\mathbf E\).

## 6. Flexibility from energy

Solve for each unit column of \(\mathbf F_0\). This gives \(\mathbf e_w(0, s) = \boldsymbol\Phi(s)\mathbf F_0\). The strain energy per unit length at \(x = 0\) is \(\tfrac12\mathbf F_0^\top\left(\int\boldsymbol\Phi^\top\mathbf K\boldsymbol\Phi\,ds\right)\mathbf F_0\), so

\[
\mathbf F_{\text{flex}} = \int\boldsymbol\Phi^\top\mathbf K\boldsymbol\Phi\,ds,\qquad \mathbf K_{\text{sec}} = \mathbf F_{\text{flex}}^{-1}.
\]

This is symmetric and positive definite by construction. For forces with no shear (\(\mathbf T\mathbf F_0 = 0\)), the first level gives zero, and the second level is the energy minimisation of the classical problem. So `classical_stiffness`, the inverse of the \([N, M_y, M_z, T]\) block of \(\mathbf F_{\text{flex}}\), is the shear-rigid Saint-Venant stiffness.

## 7. Exactness of one cubic element per wall

On a straight wall, the \(\boldsymbol\psi\) parts of the strain (\(\mathbf Z\) and \(\mathbf Y\)) are linear in \(s\). Wall equilibrium along \(x\) is \(\partial_s N_{xs} = -\partial_x N_{xx}\).

- **Level 1.** \(\partial_x N_{xx}\) is zero, so \(N_{xs}\) is constant. \(\gamma_{xs}\) is then linear (coupled to the linear \(\varepsilon_{xx}\)), so \(w_1\) is quadratic.
- **Level 0.** \(\partial_x N_{xx}\) comes from \(\mathbf Z\boldsymbol\psi_1 + \mathbf B\hat{\mathbf w}_1\), which is linear in \(s\), so \(N_{xs}\) is quadratic. \(\varepsilon_{xx}\) contains \(w_1\), which is quadratic, so \(\gamma_{xs}\) is quadratic and \(w_0\) is cubic.

Cubic Lagrange elements contain both. The integrands are polynomials of degree at most 6, which 5-point Gauss integrates exactly. Test: `test_one_cubic_element_per_wall_is_exact`.

## 8. Centres and warping constant

**Shear centre.** A force \((V_y, V_z)\) at \((y_s, z_s)\) adds the torque \(T = y_sV_z - z_sV_y\). Setting the twist \(\phi = F_{5,j}F_j = 0\) gives \(y_s = -F_{54}/F_{55}\) and \(z_s = F_{53}/F_{55}\).

**Tension centre.** An axial force at \((y_t, z_t)\) adds \(M_y = z_tN\) and \(M_z = -y_tN\). Setting \(\kappa_y = \kappa_z = 0\) gives two linear equations for \((y_t, z_t)\).

**Warping constant.** The torsional warping \(\omega\) solves the warping equation with \(\boldsymbol\psi = \mathbf e_\phi\).

- A varying twist rate \(\phi'\) gives the axial strains \(\varepsilon_{xx} = \omega\phi'\) and \(\kappa_{xx} = -r_n\phi'\). The second of these is secondary warping, from \(\mathbf Y\).
- The parts of this field equal to extension and bending (\([1, 0]\), \([z, n_z]\), \([-y, -n_y]\)) are removed by least squares in the axial wall energy. That places the pole at the shear centre for uncoupled walls.
- Then \(E\Gamma = \int\mathbf a^\top\mathbf K_{ax}\mathbf a\,ds\), with \(\mathbf K_{ax}\) the \([\varepsilon_{xx}, \kappa_{xx}]\) block of the wall law.

For an isotropic I-beam this gives \(E(t_fb^3h^2/24 + t_f^3b^3/72 + t_w^3h^3/144)\).

## References

- V. Giavotto, M. Borri, P. Mantegazza, G. Ghiringhelli et al., *Anisotropic beam theory and applications*, Computers & Structures 16 (1983) 403–413.
- J. P. Blasques, M. Stolpe, *Multi-material topology optimization of laminated composite beam cross sections*, Composite Structures 94 (2012) 3278–3289 (the BECAS formulation).
- V. Z. Vlasov, *Thin-walled elastic beams* (1961).
