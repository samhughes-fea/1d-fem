# Selective integration with a fully coupled section stiffness

The linear Timoshenko element interpolates every field linearly: \(u_x, u_y, u_z, \theta_x, \theta_y, \theta_z\). So within an element:

- \(\varepsilon_x\), \(\kappa_y\), \(\kappa_z\) and \(\phi_x\) are **constant**;
- \(\gamma_{xy} = u_y' - \theta_z\) and \(\gamma_{xz} = u_z' + \theta_y\) are **linear** in \(\xi\): a constant part plus \(\xi\) times a term that depends only on the nodal rotations.

Shear locking comes from the \(\xi\)-linear part of \(\gamma\) entering \(\int\gamma^2\). Integrating the shear block with one point keeps only the mean, \(\bar\gamma\).

## Which couplings matter

Write \(\gamma = \bar\gamma + \xi\hat\gamma\) and let \(c\) be any constant strain.

- **Shear with a constant strain.** For a shear-constant coupling \(D_{cs}\), \(\int_{-1}^{1} c\,D_{cs}(\bar\gamma + \xi\hat\gamma)\,d\xi = 2c\,D_{cs}\bar\gamma\). Any rule exact for linears, including the one-point rule, gives this same value. So the full rule and the reduced rule agree on these terms, and leaving them in the full-rule baseline cannot lock.
- **Shear with shear.** For \(D_{34}\) (\(\gamma_{xy}\)–\(\gamma_{xz}\)), \(\int\gamma_{xy}D_{34}\gamma_{xz}\) contains \(\xi^2\hat\gamma_{xy}\hat\gamma_{xz}\). The full rule integrates this to a non-zero value and the one-point rule to zero. This term locks exactly like the diagonal shear terms, so it must be integrated with the shear block.
- **Bending with bending.** \(D_{12}\) (\(\kappa_y\)–\(\kappa_z\)) multiplies two constants, so every rule gives the same value. It is grouped with the bending block for consistency.

## Scheme

`assemble_timoshenko_K0` forms

\[
\mathbf K_e = \mathbf K_{\text{full}}(\mathbf D) - \mathbf K_{\text{full}}(\mathbf D_{bb}) - \mathbf K_{\text{full}}(\mathbf D_{ss}) + \mathbf K_{\text{bend}}(\mathbf D_{bb}) + \mathbf K_{\text{shear}}(\mathbf D_{ss}),
\]

where \(\mathbf D_{bb} = \mathbf D[1{:}3, 1{:}3]\) and \(\mathbf D_{ss} = \mathbf D[3{:}5, 3{:}5]\) are the full 2 × 2 blocks, including \(D_{12}\) and \(D_{34}\). The subscripts give the Gauss rule. For a diagonal \(\mathbf D\) this is the previous diagonal-block scheme exactly.

## Checks

All in `tests/test_coupled_beam_elements.py`:

- coupled cantilevers converge to the exact solution of the same theory;
- a slender coupled cantilever (\(L/h \approx 270\)) does not lock;
- the element stiffness has exactly six zero eigenvalues.
