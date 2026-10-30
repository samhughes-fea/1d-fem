# Element frames: beams in any direction, section twist, curved and pretwisted beams

The linear Euler-Bernoulli, Timoshenko and Benscoter beams can lie in any direction,
carry a section twist about their axis, and be curved and pretwisted. Code:
`pre_processing/element_library/element_frame.py` (geometry) and
`pre_processing/element_library/framed_element.py` (the element kernels).

## Inputs

- **Node coordinates** (`grid.txt`) set each element's axis `x`, from node 1 to node 2.
- **`[twist]`** (optional last column of `element.txt`, degrees) turns the section about
  `x`, right handed. With zero twist the section's `y` axis is horizontal,
  `y = Z × x / |Z × x|`, so its `z` axis points up. For a vertical element `y = Y`.
  An element along `+X` with zero twist has the global axes, so existing jobs are
  unchanged.
- **`precurvature.txt`** (`[element_id] [k_x0] [k_y0] [k_z0]`) gives an element a constant
  reference twist rate `k_x0` and bending curvatures `k_y0`, `k_z0`, in its own axes. The
  element is then the arc of that curvature between its nodes, stress free.

## Straight elements in any direction

Each kernel runs in the element's own axes, and the results are rotated to global axes
with `T` (the node frame `R`, whose columns are the element axes, on each node's
translations and rotations; warping DOFs are scalars):

    K = Tᵀ K_local T,   M = Tᵀ M_local T,   F = Tᵀ F_local,   K_σ(U) = Tᵀ K_σ,local(T U) T

The cached Gauss-point operators are turned the same way (`B T`, `Tᵀ N`), so the result
pipelines get the element's own strains from global displacements, and the wall stress
recovery turns the element displacements to local axes first.

## Curved and pretwisted elements

The arc's frame is `R(s) = R_mid exp(s [k]×)` for arc length `s` from its middle. `R_mid`
puts the arc's chord along node 1 to node 2, with the section turned about the chord by
`[twist]`. The arc length follows from the chord length. For a helix of radius `R` and
helix angle `α`, with the default section orientation `y` points to the axis and
`k = [sin α cos α / R, 0, cos² α / R]`.

The element is `m = 8` straight sub-elements of its own type and section along the arc,
each in the frame at its middle. Their chain is assembled in global axes and its inner
nodes are condensed out:

    K = Kbb − Kbi Kii⁻¹ Kib

The mass and `K_σ` use the same static shapes, `Φ = [I; −Kii⁻¹ Kib]`. So the curved
element is exact for the beam theory of its type as `m` grows, it has no membrane or
shear locking, and it works for every section the straight element takes. Wall stresses
at `ξ` come from the sub-element there.

A twist rate also stretches the fibres away from the axis: a fibre at radius `r` lies on
a helix, and a change `φ'` of the twist rate stretches it by `k_x0 r² φ'`. The section
gets the coupling terms (`SectionAnalysis.with_pretwist`):

    D[ε, φ'] += k ∫ E r² dA,   D[κ_y, φ'] −= k ∫ E z r² dA,   D[κ_z, φ'] −= k ∫ E y r² dA
    D[φ', φ'] += k² ∫ E r⁴ dA

These are the extension-twist coupling of a pretwisted strip (Rosen 1980; Hodges 2006).
The integrals use the section's `axial_points`: the wall quadrature of a thin-walled
section, or for `section.txt` inputs the rectangle with the same `A`, `I_y` and `I_z`.

The curvature's effect on the section itself is left out: fibres on the inside of a bend
being shorter than those outside (the Winkler correction). That is a `(k h)²` effect,
small for a slender beam.

## Loads and limits

- On a mesh with any rotated or curved element, point loads must sit on nodes. They are
  global components, located by their `x`, `y`, `z` and applied once at the node.
- Distributed loads (`distributed_load.txt`) are global components per length of the
  beam axis, sampled along global `X`: each element interpolates them at points along its
  own geometry, turns them into its axes and spreads them with its work-conjugate fields
  (Hermite for Euler-Bernoulli, linear otherwise). A curved element does this on each of
  its sub-elements before condensing. Every point of the axis must lie inside the
  samples' `X` range, so a beam must not run backwards in `X` beyond them.
- Wall loads (`wall_load.txt`) are section loads, built in the element's own axes (each
  sub-element's, on a curved element) and turned to global.
- GBT, Reddy and Levinson beams take frames too. GBT mode amplitudes belong to the
  section, so they map to themselves under the frame rotation. The Reddy and Levinson
  slope DOFs are held at each node as the shear strains `γ_xy`, `γ_xz` in one frame for
  the node (that of its lowest-numbered element); each element recovers its own slopes as
  `v' = γ_xy + θ_z` and `w' = γ_xz − θ_y`, so slopes stay continuous through a kink.
- The nonlinear elements refuse rotated or curved geometry. On those elements
  `precurvature.txt` keeps its older meaning, a reference strain on a straight element.

## Checks (`tests/test_element_frames.py`)

| case | reference | result |
|---|---|---|
| Cantilever turned to any direction and twisted (EB, Timoshenko, Benscoter on a coupled channel); displacements, mass, `K_σ`, wall resultants, buckling load | the unturned beam, turned | equal to round-off |
| Quarter-circle cantilever, out-of-plane tip force: 48 straight elements, or 4 curved ones | Roark: `P R³ (π / 4EI + (3π/4 − 2) / GJ)` | within 0.2 % |
| Quarter circle, in-plane radial tip force, 4 curved elements | `π P R³ / 4EI` | within 0.5 % |
| Helical spring, 3 turns, α = 8°, 48 straight elements per turn or 8 curved | Timoshenko / Wahl: `2π n P R³ / cos α (cos² α / GJ + sin² α / EI)` | within 0.4 % |
| Pretwisted strip under tension (90° over its length) | Rosen's coupled extension and twist | within 0.1 %; it untwists |
| Curved cantilever through the job runner (`precurvature.txt`, tip point load at the node) | Roark | within 0.2 % |
| Reddy and Levinson cantilever turned and twisted (channel section) | the unturned beam, turned; slopes as shear strains | equal to round-off |
| Quarter-circle cantilever, Reddy and Levinson: 48 straight or 4 curved elements | Roark, as above | within 0.06 % |
| GBT box cantilever turned and twisted, 2 modes, through the job runner (`tests/test_frames_gbt.py`) | the unturned beam, turned; mode amplitudes unchanged | equal to round-off |
| Quarter-circle GBT box cantilever, 8 curved elements, through the job runner | 64 straight chords through the same arc points (each curved element is 8 condensed sub-elements) | equal to round-off, mode amplitudes included |
