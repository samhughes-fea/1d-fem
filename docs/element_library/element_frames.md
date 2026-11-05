# Element frames: beams in any direction, section twist, curved and pretwisted beams

The linear Euler-Bernoulli, Timoshenko and Benscoter beams can lie in any direction,
carry a section twist about their axis, and be curved and pretwisted. Code:
`src/onedfem/pre_processing/element_library/element_frame.py` (geometry) and
`src/onedfem/pre_processing/element_library/framed_element.py` (the element kernels).

## Inputs

- **Node coordinates** (`grid.txt`) set each element's axis `x`, from node 1 to node 2.
- **`[twist]`** (optional last column of `element.txt`, degrees) turns the section about
  `x`, right handed. With zero twist the section's `y` axis is horizontal,
  `y = Z × x / |Z × x|`, so its `z` axis points up. For a vertical element `y = Y`.
  An element along `+X` with zero twist has the global axes, so existing jobs are
  unchanged.
- **`precurvature.txt`** (`[element_id] [k_x0] [k_y0] [k_z0]`) gives an element a constant
  reference twist rate `k_x0` and bending curvatures `k_y0`, `k_z0`, in its own axes. The
  element is then the arc of that curvature between its nodes, stress free. An optional
  fifth column `[winkler]` (1 or 0, default 1) turns the curvature correction of the
  element's section on or off (see below).

### Blade terms

"Precurvature" in 1d-fem means the whole reference curvature `[k_x0, k_y0, k_z0]`. In
blade terms that is three separate things. In the axes of a blade job (`X` span, `Y`
chord towards the trailing edge, `Z` flap with the suction side positive, see
[JOB_INPUT_THIN_WALLED_SECTIONS.md](../conventions/JOB_INPUT_THIN_WALLED_SECTIONS.md#blade-beam-jobs-from-implicit-geometry-kernel)),
an element along the span has `y = Y` and `z = Z`, so:

| Blade term | What it is | 1d-fem quantity | How `blade_beam_job` builds it |
|---|---|---|---|
| Twist | the section's angle about the span | `[twist]` column of `element.txt` (degrees) | not used: the cuts are already twisted |
| Pretwist | the rate of twist along the span | `k_x0` (rad/m) | `precurvature.txt`, one rate per element |
| Prebend | out-of-plane (flapwise) curvature of the reference axis | `k_y0` (1/m); `k_y0 > 0` bends the axis towards `−Z`, the pressure side | node offsets, a chain of straight elements |
| Presweep | in-plane (edgewise) curvature of the reference axis | `k_z0` (1/m); `k_z0 > 0` bends the axis towards `+Y`, the trailing edge | node offsets, a chain of straight elements |

The signs follow from the frame `R(s) = R_mid exp(s [k]×)`: the axis turns as
`dx/ds = k × x`. A prebent or swept blade can be given either way: as node offsets on
straight elements (what `blade_beam_job` does), or as curved elements with `k_y0`,
`k_z0` in `precurvature.txt`. Nodes on the curve give the same geometry in both, and
the curved elements' chords then run between those nodes.

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

The curvature also acts on the section itself: fibres on the inside of the bend are
shorter than those outside. A fibre at `(y, z)` has length `g = 1 − z k_y − y k_z` per length
of the axis and strain `a·e / g`, with `a = (1, −z, −y)` the fibre strain operator of
`(ε, κ_y, κ_z)`. Its energy per length of the axis is `½ ∫ E (a·e)² / g dA`, so the
sub-elements' sections take (`section_curvature.with_curvature`):

    D[i, j] += ∫ E a_i a_j (1/g − 1) dA,    i, j ∈ (ε, κ_y, κ_z)

This is the Winkler-Bach curved beam. To first order it couples extension and bending by
`−k EI`; otherwise it is a `(k h)²` change, small for a slender beam. On a bar of depth
`h = R/2` under a tip moment the beam turns as Winkler-Bach's `M L / (E A e R)`
(`e` the centroid to neutral-axis distance), 1.4 % more than `M L / EI`
(`tests/test_winkler_curvature.py`). The shear, torsion and warping terms and the mass
keep their straight-beam values. The integrals use the section's `axial_points`: the wall
quadrature of a thin-walled section, the 2D section's quadrature, or for `section.txt`
inputs a 5 × 5 Gauss rule on the rectangle with the same `A`, `I_y` and `I_z` (the
provider still stores that rectangle on 3 × 3 points, which is exact for the pretwist
polynomials and is what pretwist keeps using). A curvature
that puts the centre of curvature inside the section is refused. `[winkler] 0` leaves the
correction out, so a curved element is again exactly its chain of straight sub-elements.

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
- The corotational Euler-Bernoulli, Timoshenko and Benscoter beams take rotated and
  twisted elements, and kinked lines of them: each element's frame is its undeformed
  frame `E_0`, so a curved, prebent or swept axis is given as nodes along it and the
  kinks between the straight elements carry its curvature, stress-free. In
  `precurvature.txt` they take only the pretwist rate `k_x0`: the element stays straight,
  with its section frame at mid-length and the pretwist coupling in its section. A
  bending curvature (`k_y0`, `k_z0`) is refused. One corotational element is one straight
  element, while a linear pretwisted element is condensed from 8 sub-elements, so a
  corotational blade needs a finer mesh for the same accuracy.
- The other nonlinear elements refuse rotated or curved geometry. On those elements
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
| Corotational cantilever turned and twisted (EB, Timoshenko, Benscoter), tip deflection up to 0.3 L (`tests/test_corotational_frames.py`) | the unturned path, turned | to 1e-8 L |
| Quarter circle of 48 chords under the tip moment `EI / R` (corotational EB, Timoshenko) | exact: a semicircle of radius `R / 2`, tip at `(0, R)` | within 0.2 % of `R` |
| Prebent corotational cantilever, `P L² / EI = 1`, prebend 0.04, 0.02, 0.01 L | the straight path (Bisshopp and Drucker) | the gap halves with the prebend |
| Bent, swept, twisted corotational line at rest | the linear stiffness of the same mesh | equal to round-off |
| Prebent, swept, pretwisted blade with corotational Timoshenko elements through the nonlinear runner (`tests/test_blade_beam_job.py`) | the linear blade at small load | within 1.3 % (finer linear sub-elements); under a large flap load the tip runs in |
